#include "SceneNodeKDTreeAccel.hpp"
#include <algorithm>
#include <cmath>
#include "SceneManager/SceneManager.hpp"

namespace Panda
{
    struct KdAccelNode
    {
    public:
        void InitLeaf(int32_t* nodeNums, int32_t np,
            std::vector<int32_t>& nodeIndices)
        {
            flags = 3;  // which means this is a leaf node.
            nNodes |= (np << 2);    // Stores the number of nodes, but we can't use the low-order 2 bits.

            // Store node ids for leaf node
            if (np == 0)
                oneNode = 0;
            else if (np == 1)
                oneNode = nodeNums[0];
            else 
            {
                nodeIndicesOffset = nodeIndices.size();
                for (int32_t i = 0; i < np; ++i)
                    nodeIndices.push_back(nodeNums[i]);
            }
        }
        void InitInterior(int32_t axis, int32_t ac, float s)
        {
            flags = axis;
            split = s;
            aboveChild != (ac << 2);
        }

        float SplitPos() const {return split;}
        int32_t NumOfNodes() const {return nNodes >> 2;}
        int32_t SplitAxis() const {return flags & 3;}
        bool IsLeaf() const {return (flags & 3) == 3;}
        int32_t AboveChild() const {return aboveChild >> 2;}

    public:
        union {
            float split;                // Interior
            int32_t oneNode;            // Leaf
            int32_t nodeIndicesOffset;  // Leaf
        };

    private:
        union {
            int32_t flags;      // Both
            int32_t nNodes;     // Leaf
            int32_t aboveChild; // Interior
        };
    };

    enum class EdgeType { Start, End };

    struct BoundEdge
    {
        BoundEdge() {}
        BoundEdge(float t, int32_t nodeNum, bool starting)
            : t(t), nodeNum(nodeNum)
        {
            type = starting? EdgeType::Start : EdgeType::End;
        }

        float t;
        int32_t nodeNum;
        EdgeType type;
    };

    struct KdToDo
    {
        const KdAccelNode* node;
        float tMin, tMax;
    };  

    SceneNodeKdTreeAccel::SceneNodeKdTreeAccel(std::vector<std::shared_ptr<BaseSceneNode>> p,
        int32_t isectCost, int32_t traversalCost, 
        float emptyBonus, int32_t maxNodes, int32_t maxDepth)
        : m_IsectCost(isectCost), m_TraversalCost(traversalCost), m_MaxNodes(maxNodes),
        m_EmptyBonus(emptyBonus), m_SceneNodes(std::move(p))
    {   // Build kd-tree for accelerator
        m_NextFreeNode = m_nAllocedNodes = 0;
        if (maxDepth <= 0)
            maxDepth = std::round(8 + 1.3f * std::log2f(m_SceneNodes.size()));

        // Compute bounds for kd-tree construction
        std::vector<Bounds3Df> nodeBounds;
        nodeBounds.reserve(m_SceneNodes.size());
        for (const std::shared_ptr<BaseSceneNode>& node : m_SceneNodes)
        {
            Bounds3Df b = node->WorldBound();
            m_Bounds = Union(m_Bounds, b);
            nodeBounds.push_back(b);
        }

        // Allocate working memory for kd-tree construction
        std::unique_ptr<BoundEdge[]> edges[3];
        for (int32_t i = 0; i < 3; ++i)
            edges[i].reset(new BoundEdge[2 * m_SceneNodes.size()]);
        std::unique_ptr<int32_t[]> nodes0(new int32_t[m_SceneNodes.size()]);
        std::unique_ptr<int32_t[]> nodes1(new int32_t[(maxDepth + 1) * m_SceneNodes.size()]);

        // Initialize _nodeNums_ for kd-tree construction
        std::unique_ptr<int32_t[]> nodeNums(new int32_t[m_SceneNodes.size()]);
        for (size_t i = 0; i < m_SceneNodes.size(); ++i)
            nodeNums[i] = i;

        // Start recursive construction of kd-tree
        BuildTree(0, m_Bounds, nodeBounds, nodeNums.get(), m_SceneNodes.size(),
            maxDepth, edges, nodes0.get(), nodes1.get());
    }

    void SceneNodeKdTreeAccel::BuildTree(int32_t nodeNum, const Bounds3Df& nodeBounds,
        const std::vector<Bounds3Df>& allNodeBounds, int32_t* nodeNums,
        int32_t nNodes, int32_t depth,
        const std::unique_ptr<BoundEdge[]> edges[3], int32_t* nodes0,
        int32_t* nodes1, int32_t badRefines)
    {
        // Get next free node from _nodes_ array
        if (m_NextFreeNode == m_nAllocedNodes)
        {
            int32_t nNewAllocNodes = (std::max)(2 * m_nAllocedNodes, 512);
            KdAccelNode* n = (KdAccelNode*)g_pMemoryManager->Allocate(sizeof(KdAccelNode) * nNewAllocNodes, PANDA_L1_CACHE_LINE_SIZE);
            if (m_nAllocedNodes > 0)
            {
                memcpy(n, m_Nodes, m_nAllocedNodes * sizeof(KdAccelNode));
                g_pMemoryManager->Free(m_Nodes, m_nAllocedNodes * sizeof(KdAccelNode));
            }
            m_Nodes = n;
            m_nAllocedNodes = nNewAllocNodes;
        }
        ++m_NextFreeNode;

        // Inilialize leaf node if termination criteria met
        if (nNodes <= m_MaxNodes || depth == 0)
        {
            m_Nodes[nodeNum].InitLeaf(nodeNums, nNodes, m_NodeIndices);
            return;
        }

        // Initialize interior node and continue recursion

        // Choose split axis position for interior node
        int32_t bestAxis = -1, bestOffset = -1;
        float bestCost = Infinity;
        float oldCost = m_IsectCost * float(nNodes);
        float totalSA = nodeBounds.SurfaceAreas();
        float invTotalSA = 1.f / totalSA;
        Vector3Df d = nodeBounds.Diagonal();

        // Choose which axis to split along
        int32_t axis = nodeBounds.MaximumExtent();
        int32_t retries = 0;
        while (bestAxis == -1 && retries < 2)
        {
            // Initialize edges for _axis_
            for (int32_t i = 0; i < nNodes; ++i)
            {
                int32_t pn = nodeNums[i];
                const Bounds3Df& bounds = allNodeBounds[pn];
                edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
                edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
            }

            // Sort _edges_ for _axis_
            std::sort(&edges[axis][0], &edges[axis][2 * nNodes],
                [](const BoundEdge& e0, const BoundEdge& e1)->bool
                {
                    if (e0.t == e1.t)
                        return (int32_t)e0.type < (int32_t)e1.type;
                    else 
                        return e0.t < e1.t;
                });

            ++retries;
            axis = (axis + 1) % 3;

            // Compute cost of all splits for _axis_ to find best
            int32_t nBelow = 0, nAbove = nNodes;
            for (int32_t i = 0; i < 2 * nNodes; ++i)
            {
                if (edges[axis][i].type == EdgeType::End) --nAbove;
                float edgeT = edges[axis][i].t;
                if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis])
                {
                    // Compute cost for split at _i_th edge

                    // Compute child surface areas for split at _edgeT_
                    int32_t otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
                    float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                         (edgeT - nodeBounds.pMin[axis]) * 
                                            (d[otherAxis0] + d[otherAxis1]));
                    float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] + 
                                         (nodeBounds.pMax[axis] - edgeT) *
                                            (d[otherAxis0] + d[otherAxis1]));
                    float pBelow = belowSA * invTotalSA;
                    float pAbove = aboveSA * invTotalSA;
                    float eb = (nAbove == 0 || nBelow == 0)? m_EmptyBonus : 0;
                    float cost = m_TraversalCost + m_IsectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

                    // Update best split if this is lowest cost so far
                    if (cost < bestCost)
                    {
                        bestCost = cost;
                        bestAxis = axis;
                        bestOffset = i;
                    }
                }

                if (edges[axis][i].type == EdgeType::Start) ++nBelow;   // The node in this iteration must become below for next iteration.
            }
            assert(nBelow == nNodes && nAbove == 0);
        }

        // Create leaf if no good splits were found
        if (bestCost > oldCost) ++badRefines;
        if ((bestCost > 4 * oldCost && nNodes < 16) || bestAxis == -1 |
            badRefines == 3)
        {
            m_Nodes[nodeNum].InitLeaf(nodeNums, nNodes, m_NodeIndices);
            return;
        }

        // Classify nodes with respect to split
        int32_t n0 = 0, n1 = 0;
        for (int32_t i = 0; i < bestOffset; ++i)
            if (edges[bestAxis][i].type == EdgeType::Start)
                nodes0[n0++] = edges[bestAxis][i].nodeNum;
        for (int32_t i = bestOffset + 1; i < 2 * nNodes; ++i)
            if (edges[bestAxis][i].type == EdgeType::End)
                nodes1[n1++] = edges[bestAxis][i].nodeNum;

        // Recursively initialize children nodes
        float tSplit = edges[bestAxis][bestOffset].t;
        Bounds3Df bounds0 = nodeBounds, bounds1 = nodeBounds;
        bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;
        BuildTree(nodeNum + 1, bounds0, allNodeBounds, nodes0, n0, depth -1,
            edges, nodes0, nodes1 + nNodes, badRefines);
        int32_t aboveChild = m_NextFreeNode;
        m_Nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
        BuildTree(aboveChild, bounds1, allNodeBounds, nodes1, n1, depth - 1,
            edges, nodes0, nodes1 + nNodes, badRefines);
    }

    bool SceneNodeKdTreeAccel::Intersect(const Ray& ray, SceneObjectSurfaceInteraction& isect)
    {
        // Compute initial parametric range of ray inside kd-tree extnt
        float tMin, tMax;
        if (!m_Bounds.IntersectP(ray, tMin, tMax))
            return false;

        // Prepare to traverse kd-tree for ray
        Vector3Df invDir({1.f / ray.d[0], 1.f / ray.d[1], 1.f / ray.d[2]});
        constexpr int32_t maxTodo = 64;
        KdToDo todo[maxTodo];
        int32_t todoPos = 0;

        // Traverse kd-tree nodes in order for ray
        bool hit = false;
        const KdAccelNode* node = &m_Nodes[0];
        while (node != nullptr)
        {
            // Bail out if we found a hit closer than the current node
            if (ray.tMax < tMin) break;

            if (!node->IsLeaf())
            {   // Process kd-tree interior node

                // Compute parametric distance along ray to split plane
                int32_t axis = node->SplitAxis();
                float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

                // Get node children pointers for ray
                const KdAccelNode* firstChild;
                const KdAccelNode* secondChild;
                bool belowFirst = (ray.o[axis] < node->SplitPos()) || (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0.f);
                if (belowFirst) 
                {
                    firstChild = node + 1;
                    secondChild = &m_Nodes[node->AboveChild()];
                }
                else 
                {
                    firstChild = &m_Nodes[node->AboveChild()];
                    secondChild = node + 1;
                }

                // Advance to next child node, possible enqueue other child
                if (tPlane > tMax || tPlane <= 0)
                    node = firstChild;
                else if (tPlane < tMin)
                    node = secondChild;
                else 
                {
                    // Enqueue _secondChild_ in _todo_ list
                    todo[todoPos].node = secondChild;
                    todo[todoPos].tMin = tPlane;
                    todo[todoPos].tMax = tMax;
                    ++todoPos;
                    node = firstChild;
                    tMax = tPlane;
                }
            }
            else 
            {   // Check for intersections inside leaf node
                int32_t nNodes = node->NumOfNodes();
                if (nNodes == 1)
                {
                    const std::shared_ptr<BaseSceneNode>& p = m_SceneNodes[node->oneNode];
                    // Check one node inside leaf node
                    if (p->Intersect(ray, isect)) hit = true;
                }
                else 
                {
                    for (int32_t i = 0; i < nNodes; ++i)
                    {
                        int32_t index = m_NodeIndices[node->nodeIndicesOffset + i];
                        const std::shared_ptr<BaseSceneNode>& p = m_SceneNodes[index];
                        // Check one node inside leaf node
                        if (p->Intersect(ray, isect)) hit = true;
                    }
                }

                // Grab next node to process from _todo_ list
                if (todoPos > 0)
                {
                    --todoPos;
                    node = todo[todoPos].node;
                    tMin = todo[todoPos].tMin;
                    tMax = todo[todoPos].tMax;
                }
                else break;
            }
        }

        return hit;
    }

    bool SceneNodeKdTreeAccel::IntersectP(const Ray& ray)
    {
        // Compute initial parametric range of ray inside kd-tree extent
        float tMin, tMax;
        if (!m_Bounds.IntersectP(ray, tMin, tMax))
            return false;

        // Prepare to traverse kd-tree for ray
        Vector3Df invDir({1.f / ray.d[0], 1.f / ray.d[1], 1.f / ray.d[2]});
        constexpr int32_t maxTodo = 64;
        KdToDo todo[maxTodo];
        int32_t todoPos = 0;
        const KdAccelNode* node = &m_Nodes[0];
        while(node != nullptr)
        {
            if (node->IsLeaf())
            {
                // Check for shadow ray intersections inside leaf node
                int32_t nNodes = node->NumOfNodes();
                if (nNodes == 1)
                {
                    const std::shared_ptr<BaseSceneNode>& p = 
                        m_SceneNodes[node->oneNode];
                    if (p->IntersectP(ray))
                        return true;
                }
                else 
                {
                    for (int32_t i = 0; i < nNodes; ++i)
                    {
                        int32_t nodeIndex = m_NodeIndices[node->nodeIndicesOffset + i];
                        const std::shared_ptr<BaseSceneNode>& p = m_SceneNodes[nodeIndex];
                        if (p->IntersectP(ray))
                            return true;
                    }
                }

                // Grab next node to process from todo list
                if (todoPos > 0)
                {
                    --todoPos;
                    node = todo[todoPos].node;
                    tMin = todo[todoPos].tMin;
                    tMax = todo[todoPos].tMax;
                }
                else 
                    break;
            }
            else 
            {   // Process kd-tree interior node
                // Compute parametric distance along ray to split plane
                int32_t axis = node->SplitAxis();
                float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

                // Get node children pointers for ray
                const KdAccelNode* firstChild;
                const KdAccelNode* secondChild;
                int32_t belowFirst = (ray.o[axis] < node->SplitPos()) || (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
                if (belowFirst) {
                    firstChild = node + 1;
                    secondChild = &m_Nodes[node->AboveChild()];
                }
                else {
                    firstChild = &m_Nodes[node->AboveChild()];
                    secondChild = node + 1;
                }

                // Advance to next child node, possibly enqueue other child
                if (tPlane > tMax || tPlane <= 0)
                    node = firstChild;
                else if (tPlane < tMin)
                    node = secondChild;
                else 
                {
                    // Enqueue _secondChild_ in _todo_ list
                    todo[todoPos].node = secondChild;
                    todo[todoPos].tMin = tPlane;
                    todo[todoPos].tMax = tMax;
                    ++todoPos;
                    node = firstChild;
                    tMax = tPlane;
                }
            }
        }

        return false;
    }
}