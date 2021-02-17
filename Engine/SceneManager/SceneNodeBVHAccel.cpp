#include "SceneNodeBVHAccel.hpp"
#include "Core/Math/Tree.hpp"
#include "Core/Parallel.hpp"


namespace Panda
{
    struct BucketInfo
    {
        int32_t count = 0;
        Bounds3Df bounds;
    };

    struct MortonNode
    {
        int32_t nodeIndex;
        uint32_t mortonCode;
    };

    struct BVHNodeInfo
    {
        BVHNodeInfo() {}
        BVHNodeInfo(size_t nodeNum, const Bounds3Df& bounds)
            : nodeNum(nodeNum), bounds(bounds),
              centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
        
        size_t nodeNum;
        Bounds3Df bounds;
        Vector3Df centroid;
    };

    struct BVHBuildNode : public TreeNode
    {
        void InitLeaf(int32_t first, int32_t n, const Bounds3Df& b)
        {
            firstNodeOffset = first;
            nodeCount = n;
            bounds = b;
            // two children
            // If the children are null, then it is a leaf.
            m_Children.push_back(std::make_shared<BVHBuildNode>());
            m_Children.push_back(std::make_shared<BVHBuildNode>());
        }

        void InitInterior (int32_t axis, 
            std::shared_ptr<BVHBuildNode> c0, std::shared_ptr<BVHBuildNode> c1)
        {
            // two children
            if (m_Children.size() == 0)
            {
                m_Children.push_back(c0);
                m_Children.push_back(c1);
            }
            else if (m_Children.size() == 1)
            {
                m_Children[0] = c0;
                m_Children.push_back(c1);
            }
            else
            {
                m_Children[0] = c0;
                m_Children[1] = c1;
            }

            bounds = Union(c0->bounds, c1->bounds);
            splitAxis = axis;
            nodeCount = 0;
        }

        void InitInterior (int32_t axis, BVHBuildNode* c0, BVHBuildNode* c1)
        {
            // two children
            if (m_Children.size() == 0)
            {
                m_Children.push_back(std::make_shared<BVHBuildNode>(c0));
                m_Children.push_back(std::make_shared<BVHBuildNode>(c1));
            }
            else if (m_Children.size() == 1)
            {
                m_Children[0] = std::make_shared<BVHBuildNode>(c0);
                m_Children.push_back(std::make_shared<BVHBuildNode>(c1));
            }
            else
            {
                m_Children[0] = std::make_shared<BVHBuildNode>(c0);
                m_Children[1] = std::make_shared<BVHBuildNode>(c1);
            }

            bounds = Union(c0->bounds, c1->bounds);
            splitAxis = axis;
            nodeCount = 0;
        }

        Bounds3Df bounds;
        int32_t splitAxis;  // The axis splitted according to 
        int32_t firstNodeOffset;    // offset in the array
        int32_t nodeCount;  // the cound of nodes contained by _buildnode_
    };

    // Node cluster
    struct LBVHTreelet
    {
        int32_t startIndex, nNodes;
        BVHBuildNode* buildNodes;
    };

    struct LinearBVHNode
    {
        Bounds3Df bounds;
        union {
            int32_t nodesOffset;         // leaf
            int32_t secondChildOffset;  // interior 
        };
        uint16_t nNodes;    // 0 -> interior node
        uint8_t axis;       // interior node: xyz
        uint8_t pad[1];     // ensure 32 byte total size
    };

    inline uint32_t LeftShift3(uint32_t x)
    {
        if ( x == (1 << 10)) --x;

        x = (x | (x << 16)) & 0b00000011000000000000000011111111;
        // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x | (x << 8)) & 0b00000011000000001111000000001111;
        // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x | (x << 4)) & 0b00000011000011000011000011000011;
        // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x | (x << 2)) & 0b00001001001001001001001001001001;
        // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        return x;
    }

    inline uint32_t EncodeMorton3(const Vector3Df& v)
    {
        return (LeftShift3(v[2]) << 2) | (LeftShift3(v[1]) << 1) || LeftShift3(v[0]);
    }

    static void RadixSort(std::vector<MortonNode>& v)
    {
        std::vector<MortonNode> temp(v.size());
        constexpr int32_t bitsPerPass = 8;
        constexpr int32_t nBits = 32;   // Actually, we only use 30 bits of the value.
                                        // Set nBits to 32 is just for convinience.
        constexpr int32_t nPasses = nBits / bitsPerPass;

        for (int32_t pass = 0; pass < nPasses; ++pass)
        {
            // Perform one pass of radix sort, sorting _bitsPerPass_ bits
            int32_t lowBit = pass * bitsPerPass;

            // Set in and out vector references for radix sort pass
            std::vector<MortonNode>& in = (pass & 1)? temp : v;
            std::vector<MortonNode>& out = (pass & 1)? v : temp;

            // Count number of zero bits in array for current radix sort bit
            constexpr int32_t nBuckets = 1 << bitsPerPass;
            int32_t bucketCount[nBuckets] = {0};
            constexpr int32_t bitMask = (1 << bitsPerPass) - 1;
            for (const MortonNode& mn : in)
            {
                int32_t bucket = (mn.mortonCode >> lowBit) & bitMask;
                ++bucketCount[bucket];
            }

            // Compute starting index in output array for each bucket
            int32_t outIndex[nBuckets];
            outIndex[0] = 0;
            for (int32_t i = 1; i < nBuckets; ++i)
                outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

            // Store sorted values in output array
            for (const MortonNode& mn : in)
            {
                int32_t bucket = (mn.mortonCode >> lowBit) & bitMask;
                out[outIndex[bucket]++] = mn;
            }
        }

    }

    SceneNodeBVHAccel::SceneNodeBVHAccel(std::vector<std::shared_ptr<BaseSceneNode>> p, 
            int32_t maxPrimsInNode = 1, SplitMethod splitMethod = SplitMethod::SAH)
            : m_MaxPrimsInNode((std::min)(maxPrimsInNode, 255)), 
            m_SceneNodes(std::move(p)), m_SplitMethod(splitMethod)
    {
        if (p.empty()) return;

        // Build BVH from _scenenodes_

        // Initialize _NodeInfo_ array for nodes
        std::vector<BVHNodeInfo> nodeInfo(m_SceneNodes.size());
        for (size_t i = 0; i < m_SceneNodes.size(); ++i)
            nodeInfo[i] = {i, m_SceneNodes[i]->WorldBound()};

        // Build BVH tree for nodes using _nodeInfo_
        int32_t totalNodes = 0;
        std::vector<std::shared_ptr<BaseSceneNode>> orderedNodes;
        orderedNodes.reserve(m_SceneNodes.size());
        BVHBuildNode* root;
        MemoryLocal memLocal(1024 * 1024);
        if (m_SplitMethod == SplitMethod::HLBVH)
            root = HLBVHBuild(memLocal, nodeInfo, totalNodes, orderedNodes);
        else 
            root = RecursiveBuild(memLocal, nodeInfo, 0, m_SceneNodes.size(), 
                    totalNodes, orderedNodes);

        m_SceneNodes.swap(orderedNodes);
        nodeInfo.resize(0);

        m_TotalNodes = totalNodes;
        m_Nodes = (LinearBVHNode*)g_pMemoryManager->Allocate(sizeof(LinearBVHNode) * totalNodes, PANDA_L1_CACHE_LINE_SIZE);
        int32_t offset = 0;
        FlattenBVHTree(root, offset);
    }

    SceneNodeBVHAccel::~SceneNodeBVHAccel()
    {
        if (m_Nodes)
            g_pMemoryManager->Free(m_Nodes, sizeof(LinearBVHNode) * m_TotalNodes);
    }

    BVHBuildNode* SceneNodeBVHAccel::RecursiveBuild(MemoryLocal& memLocal, std::vector<BVHNodeInfo>& nodeInfo,
        int32_t start, int32_t end, int32_t& totalNodes,
        std::vector<std::shared_ptr<BaseSceneNode>>& orderedNodes)
    {
        BVHBuildNode* node = memLocal.Alloc<BVHBuildNode>();
        totalNodes++;

        // Compute bounds of all primitives in BVH node
        // Contruct a big bound using all bounds we met.
        Bounds3Df bounds;
        for (int32_t i = start; i < end; ++i)
            bounds = Union(bounds, nodeInfo[i].bounds);

        int32_t nodeCount = end - start;
        if (nodeCount == 1)
        {
            // Create leaf _BVHBuildNode_
            int32_t firstNodeOffset = orderedNodes.size();
            for (int32_t i = start; i < end; ++i)
            {
                int32_t nodeNum = nodeInfo[i].nodeNum;
                orderedNodes.push_back(m_SceneNodes[nodeNum]);
            }
            node->InitLeaf(firstNodeOffset, nodeCount, bounds);
            return node;
        }
        else 
        {
            // Compute bound of node centroids, choose split dimension _dim_
            Bounds3Df centroidBounds;
            for (int32_t i = start; i < end; ++i)
                centroidBounds = Union(centroidBounds, nodeInfo[i].centroid);
                int32_t dim = centroidBounds.MaximumExtent();

            // Partition nodes into two sets and build children
            int32_t mid = (start + end) / 2;
            if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim])
            {
                // Create leaf _BVHBuildNode_
                int32_t firstNodeOffset = orderedNodes.size();
                for (int32_t i = start; i < end; ++i)
                {
                    int32_t nodeNum = nodeInfo[i].nodeNum;
                    orderedNodes.push_back(m_SceneNodes[nodeNum]);
                }
                node->InitLeaf(firstNodeOffset, nodeCount, bounds);
                return node;
            }
            else 
            {
                // Partition nodes based on _splitMethod_
                switch (m_SplitMethod)
                {
                case SplitMethod::Middle:
                {
                    // Partition nodes through node's midpoint
                    float pmid = 
                        (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                    BVHNodeInfo* midPtr = std::partition(
                        &nodeInfo[start], &nodeInfo[end - 1] + 1,
                        [dim, pmid](const BVHNodeInfo& pi) {
                            return pi.centroid[dim] < pmid;
                        });
                    mid = midPtr - &nodeInfo[0];
                    // For lots of nodes with large overlapping bounding boxes, this may fail
                    // to partition; in that case don't break and fail through
                    // to EqualCounts.
                    if (mid != start && mid != end) break;
                }
                case SplitMethod::EqualCounts:
                {
                    // Partition nodes into equally-sized subsets
                    mid = (start + end) / 2;
                    std::nth_element (&nodeInfo[start], &nodeInfo[mid],
                        &nodeInfo[end - 1] + 1,
                        [dim](const BVHNodeInfo& a, const BVHNodeInfo& b) {
                            return a.centroid[dim] < b.centroid[dim];
                        });
                    break;
                }
                break;
                case SplitMethod::SAH:
                default:
                {
                    // Partition nodes using approximate SAH
                    if (nodeCount <= 2)
                    {
                        // Partition nodes into equally-sized subsets
                        mid = (start + end) / 2;
                        std::nth_element(&nodeInfo[start], &nodeInfo[mid],
                            &nodeInfo[end - 1] + 1,
                            [dim](const BVHNodeInfo& a, const BVHNodeInfo& b)
                            {
                                return a.centroid[dim] < b.centroid[dim];
                            });
                    }
                    else 
                    {
                        // Allocate _BucketInfo_ for SAH partition buckets
                        constexpr int32_t nBuckets = 12;
                        BucketInfo buckets[nBuckets];

                        // Initailize _BucketInfo_ for SAH partition buckets
                        for (int32_t i = start; i < end; ++i)
                        {
                            int32_t b = nBuckets * centroidBounds.Offset(nodeInfo[i].centroid)[dim];
                            if (b == nBuckets) b = nBuckets - 1;
                            buckets[b].count++;
                            buckets[b].bounds = Union(buckets[b].bounds, nodeInfo[i].bounds);
                        }

                        // Compute costs for splitting after each bucket
                        float cost[nBuckets -1] = {0};
                        for (int32_t i = 0; i < nBuckets - 1; ++i)
                        {
                            Bounds3Df b0, b1;
                            int32_t count0 = 0, count1 = 0;
                            for (int32_t j = 0; j <= i; ++j)
                            {
                                b0 = Union(b0, buckets[j].bounds);
                                count0 += buckets[j].count;
                            }
                            for (int32_t j = i + 1; j < nBuckets; ++j)
                            {
                                b1 = Union(b1, buckets[j].bounds);
                                count1 += buckets[j].count;
                            }
                            cost[i] = 1 + (count0 * b0.SurfaceAreas() + count1 * b1.SurfaceAreas()) / bounds.SurfaceAreas();
                        }

                        // Find bucket to split at that minimizes SAH metric
                        float minCost = cost[0];
                        int32_t minCostSplitBucket = 0;
                        for (int32_t i = 1; i < nBuckets - 1; ++i)
                        {
                            if (cost[i] < minCost)
                            {
                                minCost = cost[i];
                                minCostSplitBucket = i;
                            }
                        }

                        // Either create leaf or split nodes at selected SAH
                        float leafCost = nodeCount;
                        if (nodeCount > m_MaxPrimsInNode || minCost < leafCost)
                        {
                            BVHNodeInfo* pmid = std::partition(
                                &nodeInfo[start], &nodeInfo[end - 1] + 1,
                                [=](const BVHNodeInfo& pi) {
                                    int32_t b = nBuckets * 
                                        centroidBounds.Offset(pi.centroid)[dim];
                                    if (b == nBuckets) b = nBuckets - 1;
                                    return b <= minCostSplitBucket;
                                });
                            mid = pmid - &nodeInfo[0];
                        }
                        else 
                        {
                            // Create leaf _BVHBuildNode_
                            int32_t firstNodeOffset = orderedNodes.size();
                            for (int32_t i = start; i < end; ++i)
                            {
                                int32_t nodeNum = nodeInfo[i].nodeNum;
                                orderedNodes.push_back(m_SceneNodes[nodeNum]);
                            }
                            node->InitLeaf(firstNodeOffset, nodeCount, bounds);
                            return node;
                        }
                    }
                }
                break;
                }

                node ->InitInterior(dim,
                    RecursiveBuild(memLocal, nodeInfo, start, mid, totalNodes, orderedNodes),
                    RecursiveBuild(memLocal, nodeInfo, mid, end, totalNodes, orderedNodes));
            }
        }

        return node;
    }

    BVHBuildNode* SceneNodeBVHAccel::HLBVHBuild(MemoryLocal& memLocal, const std::vector<BVHNodeInfo>& nodeInfo,
        int32_t& totalNodes,
        std::vector<std::shared_ptr<BaseSceneNode>>& orderedNodes)
    {
        // Compute bounding box of all primitive centroids
        // Create a bounding box which can contain all centroids.
        Bounds3Df bounds;
        for(const BVHNodeInfo& pi : nodeInfo)
            bounds = Union(bounds, pi.centroid);

        // Compute Morton indices for all nodes
        std::vector<MortonNode> mortonNodes(nodeInfo.size());
        ParallelFor([&](int32_t i){
            // Initialize _mortonNodes[i]_ for _i_th node
            constexpr int32_t mortonBits = 10;
            constexpr int32_t mortonScale = 1 << mortonBits;
            mortonNodes[i].nodeIndex = nodeInfo[i].nodeNum;
            Vector3Df centroidOffset = bounds.Offset(nodeInfo[i].centroid);
            mortonNodes[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
        }, nodeInfo.size(), 512);

        // Radix sort nodes Morton indices
        RadixSort(mortonNodes);

        // Create LBVH treelets at bottom of BVH

        // Find intervals of primitives for each treelet
        // Here we'll find sets of primitives that have the same values for the 
        // high 12 bits of their 30-bit Morton codes.
        // But why 12? That's may be just experience.
        std::vector<LBVHTreelet> treeletsToBuild;
        for (int32_t start = 0, end = 1; end <= (int32_t)mortonNodes.size(); ++end)
        {
            uint32_t mask = 0b00111111111111110000000000000000;
            if (end == (int32_t)mortonNodes.size() || 
                ((mortonNodes[start].mortonCode & mask) != 
                 (mortonNodes[end].mortonCode & mask)))
            {
                // Add entry to _treeletsToBuild_ for this treelet
                int32_t nNodes = end - start;
                // The number of nodes in a BVH is bounded by twice the number of
                // leaf nodes, which in ture is bounded by the number of nodes.
                // A full binary tree has 2^h - 1 nodes, h is the height of the tree.
                // That is to say it has 2^(h - 1) leaves, which means that 
                // 2^h - 1 < 2^h = 2 * 2^(h - 1) = 2 * leaves.
                int32_t maxBVHNodes = 2 * nNodes;
                BVHBuildNode* nodes = g_pMemoryManager->New<BVHBuildNode>(maxBVHNodes, false);
                treeletsToBuild.push_back({start, nNodes, nodes});
                start = end;
            }
        }

        // Create LBVHs for treelets in parallel
        std::atomic<int32_t> atomicTotal(0), orderedNodesOffset(0);
        orderedNodes.resize(m_SceneNodes.size());
        ParallelFor([&](int32_t i) {
            // Generate _i_th LBVH treelet
            int32_t nodesCreated = 0;
            const int32_t firstBitIndex = 29 - 12;
            LBVHTreelet& tr = treeletsToBuild[i];
            tr.buildNodes = EmitLBVH(tr.buildNodes, nodeInfo, &mortonNodes[tr.startIndex],
                                     tr.nNodes, nodesCreated, orderedNodes,
                                     orderedNodesOffset, firstBitIndex);
            atomicTotal += nodesCreated;
        }, treeletsToBuild.size());
        totalNodes = atomicTotal;

        // Create and return SAH BVH from LBVH treelets
        std::vector<BVHBuildNode*> finishedTreelets;
        finishedTreelets.reserve(treeletsToBuild.size());
        for (LBVHTreelet& treelet : treeletsToBuild)
            finishedTreelets.push_back(treelet.buildNodes);
        return BuildUpperSAH(memLocal, finishedTreelets, 0, finishedTreelets.size(), totalNodes);

    }

    BVHBuildNode* SceneNodeBVHAccel::EmitLBVH(BVHBuildNode*& buildNodes,
        const std::vector<BVHNodeInfo>& nodeInfo,
        MortonNode* mortonNodes, int32_t nNodes, int32_t& totalNodes,
        std::vector<std::shared_ptr<BaseSceneNode>>& orderedNodes,
        std::atomic<int32_t>& orderedNodesOffset, int32_t bitIndex)
    {
        if (bitIndex == -1 || nNodes < m_MaxPrimsInNode)
        {
            totalNodes++;
            BVHBuildNode* node = buildNodes++;  // I can't understand why we need '++' here!
            Bounds3Df bounds;
            int32_t firstNodeOffset = orderedNodesOffset.fetch_add(nNodes); // _orderedNodesOffset_ will return old value.
            for (int32_t i = 0; i < nNodes; ++i)
            {
                int32_t nodeIndex = mortonNodes[i].nodeIndex;
                orderedNodes[firstNodeOffset + i] = m_SceneNodes[nodeIndex];
                bounds = Union(bounds, nodeInfo[nodeIndex].bounds);
            }
            node->InitLeaf(firstNodeOffset, nNodes, bounds);
            return node;
        }
        else 
        {
            int32_t mask = 1 << bitIndex;
            // Advance to new subtree level if there's no LBVH split for this bit
            if ((mortonNodes[0].mortonCode & mask) ==
                (mortonNodes[nNodes - 1].mortonCode & mask))
                return EmitLBVH(buildNodes, nodeInfo, mortonNodes, nNodes,
                                totalNodes, orderedNodes, orderedNodesOffset, 
                                bitIndex - 1);

            // Find LBVH split point for this dimension
            int32_t searchStart = 0, searchEnd = nNodes - 1;
            while (searchStart + 1 != searchEnd)
            {
                int32_t mid = (searchStart + searchEnd) / 2;
                if ((mortonNodes[searchStart].mortonCode & mask) ==
                    (mortonNodes[mid].mortonCode & mask))
                    searchStart = mid;
                else 
                    searchEnd = mid;
            }
            int32_t splitOffset = searchEnd;

            // Create and return interior LBVH node
            totalNodes++;
            BVHBuildNode* node = buildNodes++;
            BVHBuildNode* lbvh[2] = {
                EmitLBVH(buildNodes, nodeInfo, mortonNodes, splitOffset, 
                         totalNodes, orderedNodes, orderedNodesOffset, bitIndex - 1),
                EmitLBVH(buildNodes, nodeInfo, mortonNodes, nNodes - splitOffset, 
                         totalNodes, orderedNodes, orderedNodesOffset, bitIndex - 1)
            };
            int32_t axis = bitIndex % 3;
            node->InitInterior(axis, std::make_shared<BVHBuildNode>(lbvh[0]), std::make_shared<BVHBuildNode>(lbvh[0]));
            return node;
        }
    }

    BVHBuildNode* SceneNodeBVHAccel::BuildUpperSAH(MemoryLocal& memLocal, std::vector<BVHBuildNode *>& treeletRoots,
        int32_t start, int32_t end, int32_t& totalNodes)
    {
        int32_t nTreeNodes = end - start;
        if (nTreeNodes == 1) return treeletRoots[0];
        totalNodes++;
        BVHBuildNode* treeNode = memLocal.Alloc<BVHBuildNode>();

        // Compute bounds of all nodes under this HLBVH node
        Bounds3Df bounds;
        for (int32_t i = start; i < end; ++i)
            bounds = Union(bounds, treeletRoots[i]->bounds);

        // Compute bound of HLBVH node centroids, choose split dimesnion _dim_
        Bounds3Df centroidBounds;
        for (int32_t i = start; i < end; ++i)
        {
            Vector3Df centroid =
                (treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) * 0.5f;
            centroidBounds = Union(centroidBounds, centroid);
        }
        int32_t dim = centroidBounds.MaximumExtent();

        // Allocate _BucketInfo_ for SAH partition buckets
        constexpr int32_t nBuckets = 12;
        struct BucketInfo
        {
            int32_t count = 0;
            Bounds3Df bounds;
        };
        BucketInfo buckets[nBuckets];

        // Initialize _BucketInfo_ for HLBVH SAH partition buckets
        for (int32_t i = start; i < end; ++i)
        {
            float centroid = (treeletRoots[i]->bounds.pMin[dim] + 
                              treeletRoots[i]->bounds.pMax[dim]) * 0.5f;
            int32_t b = nBuckets * ((centroid - centroidBounds.pMin[dim]) /
                                  (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
            if (b == nBuckets) b = nBuckets - 1;
            buckets[b].count++;
            buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
        }

        // Compute costs for splitting after each bucket
        float cost[nBuckets - 1];
        for (int32_t i = 0; i < nBuckets - 1; ++i)
        {
            Bounds3Df b0, b1;
            int32_t count0 = 0, count1 = 0;
            for (int32_t j = 0; j <= i; ++j)
            {
                b0 = Union(b0, buckets[j].bounds);
                count0 += buckets[j].count;
            }
            for (int32_t j = i + 1; j < nBuckets; ++j)
            {
                b1 = Union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }
            cost[i] = .125f + (count0 * b0.SurfaceAreas() + count1 * b1.SurfaceAreas()) / bounds.SurfaceAreas();
        }

        // Find bucket to split at that minimizes SAH metric
        float minCost = cost[0];
        int32_t minCostSplitBucket = 0;
        for (int32_t i = 1; i < nBuckets - 1; ++i)
        {
            if (cost[i] < minCost)
            {
                minCost = cost[i];
                minCostSplitBucket = i;
            }
        }

        // Split nodes and create interior HLBVH SAH node
        BVHBuildNode** pmid = std::partition(
            &treeletRoots[start], &treeletRoots[end - 1] + 1,
            [=](const BVHBuildNode* node)
            {
                float centroid = (node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
                int32_t b = nBuckets * 
                            ((centroid - centroidBounds.pMin[dim]) / 
                             (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                if (b == nBuckets) b = nBuckets - 1;
                return b <= minCostSplitBucket;
            });
        int32_t mid = pmid - &treeletRoots[0];
        treeNode->InitInterior(dim, this->BuildUpperSAH(memLocal, treeletRoots, start, mid, totalNodes),
            this->BuildUpperSAH(memLocal, treeletRoots, mid, end, totalNodes));
        return treeNode;
    }

    int32_t SceneNodeBVHAccel::FlattenBVHTree(BVHBuildNode* node, int32_t& offset)
    {
        LinearBVHNode* linearNode = &m_Nodes[offset];
        linearNode->bounds = node->bounds;
        int32_t myOffset = offset++;
        if (node->nodeCount > 0)
        {
            linearNode->nodesOffset = node->firstNodeOffset;
            linearNode->nNodes = node->nodeCount;
        }
        else 
        {
            // Create interior flattened BVH node
            linearNode->axis = node->splitAxis;
            linearNode->nNodes = 0;
            FlattenBVHTree(node->m_Children[0].get(), offset);
            linearNode->secondChildOffset =
                FlattenBVHTree(node->m_Children[1].get(), offset);
        }
        return myOffset;
    }

    Bounds3Df SceneNodeBVHAccel::WorldBound()
    {
        return m_Nodes? m_Nodes[0].bounds : Bounds3Df();
    }

    bool SceneNodeBVHAccel::Intersect(const Ray& ray, SceneObjectSurfaceInteraction& isect)
    {
        if (!m_Nodes)
            return false;

        bool hit = false;
        Vector3Df invDir({1.0f / ray.d[0], 1.0f / ray.d[1], 1.0f / ray.d[2]});
        int32_t dirIsNeg[3] = {invDir[0] < 0, invDir[1] < 0, invDir[2] < 0};
        // Follow ray through BVH nodes to find primitive intersections
        int32_t toVisitOffset = 0, currentNodeIndex = 0;
        int32_t nodesToVisit[64];
        while (true)
        {
            const LinearBVHNode* node = &m_Nodes[currentNodeIndex];
            // Check ray against BVH node
            if (node->bounds.IntersectP(ray, invDir, dirIsNeg))
            {
                if (node->nNodes > 0)
                {   // This is a leaf.
                    // Intersect ray with nodes in leaf BVH node
                    for (int32_t i = 0; i < node->nNodes; ++i)
                        if (m_SceneNodes[node->nodesOffset + i]->Intersect(ray, isect))
                            hit = true;
                    if (toVisitOffset == 0) break;
                    currentNodeIndex = nodesToVisit[--toVisitOffset];
                }
                else 
                {   // This is a interior node
                    // Put far BVH node on _nodesToVisit_ stack, advance to near node
                    if (dirIsNeg[node->axis])
                    {
                        nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                        currentNodeIndex = node->secondChildOffset;
                    }
                    else 
                    {
                        nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                        currentNodeIndex = currentNodeIndex + 1;
                    }
                }
            }
            else {
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            }
        }

        return hit;
    }

    bool SceneNodeBVHAccel::IntersectP(const Ray& ray)
    {
        if (!m_Nodes) return false;
        Vector3Df invDir({1.f / ray.d[0], 1.f / ray.d[1], 1.f / ray.d[2]});
        int32_t dirIsNeg[3] = {invDir[0] < 0, invDir[1] < 0, invDir[2] < 0};
        int32_t nodesToVisit[64];
        int32_t toVisitOffset = 0, currentNodeIndex = 0;
        while (true)
        {
            const LinearBVHNode* node = &m_Nodes[currentNodeIndex];
            if (node->bounds.IntersectP(ray, invDir, dirIsNeg))
            {   // Process BVH node _node_ for traversal
                if (node->nNodes > 0)
                {
                    for (int32_t i = 0; i < node->nNodes; ++i)
                        if (m_SceneNodes[node->nodesOffset + i]->IntersectP(ray))
                            return true;
                    if (toVisitOffset == 0) break;
                    currentNodeIndex = nodesToVisit[--toVisitOffset];
                }
                else 
                {
                    if (dirIsNeg[node->axis])
                    {   // second child first
                        nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                        currentNodeIndex = node->secondChildOffset;
                    }
                    else 
                    {
                        nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                        currentNodeIndex = currentNodeIndex + 1;
                    }
                }
            }
            else 
            {
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            }
        }
        return false;
    }
}