#pragma once

#include "BaseSceneNode.hpp"
#include "Core/Math/Ray.hpp"

namespace Panda
{
    struct KdAccelNode;
    struct BoundEdge;

    class SceneNodeKdTreeAccel : public BaseSceneNode
    {
    public:
        SceneNodeKdTreeAccel(std::vector<std::shared_ptr<BaseSceneNode>> p,
            int32_t isectCost = 80, int32_t traversalCost = 1, 
            float emptyBonus = 0.5f, int32_t maxNodes = 1, int32_t maxDepth = -1);
        ~SceneNodeKdTreeAccel();
        Bounds3Df WorldBound() { return m_Bounds;}
        bool Intersect(const Ray& ray, SceneObjectSurfaceInteraction& isect);
        bool IntersectP(const Ray& ray);

    private:
        void BuildTree(int32_t nodeNum, const Bounds3Df& bounds,
            const std::vector<Bounds3Df>& nodeBounds, int32_t* nodeNums,
            int32_t nNodes, int32_t depth,
            const std::unique_ptr<BoundEdge[]> edges[3], int32_t* nodes0,
            int32_t* nodes1, int32_t badRefines = 0);

        const int32_t m_IsectCost, m_TraversalCost, m_MaxNodes;
        const float m_EmptyBonus;
        std::vector<std::shared_ptr<BaseSceneNode>> m_SceneNodes;
        std::vector<int32_t> m_NodeIndices;
        KdAccelNode* m_Nodes;
        int32_t m_nAllocedNodes, m_NextFreeNode;
        Bounds3Df m_Bounds;
    };
}