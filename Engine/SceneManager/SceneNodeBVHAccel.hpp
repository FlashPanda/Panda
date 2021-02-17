#pragma once
#include "BaseSceneNode.hpp"
#include "Core/MemoryManager.hpp"
#include "Core/Math/Ray.hpp"

namespace Panda
{
    struct BVHBuildNode;
    struct BVHNodeInfo;
    struct MortonNode;
    struct LinearBVHNode;

    // Accel is a tree structure.
    class SceneNodeBVHAccel : public BaseSceneNode
    {
    public:
        // Split methods
        enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

        SceneNodeBVHAccel(std::vector<std::shared_ptr<BaseSceneNode>> p, 
            int32_t maxPrimsInNode = 1, SplitMethod splitMethod = SplitMethod::SAH);
        ~SceneNodeBVHAccel();
        Bounds3Df WorldBound();
        bool Intersect(const Ray& ray, SceneObjectSurfaceInteraction& isect);
        bool IntersectP(const Ray& ray);

    private:
        BVHBuildNode* RecursiveBuild(MemoryLocal& memLocal, std::vector<BVHNodeInfo>& nodeInfo,
            int32_t start, int32_t end, int32_t& totalNodes,
            std::vector<std::shared_ptr<BaseSceneNode>>& orderedNodes);
        BVHBuildNode* HLBVHBuild(MemoryLocal& memLocal, const std::vector<BVHNodeInfo>& nodeInfo,
            int32_t& totalNodes,
            std::vector<std::shared_ptr<BaseSceneNode>>& orderedNodes);
        BVHBuildNode* EmitLBVH(BVHBuildNode*& buildNodes,
            const std::vector<BVHNodeInfo>& nodeInfo,
            MortonNode* mortonNodes, int32_t nNodes, int32_t& totalNodes,
            std::vector<std::shared_ptr<BaseSceneNode>>& orderedNodes,
            std::atomic<int32_t>& orderedNodesOffset, int32_t bitIndex);
        BVHBuildNode* BuildUpperSAH(MemoryLocal& memLocal, std::vector<BVHBuildNode *> &treeletRoots,
            int32_t start, int32_t end, int32_t& totalNodes);
        int32_t FlattenBVHTree(BVHBuildNode* node, int32_t& offset);

    private:
        const int32_t m_MaxPrimsInNode;
        const SplitMethod m_SplitMethod;
        std::vector<std::shared_ptr<BaseSceneNode>> m_SceneNodes;
        int32_t m_TotalNodes;
        LinearBVHNode* m_Nodes = nullptr;
    };
}