#include "BaseSceneNode.hpp"

namespace Panda
{
    // I think I don't need this structure. 
    // Try not to inherit from this class first, find if it cound work well.
    class SceneAggregateNode : public BaseSceneNode
    {
    public:
        virtual Bounds3Df WorldBound();
        virtual bool Intersect(const Ray& r, SceneObjectSurfaceInteraction& interaction);
        virtual bool IntersectP(const Ray& r);
    };
}