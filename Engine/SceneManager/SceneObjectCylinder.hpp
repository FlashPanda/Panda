#pragma once
#include "Core/Math/EFloat.hpp"
#include "SceneObjectInteraction.hpp"
#include "SceneObjectShape.hpp"

namespace Panda
{
    class SceneObjectCylinder : public SceneObjectShape
    {
    public:
        SceneObjectCylinder(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
            float radius, float zMin, float zMax, float phiMax)
            : SceneObjectShape(objectToWorld, worldToObject),
              m_Radius(radius), m_zMin(zMin), m_zMax(zMax), m_PhiMax(phiMax)
            {}

    public:
        Bounds3Df ObjectBound();
        Bounds3Df WorldBound();
        bool Intersect(const Ray& ray, float& tHit, SceneObjectSurfaceInteraction& isect, bool testAlphaTexture);
        bool IntersectP(const Ray& ray, bool testAlphaTexture);
        float Area();

    private:
        float m_Radius;
        float m_zMin;
        float m_zMax;
        float m_PhiMax;
    };
}