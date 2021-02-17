#pragma once

#include <vector>
#include "Core/Math/EFloat.hpp"
#include "SceneObjectInteraction.hpp"
#include "SceneObjectShape.hpp"

namespace Panda
{
    class SceneObjectSphere : public SceneObjectShape
    {
    public:
        SceneObjectSphere(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
            float radius, float zMin, float zMax, float phiMax) : 
            SceneObjectShape(objectToWorld, worldToObject), m_Origin(Vector3Df({0.f, 0.f, 0.f}))
        {
            m_zMin = Clamp(std::min(zMin, zMax), -radius, radius);
            m_zMax = Clamp(std::max(zMin, zMax), -radius, radius);
            m_ThetaMax = std::acosf(Clamp(zMin / radius, -1.f, 1.f));
            m_ThetaMin = std::acosf(Clamp(zMax / radius, -1.f, 1.f));
            m_PhiMax = Radians(Clamp(phiMax, 0.f, 360.f));
        }
    
    public:
        Bounds3Df ObjectBound();
        Bounds3Df WorldBound();

        bool Intersect(const Ray& ray, float& tHit, SceneObjectSurfaceInteraction& isect,
            bool testAlphaTexture);
        bool IntersectP(const Ray& ray, bool testAlphaTexture);
        float Area();

    public:
        friend std::ostream& operator<<(std::ostream& out, const SceneObjectSphere& obj);

	protected:
		Vector3Df m_Origin;
		float m_Radius;
		float m_zMin, m_zMax;
		float m_ThetaMin, m_ThetaMax, m_PhiMax;
    };
}