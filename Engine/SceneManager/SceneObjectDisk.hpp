#pragma once
#include "Core/Math/EFloat.hpp"
#include "SceneObjectInteraction.hpp"
#include "SceneObjectShape.hpp"

namespace Panda
{
    class SceneObjectDisk : public SceneObjectShape
    {
    public:
		SceneObjectDisk(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
			float _height, float _radius, float _innerRadius, float _phiMax); 

    public:
        Bounds3Df ObjectBound();
        Bounds3Df WorldBound();
		bool Intersect(const Ray& ray, float& tHit, SceneObjectSurfaceInteraction& isect, bool testAlphaTexture);
		bool IntersectP(const Ray& ray, bool testAlphaTexture);

		float Area();

    private:
        float m_Height;
        float m_Radius;
        float m_InnerRadius;
        float m_PhiMax;
    };
}