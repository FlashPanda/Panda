#include "SceneObjectDisk.hpp"

namespace Panda
{
	SceneObjectDisk::SceneObjectDisk(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
		float _height, float _radius, float _innerRadius, float _phiMax)
		: SceneObjectShape(objectToWorld, worldToObject),
		  m_Height(_height),
		  m_Radius(_radius),
		  m_InnerRadius(_innerRadius),
		  m_PhiMax(_phiMax)
	{}

    Bounds3Df SceneObjectDisk::ObjectBound()
    {
        return Bounds3Df (Vector3Df({-m_Radius, -m_Radius, m_Height}),
                          Vector3Df({m_Radius, m_Radius, m_Height}));
    }

    Bounds3Df SceneObjectDisk::WorldBound()
    {
        return Bounds3Df();
    }

    bool SceneObjectDisk::Intersect(const Ray& ray, float& tHit, SceneObjectSurfaceInteraction& isect, bool testAlphaTexture)
    {
        // TO BE UNDERSTOOD
        // ProfilePhase p(Prof::ShapeIntersect);

        // Transform _Ray_ to object space
        Vector3Df oErr, dErr;
        Ray tRay(ray);
        TransformRay(tRay, *m_WorldToObject, oErr, dErr);

        // Compute plane intersection for disk

        // Reject disk intersections for rays parallel to the disk's plane
        if (ray.d.data[2] == 0.f) return false;
        float tShapeHit = (m_Height - tRay.o.data[2]) / tRay.d.data[2];
        if (tShapeHit <= 0.f || tShapeHit >= tRay.tMax) return false;

        // See if hit point is inside disk radius and $\phimax$
        Vector3Df pHit = tRay(tShapeHit);
        float dist2 = pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1];
        if (dist2 > m_Radius * m_Radius || dist2  < m_InnerRadius * m_InnerRadius)
            return false;

        // Test disk $\phi$ value against $\phimax$
        float phi = std::atan2f(pHit.data[1], pHit.data[0]);
        if (phi < 0.f) phi += 2 * PI;
        if (phi > m_PhiMax) return false;

        // Find parametric representation of dis hit
        float u = phi / m_PhiMax;
        float rHit = std::sqrtf(dist2);
        float v = (m_Radius - rHit) / (m_Radius - m_InnerRadius);
        Vector3Df dpdu({-m_PhiMax * pHit.data[1], m_PhiMax * pHit.data[0], 0.f});
        Vector3Df dpdv = Vector3Df({pHit.data[0], pHit.data[1], 0.f}) * (m_InnerRadius - m_Radius) / rHit;
        Vector3Df dndu({0.f, 0.f, 0.f}), dndv({0.f, 0.f, 0.f});

        // Refine disk intersection point
        pHit.data[2] = m_Height;

        // Compute error bounds for disk intersection
        Vector3Df pError({0.f, 0.f, 0.f});

        // Initialize _SurfaceInteraction_ from parametric information
        isect = TransformInteraction(SceneObjectSurfaceInteraction(pHit, pError, Vector2Df({u, v}),
                                                            -tRay.d, dpdu, dpdv, dndu, dndv, this), *m_ObjectToWorld);

        // Update _tHit_ for quadric intersection
        tHit = (float)tShapeHit;

        return true;
    }

    bool SceneObjectDisk::IntersectP(const Ray& ray, bool testAlphaTexture)
    {
        // TO BE UNDERSTOOD
        // ProfilePhase p(Prof::ShapeIntersect);

        // Transform _Ray_ to object space
        Vector3Df oErr, dErr;
        Ray tRay(ray);
        TransformRay(tRay, *m_WorldToObject, oErr, dErr);

        // Compute plane intersection for disk

        // Reject disk intersections for rays parallel to the disk's plane
        if (ray.d.data[2] == 0.f) return false;
        float tShapeHit = (m_Height - tRay.o.data[2]) / tRay.d.data[2];
        if (tShapeHit <= 0.f || tShapeHit >= tRay.tMax) return false;

        // See if hit point is inside disk radius and $\phimax$
        Vector3Df pHit = tRay(tShapeHit);
        float dist2 = pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1];
        if (dist2 > m_Radius * m_Radius || dist2  < m_InnerRadius * m_InnerRadius)
            return false;

        // Test disk $\phi$ value against $\phimax$
        float phi = std::atan2f(pHit.data[1], pHit.data[0]);
        if (phi < 0.f) phi += 2 * PI;
        if (phi > m_PhiMax) return false;

        return true;
    }

    float SceneObjectDisk::Area()
    {
        return m_PhiMax * 0.5 * (m_Radius * m_Radius - m_InnerRadius * m_InnerRadius);
    }

}