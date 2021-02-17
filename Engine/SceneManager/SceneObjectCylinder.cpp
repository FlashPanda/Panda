#include "SceneObjectCylinder.hpp"

namespace Panda
{
    Bounds3Df SceneObjectCylinder::ObjectBound()
    {
        return Bounds3Df(Vector3Df({-m_Radius, -m_Radius, m_zMin}), Vector3Df({m_Radius, m_Radius, m_zMax}));
    }
    
    Bounds3Df SceneObjectCylinder::WorldBound()
    {
        return Bounds3Df();
    }

    bool SceneObjectCylinder::Intersect(const Ray& ray, float& tHit, SceneObjectSurfaceInteraction& isect, bool testAlphaTexture)
    {
        // TO BE UNDERSTOOD
        // ProfilePhase p(Prof::ShapeIntersect);
        float phi;
        Vector3Df pHit;
        // Transform _Ray_ to object space
        Vector3Df oErr, dErr;
        Ray tRay(ray);
        TransformRay(tRay, *m_WorldToObject, oErr, dErr);

        // Compute quadratic cylinder coefficients
        
        // Initialize _EFloat_ ray coordinate values
        EFloat ox(tRay.o.data[0], oErr.data[0]), oy(tRay.o.data[1], oErr.data[1]), oz(tRay.o.data[2], oErr.data[2]);
        EFloat dx(tRay.d.data[0], dErr.data[0]), dy(tRay.d.data[1], dErr.data[1]), dz(tRay.d.data[2], dErr.data[2]);
        EFloat a = dx * dx + dy * dy;
        EFloat b = 2 * (dx * ox + dy * oy);
        EFloat c = ox * ox + oy * oy - EFloat(m_Radius) * EFloat(m_Radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, t0, t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > tRay.tMax || t1.LowerBound() < 0.f) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0.f)
        {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > tRay.tMax) return false;
        }

        // Compute cylinder hit point and $\phi$
        pHit = tRay((float)tShapeHit);

        // Refine cylinder intersection point
        float hitRad = std::sqrtf(pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1]);
        pHit.data[0] *= m_Radius / hitRad;
        pHit.data[1] *= m_Radius / hitRad;
        phi = std::atan2f(pHit.data[1], pHit.data[0]);
        if (phi < 0.f) phi += 2 * PI;

        // Test cylinder intersection against clipping parameters
        if (pHit.data[2] < m_zMin || pHit.data[2] > m_zMax || phi > m_PhiMax)
        {
            if (tShapeHit == t1) return false;
            tShapeHit = t1;
            if (t1.UpperBound() > ray.tMax) return false;
            // Compute cylinder hit point and $\phi$
            pHit = tRay((float)tShapeHit);

            // Refine cylinder intersection point
            float hitRad = std::sqrt(pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1]);
            pHit.data[0] *= m_Radius / hitRad;
            pHit.data[1] *= m_Radius / hitRad;
            phi = std::atan2f(pHit.data[1], pHit.data[0]);
            if (phi < 0) phi += 2 * PI;
            if (pHit.data[2] < m_zMin || pHit.data[2] > m_zMax || phi > m_PhiMax) return false;
        }

        // Find parametric representation of cylinder hit
        float u = phi / m_PhiMax;
        float v = (pHit.data[2] - m_zMin) / (m_zMax - m_zMin);

        // Compute cylinder $\dpdu$ and $\dpdv$
        Vector3Df dpdu ({-m_PhiMax * pHit.data[1], m_PhiMax * pHit.data[0], 0.f});
        Vector3Df dpdv ({0.f, 0.f, m_zMax - m_zMin});

        // Compute cylinder $\dndu$ and $\dndv$
        Vector3Df d2Pduu = -m_PhiMax * m_PhiMax * Vector3Df({pHit.data[0], pHit.data[1], 0.f});
        Vector3Df d2Pduv({0.f, 0.f, 0.f});
        Vector3Df d2Pdvv({0.f, 0.f, 0.f});

        // Compute coefficients for fundamental forms
        float E = DotProduct(dpdu, dpdu);
        float F = DotProduct(dpdu, dpdv);
        float G = DotProduct(dpdv, dpdv);
        Vector3Df N = Normalize(CrossProduct(dpdu, dpdv));
        float e = DotProduct(N, d2Pduu);
        float f = DotProduct(N, d2Pduv);
        float g = DotProduct(N, d2Pdvv);

        // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        float invEGF2 = 1.f / (E * G - F * F);
        Vector3Df dndu = (f * F - e * G) * invEGF2 * dpdu +
                            (e * F - f * E) * invEGF2 * dpdv;
        Vector3Df dndv = (g * F - f * G) * invEGF2 * dpdu +
                            (f * f - g * E) * invEGF2 * dpdv;

        // Compute error bounds for cylinder intersection
        Vector3Df pError = Gamma(3) * Abs(Vector3Df({pHit.data[0], pHit.data[1], 0}));

        // Initialize _SurfaceInteraction_ from 
        isect = TransformInteraction(SceneObjectSurfaceInteraction(pHit, pError, Vector2Df({u, v}),
                                                            -tRay.d, dpdu, dpdv, dndu, dndv, this), *m_ObjectToWorld);

        // Update _tHit_ for quadric intersection
        tHit = (float)tShapeHit;
        return true;
    }

    bool SceneObjectCylinder::IntersectP(const Ray& ray, bool testAlphaTexture)
    {
        // TO BE UNDERSTOOD
        // ProfilePhase p(Prof::ShapeIntersect);
        float phi;
        Vector3Df pHit;
        // Transform _Ray_ to object space
        Vector3Df oErr, dErr;
        Ray tRay(ray);
        TransformRay(tRay, *m_WorldToObject, oErr, dErr);

        // Compute quadratic cylinder coefficients
        
        // Initialize _EFloat_ ray coordinate values
        EFloat ox(tRay.o.data[0], oErr.data[0]), oy(tRay.o.data[1], oErr.data[1]), oz(tRay.o.data[2], oErr.data[2]);
        EFloat dx(tRay.d.data[0], dErr.data[0]), dy(tRay.d.data[1], dErr.data[1]), dz(tRay.d.data[2], dErr.data[2]);
        EFloat a = dx * dx + dy * dy;
        EFloat b = 2 * (dx * ox + dy * oy);
        EFloat c = ox * ox + oy * oy - EFloat(m_Radius) * EFloat(m_Radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, t0, t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > tRay.tMax || t1.LowerBound() < 0.f) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0.f)
        {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > tRay.tMax) return false;
        }

        // Compute cylinder hit point and $\phi$
        pHit = tRay((float)tShapeHit);

        // Refine cylinder intersection point
        float hitRad = std::sqrtf(pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1]);
        pHit.data[0] *= m_Radius / hitRad;
        pHit.data[1] *= m_Radius / hitRad;
        phi = std::atan2f(pHit.data[1], pHit.data[0]);
        if (phi < 0.f) phi += 2 * PI;

        // Test cylinder intersection against clipping parameters
        if (pHit.data[2] < m_zMin || pHit.data[2] > m_zMax || phi > m_PhiMax)
        {
            if (tShapeHit == t1) return false;
            tShapeHit = t1;
            if (t1.UpperBound() > ray.tMax) return false;
            // Compute cylinder hit point and $\phi$
            pHit = tRay((float)tShapeHit);

            // Refine cylinder intersection point
            float hitRad = std::sqrt(pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1]);
            pHit.data[0] *= m_Radius / hitRad;
            pHit.data[1] *= m_Radius / hitRad;
            phi = std::atan2f(pHit.data[1], pHit.data[0]);
            if (phi < 0) phi += 2 * PI;
            if (pHit.data[2] < m_zMin || pHit.data[2] > m_zMax || phi > m_PhiMax) return false;
        }

        return true;
    }

    float SceneObjectCylinder::Area() 
    {
        return (m_zMax - m_zMin) * m_Radius * m_PhiMax;
    }
}