#include "SceneObjectSphere.hpp"
#include "SceneObjectInteraction.hpp"

namespace Panda
{
    Bounds3Df SceneObjectSphere::ObjectBound()
    {
        return Bounds3Df (Vector3Df({-m_Radius, -m_Radius, m_zMin}),
            Vector3Df({m_Radius, m_Radius, m_zMax}));
    }
    Bounds3Df SceneObjectSphere::WorldBound()
    {
        return Bounds3Df();
    }

    bool SceneObjectSphere::Intersect(const Ray& ray, float& tHit, SceneObjectSurfaceInteraction& isect,
        bool testAlphaTexture)
    {
        // ProfilePhase p(Prof::ShapeIntersect);
        float phi;
        Vector3Df pHit;
        Ray tRay(ray);

        // Transform _tRay_ to object space
        Vector3Df oErr, dErr;
        TransformRay(tRay, *m_WorldToObject, oErr, dErr);

        /// Compute quadratic sphere coefficients
        // Initialize _EFloat_ ray coordinate values
        EFloat ox(tRay.o.data[0], oErr.data[0]), oy(tRay.o.data[1], oErr.data[1]), oz(tRay.o.data[2], oErr.data[2]);
        EFloat dx(tRay.o.data[0], oErr.data[0]), dy(tRay.o.data[1], oErr.data[1]), dz(tRay.o.data[2], oErr.data[2]);
        EFloat a(dx * dx + dy * dy + dz * dz);
        EFloat b(2 * (dx * ox + dy * oy + dz * oz));
        EFloat c(ox * ox + oy * oy + oz * oz - EFloat(m_Radius) * EFloat(m_Radius));

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, t0, t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > tRay.tMax || t1.LowerBound() <= 0.f) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0.f)
        {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > tRay.tMax) return false;
        }

        // Compute sphere hit position and $\phi$
        pHit = tRay((float)tShapeHit);

        // Refine sphere intersection point
        pHit *= m_Radius / GetLength(pHit);
        if (pHit[0] == 0.f && pHit[1] == 0.f) pHit[0] = 1e-5f * m_Radius;
        phi = std::atan2f(pHit[1], pHit[0]);
        if (phi < 0.f) phi += 2 * PI;

        // Test sphere intersection against clipping parameters
        if ((m_zMin > -m_Radius && pHit[2] < m_zMin) || (m_zMax < m_Radius && pHit[2] > m_zMax) ||
            phi > m_PhiMax)
        {
            /// Test t1
            if (tShapeHit == t1) return false;
            if (t1.UpperBound() > ray.tMax) return false;
            tShapeHit = t1;

            // Compute sphere hit positin and $\phi$
            pHit = tRay((float)tShapeHit);

            // Refine sphere intersection point
            pHit *= m_Radius / GetLength(pHit);
            if (pHit[0] == 0.f && pHit[1] == 0.f) pHit[0] = 1e-5f * m_Radius;
            phi = std::atan2f(pHit[1], pHit[0]);
            if (phi < 0.f) phi += 2 * PI;
            if ((m_zMin > -m_Radius && pHit[2] < m_zMin) || 
                (m_zMax < m_Radius && pHit[2] > m_zMax) ||
                phi > m_PhiMax)
                return false;
        }

        // Find parametric representation of sphere hit
        float u = phi / m_PhiMax;
        float theta = std::acosf(Clamp(pHit.data[2] / m_Radius, -1.f, 1.f));
        float v = (theta - m_ThetaMin) / (m_ThetaMax - m_ThetaMin);

        // Compute sphere $\dpdu$ and $\dpdv$
        float zRadius = std::sqrtf(pHit.data[0] * pHit.data[0] + pHit.data[1] * pHit.data[1]);
        float invZRadius = 1.f / zRadius;
        float cosPhi = pHit.data[0] * invZRadius;
        float sinPhi = pHit.data[1] * invZRadius;
        Vector3Df dpdu ({-m_PhiMax * pHit.data[1], m_PhiMax * pHit.data[0], 0});
        Vector3Df dpdv = (m_ThetaMax - m_ThetaMin) * 
                            Vector3Df({pHit.data[2] * cosPhi, pHit.data[2] * sinPhi, -m_Radius * std::sinf(theta)});
        
        // Compute sphere $\dndu$ and $\dndv$
        Vector3Df d2Pduu = -m_PhiMax * m_PhiMax * Vector3Df({pHit.data[0], pHit.data[1], 0});
        Vector3Df d2Pduv = (m_ThetaMax - m_ThetaMin) * pHit.data[2] * m_PhiMax * Vector3Df({-sinPhi, cosPhi, 0.f});
        Vector3Df d2Pdvv = -(m_ThetaMax - m_ThetaMin) * (m_ThetaMax - m_ThetaMin) * Vector3Df({pHit.data[0], pHit.data[1], pHit.data[2]});

        // Compute coefficients for first and second fundamental forms.
        // You can refer to 'Differential Geometry of Curves and Surfaces' ebook to 
        // try to understand the reason to do so.
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
                            (f * F - g * E) * invEGF2 * dpdv;

        // Compute error bounds for sphere intersection
        Vector3Df pError = Gamma(5) * Abs(pHit);

        // Initialize _SurfaceInteraction_ from parametric information
        isect = TransformInteraction(SceneObjectSurfaceInteraction(pHit, pError, Vector2Df({u, v}),
                                                            -tRay.d, dpdu, dpdv, dndu, dndv, this), *m_ObjectToWorld);

        // Update _tHit_ for quadric intersection
        tHit = (float)tShapeHit;
        return true;
    }
    bool SceneObjectSphere::IntersectP(const Ray& ray, bool testAlphaTexture)
    {
        // ProfilePhase p(Prof::ShapeIntersectP);
        float phi;
        Vector3Df pHit;
        // Transform _Ray_ to object scene
        Ray tRay(ray);
        Vector3Df oErr, dErr;
        TransformRay(tRay, *m_WorldToObject, oErr, dErr);

        // Compute quadratic sphere coefficients

        // Initialzie _EFloat_ ray coordinate values
        EFloat ox(tRay.o.data[0], oErr.data[0]), oy(tRay.o.data[1], oErr.data[1]), oz(tRay.o.data[2], oErr.data[2]);
        EFloat dx(tRay.d.data[0], dErr.data[0]), dy(tRay.d.data[1], dErr.data[1]), dz(tRay.d.data[2], dErr.data[2]);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2.f * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(m_Radius) * EFloat(m_Radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if(!Quadratic(a, b, c, t0, t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > tRay.tMax || t1.LowerBound() < 0.f) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0.f)
        {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > tRay.tMax) return false;
        }

        // Compute sphere hit position and $\phi$
        pHit = tRay((float)tShapeHit);

        // Refine sphere intersection point
        pHit *= m_Radius / GetLength(pHit);
        if (pHit.data[0] == 0.f && pHit.data[1] == 0.f)
            pHit.data[0] = 1e-5f * m_Radius;
        phi = std::atan2(pHit.data[1], pHit.data[0]);
        if (phi < 0.f) phi += 2 * PI;

        // Test sphere intersection against clipping parameters
        if ((m_zMin > -m_Radius && pHit.data[2] < m_zMin) || (m_zMax < m_Radius && pHit.data[2] > m_zMax) || phi > m_PhiMax)
        {
            if (tShapeHit == t1) return false;
            if (t1.UpperBound() > ray.tMax) return false;
            tShapeHit = t1;
            // Compute sphere hit poisition and $\phi$
            pHit = tRay((float)tShapeHit);

            // Refine sphere intersection point
            pHit *= m_Radius / GetLength(pHit);
            if (pHit.data[0] == 0.f && pHit.data[1] == 0.f) pHit.data[0] = 1e-5f * m_Radius;
            phi = std::atan2f(pHit.data[1], pHit.data[0]);
            if (phi < 0.f) phi += 2 * PI;
            if ((m_zMin > -m_Radius && pHit.data[2] < m_zMin) || (m_zMax < m_Radius && pHit.data[2] > m_zMax) || phi > m_PhiMax)
                return false;
        }

        return true;
    }
    float SceneObjectSphere::Area()
    {
        return m_PhiMax * m_Radius * (m_zMax - m_zMin);
    }
}