#pragma once

#include "BaseSceneObject.hpp"
#include "SceneObjectTransform.hpp"
#include "Math/EFloat.hpp"
#include "NumericalExt.hpp"
#include "Shape.hpp"
#include "SceneObjectTypeDef.hpp"
//#include "Math/Ray.hpp"

namespace Panda
{
    static Vector3Df BlossomBezier(const Vector3Df p[4], float u0, float u1, float u2)
    {
        Vector3Df a[3] = {Lerp(u0, p[0], p[1]), Lerp(u0, p[1], p[2]),
                          Lerp(u0, p[2], p[3])};
        Vector3Df b[2] = {Lerp(u1, a[0], a[1]), Lerp(u1, a[1], a[2])};
        return Lerp(u2, b[0], b[1]);
    }

    void SubdivideBezier(const Vector3Df cp[4], Vector3Df cpSplict[7])
    {
        cpSplict[0] = cp[0];
        cpSplict[1] = (cp[0] + cp[1]) / 2;
        cpSplict[2] = (cp[0] + 2 * cp[1] + cp[2]) / 4;
        cpSplict[3] = (cp[0] + 3 * cp[1] + 3 * cp[2] + cp[3]) / 8;
        cpSplict[4] = (cp[1] + 2 * cp[2] + cp[3]) / 4;
        cpSplict[5] = (cp[2] + cp[3]) / 2;
        cpSplict[6] = cp[3];
    }

    static Vector3Df EvalBezier(const Vector3Df cp[4], float u, Vector3Df* deriv = nullptr)
    {
        Vector3Df cp1[3] = {Lerp(u, cp[0], cp[1]), Lerp(u, cp[1], cp[2]),
                            Lerp(u, cp[2], cp[3])};
        Vector3Df cp2[3] = {Lerp(u, cp1[0], cp1[1]), Lerp(u, cp1[1], cp1[2])};
        if (deriv)
        {
            if (GetLengthSquare(cp2[1] - cp2[0]) > 0)
                *deriv = 3 * (cp2[1] - cp2[0]);
            else 
            {
                // FOr a cubic Bezier, if the first three control points (say) are
                // coincident, then the derivative of the curve is legitimately (0, 0, 0)
                // at u = 0. This is problematic for us, though, since we'd like to be
                // able to compute a surface normal there. In that case, just punt and 
                // take the difference between the first and last control points, which
                // ain't great. but will hopefully do.
                *deriv = cp[3] - cp[0];
            }
        }
        return Lerp(u, cp2[0], cp2[1]);
    }

    struct CurveCommon
    {
        CurveCommon(const Vector3Df c[4], float w0, float w1, CurveObjectType type, const Vector3Df* normal)
            : type(type)
        {
            width[0] = w0; width[1] = w1;
            for (int32_t i = 0; i < 4; ++i)
                cpObj[i] = c[i];
            if (normal)
            {
                n[0] = Normalize(normal[0]);
                n[1] = Normalize(normal[1]);
                normalAngle = std::acosf(Clamp(DotProduct(n[0], n[1]), 0.f, 1.f));
                invSinNormalAngle = 1.f / std::sinf(normalAngle);
            }
        }

        const CurveObjectType type;
        Vector3Df cpObj[4];
        float width[2];
        Vector3Df n[2];
        float normalAngle, invSinNormalAngle;
    };

    class SceneObjectCurve : public BaseSceneObject, public Shape
    {
	protected:
		std::shared_ptr<SceneObjectTransform> m_ObjectToWorld;
        std::shared_ptr<SceneObjectTransform> m_WorldToObject;

    public:
        SceneObjectCurve(std::shared_ptr<SceneObjectTransform> ObjectToWorld, std::shared_ptr<SceneObjectTransform> WorldToObject,
                        const CurveCommon& common, float uMin, float uMax)
                        : BaseSceneObject(SceneObjectType::kSceneObjectTypeMath),
                          Shape(),
                          m_CurveCommon(common), m_uMin(uMin), m_uMax(uMax) {}

        Bounds3Df ObjectBound()
        {
            // Compute object-space control points for curve segment, _cpOvj_
            Vector3Df cpObj[4];
            cpObj[0] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMin, m_uMin);
            cpObj[1] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMin, m_uMax);
            cpObj[2] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMax, m_uMax);
            cpObj[3] = BlossomBezier(m_CurveCommon.cpObj, m_uMax, m_uMax, m_uMax);

            Bounds3Df b = Union(Bounds3Df(cpObj[0], cpObj[1]), Bounds3Df(cpObj[2], cpObj[3]));
            float width[2] = {Lerp(m_uMin, m_CurveCommon.width[0], m_CurveCommon.width[1]),
                              Lerp(m_uMax, m_CurveCommon.width[0], m_CurveCommon.width[1])};
            return Expand(b, (std::max)(width[0], width[1] * 0.5f));
        }

        bool Intersect(const Ray& r, float& tHit, SurfaceInteraction& isect, bool testAlphaTexture)
        {
            // TO BE UNDERSTOOD
            // ProfilePhase p(Prof::CurveIntersect);

            // Transform _Ray_ to object space
            Vector3Df oErr, dErr;
            Ray ray = TransformRay(r, *m_WorldToObject, oErr, dErr);

            // Compute object-space control points for curve segment, _cpObj_
            Vector3Df cpObj[4];
            cpObj[0] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMin, m_uMin);
            cpObj[1] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMin, m_uMax);
            cpObj[2] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMax, m_uMax);
            cpObj[2] = BlossomBezier(m_CurveCommon.cpObj, m_uMax, m_uMax, m_uMax);

            // Project curve control points to plane perpendicular to ray

            // Be careful to set the "up" direction passed to equal the vector from
            // the first to the last control points. In turn, this helps orient the curve
            // to be roughly parallel to the x axis in the ray coordinate system.
            // 
            // In turn (especially for curves that are approaching stright lines.)
            // we get curve bounds with minimal extent in y, which in turn let us
            // early out more quickly in RecursiveIntersect().
            Vector3Df dx = CrossProduct(ray.d, cpObj[3] - cpObj[0]);
            if (GetLengthSquare(dx) == 0.f)
            {
                // If the ray and the vector between the first and last control
                // points are parallel, dx will be zero. Generate an arbitrary xy orientation
                // for the ray coordinate system so that intersection
                // tests can procceed in this unusual case.
                Vector3Df dy;
                CoordinateSystem(ray.d, dx, dy);
            }

            Matrix4f lookAt;
            BuildViewMatrix(lookAt, ray.o, ray.o + ray.d, dx);
            Vector3Df cp[4] = {cpObj[0], cpObj[1], cpObj[2], cpObj[3]};
            TransformPoint(cp[0], lookAt); TransformPoint(cp[1], lookAt); TransformPoint(cp[2], lookAt); TransformPoint(cp[3], lookAt);
            
            // Before going any further, see if the ray's bounding box intersects
            // the curve's bounding box. We start with the y dimension, since the y 
            // extent is generally the smallest (and is often tiny) due to our 
            // careful orientation of the ray coordinate system above.
            float maxWidth = (std::max)(Lerp(m_uMin, m_CurveCommon.width[0], m_CurveCommon.width[1]),
                                        Lerp(m_uMax, m_CurveCommon.width[0], m_CurveCommon.width[1]));
            if ((std::max)((std::max)(cp[0][1], cp[1][1]), (std::max)(cp[2][1], cp[3][1])) +
                0.5f * maxWidth < 0.f || 
                (std::min)((std::min)(cp[0][1], cp[1][1]), (std::min)(cp[2][1], cp[3][1])) - 
                0.5f * maxWidth > 0.f)
                return false;

            // Check for non-overlap in x.
            if ((std::max)((std::max)(cp[0][0], cp[1][0]), (std::max)(cp[2][0], cp[3][0])) + 
                0.5f * maxWidth < 0.f ||
                (std::min)((std::min)(cp[0][0], cp[1][0]), (std::min)(cp[2][0], cp[3][0])) - 
                0.5f * maxWidth > 0.f)
                return false;

            // Check for non-overlap in z.
            float rayLength = GetLength(ray.d);
            float zMax = rayLength * ray.tMax;
            if ((std::max)((std::max)(cp[0][2], cp[1][2]), (std::max)(cp[2][2], cp[3][2])) + 
                0.5f * maxWidth < 0.f ||
                (std::min)((std::min)(cp[0][2], cp[1][2]), (std::min)(cp[2][2], cp[3][2])) - 
                0.5f * maxWidth > zMax)
                return false;

            // Compute refinement depth for curve, _maxDepth_
            float l0 = 0.f;
            for (int32_t i = 0; i < 2; ++i)
            {
                l0 = (std::max)(l0,
                        (std::max)(
                            (std::max)(Abs(cp[i][0] - 2 * cp[i + 1][0] + cp[i + 2][0]),
                                       Abs(cp[i][1] - 2 * cp[i + 1][1] + cp[i + 2][1])),
                            Abs(cp[i][2] - 2 * cp[i + 1][2 + cp[i + 2][2]])
                        ));
            }
            float eps = (std::max)(m_CurveCommon.width[0], m_CurveCommon.width[1]) * .05f; // width / 20
            auto Log2 = [](float v) -> int32_t {
                if (v < 1.f ) return 0;
                uint32_t bits = FloatToBits(v);
                // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
                // (With an additional add so get round-to-nearest rather than
                // round down.)
                return (bits >> 23) - 127 + (bits & (1 << 2)? 1 : 0);
            };
            // Compute log base 4 by dividing log2 in half.
            int32_t r0 = Log2(PI * 6.f * l0 / (8.f * eps)) / 2;
            int32_t maxDepth = Clamp(r0, 0, 10);

            InverseMatrix(lookAt, lookAt);
            return RecursiveIntersect(ray, tHit, isect, cp, lookAt, m_uMin, m_uMax, maxDepth);
        }

        float Area()
        {
            // Compute object-space control points for curve segment, _cpObj_
            Vector3Df cpObj[4];
            cpObj[0] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMin, m_uMin);
            cpObj[1] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMin, m_uMax);
            cpObj[2] = BlossomBezier(m_CurveCommon.cpObj, m_uMin, m_uMax, m_uMax);
            cpObj[3] = BlossomBezier(m_CurveCommon.cpObj, m_uMax, m_uMax, m_uMax);

            float width0 = Lerp(m_uMin, m_CurveCommon.width[0], m_CurveCommon.width[1]);
            float width1 = Lerp(m_uMax, m_CurveCommon.width[0], m_CurveCommon.width[1]);
            float avgWidth = (width0 + width1) * 0.5f;
            float approxLength = 0.f;
            for (int32_t i = 0; i < 3; ++i)
                approxLength += GetLength(cpObj[i] - cpObj[i + 1]);
            return approxLength * avgWidth;
        }

        Interaction Sample(const Vector2Df& u, float& pdf)
        {
            std::cout << "Curve::Sample not implemented." << std::endl;
            return Interaction();
        }

    private:
        bool RecursiveIntersect(const Ray& r, float& tHit, SurfaceInteraction& isect, 
                            const Vector3Df cp[4], const Matrix4f& rayToObject,
                            float u0, float u1, int32_t depth)
        {
            Ray ray(r);
            float rayLength = GetLength(ray.d);

            if (depth > 0.f)
            {
                // Split curve segment into sub-segments and test for intersection
                Vector3Df cpSplit[7];
                SubdivideBezier(cp, cpSplit);

                // For each of the two segments, see if the ray's bounding box
                // overlaps the segment before recursively checking for 
                // intersection with it.
                bool hit = false;
                float u[3] = {u0, (u0 + u1) * .5f, u1};
                // Pointer to the 4 control points for the current segment.
                const Vector3Df* cps = cpSplit;
                for(int32_t seg = 0; seg < 2; ++seg, cps += 3)
                {
                    float maxWidth = (std::max)(Lerp(u[seg], m_CurveCommon.width[0], m_CurveCommon.width[1]),
                                                Lerp(u[seg + 1], m_CurveCommon.width[0], m_CurveCommon.width[1]));

                    // As above, check y first, since it most commonly lets us exit out early.
                    if ((std::max)((std::max)(cp[0][1], cp[1][1]), (std::max)(cp[2][1], cp[3][1])) +
                            0.5f * maxWidth < 0.f || 
                            (std::min)((std::min)(cp[0][1], cp[1][1]), (std::min)(cp[2][1], cp[3][1])) - 
                            0.5f * maxWidth > 0.f)
                        continue;

                    // Check for non-overlap in x.
                    if ((std::max)((std::max)(cp[0][0], cp[1][0]), (std::max)(cp[2][0], cp[3][0])) + 
                        0.5f * maxWidth < 0.f ||
                        (std::min)((std::min)(cp[0][0], cp[1][0]), (std::min)(cp[2][0], cp[3][0])) - 
                        0.5f * maxWidth > 0.f)
                        continue;

                    // Check for non-overlap in z.
                    float rayLength = GetLength(ray.d);
                    float zMax = rayLength * ray.tMax;
                    if ((std::max)((std::max)(cp[0][2], cp[1][2]), (std::max)(cp[2][2], cp[3][2])) + 
                        0.5f * maxWidth < 0.f ||
                        (std::min)((std::min)(cp[0][2], cp[1][2]), (std::min)(cp[2][2], cp[3][2])) - 
                        0.5f * maxWidth > zMax)
                        continue;

                    hit |= RecursiveIntersect(ray, tHit, isect, cps, rayToObject,
                                              u[seg], u[seg + 1], depth - 1);
                    if (hit && !tHit) return true;
                }
                return hit;
            }
            else 
            {
                // Intersect ray with curve segment

                // Test ray agianst segment endpoint boundaries

                // Test sample point against tangent perpendicular at curve start
                float edge = (cp[1][1] - cp[0][1]) * -cp[0][1] + cp[0][0] * (cp[0][0] - cp[1][0]);
                if (edge < 0.f) return false;

                // Test sample point against tangent perpendicular at curve end
                edge = (cp[2][1] - cp[3][1]) * -cp[3][1] + cp[3][0] * (cp[3][0] - cp[2][0]);
                if (edge < 0) return false;

                // Compute line $w$ that gives minimum distance to sample point.
                Vector2Df segmengtDirection = Vector2Df({cp[3][0], cp[3][1]}) - Vector2Df({cp[0][0], cp[0][1]});
                float denom = GetLengthSquare(segmengtDirection);
                if (denom == 0.f) return false;
                float w = DotProduct(-Vector2Df({cp[0][0], cp[0][1]}), segmengtDirection) / denom;

                // Compute $u$ coordinate of curve intersection point and _hitWidth_
                float u = Clamp(Lerp(w, u0, u1), u0, u1);
                float hitWidth = Lerp(u, m_CurveCommon.width[0], m_CurveCommon.width[1]);
                Vector3Df nHit;
                if (m_CurveCommon.type == CurveObjectType::kRibbon)
                {
                    // Scale _hitWidth_ based on ribbon orientation
                    float sin0 = std::sinf((1.f - u) * m_CurveCommon.normalAngle) *
                                 m_CurveCommon.invSinNormalAngle;
                    float sin1 = std::sinf(u * m_CurveCommon.normalAngle) *
                                 m_CurveCommon.invSinNormalAngle;
                    nHit = sin0 * m_CurveCommon.n[0] + sin1 * m_CurveCommon.n[1];
                    hitWidth *= AbsDotProduct(nHit, ray.d) / rayLength;
                }

                // Test intersection point agianst curve width
                Vector3Df dpcdw;
                Vector3Df pc = EvalBezier(cp, Clamp(w, 0, 1), &dpcdw);
                float ptCurveDist2 = pc[0] * pc[0] + pc[1] * pc[1];
                if (ptCurveDist2 > hitWidth * hitWidth * .25f) return false;
                float zMax = rayLength * ray.tMax;
                if (pc[2] < 0 || pc[2] > zMax) return false;

                // Compute $v$ coordinate of curve intersection point
                float ptCurveDist = std::sqrtf(ptCurveDist2);
                float edgeFunc = dpcdw[0] * -pc[1] + pc[0] * dpcdw[1];
                float v = (edgeFunc > 0.f)? .5f + ptCurveDist / hitWidth : .5f - ptCurveDist / hitWidth;

                // Compute hit _t_ and partial derivatives for curve intersection
                {
                    tHit = pc[2] / rayLength;
                    // COmpute error bounds for curve intersection
                    Vector3Df pError({2.f * hitWidth, 2.f * hitWidth, 2.f * hitWidth});

                    // Compute $\dpdu$ and $\dpdv$ for curve intersection
                    Vector3Df dpdu, dpdv;
                    EvalBezier(m_CurveCommon.cpObj, u, &dpdu);
                    
                    if (m_CurveCommon.type == CurveObjectType::kRibbon)
                    {
                        dpdv = Normalize(CrossProduct(nHit, dpdu)) * hitWidth;
                    }
                    else 
                    {
                        // COmpute curve $\dpdv$ for flat and cylinder curves
                        Matrix4f ObjectToRay;
                        InverseMatrix(rayToObject, ObjectToRay);
                        Vector3Df dpduPlane(dpdu);
                        TransformDirection(dpduPlane, ObjectToRay);
                        Vector3Df dpdvPlane = Normalize(Vector3Df({-dpduPlane[1], dpduPlane[0], 0.f})) * hitWidth;
                        if(m_CurveCommon.type == CurveObjectType::kCylinder)
                        {
                            // Rotate _dpdvPlane_ to give cylinderical appearence
                            float theta = Lerp(v, -90.f, 90.f);
                            Matrix4f rot;
                            MatrixRotationAxis(rot, dpduPlane, -theta);
                            TransformDirection(dpdvPlane, rot);
                        }
                        dpdv = dpdvPlane;
                        TransformDirection(dpdv, rayToObject);
                    }
                    std::shared_ptr<Shape> pShape(this);
                    isect = TransformInteraction(SurfaceInteraction(
                        ray(tHit), pError, Vector2Df({u, v}), -ray.d, dpdu, dpdv,
                        Vector3Df({0.f, 0.f, 0.f}), Vector3Df({0.f, 0.f, 0.f}), pShape), *m_ObjectToWorld);
                }
                return true;
            }
        }
        const CurveCommon m_CurveCommon;
        const float m_uMin, m_uMax;
    };
}