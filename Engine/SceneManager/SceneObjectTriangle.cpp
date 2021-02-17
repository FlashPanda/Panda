#include "SceneObjectTriangle.hpp"

namespace Panda
{
	SceneObjectTriangle::SceneObjectTriangle(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
		const std::shared_ptr<SceneObjectTriangleMesh>& mesh, int32_t triNumber)
        : SceneObjectShape(objectToWorld, worldToObject),
        m_TriangleMesh(mesh)
    {
        const std::unique_ptr<SceneObjectIndexArray>& ia = m_TriangleMesh->GetIndexArray();
        const void* pData = ia->GetData();
        m_IndexDataType = ia->GetIndexType();
        size_t sizeBytes = sizeof(int8_t);
        switch (m_IndexDataType)
        {
        case IndexDataType::kIndexDataTypeInt8:
            sizeBytes = sizeof(int8_t);
            break;
        case IndexDataType::kIndexDataTypeInt16:
            sizeBytes = sizeof(int16_t);
            break;
        case IndexDataType::kIndexDataTypeInt32:
            sizeBytes = sizeof(int32_t);
            break;
        case IndexDataType::kIndexDataTypeInt64:
            sizeBytes = sizeof(int64_t);
            break;
        default:
            assert(false);
            break;
        }
        // Each triangle has 3 indices, each index has _sizeBytes_ memory space.
        m_pIndices = reinterpret_cast<const int8_t*>(pData) + triNumber * 3 * sizeBytes;
        m_FaceIndex = m_TriangleMesh->GetFaceIndex(triNumber);
    }

    Bounds3Df SceneObjectTriangle::ObjectBound()
    {
        // Get triangle vertices in _p0_, _p1_ and _p2_
		const std::unique_ptr<SceneObjectVertexArray>& arr = m_TriangleMesh->GetVertexArray();
        const void* pData = arr->GetData();
        const Vector3Df* pVertices = reinterpret_cast<const Vector3Df*>(pData);
        // HINT: the index format maybe uint16_t in dae outputed by blender
        const int32_t* p = reinterpret_cast<const int32_t*>(m_pIndices);
        const Vector3Df& p0 = pVertices[p[0]]; const Vector3Df& p1 = pVertices[p[1]]; const Vector3Df& p2 = pVertices[p[2]];
        Vector3Df oP0 = p0; Vector3Df oP1 = p1; Vector3Df oP2 = p2;
        TransformPoint(oP0, *m_WorldToObject); TransformPoint(oP1, *m_WorldToObject); TransformPoint(oP2, *m_WorldToObject);
        return Union(Bounds3Df(oP0, oP1), oP2);
    }

    Bounds3Df SceneObjectTriangle::WorldBound()
    {
		const std::unique_ptr<SceneObjectVertexArray>& arr = m_TriangleMesh->GetVertexArray();
        const void* pData = arr->GetData();
        const Vector3Df* pVertices = reinterpret_cast<const Vector3Df*>(pData);
        // HINT: the index format maybe uint16_t in dae outputed by blender
        const int32_t* p = reinterpret_cast<const int32_t*>(m_pIndices);
        const Vector3Df& p0 = pVertices[p[0]]; const Vector3Df& p1 = pVertices[p[1]]; const Vector3Df& p2 = pVertices[p[2]];
        return Union(Bounds3Df(p0, p1), p2);
    }

    bool SceneObjectTriangle::Intersect(const Ray& r, float& tHit, SceneObjectSurfaceInteraction& isect,
        bool testAlphaTexture) 
    {
        // TO BE UNDERSTOOD
        // ProfilePhase p(Proof::TriIntersect);
		const std::unique_ptr<SceneObjectVertexArray>& arr = m_TriangleMesh->GetVertexArray();
        const void* pData = arr->GetData();
        const Vector3Df* pVertices = reinterpret_cast<const Vector3Df*>(pData);
        // HINT: the index format maybe uint16_t in dae outputed by blender
        const int32_t* p = reinterpret_cast<const int32_t*>(m_pIndices);
        const Vector3Df& p0 = pVertices[p[0]]; const Vector3Df& p1 = pVertices[p[1]]; const Vector3Df& p2 = pVertices[p[2]];

        Ray ray(r);
        // Perform ray--triangle intersection test

        // Transform triangle vertices to ray coordinate space

        // Translate vertices based on ray origin
        Vector3Df p0t = p0 - ray.o; // Gamma(1)
        Vector3Df p1t = p1 - ray.o; // Gamma(1)
        Vector3Df p2t = p2 - ray.o; // Gamma(1)
        int32_t kz = MaxAbsDimension(ray.d);
        int32_t kx = kz + 1;
        if (kx == 3) kx = 0;
        int32_t ky = kx + 1;
        if (ky == 3) ky = 0;
        Vector3Df d = Permute(ray.d, kx, ky, kz);
        p0t = Permute(p0t, kx, ky, kz);
        p1t = Permute(p1t, kx, ky, kz);
        p2t = Permute(p2t, kx, ky, kz);

        // Apply shear transformation to translate vertex positions
        // For now, only the x and y dimensions are sheared. We can
        // wait and shear the z dimension only if the ray acturally
        // intersects the triangle.
        float Sx = -d[0] / d[2]; // Gamma(1)
        float Sy = -d[1] / d[2]; // Gamma(1)
        float Sz = 1.f / d[2]; // Gamma(1)
        p0t[0] += Sx * p0t[2]; // Gamma(1)*Gamma(1) + Gamma(1)=Gamma(4)+Gamma(2)
        p0t[1] += Sy * p0t[2];
        p1t[0] += Sx * p1t[2];
        p1t[1] += Sy * p1t[2];
        p2t[0] += Sx * p2t[2];
        p2t[1] += Sy * p2t[2];

        // Compute edge function coefficients _e0_, _e1_, _e2_
        float e2 = p0t[0] * p1t[1] - p0t[1] * p1t[0]; // Gamma(4)*Gamma(4)-Gamma(4)*Gamma(4)=Gamma(10)
        float e0 = p1t[0] * p2t[1] - p1t[1] * p2t[0];
        float e1 = p2t[0] * p0t[1] - p2t[1] * p0t[0];

        // Fall back to double precision test at triangle edges.
        if (e0 == 0.f || e1 == 0.f || e2 == 0.f)
        {
            double p2txp1ty = (double)p2t.data[0] * (double)p1t.data[1];
            double p2typ1tx = (double)p2t.data[1] * (double)p1t.data[0];
            e0 = (float)(p2typ1tx - p2txp1ty);
            double p0txp2ty = (double)p0t.data[0] * (double)p2t.data[1];
            double p0typ2tx = (double)p0t.data[1] * (double)p2t.data[0];
            e1 = (float)(p0typ2tx - p0txp2ty);
            double p1txp0ty = (double)p1t.data[0] * (double)p0t.data[1];
            double p1typ0tx = (double)p1t.data[1] * (double)p0t.data[0];
            e2 = (float)(p1typ0tx - p1txp0ty);
        }

        // Perform triangle edge and determinant tests.
        if ((e0 < 0.f || e1 < 0.f || e2 < 0.f) && (e0 > 0.f || e1 > 0.f || e2 > 0.f))
            return false;
        float det = e0 + e1 + e2; // Gamma(10) + Gamma(10) +Gamma(10) = Gamma(12)
        if (det == 0.f)  return false;

        // Compute scaled hit distance to triangle and test agianst ray $t$ range
        // Save the cost of the floating-point division.
        p0t[2] *= Sz; // Gamma(2)
        p1t[2] *= Sz;
        p2t[2] *= Sz;
        float tScaled = e0 * p0t[2] + e1 * p1t[2] + e2 * p2t[2]; // Gamma(10)*Gamma(4)+Gamma(10)*Gamma(4)+Gamma(10)*Gamma(4) = Gamma(17)
        if (det < 0.f && (tScaled >= 0.f || tScaled < ray.tMax * det))
            return false;
        else if (det > 0.f && (tScaled <= 0.f || tScaled > ray.tMax * det))
            return false;

        // Compute barycentric coordinates and $t$ value for triangle intersection
        float invDet = 1.f / det; // Gamma(13)
        float b0 = e0 * invDet; // Gamma(10) / Gamma(13) = Gamma(23)
        float b1 = e1 * invDet;
        float b2 = e2 * invDet;
        float t = tScaled * invDet; // Gamma(17) * Gamma(13) = Gamma(31)

        // Ensure that computed triangle $t$ is conservatively greater thatn zero

        // Compute $\delta_z$ term for tirangle $t$ error bounds
        float maxZt = MaxComponent(Abs(Vector3Df({p0t[2], p1t[2], p2t[2]})));
        float deltaZ = Gamma(3) * maxZt;

        // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        float maxXt = MaxComponent(Abs(Vector3Df({p0t[0], p1t[0], p2t[0]})));
        float maxYt = MaxComponent(Abs(Vector3Df({p0t[1], p1t[1], p2t[1]})));
        float deltaX = Gamma(5) * (maxXt + maxZt);
        float deltaY = Gamma(5) * (maxYt + maxZt);

        // Compute $\delta_e$ term for triangle $t$ error bounds
        float deltaE = 2 * (Gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

        // Compute $\delta_t$ term for tirnagle $t$ error bounds and check _t_
        float maxE = MaxComponent(Abs(Vector3Df({e0, e1, e2})));
        float deltaT = 3 * (Gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * Abs(invDet);
        if (t <= deltaT) return false;

        // Compute triangle partial derivatives
        Vector3Df dpdu, dpdv;
        Vector2Df uv[3];
        GetUVs(uv);

        // Compute deltas for triangle partial derivatives
        Vector2Df duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
        Vector3Df dp02 = p0 - p1, dp12 = p1 - p2;
        float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
        bool degenerateUV = Abs(determinant) < 1e-8;
        if (!degenerateUV)
        {
            float invdet = 1.f / determinant;
            dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
            dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
        }
        if (degenerateUV || GetLengthSquare(CrossProduct(dpdu, dpdv)) == 0.f)
        {
            // Handle zero determinant for triangle partial derivative matrix
            Vector3Df ng = CrossProduct(p2 - p0, p1 - p0);
            if (GetLengthSquare(ng) == 0.f)
                // The triangle is actually degenerate; the intersection is bogus.
                return false;
            CoordinateSystem(Normalize(ng), dpdu, dpdv);
        }

        // Compute error bounds for triangle intersection
        float xAbsSum = (Abs(b0 * p0[0]) + Abs(b1 * p1[0]) + Abs(b2 * p2[0]));
        float yAbsSum = (Abs(b0 * p0[1]) + Abs(b1 * p1[1]) + Abs(b2 * p2[1]));
        float zAbsSum = (Abs(b0 * p0[2]) + Abs(b1 * p1[2]) + Abs(b2 * p2[2]));
        Vector3Df pError = Gamma(7) * Vector3Df({xAbsSum, yAbsSum, zAbsSum});

        // Interpolate $(u,v)$ parametric coordinates and hit point
        Vector3Df pHit = b0 * p0 + b1 * p1 + b2 * p2;
        Vector2Df uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

        // Test intersection agianst alpha texture, if present
        // if (testAlphaTexture && m_TriangleMesh->DoesAlphaMaskExist())
        // {
        //     std::shared_ptr<Shape> pShape(this);
        //     SurfaceInteraction isectLocal(pHit, Vector3Df({0.f, 0.f, 0.f}), uvHit, -ray.d, dpdu, dpdv,
        //                                   Vector3Df({0.f, 0.f, 0.f}), Vector3Df({0.f, 0.f, 0.f}), pShape);
        //     const SceneObjectTexture& pTexture = m_TriangleMesh->GetAlphaMask();
        //     pTexture->Evaluate(isectLocal);
        // }

        // Fill in _SurfaceInteraction_ from triangle hit
        isect = SceneObjectSurfaceInteraction(pHit, pError, uvHit, -ray.d, dpdu, dpdv,
                                    Vector3Df({0.f, 0.f, 0.f}), Vector3Df({0.f, 0.f, 0.f}), this, m_FaceIndex);
        
        // Override surface normal in _isect_ for triangle
        isect.normal = isect.shading.normal = Vector3Df(Normalize(CrossProduct(dp02, dp12)));
        if (m_TriangleMesh->DoesNormalExist() || m_TriangleMesh->DoesTangentExist())
        {
            // Initialize _Triangle_ shading geometry

            // Compute shading normal _ns_ for triangle
            Vector3Df ns;
            if (m_TriangleMesh->DoesNormalExist())
            {
                const std::unique_ptr<SceneObjectVertexArray>& normalArray = m_TriangleMesh->GetNormalArray();
                VertexDataType dataType = normalArray->GetDataType();
                assert(dataType == VertexDataType::kVertexDataTypeFloat3);
                const Vector3Df* pNormal = reinterpret_cast<const Vector3Df*>(normalArray->GetData());
                switch (m_IndexDataType)
                {
                case IndexDataType::kIndexDataTypeInt8:
                {
                    const int8_t* pData = reinterpret_cast<const int8_t*>(m_pIndices);
                    ns = (b0 * pNormal[pData[0]] + b1 * pNormal[pData[1]] + b2 * pNormal[pData[2]]);
                }
                    break;
                case IndexDataType::kIndexDataTypeInt16:
                {
                    const int16_t* pData = reinterpret_cast<const int16_t*>(m_pIndices);
                    ns = (b0 * pNormal[pData[0]] + b1 * pNormal[pData[1]] + b2 * pNormal[pData[2]]);
                }
                    break;
                case IndexDataType::kIndexDataTypeInt32:
                {
                    const int32_t* pData = reinterpret_cast<const int32_t*>(m_pIndices);
                    ns = (b0 * pNormal[pData[0]] + b1 * pNormal[pData[1]] + b2 * pNormal[pData[2]]);
                }
                    break;
                case IndexDataType::kIndexDataTypeInt64:
                {
                    const int64_t* pData = reinterpret_cast<const int64_t*>(m_pIndices);
                    ns = (b0 * pNormal[pData[0]] + b1 * pNormal[pData[1]] + b2 * pNormal[pData[2]]);
                }
                    break;
                default:
                    assert(false);
                    break;
                }
                if(GetLengthSquare(ns) > 0.f)
                    ns = Normalize(ns);
                else 
                    ns = isect.normal;
            }
            else
                ns = isect.normal;

            // Compute shading tangent _ss_ for triangle
            Vector3Df ss;
            if (m_TriangleMesh->DoesTangentExist())
            {
                const std::unique_ptr<SceneObjectVertexArray>& tangentArray = m_TriangleMesh->GetTangentArray();
                VertexDataType dataType = tangentArray->GetDataType();
                assert(dataType == VertexDataType::kVertexDataTypeFloat3);
                const Vector3Df* pTangent = reinterpret_cast<const Vector3Df*>(tangentArray->GetData());
                switch (m_IndexDataType)
                {
                case IndexDataType::kIndexDataTypeInt8:
                {
                    const int8_t* pData = reinterpret_cast<const int8_t*>(m_pIndices);
                    ss = (b0 * pTangent[pData[0]] + b1 * pTangent[pData[1]] + b2 * pTangent[pData[2]]);
                }
                    break;
                case IndexDataType::kIndexDataTypeInt16:
                {
                    const int16_t* pData = reinterpret_cast<const int16_t*>(m_pIndices);
                    ss = (b0 * pTangent[pData[0]] + b1 * pTangent[pData[1]] + b2 * pTangent[pData[2]]);
                }
                    break;
                case IndexDataType::kIndexDataTypeInt32:
                {
                    const int32_t* pData = reinterpret_cast<const int32_t*>(m_pIndices);
                    ss = (b0 * pTangent[pData[0]] + b1 * pTangent[pData[1]] + b2 * pTangent[pData[2]]);
                }
                    break;
                case IndexDataType::kIndexDataTypeInt64:
                {
                    const int64_t* pData = reinterpret_cast<const int64_t*>(m_pIndices);
                    ss = (b0 * pTangent[pData[0]] + b1 * pTangent[pData[1]] + b2 * pTangent[pData[2]]);
                }
                    break;
                default:
                    assert(false);
                    break;
                }

                if(GetLengthSquare(ss) > 0.f)
                    ss = Normalize(ss);
                else 
                    ss = Normalize(isect.dpdu);
            }
            else
            {
                ss = Normalize(isect.dpdu);
            }

            // Compute shading bitangent _ts_ for triangle and adjust _ss_
            Vector3Df ts = CrossProduct(ss, ns);
            if (GetLengthSquare(ts) > 0.f)
            {
                ts = Normalize(ts);
                ss = CrossProduct(ts, ns);
            }
            else 
                CoordinateSystem(ns, ss, ts);

            // Compute $\dndu$ and $\dndv$ for tirangle shading geometry
            Vector3Df dndu, dndv;
            if (m_TriangleMesh->DoesNormalExist())
            {
                // Compute deltas for triangle partial derivatives of normal
                Vector2Df duv02 = uv[0] - uv[2];
                Vector2Df duv12 = uv[1] - uv[2];
                const std::unique_ptr<SceneObjectVertexArray>& normalArray = m_TriangleMesh->GetNormalArray();
                VertexDataType dataType = normalArray->GetDataType();
                assert(dataType == VertexDataType::kVertexDataTypeFloat3);
                const Vector3Df* pNormal = reinterpret_cast<const Vector3Df*>(normalArray->GetData());
                Vector3Df dn1, dn2;
                Vector3Df normal0, normal1, normal2;
                switch (m_IndexDataType)
                {
                case IndexDataType::kIndexDataTypeInt8:
                {
                    const int8_t* pData = reinterpret_cast<const int8_t*>(m_pIndices);
                    normal0 = pNormal[pData[0]]; normal1 = pNormal[pData[1]]; normal2 = pNormal[pData[2]];
                }
                    break;
                case IndexDataType::kIndexDataTypeInt16:
                {
                    const int16_t* pData = reinterpret_cast<const int16_t*>(m_pIndices);
                    normal0 = pNormal[pData[0]]; normal1 = pNormal[pData[1]]; normal2 = pNormal[pData[2]];
                }
                    break;
                case IndexDataType::kIndexDataTypeInt32:
                {
                    const int32_t* pData = reinterpret_cast<const int32_t*>(m_pIndices);
                    normal0 = pNormal[pData[0]]; normal1 = pNormal[pData[1]]; normal2 = pNormal[pData[2]];
                }
                    break;
                case IndexDataType::kIndexDataTypeInt64:
                {
                    const int64_t* pData = reinterpret_cast<const int64_t*>(m_pIndices);
                    normal0 = pNormal[pData[0]]; normal1 = pNormal[pData[1]]; normal2 = pNormal[pData[2]];                      
                }
                    break;
                default:
                    assert(false);
                    break;
                }
                dn1 = normal0 - normal2;
                dn2 = normal1 - normal2;
                float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
                bool degenerateUV = Abs(determinant) < 1e-8;
                if (degenerateUV)
                {
                    // We can still compute dndu and dndv,, with respect to the 
                    // same aritrary coordinate system we use to compute dpdu
                    // and dpdv when this happens. It's important to do this
                    // (rather than giving up) so that ray differentials for 
                    // rays reflected from triangles with degenerate
                    // parameterizations are still reasonable.
                    Vector3Df dn = CrossProduct(normal2 - normal0, normal1 - normal0);
                    if (GetLengthSquare(dn) == 0.f)
                        dndu = dndv = Vector3Df(0.f);
                    else
                    {
                        Vector3Df dnu, dnv;
                        CoordinateSystem(dn, dnu, dnv);
                        dndu = Vector3Df(dnu);
                        dndv = Vector3Df(dnv);
                    }
                }
                else 
                {
                    float invDet = 1.f / determinant;
                    dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
                    dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
                }
            }
            else 
                dndu = dndv = Vector3Df(0.f);
            isect.SetShadingGeometry(ss, ts, dndu, dndv);
        }

        // Ensure correct orientation of the geometric normal
        if (m_TriangleMesh->DoesNormalExist())
            isect.normal = Faceforward(isect.normal, isect.shading.normal);
        tHit = t;
        return true;
    }

    bool SceneObjectTriangle::IntersectP(const Ray& r, bool testAlphaTexture)
    {
        // TO BE UNDERSTOOD
        // ProfilePhase p(Proof::TriIntersect);
        const std::unique_ptr<SceneObjectVertexArray>& arr = m_TriangleMesh->GetVertexArray();
        const void* pData = arr->GetData();
        const Vector3Df* pVertices = reinterpret_cast<const Vector3Df*>(pData);
        // HINT: the index format maybe uint16_t in dae outputed by blender
        const int32_t* p = reinterpret_cast<const int32_t*>(m_pIndices);
        const Vector3Df& p0 = pVertices[p[0]]; const Vector3Df& p1 = pVertices[p[1]]; const Vector3Df& p2 = pVertices[p[2]];

        Ray ray(r);
        // Perform ray--triangle intersection test

        // Transform triangle vertices to ray coordinate space

        // Translate vertices based on ray origin
        Vector3Df p0t = p0 - ray.o;
        Vector3Df p1t = p1 - ray.o;
        Vector3Df p2t = p2 - ray.o;
        int32_t kz = MaxAbsDimension(ray.d);
        int32_t kx = kz + 1;
        if (kx == 3) kx = 0;
        int32_t ky = kx + 1;
        if (ky == 3) ky = 0;
        Vector3Df d = Permute(ray.d, kx, ky, kz);
        p0t = Permute(p0t, kx, ky, kz);
        p1t = Permute(p1t, kx, ky, kz);
        p2t = Permute(p2t, kx, ky, kz);

        // Apply shear transformation to translate vertex positions
        // For now, only the x and y dimensions are sheared. We can
        // wait and shear the z dimension only if the ray acturally
        // intersects the triangle.
        float Sx = -d[0] / d[2];
        float Sy = -d[1] / d[2];
        float Sz = 1.f / d[2];
        p0t[0] += Sx * p0t[2];
        p0t[1] += Sy * p0t[2];
        p1t[0] += Sx * p1t[2];
        p1t[1] += Sy * p1t[2];
        p2t[0] += Sx * p2t[2];
        p2t[1] += Sy * p2t[2];

        // Compute edge function coefficients _e0_, _e1_, _e2_
        float e2 = p0t[0] * p1t[1] - p0t[1] * p1t[0];
        float e0 = p1t[0] * p2t[1] - p1t[1] * p2t[0];
        float e1 = p2t[0] * p0t[1] - p2t[1] * p0t[0];

        // Fall back to double precision test at triangle edges.
        if (e0 == 0.f || e1 == 0.f || e2 == 0.f)
        {
            double p2txp1ty = (double)p2t[0] * (double)p1t[1];
            double p2typ1tx = (double)p2t[1] * (double)p1t[0];
            e0 = (float)(p2typ1tx - p2txp1ty);
            double p0txp2ty = (double)p0t[0] * (double)p2t[1];
            double p0typ2tx = (double)p0t[1] * (double)p2t[0];
            e1 = (float)(p0typ2tx - p0txp2ty);
            double p1txp0ty = (double)p1t[0] * (double)p0t[1];
            double p1typ0tx = (double)p1t[1] * (double)p0t[0];
            e2 = (float)(p1typ0tx - p1txp0ty);
        }

        // Perform triangle edge and determinant tests.
        if ((e0 < 0.f || e1 < 0.f || e2 < 0.f) && (e0 > 0.f || e1 > 0.f || e2 > 0.f))
            return false;
        float det = e0 + e1 + e2;
        if (det == 0.f)  return false;

        // Compute scaled hit distance to triangle and test agianst ray $t$ range
        // Save the cost of the floating-point division.
        p0t[2] *= Sz;
        p1t[2] *= Sz;
        p2t[2] *= Sz;
        float tScaled = e0 * p0t[2] + e1 * p1t[2] + e2 * p2t[2];
        if (det < 0.f && (tScaled >= 0.f || tScaled < ray.tMax * det))
            return false;
        else if (det > 0.f && (tScaled <= 0.f || tScaled > ray.tMax * det))
            return false;

        // Compute barycentric coordinates and $t$ value for triangle intersection
        float invDet = 1.f / det;
        float b0 = e0 * invDet;
        float b1 = e1 * invDet;
        float b2 = e2 * invDet;
        float t = tScaled * invDet;

        // Ensure that computed triangle $t$ is conservatively greater than zero

        // Compute $\delta_z$ term for triangle $t$ error bounds
        float maxZt = MaxAbsComponent(Vector3Df({p0t[2], p1t[2], p2t[2]}));
        float deltaZ = Gamma(3) * maxZt;

        // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        float maxXt = MaxAbsComponent(Vector3Df({p0t[0], p1t[0], p2t[0]}));
        float maxYt = MaxAbsComponent(Vector3Df({p0t[1], p1t[1], p2t[1]}));
        float deltaX = Gamma(5) * (maxXt + maxZt);
        float deltaY = Gamma(5) * (maxYt + maxZt);

        // Compute $\delta_e$ term for triangle $t$ error bounds 
        float deltaE = 2 * (Gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

        // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        float maxE = MaxComponent(Abs(Vector3Df({e0, e1, e2})));
        float deltaT = 3 * (Gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * Abs(invDet);
        if (t <= deltaT) return false;

        // Test shadow ray intersection against alpha texture, if present
        if (testAlphaTexture && (m_TriangleMesh->DoesAlphaMaskExist() || m_TriangleMesh->DoesShadowAlphaMaskExist()))
        {
            // Compute triangle partial derivatives
            Vector3Df dpdu, dpdv;
            Vector2Df uv[3];
            GetUVs(uv);

            // Compute deltas for triangle partial derivatives
            Vector2Df duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
            Vector3Df dp02 = p0 - p2, dp12 = p1 - p2;
            float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            bool degenerateUV = Abs(determinant) < 1e-8;
            if(!degenerateUV)
            {
                float invdet = 1.f / determinant;
                dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
                dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
            }

            if (degenerateUV || GetLengthSquare(CrossProduct(dpdu, dpdv)) == 0.f)
            {
                // Handle zero determinant for triangle partial derivative matrix
                Vector3Df ng = CrossProduct(p2 - p0, p1 - p0);
                if (GetLengthSquare(ng) == 0.f)
                    // The triangle is acutally degenerate; the intersection is bogus.
                    return false;
                CoordinateSystem(Normalize(CrossProduct(p2 - p0, p1 - p0)), dpdu, dpdv);
            }

            // Interpolate $(u,v)$ parametric coordinates and hit point
            Vector3Df pHit = b0 * p0 + b1 * p1  + b2 * p2;
            Vector2Df uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
            SceneObjectSurfaceInteraction isectLocal(pHit, Vector3Df(0.f), uvHit, -ray.d, 
                                            dpdu, dpdv, Vector3Df(0.f), Vector3Df(0.f), this);
            // if (m_TriangleMesh->DoesAlphaMaskExist() && m_TriangleMesh->alphaMask->Evaluate(isectLocal)==0)
            //     return false;
            // if (m_TriangleMesh->DoesShadowAlphaMaskExist() && m_TriangleMesh->ShadowAlphaMask->Evaluate(isectLocal) == 0)
            //     return false;
        }

        return true;
    }

    float SceneObjectTriangle::Area() 
    {
        const std::unique_ptr<SceneObjectVertexArray>& arr = m_TriangleMesh->GetVertexArray();
        const void* pData = arr->GetData();
        const Vector3Df* pVertices = reinterpret_cast<const Vector3Df*>(pData);
        // HINT: the index format maybe uint16_t in dae outputed by blender
        const int32_t* p = reinterpret_cast<const int32_t*>(m_pIndices);
        const Vector3Df& p0 = pVertices[p[0]]; const Vector3Df& p1 = pVertices[p[1]]; const Vector3Df& p2 = pVertices[p[2]];

        return 0.5f * GetLength(CrossProduct(p2 - p0, p1 - p0));
    }
}