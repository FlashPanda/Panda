#pragma once

#include <vector>
#include "BaseSceneObject.hpp"
#include "SceneObjectTransform.hpp"
#include "Core/Math/EFloat.hpp"
#include "SceneObjectTriangleMesh.hpp"
#include "SceneObjectShape.hpp"
#include "SceneObjectInteraction.hpp"

namespace Panda
{
    class SceneObjectTriangle : public SceneObjectShape
    {
    public:
        SceneObjectTriangle(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
                            const std::shared_ptr<SceneObjectTriangleMesh>& mesh, int32_t triNumber);

        Bounds3Df ObjectBound();

        Bounds3Df WorldBound();

        bool Intersect(const Ray& r, float& tHit, SceneObjectSurfaceInteraction& isect,
            bool testAlphaTexture); 

        bool IntersectP(const Ray& r, bool testAlphaTexture);

        float Area();

    protected:
        // Get UV coordinates
        void GetUVs(Vector2Df uv[3])
        {
            if (m_TriangleMesh->DoesUVExist())
            {
                const std::unique_ptr<SceneObjectVertexArray>& va = m_TriangleMesh->GetUVArray();
                VertexDataType vaDataType = va->GetDataType();
                assert(vaDataType == VertexDataType::kVertexDataTypeFloat2);
                const Vector2Df* pUV = reinterpret_cast<const Vector2Df*>(va->GetData());
                switch(vaDataType)
                {
                    case IndexDataType::kIndexDataTypeInt8:
                    {
                        const int8_t* pData = reinterpret_cast<const int8_t*>(m_pIndices);
                        uv[0] = pUV[pData[0]];
                        uv[1] = pUV[pData[1]];
                        uv[2] = pUV[pData[2]];
                    }
                        break;
                    case IndexDataType::kIndexDataTypeInt16:
                    {
                        const int16_t* pData = reinterpret_cast<const int16_t*>(m_pIndices);
                        uv[0] = pUV[pData[0]];
                        uv[1] = pUV[pData[1]];
                        uv[2] = pUV[pData[2]];
                    }
                        break;
                    case IndexDataType::kIndexDataTypeInt32:
                    {
                        const int32_t* pData = reinterpret_cast<const int32_t*>(m_pIndices);
                        uv[0] = pUV[pData[0]];
                        uv[1] = pUV[pData[1]];
                        uv[2] = pUV[pData[2]];
                    }
                        break;
                    case IndexDataType::kIndexDataTypeInt64:
                    {
                        const int64_t* pData = reinterpret_cast<const int64_t*>(m_pIndices);
                        uv[0] = pUV[pData[0]];
                        uv[1] = pUV[pData[1]];
                        uv[2] = pUV[pData[2]];
                    }
                        break;
                    default:
						assert(false);
						break;
                }
            }
            else 
            {
                uv[0] = Vector2Df({0.f, 0.f});
                uv[1] = Vector2Df({1.f, 0.f});
                uv[2] = Vector2Df({1.f, 1.f});
            }
        }

    protected:
        std::shared_ptr<SceneObjectTriangleMesh> m_TriangleMesh;

        const void* m_pIndices;
		IndexDataType m_IndexDataType;
        int32_t m_FaceIndex;
    };

    // Yeah, we need to create triangle meshs, 'cause we need to 
    // subdivide surface.
    // Parameters:
    //      nTriangles - the count of triangles
    //      vertexIndices - the indices array of vertices
    //      nVertices - count of vertices
    //      p - the vertices array
    //      s - the tangent array
    //      n - the normal array
    //      uv - the uv array
    //      faceIndices - the face indices array
    std::vector<std::shared_ptr<SceneObjectShape>> CreateTriangleMesh(
        const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject,
        int32_t nTriangles, const int32_t* vertexIndices, int32_t nVertices,
        const Vector3Df* p, const Vector3Df* s, const Vector3Df* n, const Vector2Df* uv,
		const std::shared_ptr<SceneObjectTexture>& alphaMask, 
		const std::shared_ptr<SceneObjectTexture>& shadowAlphaMask,
        const int32_t* faceIndices = nullptr)
    {
		std::shared_ptr<SceneObjectTriangleMesh> mesh = std::make_shared<SceneObjectTriangleMesh>(
			objectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,
			alphaMask, shadowAlphaMask, faceIndices);

		std::vector<std::shared_ptr<SceneObjectShape>> tris;
		tris.reserve(nTriangles);
		for (int32_t i = 0; i < nTriangles; ++i)
			tris.push_back(std::make_shared<SceneObjectTriangle>(objectToWorld, worldToObject, mesh, i));

		return tris;
    }
}