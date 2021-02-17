#include "SceneObjectTriangleMesh.hpp"

namespace Panda
{
    SceneObjectTriangleMesh::SceneObjectTriangleMesh(const std::shared_ptr<Matrix4f>& ObjectToWorld, int32_t nTriangles,
            const int32_t* vertexIndices, int32_t nVertices, const Vector3Df* P,
            const Vector3Df* S, const Vector3Df* N, const Vector2Df* UV,
            const std::shared_ptr<SceneObjectTexture>& alphamask,
            const std::shared_ptr<SceneObjectTexture>& shadowAlphaMask,
            const int32_t* faceIndices)
        : BaseSceneObject(SceneObjectType::kSceneObjectTypeMesh),
          m_VerticesCount(nVertices),
          m_TrianglesCount(nTriangles),
          m_AlphaMask(alphamask),
          m_ShadowAlphaMask(shadowAlphaMask)
    {
        assert(P);
        Vector3Df* p = reinterpret_cast<Vector3Df*>(g_pMemoryManager->Allocate(nVertices * sizeof(Vector3Df)));
        memcpy(p, P, sizeof(Vector3Df) * nVertices);
        for (int32_t i = 0; i < nVertices; ++i)
            TransformPoint(p[i], *ObjectToWorld);
        m_VertexArray = std::make_unique<SceneObjectVertexArray>("", 0, VertexDataType::kVertexDataTypeFloat3, p, nVertices);

        assert(vertexIndices);
        int32_t* indices = reinterpret_cast<int32_t*>(g_pMemoryManager->Allocate(nTriangles * 3 * sizeof(int32_t)));
        memcpy(indices, vertexIndices, nTriangles * 3 * sizeof(int32_t));
        m_IndexArray = std::make_unique < SceneObjectIndexArray>(0, 0, IndexDataType::kIndexDataTypeInt32, indices, nTriangles * 3);

        if (N)
        {
            Vector3Df* n = reinterpret_cast<Vector3Df*>(g_pMemoryManager->Allocate(nVertices * sizeof(Vector3Df)));
            memcpy(n, N, sizeof(Vector3Df) * nVertices);
            for (int32_t i = 0; i < nVertices; ++i)
                TransformPoint(n[i], *ObjectToWorld);
            m_NormalArray = std::make_unique<SceneObjectVertexArray>("", 0, VertexDataType::kVertexDataTypeFloat3, n, nVertices);
        }
        
        if (S)
        {
            Vector3Df* s = reinterpret_cast<Vector3Df*>(g_pMemoryManager->Allocate(nVertices * sizeof(Vector3Df)));
            memcpy(s, S, sizeof(Vector3Df) * nVertices);
            for (int32_t i = 0; i < nVertices; ++i)
                TransformPoint(s[i], *ObjectToWorld);
            m_TangentArray = std::make_unique<SceneObjectVertexArray>("", 0, VertexDataType::kVertexDataTypeFloat3, s, nVertices);
        }

        if (UV)
        {
            Vector2Df* uv = reinterpret_cast<Vector2Df*>(g_pMemoryManager->Allocate(nVertices * sizeof(Vector2Df)));
            memcpy(uv, UV, sizeof(Vector2Df) * nVertices);
            m_UVArray = std::make_unique<SceneObjectVertexArray>("", 0, VertexDataType::kVertexDataTypeFloat2, uv, nVertices);
        }

        if (faceIndices)
            m_FaceIndices = std::vector<int32_t>(faceIndices, faceIndices + nTriangles);
    }
}