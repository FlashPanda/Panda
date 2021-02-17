#pragma once

#include <vector>
#include "SceneObjectIndexArray.hpp"
#include "SceneObjectVertexArray.hpp"
#include "Core/Math/PandaMath.hpp"
#include "BaseSceneObject.hpp"
#include "Core/MemoryManager.hpp"

namespace Panda
{
    class SceneObjectTriangleMesh : public BaseSceneObject
    {
    public:
        SceneObjectTriangleMesh(const std::shared_ptr<Matrix4f>& ObjectToWorld, int32_t nTriangles,
            const int32_t* vertexIndices, int32_t nVertices, const Vector3Df* P,
            const Vector3Df* S, const Vector3Df* N, const Vector2Df* UV,
            const std::shared_ptr<SceneObjectTexture>& alphamask,
            const std::shared_ptr<SceneObjectTexture>& shadowAlphaMask,
            const int32_t* faceIndices);

        SceneObjectTriangleMesh(SceneObjectTriangleMesh&& mesh)
            : BaseSceneObject(SceneObjectType::kSceneObjectTypeMesh),
              m_PrimitiveType(mesh.m_PrimitiveType),
              m_VertexArray(std::move(mesh.m_VertexArray)),
              m_IndexArray(std::move(mesh.m_IndexArray)),
              m_TrianglesCount(mesh.m_TrianglesCount),
              m_VerticesCount(mesh.m_VerticesCount),
              m_FaceIndices(std::move(m_FaceIndices)),
              m_NormalArray(std::move(m_NormalArray)),
              m_TangentArray(std::move(mesh.m_TangentArray)),
              m_UVArray(std::move(mesh.m_UVArray)),
              m_AlphaMask(std::move(mesh.m_AlphaMask)),
              m_ShadowAlphaMask(std::move(mesh.m_ShadowAlphaMask))
        {}

    public:
		std::unique_ptr<SceneObjectVertexArray>& GetVertexArray()
		{
			return m_VertexArray;
		}

		std::unique_ptr<SceneObjectIndexArray>& GetIndexArray()
		{
			return m_IndexArray;
		}

        int32_t GetTriangleCount() const 
        {
            return m_TrianglesCount;
        }

        int32_t GetVertexCount() const 
        {
            return m_VerticesCount;
        }

        const std::vector<int32_t>& GetFaceIndices()
        {
            return m_FaceIndices;
        }

		int32_t GetFaceIndex(int32_t index = 0)
		{
			assert(m_FaceIndices.size() > index);
			return m_FaceIndices[index];
		}
        
        PrimitiveType GetPrimitiveType() {return m_PrimitiveType;}

        // some operations on $m_UVArray$
        bool DoesUVExist()
        {
            return m_UVArray->GetDataSize() > 0;
        }
        const std::unique_ptr<SceneObjectVertexArray>& GetUVArray()
        {
            return m_UVArray;
        }

        // some operations on $m_AlphaMask$
        bool DoesAlphaMaskExist()
        {
            return (!!m_AlphaMask);
        }
        const std::shared_ptr<SceneObjectTexture>& GetAlphaMask()
        {
            return m_AlphaMask;
        }

        // Some operations on $m_ShadowAlphaMask$
        bool DoesShadowAlphaMaskExist()
        {
            return (!!m_ShadowAlphaMask);
        }
        const std::shared_ptr<SceneObjectTexture>& GetShadowAlphaMask()
        {
            return m_ShadowAlphaMask;
        }

        // Some operations on $m_NormalArray$
        bool DoesNormalExist()
        {
            return m_NormalArray->GetDataSize() > 0;
        }
        const std::unique_ptr<SceneObjectVertexArray>& GetNormalArray()
        {
            return m_NormalArray;
        }

        // Some operations on $m_TrangentArray$
        bool DoesTangentExist()
        {
            return m_TangentArray->GetDataSize() > 0;
        }
        const std::unique_ptr<SceneObjectVertexArray>& GetTangentArray()
        {
            return m_TangentArray;
        }

    protected:
        std::shared_ptr<Matrix4f> m_ObjectToWorld;

        PrimitiveType m_PrimitiveType;
        std::unique_ptr<SceneObjectVertexArray> m_VertexArray; // p
        std::unique_ptr<SceneObjectIndexArray>  m_IndexArray; // vertexIndices
        int32_t m_TrianglesCount; // nTriangles
        int32_t m_VerticesCount; // nVertices
        std::vector<int32_t> m_FaceIndices; // faceIndces
        std::unique_ptr<SceneObjectVertexArray> m_NormalArray; // n
        std::unique_ptr<SceneObjectVertexArray> m_TangentArray; // s
        std::unique_ptr<SceneObjectVertexArray> m_UVArray; // uv
        const std::shared_ptr<SceneObjectTexture> m_AlphaMask;
        const std::shared_ptr<SceneObjectTexture> m_ShadowAlphaMask;
    };
}