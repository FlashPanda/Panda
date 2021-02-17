#pragma once
#include "BaseSceneObject.hpp"
#include "SceneObjectMesh.hpp"

namespace Panda
{
    class SceneObjectGeometry : public BaseSceneObject
    {
        protected:
            std::vector<std::shared_ptr<SceneObjectMesh>> m_Mesh;
            bool        m_IsVisible;
            bool        m_IsShadow;
            bool        m_IsMotionBlur;
            SceneObjectCollisionType m_CollisionType;
            float       m_CollisionParameters[10];

        public:
            SceneObjectGeometry() : BaseSceneObject(SceneObjectType::kSceneObjectTypeGeometry), m_CollisionType(SceneObjectCollisionType::kSceneObjectCollisionTypeNone) {}

            void SetVisibility(bool visible) {m_IsVisible = visible;}
            const bool Visible() {return m_IsVisible;}
            void SetIfCastShadow(bool shadow) {m_IsShadow = shadow;}
            const bool CastShadow() {return m_IsShadow;}
            void SetIfMotionBlur(bool motionBlur) {m_IsMotionBlur = motionBlur;}
            const bool MotionBlur() {return m_IsMotionBlur;}
            void SetCollisionType(SceneObjectCollisionType collisionType) {m_CollisionType = collisionType;}
            const SceneObjectCollisionType CollisionType() const {return m_CollisionType;}
            void SetCollisionParameters(const float* param, int32_t count)
            {
                assert(count > 0 && count < 10);
                memcpy(m_CollisionParameters, param, sizeof(float) * count);
            }
            const float* CollisionParameters() const { return m_CollisionParameters;}

            void AddMesh(std::shared_ptr<SceneObjectMesh>& mesh) {m_Mesh.push_back(std::move(mesh));}
			uint32_t GetMeshCount() { return m_Mesh.size(); }
            const std::weak_ptr<SceneObjectMesh> GetMesh(uint32_t index = 0) {return m_Mesh.empty()? nullptr : m_Mesh[index];}
            const std::weak_ptr<SceneObjectMesh> GetMeshLOD(size_t lod) {return lod < m_Mesh.size()? m_Mesh[lod] : nullptr;}

            friend std::ostream& operator<<(std::ostream& out, const SceneObjectGeometry& obj);
    };
}