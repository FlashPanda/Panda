#pragma once
#include "BaseSceneNode.hpp"
#include "Core/Medium.hpp"

namespace Panda
{
    class SceneGeometryNode : public BaseSceneNode
    {
        protected:
            bool    m_IsVisible;
            bool    m_IsShadow;
            bool    m_IsMotionBlur;
            std::vector<std::string> m_Materials;
            void*   m_pRigidBody = nullptr;
            MediumInterface m_MediumInterface;

        protected:
            virtual void Dump(std::ostream& out)
            {
                BaseSceneNode::Dump(out);
                out << "Visible: " << m_IsVisible << std::endl;
                out << "Shadow: " << m_IsShadow << std::endl;
                out << "Motion Blur: " << m_IsMotionBlur << std::endl;
                out << "Material(s): " << std::endl;
                for (auto material : m_Materials)
                {
                    out << material << std::endl;
                }
            }

        public:
            using BaseSceneNode::BaseSceneNode;

            void SetVisibility(bool visible) {m_IsVisible = visible;}
            const bool Visible() {return m_IsVisible;}
            void SetIfCastShadow(bool shadow) {m_IsShadow = shadow;}
            const bool CastShadow() {return m_IsShadow;}
            void SetIfMotionBlur(bool motionBlur) {m_IsMotionBlur = motionBlur;}
            const bool MotionBlur() {return m_IsMotionBlur;}
            void AddMaterialRef(const std::string& key) {m_Materials.push_back(key);}
            void AddMaterialRef(std::string&& key) {m_Materials.push_back(std::move(key));}
            virtual std::string GetMaterialRef(const size_t index)
            {
                if (index < m_Materials.size())
                    return m_Materials[index];
                else 
                    return std::string("default");
            }
            void SetMediumInterface(const MediumInterface& medium)
            {
                m_MediumInterface = medium;
            }
            const MediumInterface& GetMediumInterface()
            {
                return m_MediumInterface;
            }

            void LinkRigidBody(void* pRigidBody)
            {
                m_pRigidBody = pRigidBody;
            }

            void* UnlinkRigidBody()
            {
                void* pRigidBody = m_pRigidBody;
                m_pRigidBody = nullptr;

                return pRigidBody;
            }

            void* RigidBody() {return m_pRigidBody;}

    public:
        virtual Bounds3Df WorldBound();
        virtual bool Intersect(const Ray& r, SceneObjectSurfaceInteraction& interaction);
        virtual bool IntersectP(const Ray& r);
        //void ComputeScatteringFunctions(SceneObjectSurfaceInteraction& interaction);
    };
}