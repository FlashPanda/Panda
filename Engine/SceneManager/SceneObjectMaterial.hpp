#pragma once
#include <string>
#include "Core/Math/PandaMath.hpp"
#include "BaseSceneObject.hpp"
#include "SceneObjectTexture.hpp"
#include "SceneObjectTypeDef.hpp"

namespace Panda
{
    // PBR material
    class SceneObjectMaterial : public BaseSceneObject
    {
        protected:
            std::string         m_Name;
            Color               m_BaseColor;
            Parameter           m_Metallic;
            Parameter           m_Roughness;
            Normal              m_Normal;
            Color               m_Specular;
            Parameter           m_SpecularPower;
            Parameter           m_AmbientOcclusion;
            Color               m_Opacity;
            Color               m_Transparency;
            Color               m_Emission;

        public:
            SceneObjectMaterial() : BaseSceneObject(SceneObjectType::kSceneObjectTypeMaterial),
				m_Name(""), m_BaseColor(Vector4Df(1.0f)), m_Metallic(0.0f), m_Roughness(0.0f), m_Normal(Vector3Df({ 0.0f, 0.0f, 1.0f })),
                m_Specular(0.0f), m_SpecularPower(1.0f), m_AmbientOcclusion(1.0f), m_Opacity(1.0f), m_Transparency(0.0f), m_Emission(0.0f) {}
            SceneObjectMaterial(const char* name) : BaseSceneObject(SceneObjectType::kSceneObjectTypeMaterial), m_Name(name) {}
            SceneObjectMaterial(const std::string& name) : BaseSceneObject(SceneObjectType::kSceneObjectTypeMaterial), m_Name(name) {}
            SceneObjectMaterial(std::string&& name) : BaseSceneObject(SceneObjectType::kSceneObjectTypeMaterial), m_Name(std::move(name)) {}

            const std::string& GetName() const {return m_Name;}
            const Color& GetBaseColor() const {return m_BaseColor;}
            const Color& GetSpecularColor() const {return m_Specular;}
            const Parameter& GetSpecularPower() const {return m_SpecularPower;}
            const Normal& GetNormal() const { return m_Normal; }
            void SetName(const std::string& name) {m_Name = name;}
            void SetName(std::string&& name) {m_Name = std::move(name);}
            void SetColor(const std::string& attrib, const Vector4Df& color)
            {
                if (attrib == "diffuse")
                {
                    m_BaseColor = Color(color);
                }
                else if (attrib == "specular")
                {
                    m_Specular = Color(color);
                }
                else if (attrib == "emission")
                {
                    m_Emission = Color(color);
                }
                else if (attrib == "opacity")
                {
                    m_Opacity = Color(color);
                }
                else if (attrib == "transparency")
                {
                    m_Opacity = Color(color);
                }
            }
            void SetParam(const std::string& attrib, float param)
            {
                if (attrib == "specular_power")
                {
                    m_SpecularPower = Parameter(param);
                }
            }

            void SetTexture(const std::string& attrib, const std::string& textureName)
            {
                if (attrib == "diffuse")
                {
                    m_BaseColor = std::make_shared<SceneObjectTexture>(textureName);
                }
                else if (attrib == "specular")
                {
                    m_Specular = std::make_shared<SceneObjectTexture>(textureName);
                }
                else if (attrib == "specular_power")
                {
                    m_SpecularPower = std::make_shared<SceneObjectTexture>(textureName);
                }
                else if (attrib == "emission")
                {
                    m_Emission = std::make_shared<SceneObjectTexture>(textureName);
                }
                else if (attrib == "opacity")
                {
                    m_Opacity = std::make_shared<SceneObjectTexture>(textureName);
                }
                else if (attrib == "transparency")
                {
                    m_Transparency = std::make_shared<SceneObjectTexture>(textureName);
                }
                else if (attrib == "normal")
                {
                    m_Normal = std::make_shared<SceneObjectTexture>(textureName);
                }
            }

            void SetTexture(const std::string& attrib, const std::shared_ptr<SceneObjectTexture>& texture)
            {
                if (attrib == "diffuse")
                {
                    m_BaseColor = texture;
                }
                else if (attrib == "specular")
                {
                    m_Specular = texture;
                }
                else if (attrib == "specular_power")
                {
                    m_SpecularPower = texture;
                }
                else if (attrib == "opacity")
                {
                    m_Opacity = texture;
                }
                else if (attrib == "transparency")
                {
                    m_Transparency = texture;
                }
                else if (attrib == "normal")
                {
                    m_Normal = texture;
                }
            }

            void LoadTexture()
            {
                if (m_BaseColor.ValueMap)
                {
                    m_BaseColor.ValueMap->LoadTexture();
                }
            }

            friend std::ostream& operator<<(std::ostream& out, const SceneObjectMaterial& obj);
    };
}