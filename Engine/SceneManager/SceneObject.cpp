#include "SceneObject.hpp"

namespace Panda
{
    std::ostream& operator<<(std::ostream& out, SceneObjectType type)
    {
        int32_t n = static_cast<int32_t> (type);
        n = to_endian_net<int32_t>(n);
        char* c = reinterpret_cast<char*>(&n);

        for (size_t i = 0; i < sizeof (int32_t); ++i)
        {
            out << *c++;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, IndexDataType type)
    {
        int32_t n = static_cast<int32_t>(type);
        n = to_endian_net<int32_t>(n);
        char* c = reinterpret_cast<char*>(&n);

        for (size_t i = 0; i < sizeof(int32_t); ++i)
        {
            out << *c++;
        }
        return out;
    }

    std::ostream& operator<<(std::ostream& out, VertexDataType type)
    {
        int32_t n = static_cast<int32_t>(type);
        n = to_endian_net<int32_t> (n);
        char* c = reinterpret_cast<char*>(&n);

        for (size_t i = 0; i < sizeof(int32_t); ++i)
        {
            out << *c++;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream&out, PrimitiveType type)
    {
        int32_t n = static_cast<int32_t>(type);
        n = to_endian_net<int32_t>(n);
        char* c = reinterpret_cast<char*>(&n);

        for(size_t i = 0; i < sizeof(int32_t); ++i)
        {
            out << *c++;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, CurveType type)
    {
        int32_t n = static_cast<int32_t>(type);
        n = to_endian_net<int32_t>(n);
        char* c = reinterpret_cast<char*>(&n);
        for(size_t i = 0; i < sizeof(int32_t); ++i)
        {
            out << *c++;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const BaseSceneObject& obj)
    {
        out << "SceneObject" << std::endl;
        out << "------------" << std::endl;
        out << "GUID: " << obj.m_Guid << std::endl;
        out << "Type: " << obj.m_Type << std::endl;

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectVertexArray& obj)
    {
        out << "Attribute: " << obj.m_Attribute << std::endl;
        out << "Morph Target Index: 0x" << obj.m_MorphTargetIndex << std::endl;
        out << "Data Type: " << obj.m_DataType << std::endl;
        out << "Data Size: 0x" << obj.m_DataSize << std::endl;
        out << "Data: ";
        for (size_t i = 0; i < obj.m_DataSize; ++i)
        {
            out << *(reinterpret_cast<const float*>(obj.m_pData) + i) << ' ';
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectIndexArray& obj)
    {
        out << "Material Index: 0x" << obj.m_MaterialIndex << std::endl;
        out << "Restart Index: 0x" << obj.m_RestartIndex << std::endl;
        out << "Data Type: " << obj.m_DataType << std::endl;
        out << "Data Size: 0x" << obj.m_DataSize << std::endl;
        out << "Data: ";

        for (size_t i = 0; i < obj.m_DataSize; ++i)
        {
            switch(obj.m_DataType)
            {
                case IndexDataType::kIndexDataTypeInt8:
                    out << "0x" << *(reinterpret_cast<const uint8_t*>(obj.m_pData) + i) << ' ';
                    break;
                case IndexDataType::kIndexDataTypeInt16:
                    out << "0x" << *(reinterpret_cast<const uint16_t*>(obj.m_pData) + i) << ' ';
                    break;
                case IndexDataType::kIndexDataTypeInt32:
                    out << "0x" << *(reinterpret_cast<const uint32_t*>(obj.m_pData) + i) << ' ';
                    break;
                case IndexDataType::kIndexDataTypeInt64:
                    out << "0x" << *(reinterpret_cast<const uint64_t*>(obj.m_pData) + i) << ' ';
                    break;
                default:
                    break;
            }
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectMesh& obj)
    {
        out << static_cast<const BaseSceneObject&>(obj) << std::endl;
        out << "Primitive Type: " << obj.m_PrimitiveType << std::endl;
        out << "This mesh contains 0x" << obj.m_VertexArray.size() << " vertex properties." << std::endl;
        for (size_t i = 0; i < obj.m_VertexArray.size(); ++i)
        {
            out << obj.m_VertexArray[i] << std::endl;
        }
        out << "This mesh contains 0x" << obj.m_IndexArray.size() << " index arrays." << std::endl;
        for(size_t i = 0; i < obj.m_IndexArray.size(); ++i)
        {
            out << obj.m_IndexArray[i] << std::endl;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectTexture& obj)
    {
        out << static_cast<const BaseSceneObject&>(obj) << std::endl;
        out << "Coord Index: " << obj.m_TexCoordIndex << std::endl;
        out << "Name: " << obj.m_TexCoordIndex << std::endl;
        if (obj.m_pImage)
            out << "Image: " << *obj.m_pImage << std::endl;
        
        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectMaterial& obj)
    {
        out << static_cast<const BaseSceneObject&>(obj) << std::endl;
        out << "Name: " << obj.m_Name << std::endl;
        out << "Albedo: " << obj.m_BaseColor << std::endl;
        out << "Metallic: " << obj.m_Metallic << std::endl;
        out << "Roughness: " << obj.m_Roughness << std::endl;
        out << "Normal: " << obj.m_Normal << std::endl;
        out << "Specular: " << obj.m_Specular << std::endl;
        out << "Ambient Occlusion: " << obj.m_AmbientOcclusion << std::endl;

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectGeometry& obj)
    {
        auto count = obj.m_Mesh.size();
        for (decltype(count) i = 0; i < count; ++i)
        {
            out << "Mesh: " << *obj.m_Mesh[i] << std::endl;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectLight& obj)
    {
        out << static_cast<const BaseSceneObject&>(obj) << std::endl;
        out << "Color: " << obj.m_LightColor << std::endl;
        out << "Intensity: " << obj.m_Intensity << std::endl;
        out << "Cast Shadows: " << obj.m_IsCastShadows << std::endl;

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectPointLight& obj)
    {
        out << static_cast<const SceneObjectLight&>(obj) << std::endl;
        out << "Light Type: Point" << std::endl;

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectSpotLight& obj)
    {
        out << static_cast<const SceneObjectLight&>(obj) << std::endl;
        out << "Light Type: Spot" << std::endl;
        out << "Cos of Inner Angle: " << obj.m_CosOfInnerAngle << std::endl;
        out << "Cos of Outer Angle: " << obj.m_CosOfOuterAngle << std::endl;
        return out;
    }

	std::ostream& operator<<(std::ostream& out, const SceneObjectInfiniteLight& obj)
	{
		out << static_cast<const SceneObjectLight&>(obj) << std::endl;
		out << "Light Type: Infinite" << std::endl;

		return out;
	}

	std::ostream& operator<<(std::ostream& out, const SceneObjectAreaLight& obj)
	{
		out << static_cast<const SceneObjectLight&>(obj) << std::endl;
		out << "Light Type: Area" << std::endl;

		return out;
	}

    std::ostream& operator<<(std::ostream& out, const SceneObjectCamera& obj)
    {
        out << static_cast<const BaseSceneObject&>(obj) << std::endl;
        out << "Aspect: " << obj.m_Aspect << std::endl;
        out << "Near Clip Distance: " << obj.m_NearClipDistance << std::endl;
        out << "Far Clip Disance: " << obj.m_FarClipDistance << std::endl;

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectOrthogonalCamera& obj)
    {
        out << static_cast<const SceneObjectCamera&>(obj) << std::endl;
        out << "Camera Type: Orthogonal" << std::endl;
        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectPerspectiveCamera& obj)
    {
        out << static_cast<const SceneObjectCamera&>(obj) << std::endl;
        out << "Camera Type: Perspective" << std::endl;
        out << "FOV: " << obj.m_Fov << std::endl;

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const SceneObjectTransform& obj)
    {
		if (obj.m_Type == SceneObjectType::kSceneObjectTypeTranslate)
			out << "------kSceneObjectTypeTranslate--------" << std::endl;
		else if (obj.m_Type == SceneObjectType::kSceneObjectTypeRotate)
			out << "------kSceneObjectTypeRotate------" << std::endl;
		else if (obj.m_Type == SceneObjectType::kSceneObjectTypeScale)
			out << "------kSceneObjectTypeScale------" << std::endl;
		else
			out << "------kSceneObjectTypeTransform------" << std::endl;
        out << "Transform Matrix: " << obj.m_Matrix << std::endl;
        out << "Is Object Local: " << obj.m_IsSceneObjectOnly << std::endl;

        return out;
    }

    float DefaultAttenFunc (float intensity, float distance)
    {
        return intensity / std::pow(1 + distance, 2.0f);
    }
}