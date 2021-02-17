#include "OGEX.hpp"

namespace Panda
{
    void OgexParser::ConvertOddlStructureToSceneNode(const ODDL::Structure& structure, std::shared_ptr<BaseSceneNode>& baseNode, Scene& scene)
    {
        std::shared_ptr<BaseSceneNode> node;

        switch(structure.GetStructureType())
        {
            case OGEX::kStructureMetric:
            {
                const OGEX::MetricStructure& _structure = dynamic_cast<const OGEX::MetricStructure&>(structure);
                auto _key = _structure.GetMetricKey();
                const ODDL::Structure* subStructure = _structure.GetFirstCoreSubnode();
                if (_key == "up")
                {
                    const ODDL::DataStructure<ODDL::StringDataType>* dataStructure = static_cast<const ODDL::DataStructure<ODDL::StringDataType> *>(subStructure);
                    auto axisName = dataStructure->GetDataElement(0);
                    if (axisName == "y")
                        m_UpIsYAxis = true;
                    else
                        m_UpIsYAxis = false;
                }
                
                return;
            }
            case OGEX::kStructureNode:
            {
                node = std::make_shared<SceneEmptyNode>(structure.GetStructureName());
                break;
            }
            case OGEX::kStructureBoneNode:
            {
                auto _node = std::make_shared<SceneBoneNode>(structure.GetStructureName());
                std::string _key = structure.GetStructureName();
                scene.BoneNodes.emplace(_key, _node);
                node = _node;
                break;
            }
            case OGEX::kStructureGeometryNode:
            {
                std::string _key = structure.GetStructureName();
                auto _node = std::make_shared<SceneGeometryNode>(_key);
                const OGEX::GeometryNodeStructure& _structure = dynamic_cast<const OGEX::GeometryNodeStructure&>(structure);

                _node->SetVisibility(_structure.GetVisibleFlag());
                _node->SetIfCastShadow(_structure.GetShadowFlag());
                _node->SetIfMotionBlur(_structure.GetMotionBlurFlag());

                // ref scene objects
                _key = _structure.GetObjectStructure()->GetStructureName();
                _node->AddSceneObjectRef(_key);

                // ref materials
                auto materials = _structure.GetMaterialStructureArray();
                auto materials_count = materials.GetElementCount();
                for (auto i = 0; i < materials_count; ++i)
                {
                    auto material = materials[i];
                    std::string name = material->GetStructureName();
                    _node->AddMaterialRef(name);
                }

                std::string name = _structure.GetNodeName();
                scene.LUTNameGeometryNode.emplace(name, _node);
                scene.GeometryNodes.emplace(_key, _node);

                node = _node;
                break;
            }
            case OGEX::kStructureLightNode:
            {
                auto _node = std::make_shared<SceneLightNode>(structure.GetStructureName());
                const OGEX::LightNodeStructure& _structure = dynamic_cast<const OGEX::LightNodeStructure&>(structure);

                _node->SetIfCastShadow(_structure.GetShadowFlag());

                // ref scene objects
                std::string _key = _structure.GetObjectStructure()->GetStructureName();
                _node->AddSceneObjectRef(_key);

                scene.LightNodes.emplace(_key, _node);

                node = _node;
                break;
            }
            case OGEX::kStructureCameraNode:
            {
                auto _node = std::make_shared<SceneCameraNode>(structure.GetStructureName());
                const OGEX::CameraNodeStructure& _structure = dynamic_cast<const OGEX::CameraNodeStructure&>(structure);

                // ref scene objects
                std::string _key = _structure.GetObjectStructure()->GetStructureName();
                _node->AddSceneObjectRef(_key);

                scene.CameraNodes.emplace(_key, _node);

                node = _node;
                break;
            }
            case OGEX::kStructureGeometryObject:
            {
                const OGEX::GeometryObjectStructure& _structure = dynamic_cast<const OGEX::GeometryObjectStructure&>(structure);
                std::string _key = _structure.GetStructureName();
                auto _object = std::make_shared<SceneObjectGeometry>();

                // properties
                _object->SetVisibility(_structure.GetVisibleFlag());
                _object->SetIfCastShadow(_structure.GetShadowFlag());
                _object->SetIfMotionBlur(_structure.GetMotionBlurFlag());

                // extensions
                //// collision shape
                ODDL::Structure* extension = _structure.GetFirstExtensionSubnode();
                while (extension)
                {
                    const OGEX::ExtensionStructure* _extension = dynamic_cast<const OGEX::ExtensionStructure*>(extension);
                    auto _appid = _extension->GetApplicationString();
                    if (_appid == "MyGameEngine")
                    {
                        auto _type = _extension->GetTypeString();
                        if (_type == "collision")
                        {
                            const ODDL::Structure* subStructure = _extension->GetFirstCoreSubnode();
                            const ODDL::DataStructure<ODDL::StringDataType>* dataStructure1 = static_cast<const ODDL::DataStructure<ODDL::StringDataType>*>(subStructure);
                            auto collisionType = dataStructure1->GetDataElement(0);

                            subStructure = _extension->GetLastSubnode();
                            const ODDL::DataStructure<ODDL::FloatDataType>* dataStructure2 = static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(subStructure);
                            auto elementCount = dataStructure2->GetDataElementCount();
                            float* _data = (float*)&dataStructure2->GetDataElement(0);
                            if (collisionType == "plane")
                            {
                                _object->SetCollisionType(SceneObjectCollisionType::kSceneObjectCollisionTypePlane);
                                _object->SetCollisionParameters(_data, elementCount);
                            }
                            else if (collisionType == "sphere")
                            {
                                _object->SetCollisionType(SceneObjectCollisionType::kSceneObjectCollisionTypeSphere);
                                _object->SetCollisionParameters(_data, elementCount);
                            }
                            else if (collisionType == "box")
                            {
                                _object->SetCollisionType(SceneObjectCollisionType::kSceneObjectCollisionTypeBox);
                                _object->SetCollisionParameters(_data, elementCount);
                            }
                            break;
                        }
                    }

                    extension = extension->Next();
                }

                // meshs
                const ODDL::Map<OGEX::MeshStructure>* _meshs = _structure.GetMeshMap();
                int32_t _count = _meshs->GetElementCount();
                for (int32_t i = 0; i < _count; ++i)
                {
                    const OGEX::MeshStructure* _mesh = (*_meshs)[i];
                    std::shared_ptr<SceneObjectMesh> mesh(new SceneObjectMesh());
                    const std::string _primitive_type = static_cast<const char*>(_mesh->GetMeshPrimitive());
                    if (_primitive_type == "points")
                        mesh->SetPrimitiveType(PrimitiveType::kPrimitiveTypePointList);
                    else if (_primitive_type == "lines")
                        mesh->SetPrimitiveType(PrimitiveType::kPrimitiveTypeLineList);
                    else if (_primitive_type == "line_strip")
                        mesh->SetPrimitiveType(PrimitiveType::kPrimitiveTypeLineList);
                    else if (_primitive_type == "triangles")
                        mesh->SetPrimitiveType(PrimitiveType::kPrimitiveTypeTriList);
                    else if (_primitive_type == "triangle_strip")
                        mesh->SetPrimitiveType(PrimitiveType::kPrimitiveTypeTriStrip);
                    else if (_primitive_type == "quads")
                        mesh->SetPrimitiveType(PrimitiveType::kPrimitiveTypeQuadList);
                    else {
                        // not supported
                        mesh.reset();
                    }

                    if (mesh)
                    {
                        const ODDL::Structure* sub_structure = _mesh->GetFirstSubnode();
                        while (sub_structure)
                        {
                            switch (sub_structure->GetStructureType())
                            {
                                case OGEX::kStructureVertexArray:
                                {
                                    const OGEX::VertexArrayStructure* _v = dynamic_cast<const OGEX::VertexArrayStructure*>(sub_structure);
                                    const char* attr = _v->GetArrayAttrib();
                                    auto morph_index = _v->GetMorphIndex();

                                    const ODDL::Structure* _data_structure = _v->GetFirstCoreSubnode();
                                    const ODDL::DataStructure<FloatDataType>* dataStructure = dynamic_cast<const ODDL::DataStructure<FloatDataType>*>(_data_structure);

                                    auto arraySize = dataStructure->GetArraySize();
                                    auto elementCount = dataStructure->GetDataElementCount();
                                    const void*_data = &dataStructure->GetDataElement(0);
                                    void* data = new float[elementCount];
                                    size_t buf_size = sizeof(float) * elementCount;
                                    memcpy(data, _data, buf_size);
                                    VertexDataType vertexDataType;
                                    switch(arraySize)
                                    {
                                        case 1:
                                            vertexDataType = VertexDataType::kVertexDataTypeFloat1;
                                            break;
                                        case 2:
                                            vertexDataType = VertexDataType::kVertexDataTypeFloat2;
                                            break;
                                        case 3:
                                            vertexDataType = VertexDataType::kVertexDataTypeFloat3;
                                            break;
                                        case 4:
                                            vertexDataType = VertexDataType::kVertexDataTypeFloat4;
                                            break;
                                        default:
                                            continue;
                                    }
                                    SceneObjectVertexArray& _v_array = *new SceneObjectVertexArray(attr, morph_index, vertexDataType, data, elementCount);
                                    mesh->AddVertexArray(std::move(_v_array));

                                    break;
                                }

                                case OGEX::kStructureIndexArray:
                                {
                                    const OGEX::IndexArrayStructure* _i = dynamic_cast<const OGEX::IndexArrayStructure*>(sub_structure);
                                    auto material_index = _i->GetMaterialIndex();
                                    auto restart_index = _i->GetRestartIndex();
                                    const ODDL::Structure* _data_structure = _i->GetFirstCoreSubnode();
                                    ODDL::StructureType type = _data_structure->GetStructureType();
                                    int32_t elementCount = 0;
                                    const void* _data = nullptr;
                                    IndexDataType index_type = IndexDataType::kIndexDataTypeInt16;
                                    switch (type)
                                    {
                                        case ODDL::kDataUnsignedInt8:
                                        {
                                            index_type = IndexDataType::kIndexDataTypeInt8;
                                            const ODDL::DataStructure<UnsignedInt8DataType>* dataStructure = dynamic_cast<const ODDL::DataStructure<UnsignedInt8DataType>*>(_data_structure);
                                            elementCount = dataStructure->GetDataElementCount();
                                            _data = &dataStructure->GetDataElement(0);
                                            break;
                                        }
                                        case ODDL::kDataUnsignedInt16:
                                        {
                                            index_type = IndexDataType::kIndexDataTypeInt16;
                                            const ODDL::DataStructure<UnsignedInt16DataType>* dataStructure = dynamic_cast<const ODDL::DataStructure<UnsignedInt16DataType>*>(_data_structure);
                                            elementCount = dataStructure->GetDataElementCount();
                                            _data = &dataStructure->GetDataElement(0);
                                            break;
                                        }
                                        case ODDL::kDataUnsignedInt32:
                                        {
                                            index_type = IndexDataType::kIndexDataTypeInt32;
                                            const ODDL::DataStructure<UnsignedInt32DataType>* dataStructure = dynamic_cast<const ODDL::DataStructure<UnsignedInt32DataType>*>(_data_structure);
                                            elementCount = dataStructure->GetDataElementCount();
                                            _data = &dataStructure->GetDataElement(0);
                                            break;
                                        }
                                        case ODDL::kDataUnsignedInt64:
                                        {
                                            index_type = IndexDataType::kIndexDataTypeInt64;
                                            const ODDL::DataStructure<UnsignedInt64DataType>* dataStructure = dynamic_cast<const ODDL::DataStructure<UnsignedInt64DataType>*>(_data_structure);
                                            elementCount = dataStructure->GetDataElementCount();
                                            _data = &dataStructure->GetDataElement(0);
                                            break;
                                        }
                                        default:
                                            break;
                                    }

                                    int32_t data_size = 0;
                                    switch (index_type)
                                    {
                                        case IndexDataType::kIndexDataTypeInt8:
                                            data_size = 1;
                                            break;
                                        case IndexDataType::kIndexDataTypeInt16:
                                            data_size = 2;
                                            break;
                                        case IndexDataType::kIndexDataTypeInt32:
                                            data_size = 4;
                                            break;
                                        case IndexDataType::kIndexDataTypeInt64:
                                            data_size = 8;
                                            break;
                                        default:
                                            break;
                                    }

                                    size_t buf_size = elementCount * data_size;
                                    void* data = new uint8_t[buf_size];
                                    memcpy(data, _data, buf_size);
                                    SceneObjectIndexArray& _i_array = *new SceneObjectIndexArray(material_index, restart_index, index_type, data, elementCount);
                                    mesh->AddIndexArray(std::move(_i_array));

                                    break;
                                }

                                default:
                                    // ignore it
                                    break;
                            }

                            sub_structure = sub_structure->Next();
                        }

                        _object->AddMesh(mesh);
                    }
                }

                scene.Geometries[_key] = _object;
                return;
            }
            case OGEX::kStructureTransform:
            {
                int32_t index, count;
                const OGEX::TransformStructure& _structure = dynamic_cast<const OGEX::TransformStructure&>(structure);
                bool object_flag = _structure.GetObjectFlag();
                Matrix4f matrix;
                std::shared_ptr<SceneObjectTransform> transform;

                auto _key = _structure.GetStructureName();
                count = _structure.GetTransformCount();
                for (index = 0; index < count; ++index)
                {
                    const float* data = _structure.GetTransform(index);
                    matrix = data;
                    if (!m_UpIsYAxis)
                    {
                        // commented out due to camera is in same coordinations
                        // so no need to exchange.
                        // TODO: EXCHANGE y and z
                        // ExchangeYandZ(matrix)
                    }
                    transform = std::make_shared<SceneObjectTransform>(matrix, object_flag);
                    baseNode->AppendTransform(_key, std::move(transform));
                    //baseNode->AppendTransform(transform);
                }
                return;
            }
            case OGEX::kStructureTranslation:
            {
                const OGEX::TranslationStructure& _structure = dynamic_cast<const OGEX::TranslationStructure&>(structure);
                bool objectFlag = _structure.GetObjectFlag();
                std::shared_ptr<SceneObjectTranslation> translation;

                auto kind = _structure.GetTranslationKind();
                auto data = _structure.GetTranslation();
                if (kind == "xyz")
                {
                    translation = std::make_shared<SceneObjectTranslation>(data[0], data[1], data[2]);
                }
                else 
                {
                    translation = std::make_shared<SceneObjectTranslation>(kind[0], data[0]);
                }
                auto _key = _structure.GetStructureName();
                baseNode->AppendTransform(_key, std::move(translation));
                return;
            }
            case OGEX::kStructureRotation:
            {
                const OGEX::RotationStructure& _structure = dynamic_cast<const OGEX::RotationStructure&>(structure);
                bool object_flag = _structure.GetObjectFlag();
                std::shared_ptr<SceneObjectRotation> rotation;

                auto kind = _structure.GetRotationKind();
                auto data = _structure.GetRotation();
                if (kind == "x")
                {
                    rotation = std::make_shared<SceneObjectRotation>('x', data[0], object_flag);
                }
                else if (kind == "y")
                {
                    rotation = std::make_shared<SceneObjectRotation>('y', data[0], object_flag);
                }
                else if (kind == "z")
                {
                    rotation = std::make_shared<SceneObjectRotation>('z', data[0], object_flag);
                }
                else if (kind == "axis")
                {
                    rotation = std::make_shared<SceneObjectRotation>(Vector3Df({ data[0], data[1], data[2] }), data[3], object_flag);
                }
                else if (kind == "quaternion")
                {
                    rotation = std::make_shared<SceneObjectRotation>(Quaternion({ data[0], data[1], data[2], data[3] }), object_flag);
                }

                auto _key = _structure.GetStructureName();
                baseNode->AppendTransform(_key, std::move(rotation));
            }
            return;
            case OGEX::kStructureScale:
            {
                const OGEX::ScaleStructure& _structure = dynamic_cast<const OGEX::ScaleStructure&>(structure);
                bool object_flag = _structure.GetObjectFlag();
                std::shared_ptr<SceneObjectScale> scale;

                auto kind = _structure.GetScaleKind();
                auto data = _structure.GetScale();
                if (kind == "x")
                {
                    scale = std::make_shared<SceneObjectScale>('x', data[0], object_flag);
                }
                else if (kind == "y")
                {
                    scale = std::make_shared<SceneObjectScale>('y', data[0], object_flag);
                }
                else if (kind == "z")
                {
                    scale = std::make_shared<SceneObjectScale>('z', data[0], object_flag);
                }
                else if (kind == "xyz")
                {
                    scale = std::make_shared<SceneObjectScale>(data[0], data[1], data[2], object_flag);
                }
                auto _key = _structure.GetStructureName();
                baseNode->AppendTransform(_key, std::move(scale));
            }
            return;
            case OGEX::kStructureMaterial:
            {
                const OGEX::MaterialStructure& _structure = dynamic_cast<const OGEX::MaterialStructure&>(structure);
                std::string material_name;
                const char* _name = _structure.GetMaterialName();
                if (_name)
                {
                    material_name = _name;
                }
                std::string _key = _structure.GetStructureName();
                auto material = std::make_shared<SceneObjectMaterial>();
                material->SetName(material_name);
                
                const ODDL::Structure* _sub_structure = _structure.GetFirstCoreSubnode();
                while (_sub_structure)
                {
                    std::string attrib, textureName;
                    Vector4Df color;
                    float param;
                    switch(_sub_structure->GetStructureType())
                    {
                        case OGEX::kStructureColor:
                        {
                            attrib = dynamic_cast<const OGEX::ColorStructure*>(_sub_structure)->GetAttribString();
                            color = dynamic_cast<const OGEX::ColorStructure*>(_sub_structure)->GetColor();
                            material->SetColor(attrib, color);
                            break;
                        }
                        case OGEX::kStructureParam:
                        {
                            attrib = dynamic_cast<const OGEX::ParamStructure*>(_sub_structure)->GetAttribString();
                            param = dynamic_cast<const OGEX::ParamStructure*>(_sub_structure)->GetParam();
                            material->SetParam(attrib, param);
                            break;
                        }
                        case OGEX::kStructureTexture:
                        {
                            attrib = dynamic_cast<const OGEX::TextureStructure*>(_sub_structure)->GetAttribString();
                            textureName = dynamic_cast<const OGEX::TextureStructure*>(_sub_structure)->GetTextureName();
                            material->SetTexture(attrib, textureName);
                            break;
                        }
                        default:
                            break;
                    }


                    _sub_structure = _sub_structure->Next();
                }
                
                scene.Materials[_key] = material;
                return;
            }
            case OGEX::kStructureLightObject:
            {
                const OGEX::LightObjectStructure& _structure = dynamic_cast<const OGEX::LightObjectStructure&>(structure);
                const char* typeStr = _structure.GetTypeString();
                const bool isShadow = _structure.GetShadowFlag();
                std::string _key = _structure.GetStructureName();
                std::shared_ptr<SceneObjectLight> light;

				if (!strncmp(typeStr, "infinite", 8))
				{
					light = std::make_shared<SceneObjectInfiniteLight>();
				}
				else if (!strncmp(typeStr, "point", 5))
				{
					light = std::make_shared<SceneObjectPointLight>();
				}
				else if (!strncmp(typeStr, "spot", 4))
				{
					light = std::make_shared<SceneObjectSpotLight>();
				}
				else if (!strncmp(typeStr, "area", 4))
				{
					light = std::make_shared<SceneObjectAreaLight>();
				}
				else
					assert(0);
                light->SetIfCastShadow(isShadow);

                const ODDL::Structure* subStructure = _structure.GetFirstCoreSubnode();
                while (subStructure)
                {
                    std::string attrib, textureName;
                    Vector4Df color;
                    float param;
                    switch(subStructure->GetStructureType())
                    {
                        case OGEX::kStructureColor:
                        {
                            attrib = dynamic_cast<const OGEX::ColorStructure*>(subStructure)->GetAttribString();
                            color = dynamic_cast<const OGEX::ColorStructure*>(subStructure)->GetColor();
                            light->SetColor(attrib, color);
                            break;
                        }
                        case OGEX::kStructureParam:
                        {
                            attrib = dynamic_cast<const OGEX::ParamStructure*>(subStructure)->GetAttribString();
                            param = dynamic_cast<const OGEX::ParamStructure*>(subStructure)->GetParam();
                            light->SetParam(attrib, param);
                            break;
                        }
                        case OGEX::kStructureTexture:
                        {
                            attrib = dynamic_cast<const OGEX::TextureStructure*>(subStructure)->GetAttribString();
                            textureName = dynamic_cast<const OGEX::TextureStructure*>(subStructure)->GetTextureName();
                            light->SetTexture(attrib, textureName);
                            break;
                        }
                        case OGEX::kStructureAtten:
                        {
                            auto atten = dynamic_cast<const OGEX::AttenStructure*>(subStructure);
                            AttenCurve curve;
                            if (atten->GetCurveType() == "linear")
                            {
                                curve.type = AttenCurveType::kAttenLinear;
                                curve.u.linearParams.beginAtten = atten->GetBeginParam();
                                curve.u.linearParams.endAtten = atten->GetEndParam();
                            }
                            else if (atten->GetCurveType() == "smooth")
                            {
                                curve.type = AttenCurveType::kAttenSmooth;
                                curve.u.smoothParams.beginAtten = atten->GetBeginParam();
                                curve.u.smoothParams.endAtten = atten->GetEndParam();
                            }
                            else if (atten->GetCurveType() == "inverse")
                            {
                                curve.type = AttenCurveType::kAttenInverse;
                                curve.u.inverseParams.scale = atten->GetScaleParam();
                                curve.u.inverseParams.offset = atten->GetOffsetParam();
                                curve.u.inverseParams.kl = atten->GetLinearParam();
                                curve.u.inverseParams.kc = atten->GetConstantParam();
                            }
                            else if (atten->GetCurveType() == "inverse_square")
                            {
                                curve.type = AttenCurveType::kAttenInverseSquare;
                                curve.u.inverseSquareParams.scale = atten->GetScaleParam();
                                curve.u.inverseSquareParams.offset = atten->GetOffsetParam();
                                curve.u.inverseSquareParams.kq = atten->GetQuadraticParam();
                                curve.u.inverseSquareParams.kl = atten->GetLinearParam();
                                curve.u.inverseSquareParams.kc = atten->GetConstantParam();
                            }

                            if (atten->GetAttenKind() == "angle")
                            {
                                auto _light = dynamic_pointer_cast<SceneObjectSpotLight>(light);
                                _light->SetAngleAttenuation(curve);
                            }
                            else if (atten->GetAttenKind() == "cos_angle")
                            {
                                // TODO: mark the angle in cos value instead of rad
                                auto _light = dynamic_pointer_cast<SceneObjectSpotLight>(light);
                                _light->SetAngleAttenuation(curve);
                            }
                            else 
                            {
                                light->SetDistanceAttenuation(curve);
                            }
                            
                            break;
                        }
                        default:
                            break;
                    }

                    subStructure = subStructure->Next();
                }

				// extensions
				ODDL::Structure* extension = _structure.GetFirstExtensionSubnode();
				while (extension)
				{
					const OGEX::ExtensionStructure* _extension = dynamic_cast<const OGEX::ExtensionStructure*>(extension);
					auto _appid = _extension->GetApplicationString();
					if (_appid == "PandaGameEngine")
					{
						auto _type = _extension->GetTypeString();
						if (_type == "area_light")
						{
							const ODDL::Structure* subStructure = _extension->GetFirstCoreSubnode();
							const ODDL::DataStructure<ODDL::FloatDataType>* dataStructure1 = static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(subStructure);
							auto elementCount = dataStructure1->GetDataElementCount();
							assert(elementCount == 2);
							auto width = dataStructure1->GetDataElement(0);
							auto height = dataStructure1->GetDataElement(1);

							auto _light = std::dynamic_pointer_cast<SceneObjectAreaLight>(light);
							_light->SetDimension( width, height );
						}
					}
					extension = extension->Next();
				}

                scene.Lights[_key] = light;
                return;
            }
            case OGEX::kStructureCameraObject:
            {
                const OGEX::CameraObjectStructure& _structure = dynamic_cast<const OGEX::CameraObjectStructure&>(structure);
                std::string _key = _structure.GetStructureName();
                auto camera = std::make_shared<SceneObjectPerspectiveCamera>(_key);

                const ODDL::Structure* subStructure = _structure.GetFirstCoreSubnode();
                while(subStructure)
                {
                    std::string attrib, textureName;
                    Vector4Df color;
                    float param;
                    switch(subStructure->GetStructureType())
                    {
                        case OGEX::kStructureColor:
                        {
                            attrib = dynamic_cast<const OGEX::ColorStructure*>(subStructure)->GetAttribString();
                            color = dynamic_cast<const OGEX::ColorStructure*>(subStructure)->GetColor();
                            camera->SetColor(attrib, color);
                            break;
                        }
                        case OGEX::kStructureParam:
                        {
                            attrib = dynamic_cast<const OGEX::ParamStructure*>(subStructure)->GetAttribString();
                            param = dynamic_cast<const OGEX::ParamStructure*>(subStructure)->GetParam();
                            camera->SetParam(attrib, param);
                            break;
                        }
                        case OGEX::kStructureTexture:
                        {
                            attrib = dynamic_cast<const OGEX::TextureStructure*>(subStructure)->GetAttribString();
                            textureName = dynamic_cast<const OGEX::TextureStructure*>(subStructure)->GetTextureName();
                            camera->SetTexture(attrib, textureName);
                            break;
                        }
                        default:
                            break;
                    }

                    subStructure = subStructure->Next();
                }
                scene.Cameras[_key] = camera;
                return;
            }
            case OGEX::kStructureAnimation:
            {
                //const OGEX::AnimationStructure& _structure = dynamic_cast<const OGEX::AnimationStructure&>(structure);
                //auto clipIndex = _structure.GetClipIndex();
                //std::shared_ptr<SceneObjectAnimationClip> clip = std::make_shared<SceneObjectAnimationClip>(clipIndex);

                //const ODDL::Structure* _sub_structure = _structure.GetFirstCoreSubnode();
                //while(_sub_structure)
                //{
                //    switch(_sub_structure->GetStructureType())
                //    {
                //        case OGEX::kStructureTrack:
                //        {
                //            const OGEX::TrackStructure& trackStructure = dynamic_cast<const OGEX::TrackStructure&>(*_sub_structure);
                //            const OGEX::TimeStructure& timeStructure = dynamic_cast<const OGEX::TimeStructure&>(*trackStructure.GetTimeStructure());
                //            const OGEX::ValueStructure& valueStructure = dynamic_cast<const OGEX::ValueStructure&>(*trackStructure.GetValueStructure());
                //            auto ref = trackStructure.GetTargetRef();
                //            std::string _key(*ref.GetNameArray());
                //            std::shared_ptr<SceneObjectTransform> trans;
                //            trans = baseNode->GetTransform(_key);
                //            std::shared_ptr<SceneObjectTrack> track;

                //            auto timeKeyValue = timeStructure.GetKeyValueStructure();
                //            auto timeKeyDataCount = timeStructure.GetKeyDataElementCount();
                //            auto dataStructure = static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(timeKeyValue->GetFirstCoreSubnode());
                //            auto timeArraySize = dataStructure->GetArraySize();
                //            const float* timeKnots = &dataStructure->GetDataElement(0);
                //            assert(timeArraySize == 0);

                //            auto valueKeyValue = valueStructure.GetKeyValueStructure();
                //            auto valueKeyDataCount = valueStructure.GetKeyDataElementCount();
                //            dataStructure = static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(valueKeyValue->GetFirstCoreSubnode());
                //            auto valueArraySize = dataStructure->GetArraySize();
                //            const float* valueKnots = &dataStructure->GetDataElement(0);
                //            std::shared_ptr<CurveBase> timeCurve;
                //            std::shared_ptr<CurveBase> valueCurve;
                //            SceneObjectTrackType type = SceneObjectTrackType::kScalar;
                //            if(timeStructure.GetCurveType() == "bezier")
                //            {
                //                auto keyIncomingControl = timeStructure.GetKeyControlStructure(0);
                //                auto keyOutgoingControl = timeStructure.GetKeyControlStructure(1);
                //                dataStructure = 
                //                    static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(keyIncomingControl->GetFirstCoreSubnode());
                //                const float* inCp = &dataStructure->GetDataElement(0);
                //                dataStructure =
                //                    static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(keyOutgoingControl->GetFirstCoreSubnode());
                //                const float* outCp = &dataStructure->GetDataElement(0);
                //                timeCurve = std::make_shared<Bezier<float, float>>(timeKnots, inCp, outCp, timeKeyDataCount);
                //            }
                //            else {
                //                timeCurve = std::make_shared<Linear<float, float>>(timeKnots, timeKeyDataCount);
                //            }
                //            if (valueStructure.GetCurveType() == "bezier")
                //            {
                //                auto keyIncomingControl = valueStructure.GetKeyControlStructure(0);
                //                auto keyOutgoingControl = valueStructure.GetKeyControlStructure(1);
                //                
                //                dataStructure = 
                //                    static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(keyIncomingControl->GetFirstCoreSubnode());
                //                const float* inCp = &dataStructure->GetDataElement(0);
                //                dataStructure = 
                //                    static_cast<const ODDL::DataStructure<ODDL::FloatDataType>*>(keyOutgoingControl->GetFirstCoreSubnode());
                //                const float* outCp = &dataStructure->GetDataElement(0);
                //                
                //                switch(valueArraySize)
                //                {
                //                    case 0:
                //                    case 1:
                //                    {
                //                        valueCurve = std::make_shared<Bezier<float, float>>(
                //                            valueKnots,
                //                            inCp,
                //                            outCp,
                //                            valueKeyDataCount
                //                        );
                //                        type = SceneObjectTrackType::kScalar;
                //                        break;
                //                    }
                //                    case 3:
                //                    {
                //                        valueCurve = std::make_shared<Bezier<Vector3Df, Vector3Df>>(
                //                            reinterpret_cast<const Vector3Df*>(valueKnots),
                //                            reinterpret_cast<const Vector3Df*>(inCp),
                //                            reinterpret_cast<const Vector3Df*>(outCp),
                //                            valueKeyDataCount
                //                        );
                //                        type = SceneObjectTrackType::kVector3;
                //                        break;
                //                    }
                //                    case 4:
                //                    {
                //                        valueCurve = std::make_shared<Bezier<Vector4Df, float>>(
                //                            reinterpret_cast<const Vector4Df*>(valueKnots),
                //                            reinterpret_cast<const Vector4Df*>(inCp),
                //                            reinterpret_cast<const Vector4Df*>(outCp),
                //                            valueKeyDataCount
                //                        );
                //                        type = SceneObjectTrackType::kQuaternion;
                //                        break;
                //                    }
                //                    case 16:
                //                    {
                //                        valueCurve = std::make_shared<Bezier<Matrix4f, float>>(
                //                            reinterpret_cast<const Matrix4f*>(valueKnots),
                //                            reinterpret_cast<const Matrix4f*>(inCp),
                //                            reinterpret_cast<const Matrix4f*>(outCp),
                //                            valueKeyDataCount
                //                        );
                //                        break;
                //                    }
                //                    default:
                //                        assert(0);
                //                }
                //            }
                //            else // default to linear
                //            {
                //                switch(valueArraySize)
                //                {
                //                    case 0:
                //                    case 1:
                //                        {
                //                            valueCurve = std::make_shared<Linear<float, float>>(
                //                                valueKnots,
                //                                valueKeyDataCount
                //                            );
                //                            type = SceneObjectTrackType::kScalar;
                //                            break;
                //                        }
                //                    case 3:
                //                    {
                //                        valueCurve = std::make_shared<Linear<Vector3Df, Vector3Df>>(
                //                            reinterpret_cast<const Vector3Df*>(valueKnots),
                //                            valueKeyDataCount);
                //                        type = SceneObjectTrackType::kVector3;
                //                        break;
                //                    }
                //                    case 4:
                //                    {
                //                        valueCurve = std::make_shared<Linear<Vector4Df, float>>(
                //                            reinterpret_cast<const Vector4Df*>(valueKnots),
                //                            valueKeyDataCount);
                //                        type = SceneObjectTrackType::kQuaternion;
                //                        break;
                //                    }
                //                    case 16:
                //                    {
                //                        valueCurve = std::make_shared<Linear<Matrix4f, float>>(
                //                            reinterpret_cast<const Matrix4f*>(valueKnots),
                //                            valueKeyDataCount
                //                        );
                //                        type = SceneObjectTrackType::kMatrix;
                //                        break;
                //                    }
                //                    default:
                //                        assert(0);
                //                }
                //            }
                //            
                //            track = std::make_shared<SceneObjectTrack>(trans, timeCurve, valueCurve, type);
                //            clip->AddTrack(track);
                //        }
                //        default:
                //            ;
                //    }

                //    _sub_structure = _sub_structure->Next();
                //}

                //baseNode->AttachAnimationClip(clipIndex, clip);

                //// register the node to animatable node LUT
                //scene.AnimatableNodes.push_back(baseNode);
            }
            return;
            default:
                // just ignore it and finish
                return;
        }

        const ODDL::Structure* sub_structure = structure.GetFirstSubnode();
        while (sub_structure)
        {
            ConvertOddlStructureToSceneNode(*sub_structure, node, scene);

            sub_structure = sub_structure->Next();
        }

        baseNode->AppendChild(std::move(node));
    }

    std::unique_ptr<Scene> OgexParser::Parse(const std::string& buf)
    {
        std::unique_ptr<Scene> pScene(new Scene("Root"));
        OGEX::OpenGexDataDescription openGexDataDescription;

        ODDL::DataResult result = openGexDataDescription.ProcessText(buf.c_str());
        if (result == ODDL::kDataOkay)
        {
            const ODDL::Structure* structure = openGexDataDescription.GetRootStructure()->GetFirstSubnode();
            while (structure)
            {
                ConvertOddlStructureToSceneNode(*structure, pScene->SceneGraph, *pScene);

                structure = structure->Next();
            }
        }

        return pScene;
    }
}