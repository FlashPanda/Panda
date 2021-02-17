#include <iostream>
#include <fstream>
#include "OpenGLGraphicsManager.hpp"
#include "AssetLoader.hpp"
#include "Interface/IApplication.hpp"
#include "Utility.hpp"
#include "SceneManager.hpp"

extern struct gladGLversionStruct GLVersion;

namespace Panda
{
    int OpenGLGraphicsManager::Initialize()
    {
        int result;

        result = GraphicsManager::Initialize();

        if (result)
            return result;

        result = gladLoadGL();
        if (!result) {
            std::cerr << "OpenGL load failed!" << std::endl;
            result = -1; 
        } else {
            result = 0;
            std::cout << "OpenGL Version " << GLVersion.major << "." << GLVersion.minor << " loaded" << std::endl;

            if (GLAD_GL_VERSION_3_3) {
                // Set the depth buffer to be entirely cleared to 1.0 values.
                glClearDepth(1.0f);

                // Enable depth testing.
                glEnable(GL_DEPTH_TEST);

                // Set the polygon winding to front facing for the right handed system.
                glFrontFace(GL_CCW);
				//glFrontFace(GL_CW);

                // Enable back face culling.
                glEnable(GL_CULL_FACE);
                glCullFace(GL_BACK);

                glEnable(GL_PROGRAM_POINT_SIZE);
            }

            auto config = g_pApp->GetConfiguration();
            glViewport(0, 0, config.screenWidth, config.screenHeight);
        }

        return result;
    }

    void OpenGLGraphicsManager::Finalize()
    {
        GraphicsManager::Finalize();
    }

    void OpenGLGraphicsManager::Clear()
    {
        GraphicsManager::Clear();

        // Set the color to clear the screen to.
        glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
        // Clear the screen and depth buffer.
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

    void OpenGLGraphicsManager::Draw()
    {
        GraphicsManager::Draw();

        // Render the model using the color shader.
        //RenderBuffers();

        glFlush();
    }

    bool OpenGLGraphicsManager::SetPerFrameShaderParameters(const DrawFrameContext& context)
    {
        bool result;
        unsigned int location;

        // Set the world matrix in the vertex shader.
        result = SetShaderParameter("worldMatrix", context.WorldMatrix);
        if (!result) return result;

        // Set the view matrix in the vertex shader.
        result = SetShaderParameter("viewMatrix", context.ViewMatrix);
        if (!result) return result;

        // Set the projection matrix in the vertex shader.
        result = SetShaderParameter("projectionMatrix", context.ProjectionMatrix);
        if (!result) return result;

        // Set ambient color
        result = SetShaderParameter("ambientColor", context.AmbientColor);
        if (!result) return result;

        // Set number of lights
        result = SetShaderParameter("numLights", (int32_t)context.Lights.size());
        if (!result) return result;
        
        // Set lighting parameters for PS shader
        for (int32_t i = 0; i < context.Lights.size(); ++i)
        {
            char paramName[128] = {0};
            sprintf(paramName, "allLights[%d].lightPosition", i);
            result = SetShaderParameter(paramName, context.Lights[i].LightPosition);
            if (!result) return result;

            sprintf(paramName, "allLights[%d].lightColor", i);
            result = SetShaderParameter(paramName, context.Lights[i].LightColor);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].lightSize", i);
			result = SetShaderParameter(paramName, context.Lights[i].LightSize);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].lightType", i);
			result = SetShaderParameter(paramName, context.Lights[i].Type);
            if (!result) return result;

            sprintf(paramName, "allLights[%d].lightDirection", i);
            result = SetShaderParameter(paramName, context.Lights[i].LightDirection);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].constantAtten", i);
			result = SetShaderParameter(paramName, context.Lights[i].ConstantAtten);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].linearAtten", i);
			result = SetShaderParameter(paramName, context.Lights[i].LinearAtten);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].quadraticAtten", i);
			result = SetShaderParameter(paramName, context.Lights[i].QuadraticAtten);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].cosOfInnerAngle", i);
			result = SetShaderParameter(paramName, context.Lights[i].CosOfInnerAngle);
            if (!result) return result;

			sprintf(paramName, "allLights[%d].cosOfOuterAngle", i);
			result = SetShaderParameter(paramName, context.Lights[i].CosOfOuterAngle);
            if (!result) return result;
        }

        return true;
    }

    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const Matrix4f& param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if (location == -1)
            return false;
        glUniformMatrix4fv(location, 1, false, param);

        return true;
    }

    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const Vector4Df& param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if(location == -1)
        {
            return false;
        }
        glUniform4fv(location, 1, param);

        return true;   
    }

    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const Vector3Df& param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if(location == -1)
        {
            return false;
        }
        glUniform3fv(location, 1, param);

        return true;   
    }

    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const Vector2Df& param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if(location == -1)
        {
            return false;
        }
        glUniform2fv(location, 1, param);

        return true;   
    }

    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const float param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if(location == -1)
        {
            return false;
        }
        glUniform1f(location, param);

        return true;
    }



    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const int param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if(location == -1)
        {
            return false;
        }
        
        glUniform1i(location, param);

        return true;   
    }

    bool OpenGLGraphicsManager::SetShaderParameter(const char* paramName, const bool param)
    {
        unsigned int location;

        location = glGetUniformLocation(m_CurrentShader, paramName);
        if(location == -1)
        {
            return false;
        }
        glUniform1f(location, param);

        return true;   
    }

    void OpenGLGraphicsManager::InitializeBuffers(const Scene& scene)
    {
        // Geometries
        for (auto _it : scene.GeometryNodes)
        {
			auto pGeometryNode = _it.second;
			if (pGeometryNode && pGeometryNode->Visible())
			{
				std::string str = pGeometryNode->GetSceneObjectRef();
				auto pGeometry = scene.GetGeometry(pGeometryNode->GetSceneObjectRef());
				assert(pGeometry);

				// Append parents' transforms
				BaseSceneNode* pParent = pGeometryNode->GetParent();
				std::string parentSid = pParent->GetSid();
				while (parentSid != "Root")
				{
					auto iter = scene.GeometryNodes.find(parentSid);
					if (iter == scene.GeometryNodes.end())
						break;
					auto parentNode = iter->second;

					std::shared_ptr<SceneObjectTransform> _transform;
					_transform = std::make_shared<SceneObjectTransform>(*parentNode->GetCalculatedTransform(), false);
					pGeometryNode->AppendTransform("", std::move(_transform));

					pParent = pParent->GetParent();
					parentSid = pParent->GetSid();
				}
				
				uint32_t meshCount = pGeometry->GetMeshCount();
				for (uint32_t index = 0; index < meshCount; ++index)
				{
					auto pMesh = pGeometry->GetMesh(index).lock();
					if (!pMesh) continue;

					// Set the number of vertex properties.
					auto vertexPropertiesCount = pMesh->GetVertexPropertiesCount();

					// Allocate an OpenGL vertex array object.
					GLuint vao;
					glGenVertexArrays(1, &vao);

					// Bind the vertex array object to store all the buffers and vertex attributes we create here.
					glBindVertexArray(vao);

					GLuint buffer_id;
					for (size_t i = 0; i < vertexPropertiesCount; ++i)
					{
						const SceneObjectVertexArray& v_property_array = pMesh->GetVertexPropertyArray(i);
						auto v_property_array_data_size = v_property_array.GetDataSize();
						auto v_property_array_data = v_property_array.GetData();

						// Generated an ID for the vertex buffer
						glGenBuffers(1, &buffer_id);

						// Bind the vertex buffer and load the vertex (position and color) data into the vertex buffer
						glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
						glBufferData(GL_ARRAY_BUFFER, v_property_array_data_size, v_property_array_data, GL_STATIC_DRAW);

						glEnableVertexAttribArray(i);

						glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
						switch (v_property_array.GetDataType())
						{
						case VertexDataType::kVertexDataTypeFloat1:
							glVertexAttribPointer(i, 1, GL_FLOAT, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeFloat2:
							glVertexAttribPointer(i, 2, GL_FLOAT, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeFloat3:
							glVertexAttribPointer(i, 3, GL_FLOAT, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeFloat4:
							glVertexAttribPointer(i, 4, GL_FLOAT, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeDouble1:
							glVertexAttribPointer(i, 1, GL_DOUBLE, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeDouble2:
							glVertexAttribPointer(i, 2, GL_DOUBLE, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeDouble3:
							glVertexAttribPointer(i, 3, GL_DOUBLE, false, 0, 0);
							break;
						case VertexDataType::kVertexDataTypeDouble4:
							glVertexAttribPointer(i, 4, GL_DOUBLE, false, 0, 0);
							break;
						default:
							assert(0);
							break;
						}

						m_Buffers.push_back(buffer_id);
					}


					GLenum mode;
					switch (pMesh->GetPrimitiveType())
					{
					case PrimitiveType::kPrimitiveTypePointList:
						mode = GL_POINTS;
						break;
					case PrimitiveType::kPrimitiveTypeLineList:
						mode = GL_LINES;
						break;
					case PrimitiveType::kPrimitiveTypeLineStrip:
						mode = GL_LINE_STRIP;
						break;
					case PrimitiveType::kPrimitiveTypeTriList:
						mode = GL_TRIANGLES;
						break;
					case PrimitiveType::kPrimitiveTypeTriStrip:
						mode = GL_TRIANGLE_STRIP;
						break;
					case PrimitiveType::kPrimitiveTypeTriFan:
						mode = GL_TRIANGLE_FAN;
						break;
					default:
						// ignore
						continue;
					}

					size_t indexGroupCount = pMesh->GetIndexGroupCount();
					for (size_t i = 0; i < indexGroupCount; ++i)
					{
						// Generate an ID for the index buffer
						glGenBuffers(1, &buffer_id);

						const SceneObjectIndexArray& indexArray = pMesh->GetIndexArray(i);
						size_t indexArraySize = indexArray.GetDataSize();
						const void* indexArrayData = indexArray.GetData();

						// Bind the index buffer and load the index data into it.
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_id);
						glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexArraySize, indexArrayData, GL_STATIC_DRAW);

						// Set the number of indices in the index array
						GLsizei indexCount = static_cast<GLsizei>(indexArray.GetIndexCount());

						GLenum type;
						switch (indexArray.GetIndexType())
						{
						case IndexDataType::kIndexDataTypeInt8:
							type = GL_UNSIGNED_BYTE;
							break;
						case IndexDataType::kIndexDataTypeInt16:
							type = GL_UNSIGNED_SHORT;
							break;
						case IndexDataType::kIndexDataTypeInt32:
							type = GL_UNSIGNED_INT;
							break;
						default:
							// not supported by OpenGL
							std::cerr << "Error: unsupported index type " << indexArray << std::endl;
							std::cerr << "Mesh: " << *pMesh << std::endl;
							std::cerr << "Geometry: " << *pGeometry << std::endl;
							continue;
						}

						m_Buffers.push_back(buffer_id);

						size_t materialIndex = indexArray.GetMaterialIndex();
						std::string materialKey = pGeometryNode->GetMaterialRef(materialIndex);
						auto material = scene.GetMaterial(materialKey);
						if (material)
						{
							auto color = material->GetBaseColor();
							if (color.ValueMap)
							{
								Image texture = color.ValueMap->GetTextureImage();
								auto it = m_TextureIndex.find(materialKey);
								if (it == m_TextureIndex.end())
								{
									GLuint textureID;
									glGenTextures(1, &textureID);
									glActiveTexture(GL_TEXTURE0 + textureID);
									glBindTexture(GL_TEXTURE_2D, textureID);
									if (texture.BitCount == 24)
									{
										glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture.Width, texture.Height,
											0, GL_RGB, GL_UNSIGNED_BYTE, texture.Data);
									}
									else
									{
										glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture.Width, texture.Height,
											0, GL_RGBA, GL_UNSIGNED_BYTE, texture.Data);
									}
									glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
									glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
									glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
									glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

									m_TextureIndex[color.ValueMap->GetName()] = textureID;
									m_Textures.push_back(textureID);
								}
							}
						}

						auto dbc = std::make_shared<OpenGLDrawBatchContext>();
						dbc->vao = vao;
						dbc->mode = mode;
						dbc->type = type;
						dbc->count = indexCount;
						dbc->node = pGeometryNode;
						dbc->material = material;
						
						m_Frames[m_FrameIndex].batchContexts.push_back(dbc);
					}
				}
            }
        }

        return;
    }

    void OpenGLGraphicsManager::ClearBuffers()
    {
#ifdef DEBUG
        //ClearDebugBuffers();
#endif

        for (int32_t i = 0; i < k_FrameCount; ++i)
        {
            auto& batchContexts = m_Frames[i].batchContexts;

            for (auto dbc : batchContexts)
                glDeleteVertexArrays(1, &std::dynamic_pointer_cast<OpenGLDrawBatchContext>(dbc)->vao);
            
            batchContexts.clear();
        }

        for (auto buf : m_Buffers) {
            glDeleteBuffers(1, &buf);
        }

        for (auto texture : m_Textures) {
            glDeleteTextures(1, &texture);
        }

        m_Buffers.clear();
        m_Textures.clear();
    }

    void OpenGLGraphicsManager::UseShaderProgram(const intptr_t shaderProgram)
    {
        m_CurrentShader = static_cast<GLuint>(shaderProgram);

        // Set the color shader as the current shader program and set the matrices that it will use for rendering.
        glUseProgram(m_CurrentShader);
    }

    void OpenGLGraphicsManager::SetPerFrameConstants(const DrawFrameContext& context)
    {
        bool result = SetPerFrameShaderParameters(context);
        assert(result);
    }

    void OpenGLGraphicsManager::DrawBatch(const DrawBatchContext& context)
    {
        const OpenGLDrawBatchContext& dbc = dynamic_cast<const OpenGLDrawBatchContext&>(context);

		Matrix4f trans;// = *dbc.node ->GetCalculatedTransform();

		trans = *dbc.node->GetCalculatedTransform();
        bool result = SetShaderParameter("modelMatrix", trans);
        assert(result);

        glBindVertexArray(dbc.vao);

        result = SetShaderParameter("usingDiffuseMap", false);
        assert(result);

        if (dbc.material)
        {
            Color color = dbc.material->GetBaseColor();
            if (color.ValueMap)
            {
                result = SetShaderParameter("diffuseMap", m_TextureIndex[color.ValueMap->GetName()]);
                assert(result);
                // set this to tell shader to use texture
                result = SetShaderParameter("usingDiffuseMap", true);
                assert(result);
            }
            else
            {
                result = SetShaderParameter("diffuseColor", Vector3Df({color.Value[0], color.Value[1], color.Value[2]}));
                assert(result);
            }
            
            color = dbc.material->GetSpecularColor();
            result = SetShaderParameter("specularColor", Vector3Df({color.Value[0], color.Value[1], color.Value[2]}));
            assert(result);

            Parameter param = dbc.material->GetSpecularPower();
            result = SetShaderParameter("specularPower", param.Value);
            assert(result);
        }

        glDrawElements(dbc.mode, dbc.count, dbc.type, 0x00);
    }

    void OpenGLGraphicsManager::DrawBatchDepthOnly(const DrawBatchContext& context)
    {
        const OpenGLDrawBatchContext& dbc = dynamic_cast<const OpenGLDrawBatchContext&>(context);

        bool result = SetShaderParameter("modelMatrix", dbc.trans);
        assert(result);

        glBindVertexArray(dbc.vao);

        glDrawElements(dbc.mode, dbc.count, dbc.type, 0x00);
    }

#ifdef DEBUG
    void OpenGLGraphicsManager::ClearDebugBuffers()
    {
        for (auto dbc : m_DebugDrawBatchContext) {
            glDeleteVertexArrays(1, &dbc.vao);
        }

        m_DebugDrawBatchContext.clear();

        for (auto buf : m_DebugBuffers) {
            glDeleteBuffers(1, &buf);
        }

        m_DebugBuffers.clear();
    }

    void OpenGLGraphicsManager::DrawLine(const PointList& vertices, const Matrix4f& trans, const Vector3Df& color)
    {
        auto count = vertices.size();
        GLfloat* _vertices = new GLfloat[3 * count];

        for (auto i = 0; i < count; i++)
        {
            _vertices[3 * i] = vertices[i]->data[0];
            _vertices[3 * i + 1] = vertices[i]->data[1];
            _vertices[3 * i + 2] = vertices[i]->data[2];
        }

        GLuint vao;
        glGenVertexArrays(1, &vao);

        // Bind the vertex array object to store all the buffers and vertex attributes we create here.
        glBindVertexArray(vao);

        GLuint buffer_id;

        // Generate an ID for the vertex buffer.
        glGenBuffers(1, &buffer_id);

        // Bind the vertex buffer and load the vertex (position and color) data into the vertex buffer.
        glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * count, _vertices, GL_STATIC_DRAW);

        delete[] _vertices;

        glEnableVertexAttribArray(0);

        glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);

        m_DebugBuffers.push_back(buffer_id);

        DebugDrawBatchContext& dbc = *(new DebugDrawBatchContext);
        dbc.vao     = vao;
        dbc.mode    = GL_LINES;
        dbc.count   = static_cast<GLsizei>(count);
        dbc.color   = color;
        dbc.trans   = trans;

        m_DebugDrawBatchContext.push_back(std::move(dbc));
    }

    void OpenGLGraphicsManager::DrawLine(const PointList& vertices, const Vector3Df& color)
    {
        Matrix4f trans;
        trans.SetIdentity();

        DrawLine(vertices, trans, color);
    }

    void OpenGLGraphicsManager::DrawLine(const Point& from, const Point& to, const Vector3Df& color)
    {
        PointList point_list;
        point_list.push_back(std::make_shared<Point>(from));
        point_list.push_back(std::make_shared<Point>(to));

        DrawLine(point_list, color);
    }
#endif    
}
