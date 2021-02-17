#pragma once

#include <unordered_map>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "GraphicsManager.hpp"
#include "Math/PandaMath.hpp"
#include "glad/glad.h"
#include "SceneManager.hpp"
#include "Interface/IApplication.hpp"

namespace Panda {
    class OpenGLGraphicsManager : public GraphicsManager
    {
        public:
            int Initialize();
            void Finalize() final;

            void Clear() final;
            void Draw() final;

            void UseShaderProgram(const intptr_t shaderProgram) final;
            void SetPerFrameConstants(const DrawFrameContext& context) final;
            void DrawBatch(const DrawBatchContext& context) final;
            void DrawBatchDepthOnly(const DrawBatchContext& context) final;

            #ifdef DEBUG
            void DrawLine(const Point& from, const Point& to, const Vector3Df& color) final;
            void DrawLine(const PointList& vertices, const Vector3Df& color) final;
            void DrawLine(const PointList& vertices, const Matrix4f& trans, const Vector3Df& color) final;
            void ClearDebugBuffers();
            #endif

            void InitializeBuffers(const Scene& scene) final;
            void ClearBuffers() final;

        protected:
            bool SetShaderParameter(const char* paramName, const Matrix4f& param);
            bool SetShaderParameter(const char* paramName, const Vector4Df& param);
            bool SetShaderParameter(const char* paramName, const Vector3Df& param);
            bool SetShaderParameter(const char* paramName, const Vector2Df& param);
            bool SetShaderParameter(const char* paramName, const float param);
            bool SetShaderParameter(const char* paramName, const int param);
            bool SetShaderParameter(const char* paramName, const bool param);
            bool SetPerFrameShaderParameters(const DrawFrameContext& context);

        private:
            GLuint m_CurrentShader;

            struct OpenGLDrawBatchContext : public DrawBatchContext
            {
                GLuint vao;
                GLenum mode;
                GLenum type;
                GLsizei count;
            };
#ifdef DEBUG
            struct DebugDrawBatchContext
            {
                GLuint vao;
                GLenum mode;
                GLsizei count;
                Vector3Df color;
                Matrix4f trans;
            };
#endif
            std::vector<GLuint> m_Buffers;
            std::vector<GLuint> m_Textures;
            std::map<std::string, GLint> m_TextureIndex;

#ifdef DEBUG
            std::vector<DebugDrawBatchContext> m_DebugDrawBatchContext;
            std::vector<GLuint> m_DebugBuffers;
#endif
    };
}