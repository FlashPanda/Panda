#pragma once
#include "ShaderModule.hpp"
#include "glad/glad.h"

namespace Panda
{
    class OpenGLShaderModule : public ShaderModule
    {
        public:
            OpenGLShaderModule() = default;
            ~OpenGLShaderModule() = default;

            virtual int Initialize() final;
            virtual void Finalize() final;

            virtual void Tick() final;
            
            virtual bool InitializeShaders() final;
            virtual void ClearShaders() final;
    };
}