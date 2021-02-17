#pragma once
#include "OpenGL/OpenGLGraphicsManager.hpp"
#include "OpenGL/OpenGLShaderModule.hpp"

namespace Panda
{
    GraphicsManager* g_pGraphicsManager = static_cast<GraphicsManager*>(new OpenGLGraphicsManager);
    IShaderModule* g_pShaderModule = static_cast<IShaderModule*>(new OpenGLShaderModule);
}