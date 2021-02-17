#pragma once
#include "IRuntimeModule.hpp"
#include "Core/portable.hpp"
#include <memory>

namespace Panda
{
    ENUM(DefaultShaderIndex)
    {
        ShadowMap = "SHMP"_i32,
        Forward = "FRWD"_i32,
        Differed = "DIFR"_i32,
        Debug = "DEBG"_i32,
        Copy = "BTBL"_i32
    };

    Interface IShaderModule : implements IRuntimeModule
    {
        public:
            virtual ~IShaderModule() = default;

            virtual bool InitializeShaders() = 0;
            virtual void ClearShaders() = 0;

            virtual intptr_t GetDefaultShaderProgram(DefaultShaderIndex index) = 0;
    };

    extern IShaderModule* g_pShaderModule;
}