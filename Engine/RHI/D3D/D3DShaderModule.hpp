#include <d3d12.h>
#include "RHI/ShaderModule.hpp"

namespace Panda
{
    struct D3DShaderProgram
    {
        D3D12_SHADER_BYTECODE vertexShaderByteCode;
        D3D12_SHADER_BYTECODE pixelShaderByteCode;
    };

    class D3DShaderModule : public ShaderModule
    {
        public:
            D3DShaderModule() = default;
            ~D3DShaderModule() = default;

            virtual int Initialize() final;
            virtual void Finalize() final;

            virtual void Tick() final;

            virtual bool InitializeShaders() final;
            virtual void ClearShaders() final;
    };
}