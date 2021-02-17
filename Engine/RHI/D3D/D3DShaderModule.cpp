#include "D3DShaderModule.hpp"
#include "Core/AssetLoader.hpp"

namespace Panda
{
    const char VS_SHADER_SOURCE_FILE[] = "Shaders/basic_vs.cso";
    const char PS_SHADER_SOURCE_FILE[] = "Shaders/basic_ps.cso";
    #ifdef DEBUG
    const char DEBUG_VS_SHADER_SOURCE_FILE[] = "Shaders/debug_vs.cso";
    const char DEBUG_PS_SHADER_SOURCE_FILE[] = "Shaders/debug_ps.cso";
    #endif

    int D3DShaderModule::Initialize()
    {
        return InitializeShaders() == false;
    }

    void D3DShaderModule::Finalize()
    {
        ClearShaders();
    }

    void D3DShaderModule::Tick()
    {

    }

    bool D3DShaderModule::InitializeShaders()
    {
        HRESULT hr = S_OK;
        const char* vsFilename = VS_SHADER_SOURCE_FILE;
        const char* fsFilename = PS_SHADER_SOURCE_FILE;
        #ifdef DEBUG
        const char* debugVsFilename = DEBUG_VS_SHADER_SOURCE_FILE;
        const char* debugFsFilename = DEBUG_PS_SHADER_SOURCE_FILE;
        #endif

        // load the shaders
        // forward rendering shader
        Buffer vertexShader = g_pAssetLoader->SyncOpenAndReadBinary(vsFilename);
        Buffer pixelShader = g_pAssetLoader->SyncOpenAndReadBinary(fsFilename);

        D3DShaderProgram& shaderProgram = *(new D3DShaderProgram);
        shaderProgram.vertexShaderByteCode.pShaderBytecode = vertexShader.GetData();
        shaderProgram.vertexShaderByteCode.BytecodeLength = vertexShader.GetDataSize();

        shaderProgram.pixelShaderByteCode.pShaderBytecode = pixelShader.GetData();
        shaderProgram.pixelShaderByteCode.BytecodeLength = pixelShader.GetDataSize();

        m_DefaultShaders[DefaultShaderIndex::Forward] = reinterpret_cast<intptr_t>(&shaderProgram);

        #ifdef DEBUG
        // debug shader
        D3DShaderProgram& d_shaderProgram = *(new D3DShaderProgram);
        vertexShader = g_pAssetLoader->SyncOpenAndReadBinary(debugVsFilename);
        pixelShader = g_pAssetLoader->SyncOpenAndReadBinary(debugFsFilename);

		d_shaderProgram.vertexShaderByteCode.pShaderBytecode = vertexShader.GetData();
		d_shaderProgram.vertexShaderByteCode.BytecodeLength = vertexShader.GetDataSize();

		d_shaderProgram.pixelShaderByteCode.pShaderBytecode = pixelShader.GetData();
		d_shaderProgram.pixelShaderByteCode.BytecodeLength = pixelShader.GetDataSize();

        m_DefaultShaders[DefaultShaderIndex::Debug] = reinterpret_cast<intptr_t>(&d_shaderProgram);
        #endif

        return hr == S_OK;
    }

    void D3DShaderModule::ClearShaders()
    {}
}