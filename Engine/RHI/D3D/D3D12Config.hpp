#include "RHI/D3D/D3D12GraphicsManager.hpp"
#include "RHI/D3D/D3DShaderModule.hpp"

namespace Panda
{
    GraphicsManager* g_pGraphicsManager = static_cast<GraphicsManager*>(new D3D12GraphicsManager);
    IShaderModule* g_pShaderModule = static_cast<IShaderModule*>(new D3DShaderModule);
}