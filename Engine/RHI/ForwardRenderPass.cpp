#include "ForwardRenderPass.hpp"
#include "GraphicsManager.hpp"
#include "Core/Interface/IShaderModule.hpp"

namespace Panda
{
    void ForwardRenderPass::Draw(Frame& frame)
    {
        auto shaderProgram = g_pShaderModule->GetDefaultShaderProgram(DefaultShaderIndex::Forward);

        // Set the color shader as the current shader program and set the matrices that it will use for rendering.
        g_pGraphicsManager->UseShaderProgram(shaderProgram);

        g_pGraphicsManager->SetPerFrameConstants(frame.frameContext);

        for (auto dbc : frame.batchContexts)
        {
            g_pGraphicsManager->DrawBatch(*dbc);
        }
    }

}