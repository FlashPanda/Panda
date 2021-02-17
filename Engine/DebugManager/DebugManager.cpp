#include <iostream>
#include "DebugManager.hpp"
#include "RHI/GraphicsManager.hpp"
#include "Core/Interface/IGameLogic.hpp"

namespace Panda
{
    int DebugManager::Initialize()
    {
        return 0;
    }

    void DebugManager::Finalize()
    {

    }

    void DebugManager::Tick()
    {

    }

    void DebugManager::ToggleDebugInfo()
    {
        m_IsDrawDebugInfo = !m_IsDrawDebugInfo;
    }

    void DebugManager::DrawDebugInfo()
    {
        DrawGrid();
        DrawAxis();
    }

    void DebugManager::DrawAxis()
    {
        // x - axis
		Vector3Df from({ -1000.0f, 0.0f, 0.0f });
		Vector3Df to({ 1000.0f, 0.0f, 0.0f });
		Vector3Df color({ 1.0f, 0.0f, 0.0f });
        g_pGraphicsManager->DrawLine(from, to, color);

        // y - axis
        from = {0.0f, -1000.0f, 0.0f};
        to = {0.0f, 1000.0f, 0.0f};
        color = {0.0f, 1.0f, 0.0f};
        g_pGraphicsManager->DrawLine(from, to, color);

        // z - axis
        from = {0.0f, 0.0f, -1000.0f};
        to = {0.0f, 0.0f, 1000.0f};
        color = {0.0f, 0.0f, 1.0f};
        g_pGraphicsManager->DrawLine(from, to, color);
    }

    void DebugManager::DrawGrid()
    {
    }
}