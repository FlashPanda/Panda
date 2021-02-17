#pragma once

#include "Core/Interface/IRuntimeModule.hpp"

namespace Panda
{
    class DebugManager : implements IRuntimeModule
    {
        public:
            int Initialize();
            void Finalize();
            void Tick();

            void ToggleDebugInfo();

            void DrawDebugInfo();

        protected:
            void DrawAxis();
            void DrawGrid();

            bool m_IsDrawDebugInfo = false;
    };

    extern DebugManager* g_pDebugManager;
}