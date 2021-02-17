#pragma once
#include "Interface/IGameLogic.hpp"
#include "Math/PandaMath.hpp"
#include "portable.hpp"

namespace Panda
{
    class EditorLogic : implements IGameLogic
    {
        // overrides
        int Initialize() final;
        void Finalize() final;
        void Tick() final;

        void OnLeftKeyDown() final;
        void OnRightKeyDown() final;
        void OnUpKeyDown() final;
        void OnDownKeyDown() final;
		void OnCharKeyDown(uint32_t keyCode) final;
		void OnCharKeyUp(uint32_t keyCode) final;

        #ifdef DEBUG
        void DrawDebugInfo() final;
        #endif
    };
}