#pragma once
#include "IRuntimeModule.hpp"
#include "Core/KeyCode.hpp"

namespace Panda {
    class IGameLogic : implements IRuntimeModule
    {
    public:
        virtual int Initialize() = 0;
        virtual void Finalize() = 0;
        virtual void Tick() = 0;

        virtual void OnUpKeyDown() {};
        virtual void OnUpKeyUp() {};
        virtual void OnUpKey() {};

        virtual void OnDownKeyDown() {};
        virtual void OnDownKeyUp() {};
        virtual void OnDownKey() {};

        virtual void OnLeftKeyDown() {};
        virtual void OnLeftKeyUp() {};
        virtual void OnLeftKey() {};

        virtual void OnRightKeyDown() {};
        virtual void OnRightKeyUp() {};
        virtual void OnRightKey() {};

		virtual void OnCharKeyDown(uint32_t keyCode) {}
		virtual void OnCharKeyUp(uint32_t keyCode) {}

        virtual void OnButton1Down() {};
        virtual void OnButton1Up() {};
    };

    extern IGameLogic* g_pGameLogic;
}