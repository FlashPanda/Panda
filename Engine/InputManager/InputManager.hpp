#pragma once
#include "Core/Interface/IRuntimeModule.hpp"
#include "Core/portable.hpp"
#include "Core/KeyCode.hpp"

namespace Panda
{
    class InputManager : implements IRuntimeModule
    {
        public:
            virtual int Initialize();
            virtual void Finalize();
            virtual void Tick();

            void UpArrowKeyDown();
            void UpArrowKeyUp();
            void DownArrowKeyDown();
            void DownArrowKeyUp();
            void LeftArrowKeyDown();
            void LeftArrowKeyUp();
            void RightArrowKeyDown();
            void RightArrowKeyUp();
			void CharKeyDown(uint32_t keyCode);
			void CharKeyUp(uint32_t keyCode);

            void ResetKeyDown();
            void ResetKeyUp();

            #ifdef DEBUG
            void DebugKeyDown();
            void DebugKeyUp();
            #endif

        protected:
            bool m_IsUpKeyPressed = false;
            bool m_IsDownKeyPressed = false;
            bool m_IsLeftKeyPressed = false;
            bool m_IsRightKeyPressed = false;
    };

    extern InputManager* g_pInputManager;
}