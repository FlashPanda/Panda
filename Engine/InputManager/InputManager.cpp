#include <iostream>
#include "InputManager.hpp"
#include "RHI/GraphicsManager.hpp"
#include "SceneManager/SceneManager.hpp"
#include "Core/Interface/IGameLogic.hpp"
#include "DebugManager/DebugManager.hpp"
#include "Core/Math/PandaMath.hpp"

namespace Panda
{
    int InputManager::Initialize()
    {
        return 0;
    }

    void InputManager::Finalize()
    {

    }

    void InputManager::Tick()
    {

    }

    void InputManager::UpArrowKeyDown()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Up Arrow Key Donw!" << std::endl;
        #endif
        g_pGameLogic->OnUpKeyDown();
        if (!m_IsUpKeyPressed)
        {
            g_pGameLogic->OnUpKey();
            m_IsUpKeyPressed = true;
        }
    }

    void InputManager::UpArrowKeyUp()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Up Arrow Key Up!" << std::endl;
        #endif
        g_pGameLogic->OnUpKeyUp();
        m_IsUpKeyPressed = false;
    }

    void InputManager::DownArrowKeyDown()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Down Arrow Key Down!" << std::endl;
        #endif
        g_pGameLogic->OnDownKeyDown();
        if (!m_IsDownKeyPressed)
        {
            g_pGameLogic->OnDownKey();
            m_IsDownKeyPressed = true;
        }
    }

    void InputManager::DownArrowKeyUp()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Down Arrow Key Up!" << std::endl;
        #endif
        g_pGameLogic->OnDownKeyUp();
        m_IsDownKeyPressed = false;
    }

    void InputManager::LeftArrowKeyDown()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Left Arrow Key Down!" << std::endl;
        #endif
        g_pGameLogic->OnLeftKeyDown();
        if (!m_IsLeftKeyPressed)
        {
            g_pGameLogic->OnLeftKey();
            m_IsLeftKeyPressed = true;
        }
    }

    void InputManager::LeftArrowKeyUp()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Left Arrow Key Up!" << std::endl;
        #endif
        g_pGameLogic->OnLeftKeyUp();
        m_IsLeftKeyPressed = false;
    }

    void InputManager::RightArrowKeyDown()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Right Arrow Key Down!" << std::endl;
        #endif
        g_pGameLogic->OnRightKeyDown();
        if (!m_IsRightKeyPressed)
        {
            g_pGameLogic->OnRightKey();
            m_IsRightKeyPressed = true;
        }
    }

    void InputManager::RightArrowKeyUp()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Right Arrow Key Up!" << std::endl;
        #endif
        g_pGameLogic->OnRightKeyUp();
        m_IsRightKeyPressed = false;
    }

    void InputManager::ResetKeyDown()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Reset Key Down!" << std::endl;
        #endif
        g_pSceneManager->ResetScene();
    }

	void InputManager::CharKeyDown(uint32_t keyCode)
	{
		g_pGameLogic->OnCharKeyDown(keyCode);
	}

	void InputManager::CharKeyUp(uint32_t keyCode)
	{
		g_pGameLogic->OnCharKeyUp(keyCode);
	}

    void InputManager::ResetKeyUp()
    {
        #ifdef DEBUG
        std::cerr << "[InputManager] Reset Key Up!" << std::endl;
        #endif
    }

    #ifdef DEBUG
    void InputManager::DebugKeyDown()
    {
        std::cerr << "[InputManager] Debug Key Down!" << std::endl;
        g_pDebugManager->ToggleDebugInfo();
    }

    void InputManager::DebugKeyUp()
    {
        std::cerr << "[InputManager] Debug Key Up!" << std::endl;
    }
    #endif
}