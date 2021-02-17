#include "BaseApplication.hpp"
#include "Core/Parallel.hpp"

namespace Panda
{
	bool BaseApplication::m_Quit = false;

	BaseApplication::BaseApplication(GfxConfiguration& cfg)
		: m_Config(cfg)
	{
		
	}

	int BaseApplication::Initialize() {
		int ret = 0;
		std::cout << m_Config;

		if((ret = g_pMemoryManager->Initialize()) != 0)
		{
			std::cerr << "Fialed. err = " << ret;
			return ret;
		}

		ParallelInit();

		if ((ret = g_pAssetLoader->Initialize() != 0) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}

		if ((ret = g_pSceneManager->Initialize()) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}

		if ((ret = g_pGraphicsManager->Initialize()) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}

		if ((ret = g_pShaderModule->Initialize()) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}

		if ((ret = g_pInputManager->Initialize()) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}

		if ((ret = g_pGameLogic->Initialize()) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}

		#ifdef DEBUG
		if ((ret = g_pDebugManager->Initialize()) != 0)
		{
			std::cerr << "Failed. err = " << ret;
			return ret;
		}
		#endif

		return ret;
	}


	void BaseApplication::Finalize() {
		#ifdef DEBUG
		g_pDebugManager->Finalize();
		#endif
		g_pGameLogic->Finalize();
		g_pInputManager->Finalize();
		g_pShaderModule->Finalize();
		g_pGraphicsManager->Finalize();
		g_pSceneManager->Finalize();
		g_pAssetLoader->Finalize();
		ParallelCleanup();
		g_pMemoryManager->Finalize();
	}

	void BaseApplication::Tick() {
		g_pMemoryManager->Tick();
		g_pAssetLoader->Tick();
		g_pSceneManager->Tick();
		g_pInputManager->Tick();
		g_pShaderModule->Tick();
		g_pGameLogic->Tick();
		g_pGraphicsManager->Tick();
		#ifdef DEBUG
		g_pDebugManager->Tick();
		#endif
	}

	void BaseApplication::SetCommandLineParameters(int argc, char** argv)
	{
		m_ArgC = argc;
		m_ppArgV = argv;
	}

	bool BaseApplication::IsQuit() {
		return m_Quit;
	}

}

