#pragma once
#include <windows.h>
#include <windowsx.h>
#include "BaseApplication.hpp"

namespace Panda
{
	class WindowsApplication : public BaseApplication
	{
	public:
		WindowsApplication(GfxConfiguration& config)
			: BaseApplication(config) {}
			
		virtual int Initialize();
		virtual void Finalize();
		virtual void Tick();
		
		// 函数回调
		static LRESULT CALLBACK WindowProc (HWND hWnd,
											UINT message,
											WPARAM wParam,
											LPARAM lParam);
											
		HWND GetMainWindow() {return m_hWnd;}

	protected:
		void CreateMainWindow();
	
	protected:
		HWND m_hWnd;
		HDC m_hDC;
	};
}