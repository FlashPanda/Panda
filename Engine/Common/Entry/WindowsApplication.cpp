#include <tchar.h>
#include "WindowsApplication.hpp"
#include "Core/KeyCode.hpp"

namespace Panda
{
    void WindowsApplication::CreateMainWindow()
    {
        // 获取控制台的句柄
        HINSTANCE hInstance = GetModuleHandle(NULL);

        // 窗口类
        WNDCLASSEX wc;

        // 清空窗口类对象
        ZeroMemory(&wc, sizeof(WNDCLASSEX));

        // 填充窗口类
        wc.cbSize = sizeof(WNDCLASSEX);
        wc.style = CS_HREDRAW | CS_VREDRAW;
        wc.lpfnWndProc = WindowProc;
        wc.hInstance = hInstance;
        wc.hCursor = LoadCursor(NULL, IDC_ARROW);
        wc.hbrBackground = (HBRUSH)COLOR_WINDOW;
        wc.lpszClassName = _T("Panda");

        // 注册窗口类
        RegisterClassEx(&wc);

        // 创建窗口
        int height_adjust = (GetSystemMetrics(SM_CYFRAME) + GetSystemMetrics(SM_CYCAPTION) + GetSystemMetrics(SM_CXPADDEDBORDER));
        int width_adjust = (GetSystemMetrics(SM_CXFRAME) + GetSystemMetrics(SM_CXPADDEDBORDER));
        m_hWnd = CreateWindowEx(0,
                            "Panda",      					// 窗口类名
                            m_Config.appName,             	// 窗口名
                            WS_OVERLAPPEDWINDOW,              // 窗口风格
                            CW_USEDEFAULT,                    // 窗口的x坐标
                            CW_USEDEFAULT,                    // 窗口的y坐标
                            m_Config.screenWidth + width_adjust,             // 窗口的宽度
                            m_Config.screenHeight + height_adjust,            // 窗口的高度
                            NULL,                             // 父窗口
                            NULL,                             // 菜单
                            hInstance,                        // 应用句柄
                            NULL);                            // 不使用多窗口

        // 显示窗口
        ShowWindow(m_hWnd, SW_SHOW);
    }

    int WindowsApplication::Initialize()
    {
        int result;

        CreateMainWindow();

        m_hDC = GetDC(m_hWnd);

        // first call base class initialization
        result = BaseApplication::Initialize();

        if (result != 0)
            exit(result);

        return result;
    }

    void WindowsApplication::Finalize()
    {
        ReleaseDC(m_hWnd, m_hDC);

        BaseApplication::Finalize();
    }

    void Panda::WindowsApplication::Tick()
    {
        BaseApplication::Tick();

        // Windows消息结构
        MSG msg;

        // 使用PeekMessage代替GetMessage的目的是不阻塞线程
        // 除非是引擎正在启动的时候
        if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);

            DispatchMessage(&msg); 
        }
		//else
		//{
		//	g_pApp->OnDraw();
		//}
    }

    // 消息处理
    LRESULT CALLBACK WindowsApplication::WindowProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
    {
		WindowsApplication* pThis;
        if (message == WM_NCCREATE)
        {
            pThis = static_cast<WindowsApplication*>(reinterpret_cast<CREATESTRUCT*>(lParam)->lpCreateParams);

            SetLastError(0);
            if (!SetWindowLongPtr(hWnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(pThis)))
            {
                if (GetLastError() != 0)
                    return FALSE;
            }
        }
        else
        {
            pThis = reinterpret_cast<WindowsApplication*>(GetWindowLongPtr(hWnd, GWLP_USERDATA));
        }

        // sort through and find what code to run for the message given
        switch(message)
        {
            //case WM_PAINT:
            //    {
            //        g_pApp->OnDraw();
            //    }
            //    break;
            case WM_KEYDOWN:
                {
                    switch(wParam)
                    {
                        case VK_LEFT:
                            g_pInputManager->LeftArrowKeyDown();
                            break;
                        case VK_RIGHT:
                            g_pInputManager->RightArrowKeyDown();
                            break;
                        case VK_UP:
                            g_pInputManager->UpArrowKeyDown();
                            break;
                        case VK_DOWN:
                            g_pInputManager->DownArrowKeyDown();
                            break;
                        case 0x52: // R key
                            g_pInputManager->ResetKeyDown();
                            break;
                        case 0x41:  // A key
                            g_pInputManager->CharKeyDown(KEY_CODE_A);
                            break;
                        case 0x44: // D key
                            g_pInputManager->CharKeyDown(KEY_CODE_D);
                            break;
                        case 0x57: // W key
                            g_pInputManager->CharKeyDown(KEY_CODE_W);
                            break;
                        case 0x53: // S key
                            g_pInputManager->CharKeyDown(KEY_CODE_S);
                            break;
						case 0x45: // E key
							g_pInputManager->CharKeyDown(KEY_CODE_E);
							break;
						case 0x51: // Q key
							g_pInputManager->CharKeyDown(KEY_CODE_Q);
							break;
						case 0x31:	// 1 key
							g_pInputManager->CharKeyDown(KEY_CODE_1);
                        default:
                            break;
                    }
                } 
                break;

            case WM_KEYUP:
            {
                switch(wParam)
                {
                    case VK_LEFT:
                        g_pInputManager->LeftArrowKeyUp();
                        break;
                    case VK_RIGHT:
                        g_pInputManager->RightArrowKeyUp();
                        break;
                    case VK_UP:
                        g_pInputManager->UpArrowKeyUp();
                        break;
                    case VK_DOWN:
                        g_pInputManager->DownArrowKeyUp();
                        break;
                    case 0x52: // R key
                        g_pInputManager->ResetKeyUp();
                        break;
                    case 0x41:  // A key
                        g_pInputManager->CharKeyUp(KEY_CODE_A);
                        break;
                    case 0x44: // D key
                        g_pInputManager->CharKeyUp(KEY_CODE_D);
                        break;
                    case 0x57: // W key
                        g_pInputManager->CharKeyUp(KEY_CODE_W);
                        break;
                    case 0x53: // S key
                        g_pInputManager->CharKeyUp(KEY_CODE_S);
                        break;
					case 0x45: // E key
						g_pInputManager->CharKeyUp(KEY_CODE_E);
						break;
					case 0x51: // Q key
						g_pInputManager->CharKeyUp(KEY_CODE_Q);
						break;
                    default:
                        break;
                }
                break;
            }

                // this message is read when the window is closed
            case WM_DESTROY:
                {
                    // close the application entirely
                    PostQuitMessage(0);
                    m_Quit = true;
                }
        }

        // Handle any messages the switch statement didn't
        return DefWindowProc (hWnd, message, wParam, lParam);
    }
}
