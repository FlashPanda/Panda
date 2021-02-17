// #include <objbase.h>
// #include <d3dcompiler.h>
// #include "D3D11GraphicsManager.hpp"
// #include "WindowsApplication.hpp"

// namespace Panda
// {
// 	extern IApplication* g_pApp;
	
// 	template<class T>
// 	inline void SafeRelease(T **ppInterfaceToRelease)
// 	{
// 		if (*ppInterfaceToRelease != nullptr)
// 		{
// 			(*ppInterfaceToRelease)->Release();

// 			(*ppInterfaceToRelease) = nullptr;
// 		}
// 	}
// }

// void Panda::D3D11GraphicsManager::CreateRenderTarget() {
//     HRESULT hr;
//     ID3D11Texture2D *pBackBuffer;

//     // Get a pointer to the back buffer
//     m_pSwapchain->GetBuffer( 0, __uuidof( ID3D11Texture2D ),
//                                  ( LPVOID* )&pBackBuffer );

//     // Create a render-target view
//     m_pDev->CreateRenderTargetView( pBackBuffer, NULL,
//                                           &g_pRTView );
//     pBackBuffer->Release();

//     // Bind the view
//     m_pDevcon->OMSetRenderTargets( 1, &g_pRTView, NULL );
// }

// // this function prepare graphic resources for use
// HRESULT CreateGraphicsResources(HWND hWnd)
// {
//     HRESULT hr = S_OK;
//     if (g_pSwapchain == nullptr)
//     {
//         // create a struct to hold information about the swap chain
//         DXGI_SWAP_CHAIN_DESC scd;

//         // clear out the struct for use
//         ZeroMemory(&scd, sizeof(DXGI_SWAP_CHAIN_DESC));

//         // fill the swap chain description struct
//         scd.BufferDesc.Width = SCREEN_WIDTH;
//         scd.BufferDesc.Height = SCREEN_HEIGHT;
//         scd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;     // use 32-bit color
//         scd.BufferDesc.RefreshRate.Numerator = 60;
//         scd.BufferDesc.RefreshRate.Denominator = 1;
// 		scd.SampleDesc.Count = 4;
		
// 		scd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
// 		scd.BufferCount = 1;
// 		scd.OutputWindow = hWnd;
// 		scd.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
// 		scd.Windowed = TRUE;
// 		scd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;

//         const D3D_FEATURE_LEVEL featurelevels[] = { //D3D_FEATURE_LEVEL_11_1,
//                                                     D3D_FEATURE_LEVEL_11_0,
//                                                     D3D_FEATURE_LEVEL_10_1,
//                                                     D3D_FEATURE_LEVEL_10_0,
//                                                     D3D_FEATURE_LEVEL_9_3,
//                                                     D3D_FEATURE_LEVEL_9_2,
//                                                     D3D_FEATURE_LEVEL_9_1
// 													};
//         D3D_FEATURE_LEVEL featurelevelsupported;

//         // create a device, device context and swap chain using the information in the scd struct
//         hr = D3D11CreateDeviceAndSwapChain(NULL,
//                                       D3D_DRIVER_TYPE_HARDWARE,
//                                       NULL,
//                                       0,
//                                       featurelevels,
//                                       _countof(featurelevels),
//                                       D3D11_SDK_VERSION,
//                                       &scd,
//                                       &g_pSwapchain,
//                                       &g_pDev,
//                                       &featurelevelsupported,
//                                       &g_pDevcon);

//         if (hr == E_INVALIDARG) {
//             hr = D3D11CreateDeviceAndSwapChain(NULL,
//                                      D3D_DRIVER_TYPE_HARDWARE,
//                                       NULL,
//                                       0,
//                                       &featurelevelsupported,
//                                       1,
//                                       D3D11_SDK_VERSION,
//                                       &scd,
//                                      &g_pSwapchain,
//                                       &g_pDev,
//                                       NULL,
//                                       &g_pDevcon);
//         }
		
//         if (hr == S_OK) {
//             CreateRenderTarget();
//         }
//     }
//     return hr;
// }

// int Panda::D3D11GraphicsManager::Initialize()
// {
// 	int result = 0;
	
// 	result = static_cast<int>(CreateGraphicsResources());
	
// 	return result;
// }

// void Panda::D3D11GraphicsManager::Tick()
// {
// }

// void Panda::D3D11GraphicsManager::Finalize()
// {
// 	SafeRelease(&m_pLayout);
//     SafeRelease(&m_pVS);
//     SafeRelease(&m_pPS);
//     SafeRelease(&m_pVBuffer);
//     SafeRelease(&m_pSwapchain);
//     SafeRelease(&m_pRTView);
//     SafeRelease(&m_pDev);
//     SafeRelease(&m_pDevcon);
// }
