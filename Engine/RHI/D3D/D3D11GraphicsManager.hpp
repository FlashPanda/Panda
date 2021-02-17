// #pragma once
// #include <stdint.h>
// #include <d3d11.h>
// #include <d3d11_1.h>
// #include "GraphicsManager.hpp"

// namespace Panda {
//     class D3D11GraphicsManager : public GraphicsManager
//     {
//     public:
//        	virtual int Initialize();
// 	    virtual void Finalize();

// 	    virtual void Tick();

//     private:
//         HRESULT CreateRenderTarget();
//         HRESULT CreateGraphicsResources();

//     private:
//         IDXGISwapChain          *m_pSwapchain = nullptr;              // the pointer to the swap chain interface
// 		ID3D11Device            *m_pDev       = nullptr;              // the pointer to our Direct3D device interface
// 		ID3D11DeviceContext     *m_pDevcon    = nullptr;              // the pointer to our Direct3D device context

// 		ID3D11RenderTargetView  *m_pRTView    = nullptr;
		
// 		ID3D11Buffer            *m_pVBuffer   = nullptr;              // Vertex Buffer
//     };
// }