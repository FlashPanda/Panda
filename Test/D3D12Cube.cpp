#include <windows.h>
#include <windowsx.h>
#include <stdio.h>
#include <tchar.h>
#include <stdint.h>
#include <array>
#include <vector>
#include <xstring>
#include <DirectXMath.h>
#include <DirectXPackedVector.h>
#include <DirectXColors.h>
#include <DirectXCollision.h>
#include <d3d12.h>

#include <DXGI1_4.h>
#include <d3dcompiler.h>
#include <wrl/client.h>
#include <exception>

#include "Utility.hpp"
#include "Math/PandaMath.hpp"
#include "d3dx12.h"

using namespace DirectX;
using namespace DirectX::PackedVector;

namespace Panda
{
    Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipZeroToOne;
}

using namespace Panda;

const uint32_t nScreenWidth = 960;
const uint32_t nScreenHeight = 480;
const uint32_t nFrameCount = 2;

HWND g_hWnd; // window handle
IDXGIFactory4*  g_pFactory = nullptr; // factory
ID3D12Device* g_pDevice = nullptr;  // Device
IDXGISwapChain* g_pSwapChain = nullptr; // swap chain

/// Heaps
ID3D12DescriptorHeap*   g_pRtvHeap = nullptr;   // render target heap
ID3D12DescriptorHeap*   g_pDsvHeap = nullptr;   // depth / stencil heap
ID3D12DescriptorHeap*   g_pCbvHeap = nullptr;   // const buffer heap

/// resources
ID3D12Resource*     g_pRenderTargets[nFrameCount] = {nullptr, nullptr}; // render target view
ID3D12Resource*     g_pDefaultVertexBuffer = nullptr; // vertex buffer
ID3D12Resource*     g_pUploadVertexBuffer = nullptr; // vertex buffer
ID3D12Resource*     g_pDefaultIndexBuffer = nullptr; // index buffer
ID3D12Resource*     g_pUploadIndexBuffer = nullptr; // index buffer
ID3D12Resource*     g_pDepthStencilBuffer = nullptr; // depth stencil buffer
ID3D12Resource*     g_pConstantBuffer = nullptr; // constant buffer (upload)

/// views
D3D12_VERTEX_BUFFER_VIEW g_VertexBufferView;
D3D12_INDEX_BUFFER_VIEW g_IndexBufferView;

/// commands
ID3D12CommandQueue*         g_pCommandQueue = nullptr;  // command queue
ID3D12CommandAllocator*     g_pCommandAllocator = nullptr;  // command allocator
ID3D12GraphicsCommandList*  g_pCommandList = nullptr;   // command list

/// root signature
ID3D12RootSignature*    g_pRootSignature = nullptr; // root signature

/// pso
ID3D12PipelineState*    g_pPso;     // pipeline state object

/// sizes
uint32_t    g_nRtvDescriptorSize;
uint32_t    g_nCbvDescriptorSize;
uint32_t    g_nConstantBufferSize;

/// synchronization object
uint32_t        g_nFrameIndex = 0;
ID3D12Fence*    g_pFence = nullptr;
uint32_t        g_nFenceValue = 0;

/// screen
D3D12_VIEWPORT  g_ViewPort = {0.0f, 0.0f, static_cast<float>(nScreenWidth), static_cast<float>(nScreenHeight)};
D3D12_RECT      g_ScissorPort = {0, 0, nScreenWidth, nScreenHeight};

// vertex struct
struct Vertex
{
    Vector3Df Pos;
    Vector4Df Color;
};

struct ObjectConstants
{
	Matrix4f WorldViewProj;

    ObjectConstants()
    {
        WorldViewProj.SetIdentity();
    }
};

std::array<Vertex, 8> vertices = 
{
	Vertex({Vector3Df({-1.0f, -1.0f, -1.0f}), Vector4Df({1.000000000f, 1.000000000f, 1.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({-1.0f, +1.0f, -1.0f}), Vector4Df({0.000000000f, 0.000000000f, 0.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({+1.0f, +1.0f, -1.0f}), Vector4Df({1.000000000f, 0.000000000f, 0.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({+1.0f, -1.0f, -1.0f}), Vector4Df({0.000000000f, 0.501960814f, 0.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({-1.0f, -1.0f, +1.0f}), Vector4Df({0.000000000f, 0.000000000f, 1.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({-1.0f, +1.0f, +1.0f}), Vector4Df({1.000000000f, 1.000000000f, 0.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({+1.0f, +1.0f, +1.0f}), Vector4Df({0.000000000f, 1.000000000f, 1.000000000f, 1.000000000f})}),
	Vertex({Vector3Df({+1.0f, -1.0f, +1.0f}), Vector4Df({1.000000000f, 0.000000000f, 1.000000000f, 1.000000000f})})
};

std::array<uint16_t, 36> indices =
{
    // front face
    0, 1, 2,
    0, 2, 3,

    // back face
    4, 6, 5,
    4, 7, 6,

    // left face
    4, 5, 1,
    4, 1, 0,

    // right face
    3, 2, 6,
    3, 6, 7,

    // top face
    1, 5, 6,
    1, 6, 2,

    // bottom face
    4, 0, 3,
    4, 3, 7
};

// input layout
D3D12_INPUT_ELEMENT_DESC g_ieds[2] =
{
	{"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
	{"COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0}
};

// shaders
D3D12_SHADER_BYTECODE g_VertexShaderByteCode;
D3D12_SHADER_BYTECODE g_PixelShaderByteCode;

// D3D initialization
void CreateDevice();
void CreateCommandObjects();
void CreateSwapChain();
void CreateDescriptorHeaps();
void CreateRenderTarget();
void CreateDepthStencilBuffer();
void CreateConstantBuffer();
void CreateVertexBuffer();
void CreateIndexBuffer();
void CreateRootSignature();
void InitializeShaders();
void CreatePSOs();
void PopulateCommandList();
void Update();
void Draw();
void FlushCommandQueue();
void InitDirect3D();
ID3DBlob* CompileShader(const std::wstring& filename,
	const D3D_SHADER_MACRO* defines,
	const std::string& entrypoint,
	const std::string& target);

// The Main window
// the WindowProc function prototype
LRESULT CALLBACK WindowProc(HWND hWnd,
                         UINT message,
                         WPARAM wParam,
                         LPARAM lParam);

// the entry point for any Windows program
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPTSTR lpCmdLine, int nCmdShow)
{
    // this struct holds information for the window class
    WNDCLASSEX wc;

    // clear out the window class for use
    ZeroMemory(&wc, sizeof(WNDCLASSEX));

    // fill in the struct with the needed information
    wc.cbSize = sizeof(WNDCLASSEX);
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)COLOR_WINDOW;
    wc.lpszClassName = _T("WindowClass1");

    // register the window class
    RegisterClassEx(&wc);

    // create the window and use the result as the handle
	g_hWnd = CreateWindowEx(0,
                          _T("WindowClass1"),                   // name of the window class
                          _T("Hello, Engine![Direct 3D]"),      // title of the window
                          WS_OVERLAPPEDWINDOW,                  // window style
                          100,                                  // x-position of the window
                          100,                                  // y-position of the window
                          nScreenWidth,                         // width of the window
                          nScreenHeight,                        // height of the window
                          NULL,                                 // we have no parent window, NULL
                          NULL,                                 // we aren't using menus, NULL
                          hInstance,                            // application handle
                          NULL);                                // used with multiple windows, NULL

    // display the window on the screen
    ShowWindow(g_hWnd, nCmdShow);
	UpdateWindow(g_hWnd);

    // Do the D3D initialization here
    InitDirect3D();

    // this struct holds Windows event messages
    MSG msg = {0};

    // wait for the next message in the queue, store the result in 'msg'
	while(msg.message != WM_QUIT)
	{
		// If there are Window messages then process them.
		if(PeekMessage( &msg, 0, 0, 0, PM_REMOVE ))
		{
            TranslateMessage( &msg );
            DispatchMessage( &msg );
		}
		// Otherwise, do animation/game stuff.
		else
        {	
			Update();
            PopulateCommandList();
			Draw();
			FlushCommandQueue();
        }
    }

    // return this part of the WM_QUIT message to Windows
    return msg.wParam;
}

// this is the main message handler for the program
LRESULT CALLBACK WindowProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    LRESULT result = 0;
    bool wasHandled = false;

    // sort through and find what code to run for the message given
    switch(message)
    {
    case WM_DESTROY:
        PostQuitMessage(0);
        wasHandled = true;
        break;

    case WM_DISPLAYCHANGE:
        InvalidateRect(hWnd, nullptr, false);
        wasHandled = true;
        break;
    }

    // Handle any messages the switch statement didn't
    return DefWindowProc (hWnd, message, wParam, lParam);
}

void CreateDevice()
{
    if (g_pSwapChain == nullptr)
    {
        // 1、the factory
        if (FAILED(CreateDXGIFactory1(IID_PPV_ARGS(&g_pFactory))))
        {
            return;
        }

        // 2、the device
        HRESULT hardwareResult = D3D12CreateDevice(
            nullptr,    // default adapter
            D3D_FEATURE_LEVEL_11_0,
            IID_PPV_ARGS(&g_pDevice)
        );
        if (FAILED(hardwareResult))
        {
            IDXGIAdapter* pAdapter;
            g_pFactory->EnumWarpAdapter(IID_PPV_ARGS(&pAdapter));

            if (FAILED(D3D12CreateDevice(
                pAdapter,
                D3D_FEATURE_LEVEL_11_0,
                IID_PPV_ARGS(&g_pDevice)
            )))
            {
                return;
            }
        }

		g_pDevice->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&g_pFence));
    }

    return;
}

void CreateCommandObjects()
{
    // Command Queue
    D3D12_COMMAND_QUEUE_DESC queueDesc = {};
    queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
    queueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
    g_pDevice->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&g_pCommandQueue));

	g_pDevice->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&g_pCommandAllocator));

	// Command List
	g_pDevice->CreateCommandList(
		0,
		D3D12_COMMAND_LIST_TYPE_DIRECT,
		g_pCommandAllocator,
		nullptr,
		IID_PPV_ARGS(&g_pCommandList)
	);

	// We must close the command list.
	g_pCommandList->Close();
}

void CreateSwapChain()
{
    DXGI_SWAP_CHAIN_DESC sd;
    sd.BufferDesc.Width = nScreenWidth;
    sd.BufferDesc.Height = nScreenHeight;
    sd.BufferDesc.RefreshRate.Numerator = 60;
    sd.BufferDesc.RefreshRate.Denominator = 1;
    sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
    sd.BufferDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED;
    sd.BufferDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;
    sd.SampleDesc.Count = 1;
    sd.SampleDesc.Quality = 0;
    sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
    sd.BufferCount = nFrameCount;
    sd.OutputWindow = g_hWnd;
    sd.Windowed = true;
    sd.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
    sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;

    // Create swap chain
    // Note: Swap chain uses queue to perform flush
    g_pFactory->CreateSwapChain(
        g_pCommandQueue,
        &sd,
        &g_pSwapChain);
}

void CreateDescriptorHeaps()
{
    // create renter target view descriptor heap
    D3D12_DESCRIPTOR_HEAP_DESC rtvHeapDesc = {};
    rtvHeapDesc.NumDescriptors = nFrameCount;
    rtvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
    rtvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
	rtvHeapDesc.NodeMask = 0;
    g_pDevice->CreateDescriptorHeap(&rtvHeapDesc, IID_PPV_ARGS(&g_pRtvHeap));

    g_nRtvDescriptorSize = g_pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);

    // create depth stencil view descriptor heap
    D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
    dsvHeapDesc.NumDescriptors = 1;
    dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
    dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
	dsvHeapDesc.NodeMask = 0;
    g_pDevice->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&g_pDsvHeap));

    // create const buffer view descriptor heap
    D3D12_DESCRIPTOR_HEAP_DESC cbvHeapDesc = {};
	cbvHeapDesc.NumDescriptors = 1;
    cbvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    cbvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE; // we want shaders to access the buffer
	cbvHeapDesc.NodeMask = 0;
    g_pDevice->CreateDescriptorHeap(&cbvHeapDesc, IID_PPV_ARGS(&g_pCbvHeap));

    g_nCbvDescriptorSize = g_pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
}

void CreateRenderTarget()
{
    D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = g_pRtvHeap->GetCPUDescriptorHandleForHeapStart();

    // Create a RTV for each frame
    for (size_t i = 0; i < nFrameCount; ++i)
    {
        if (FAILED(g_pSwapChain->GetBuffer(i, IID_PPV_ARGS(&g_pRenderTargets[i])))){
            return;
        }
        g_pDevice->CreateRenderTargetView(g_pRenderTargets[i], nullptr, rtvHandle);
        rtvHandle.ptr += g_nRtvDescriptorSize;
    }

    return;
}

void CreateDepthStencilBuffer()
{
	D3D12_RESOURCE_DESC depthStencilDesc;
	depthStencilDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
	depthStencilDesc.Alignment = 0;
	depthStencilDesc.Width = nScreenWidth;
	depthStencilDesc.Height = nScreenHeight;
	depthStencilDesc.DepthOrArraySize = 1;
	depthStencilDesc.MipLevels = 1;
	depthStencilDesc.Format = DXGI_FORMAT_R24G8_TYPELESS;
	depthStencilDesc.SampleDesc.Count = 1;
	depthStencilDesc.SampleDesc.Quality = 0;
	depthStencilDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
	depthStencilDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;

	D3D12_CLEAR_VALUE optClear;
	optClear.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
	optClear.DepthStencil.Depth = 1.0f;
	optClear.DepthStencil.Stencil = 0;
	g_pDevice->CreateCommittedResource(
		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
		D3D12_HEAP_FLAG_NONE,
		&depthStencilDesc,
		D3D12_RESOURCE_STATE_COMMON,
		&optClear,
		IID_PPV_ARGS(&g_pDepthStencilBuffer));

	// Create descriptor to mip level 0 of entire resource using the format of the resource.
	D3D12_DEPTH_STENCIL_VIEW_DESC dsvDesc;
	dsvDesc.Flags = D3D12_DSV_FLAG_NONE;
	dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2D;
	dsvDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
	dsvDesc.Texture2D.MipSlice = 0;
	g_pDevice->CreateDepthStencilView(g_pDepthStencilBuffer, &dsvDesc, g_pDsvHeap->GetCPUDescriptorHandleForHeapStart());

	// Transition the resource from its initial state to be used as a depth buffer.
	g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pDepthStencilBuffer,
		D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_STATE_DEPTH_WRITE));
}

void CreateConstantBuffer()
{
    // The heap on which constant buffer will create
    g_nConstantBufferSize = ((sizeof(ObjectConstants) + 255) & ~255);

    if (FAILED(g_pDevice->CreateCommittedResource(
        &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
        D3D12_HEAP_FLAG_NONE,
        &CD3DX12_RESOURCE_DESC::Buffer(g_nConstantBufferSize),
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(&g_pConstantBuffer)
    )))
    {
        return;
    }

    // create const buffer view
	D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc = {};
	cbvDesc.BufferLocation = g_pConstantBuffer->GetGPUVirtualAddress();
	cbvDesc.SizeInBytes = g_nConstantBufferSize;

	g_pDevice->CreateConstantBufferView(&cbvDesc, g_pCbvHeap->GetCPUDescriptorHandleForHeapStart());

    return;
}

void CreateVertexBuffer()
{
	g_pCommandList->Reset(g_pCommandAllocator, nullptr);

    // create default vertex buffer
    uint32_t verticesByteSize = vertices.size() * sizeof(Vertex);

    if (FAILED(g_pDevice->CreateCommittedResource(
        &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
        D3D12_HEAP_FLAG_NONE,
        &CD3DX12_RESOURCE_DESC::Buffer(verticesByteSize),
		D3D12_RESOURCE_STATE_COMMON,
        nullptr,
        IID_PPV_ARGS(&g_pDefaultVertexBuffer)
    )))
    {
        return;
    }

    if (FAILED(g_pDevice->CreateCommittedResource(
		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
        D3D12_HEAP_FLAG_NONE,
		&CD3DX12_RESOURCE_DESC::Buffer(verticesByteSize),
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(&g_pUploadVertexBuffer)
    )))
    {
        return;
    }

    // copy the data
    D3D12_SUBRESOURCE_DATA vertexData = {};
    vertexData.pData = vertices.data();
	vertexData.RowPitch = verticesByteSize;
	vertexData.SlicePitch = verticesByteSize;

	g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pDefaultVertexBuffer,
		D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_STATE_COPY_DEST));

    UpdateSubresources<1>(g_pCommandList, g_pDefaultVertexBuffer, g_pUploadVertexBuffer, 0, 0, 1, &vertexData);

    g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pDefaultVertexBuffer,
		D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_GENERIC_READ));

    // Initialize the vertex buffer view
    g_VertexBufferView.BufferLocation = g_pDefaultVertexBuffer->GetGPUVirtualAddress();
    g_VertexBufferView.StrideInBytes = sizeof(Vertex);
    g_VertexBufferView.SizeInBytes = vertices.size() * sizeof(Vertex);

	g_pCommandList->Close();

	// execute the command list
	ID3D12CommandList* ppCommandLists[] = { g_pCommandList };
	g_pCommandQueue->ExecuteCommandLists(_countof(ppCommandLists), ppCommandLists);

	FlushCommandQueue();
}

void CreateIndexBuffer()
{
	g_pCommandList->Reset(g_pCommandAllocator, nullptr);

    /// create default index buffer
	uint32_t indicesByteSize = indices.size() * sizeof(uint16_t);

	if (FAILED(g_pDevice->CreateCommittedResource(
		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
		D3D12_HEAP_FLAG_NONE,
		&CD3DX12_RESOURCE_DESC::Buffer(indicesByteSize),
		D3D12_RESOURCE_STATE_COMMON,
		nullptr,
		IID_PPV_ARGS(&g_pDefaultIndexBuffer)
	)))
	{
		return;
	}

	if (FAILED(g_pDevice->CreateCommittedResource(
		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
		D3D12_HEAP_FLAG_NONE,
		&CD3DX12_RESOURCE_DESC::Buffer(indicesByteSize),
		D3D12_RESOURCE_STATE_GENERIC_READ,
		nullptr,
		IID_PPV_ARGS(&g_pUploadIndexBuffer)
	)))
	{
		return;
	}

    // create upload index buffer
    D3D12_SUBRESOURCE_DATA  indexData = {};
    indexData.pData = indices.data();
	indexData.RowPitch = indicesByteSize;
	indexData.SlicePitch = indicesByteSize;

	g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pDefaultIndexBuffer,
		D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_STATE_COPY_DEST));

	UpdateSubresources<1>(g_pCommandList, g_pDefaultIndexBuffer, g_pUploadIndexBuffer, 0, 0, 1, &indexData);

	g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pDefaultIndexBuffer,
		D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_GENERIC_READ));

    // initialize the index buffer view
    g_IndexBufferView.BufferLocation = g_pDefaultIndexBuffer->GetGPUVirtualAddress();
    g_IndexBufferView.Format = DXGI_FORMAT_R16_UINT;
    g_IndexBufferView.SizeInBytes = indicesByteSize;

	g_pCommandList->Close();

	// execute the command list
	ID3D12CommandList* ppCommandLists[] = { g_pCommandList };
	g_pCommandQueue->ExecuteCommandLists(_countof(ppCommandLists), ppCommandLists);

	FlushCommandQueue();

    return;
}

void CreateRootSignature()
{
	// Shader programs typically require resources as input (constant buffers,
	// textures, samplers).  The root signature defines the resources the shader
	// programs expect.  If we think of the shader programs as a function, and
	// the input resources as function parameters, then the root signature can be
	// thought of as defining the function signature.

	CD3DX12_ROOT_PARAMETER slotRootParameter[1];

	// Create a single descriptor table of CBVs.
	CD3DX12_DESCRIPTOR_RANGE cbvTable;
	cbvTable.Init(D3D12_DESCRIPTOR_RANGE_TYPE_CBV, 1, 0);
	slotRootParameter[0].InitAsDescriptorTable(1, &cbvTable);

	// A root signature is an array of root parameters.
	CD3DX12_ROOT_SIGNATURE_DESC rootSigDesc(1, slotRootParameter, 0, nullptr,
		D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);

    ID3DBlob* signature = nullptr;
    ID3DBlob* error = nullptr;
    if (FAILED(D3D12SerializeRootSignature(&rootSigDesc, D3D_ROOT_SIGNATURE_VERSION_1, &signature, &error)))
    {
        return;
    }

    if (FAILED(g_pDevice->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&g_pRootSignature))))
    {
        return;
    }

    //SafeRelease(&signature);
    //SafeRelease(&error);

    return;
}

void InitializeShaders()
{
	ID3DBlob* pVS = CompileShader(L"color.hlsl", nullptr, "VS", "vs_5_0");
	ID3DBlob* pPS = CompileShader(L"color.hlsl", nullptr, "PS", "ps_5_0");

	g_VertexShaderByteCode.pShaderBytecode = reinterpret_cast<BYTE*>(pVS->GetBufferPointer());
    g_VertexShaderByteCode.BytecodeLength = pVS->GetBufferSize();

    g_PixelShaderByteCode.pShaderBytecode = reinterpret_cast<BYTE*>(pPS->GetBufferPointer());
    g_PixelShaderByteCode.BytecodeLength = pPS->GetBufferSize();

}

void CreatePSOs()
{
	D3D12_GRAPHICS_PIPELINE_STATE_DESC psoDesc;
	ZeroMemory(&psoDesc, sizeof(D3D12_GRAPHICS_PIPELINE_STATE_DESC));
	psoDesc.InputLayout = { g_ieds, _countof(g_ieds) };
	psoDesc.pRootSignature = g_pRootSignature;
	psoDesc.VS = g_VertexShaderByteCode;
	psoDesc.PS = g_PixelShaderByteCode;
	psoDesc.RasterizerState = CD3DX12_RASTERIZER_DESC(D3D12_DEFAULT);
	psoDesc.BlendState = CD3DX12_BLEND_DESC(D3D12_DEFAULT);
	psoDesc.DepthStencilState = CD3DX12_DEPTH_STENCIL_DESC(D3D12_DEFAULT);
	psoDesc.SampleMask = UINT_MAX;
	psoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
	psoDesc.NumRenderTargets = 1;
	psoDesc.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
	psoDesc.SampleDesc.Count = 1;
	psoDesc.SampleDesc.Quality = 0;
	psoDesc.DSVFormat = DXGI_FORMAT_D24_UNORM_S8_UINT;

    if (FAILED(g_pDevice->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pPso))))
    {
        return;
    }

	return;
}

void PopulateCommandList()
{
    // command list allocators can only be reset when the associated
    // command lists have finished execution on GPU; apps should use
    // fences to determine GPU execution progress.
    if (g_pCommandAllocator->Reset())
    {
        return;
    }

    // however, when ExecuteCommandList() is called on a particular command
    // list, that command list can then be reset at any time and must be before
    // re-recording
    if (g_pCommandList->Reset(g_pCommandAllocator, g_pPso))
    {
        return;
    }

    g_pCommandList->RSSetViewports(1, &g_ViewPort);
    g_pCommandList->RSSetScissorRects(1, &g_ScissorPort);

    // Indicate a state transition on the resource usage.
	g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pRenderTargets[g_nFrameIndex],
		D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET));

    // clear the back buffer and depth buffer
    const FLOAT clearColor[] = {0.0f, 0.1f, 0.2f, 1.0f};
	g_pCommandList->ClearRenderTargetView(CD3DX12_CPU_DESCRIPTOR_HANDLE(
		g_pRtvHeap->GetCPUDescriptorHandleForHeapStart(),
		g_nFrameIndex,
		g_nRtvDescriptorSize), clearColor, 0, nullptr);
    g_pCommandList->ClearDepthStencilView(g_pDsvHeap->GetCPUDescriptorHandleForHeapStart(), D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);

	g_pCommandList->OMSetRenderTargets(1, &CD3DX12_CPU_DESCRIPTOR_HANDLE(
		g_pRtvHeap->GetCPUDescriptorHandleForHeapStart(),
		g_nFrameIndex,
		g_nRtvDescriptorSize), true, &g_pDsvHeap->GetCPUDescriptorHandleForHeapStart());

    ID3D12DescriptorHeap* ppHeaps[] = {g_pCbvHeap};
    g_pCommandList->SetDescriptorHeaps(_countof(ppHeaps), ppHeaps);

    // Set necessary state.
    g_pCommandList->SetGraphicsRootSignature(g_pRootSignature);

    g_pCommandList->IASetVertexBuffers(0, 1, &g_VertexBufferView);
    g_pCommandList->IASetIndexBuffer(&g_IndexBufferView);
    g_pCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    g_pCommandList->SetGraphicsRootDescriptorTable(0, g_pCbvHeap->GetGPUDescriptorHandleForHeapStart());

    g_pCommandList->DrawIndexedInstanced(
        indices.size(), 1, 0, 0, 0);

	g_pCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(g_pRenderTargets[g_nFrameIndex],
		D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PRESENT));

    return;
}

void Update()
{
	Vector3Df cameraPos({ 5.0f, 5.0f, 5.0f });
	Vector3Df target({ 0.0f, 0.0f, 0.0f });
	Vector3Df up({ 0.0f, 1.0f, 0.0f });

	Matrix4f world;
	world.SetIdentity();
	Matrix4f view;
	BuildViewMatrixLH(view, cameraPos, target, up);
	Matrix4f proj;
	BuildPerspectiveFovLHMatrix(proj, PI / 4, (float)nScreenWidth / (float)nScreenHeight, 1.0f, 1000.0f);

	Matrix4f combine;
	combine = proj * view * world;
	Matrix4f aCombine;
	aCombine = world * view * proj;
	TransposeMatrix(aCombine, aCombine);

	BYTE* mappedData = nullptr;
	g_pConstantBuffer->Map(0, nullptr, reinterpret_cast<void**>(&mappedData));
	memcpy(&mappedData[0], aCombine, sizeof(ObjectConstants));
	g_pConstantBuffer->Unmap(0, nullptr);

	//////////////////////////////////////////////////////////////////////////////////

	// Convert Spherical to Cartesian coordinates.
	//float x = mRadius * sinf(mPhi)*cosf(mTheta);
	//float z = mRadius * sinf(mPhi)*sinf(mTheta);
	//float y = mRadius * cosf(mPhi);

	//XMFLOAT4X4 mWorld(
	//	1.0f, 0.0f, 0.0f, 0.0f,
	//	0.0f, 1.0f, 0.0f, 0.0f,
	//	0.0f, 0.0f, 1.0f, 0.0f,
	//	0.0f, 0.0f, 0.0f, 1.0f);
	//XMFLOAT4X4 mView(
	//	1.0f, 0.0f, 0.0f, 0.0f,
	//	0.0f, 1.0f, 0.0f, 0.0f,
	//	0.0f, 0.0f, 1.0f, 0.0f,
	//	0.0f, 0.0f, 0.0f, 1.0f
	//);
	//XMFLOAT4X4 mProj(
	//	1.0f, 0.0f, 0.0f, 0.0f,
	//	0.0f, 1.0f, 0.0f, 0.0f,
	//	0.0f, 0.0f, 1.0f, 0.0f,
	//	0.0f, 0.0f, 0.0f, 1.0f
	//);
	//// Build the view matrix.
	//XMVECTOR pos = XMVectorSet(0.0f, 0.0f, 5.0f, 1.0f);
	//XMVECTOR target = XMVectorZero();
	//XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);

	////mProj = XMMatrixPerspectiveFovLH(0.25f*PI, static_cast<float>(nScreenWidth) / nScreenHeight, 1.0f, 1000.0f);

	//XMMATRIX world = XMLoadFloat4x4(&mWorld);
	//XMMATRIX proj = XMMatrixPerspectiveFovLH(0.25f*PI, static_cast<float>(nScreenWidth) / nScreenHeight, 1.0f, 1000.0f);
	//XMMATRIX view = XMMatrixLookAtLH(pos, target, up);
	//XMStoreFloat4x4(&mView, view);
	//XMMATRIX worldViewProj = world * view*proj;


	//// Update the constant buffer with the latest worldViewProj matrix.
	//ObjectConstants objConstants;
	//XMStoreFloat4x4(objConstants.WorldViewProj.GetAddressOf(), XMMatrixTranspose(worldViewProj));

	//uint8_t* data;
	//g_pConstantBuffer->Map(0, nullptr, reinterpret_cast<void**>(&data));
	//memcpy(data, &objConstants, sizeof(ObjectConstants));
	//g_pConstantBuffer->Unmap(0, nullptr);
}

void Draw()
{
	g_pCommandList->Close();

    // execute the command list
    ID3D12CommandList* ppCommandLists[] = {g_pCommandList};
    g_pCommandQueue->ExecuteCommandLists(_countof(ppCommandLists), ppCommandLists);

    // swap the back buffer and the front buffer
    g_pSwapChain->Present(0, 0);

	g_nFrameIndex = (g_nFrameIndex + 1) % nFrameCount;
}

void FlushCommandQueue()
{
	// Advance the fence value to mark commands up to this fence point.
    g_nFenceValue++;

    // Add an instruction to the command queue to set a new fence point.  Because we 
	// are on the GPU timeline, the new fence point won't be set until the GPU finishes
	// processing all the commands prior to this Signal().
    g_pCommandQueue->Signal(g_pFence, g_nFenceValue);

	// Wait until the GPU has completed commands up to this fence point.
    if(g_pFence->GetCompletedValue() < g_nFenceValue)
	{
		HANDLE eventHandle = CreateEventEx(nullptr, false, false, EVENT_ALL_ACCESS);

        // Fire event when GPU hits current fence.  
        g_pFence->SetEventOnCompletion(g_nFenceValue, eventHandle);

        // Wait until the GPU hits current fence event is fired.
		WaitForSingleObject(eventHandle, INFINITE);
        CloseHandle(eventHandle);
	}
}

void InitDirect3D()
{
    CreateDevice();
    CreateCommandObjects();
    CreateSwapChain();
    CreateDescriptorHeaps();
    CreateRenderTarget();
	CreateDepthStencilBuffer();
    CreateConstantBuffer();
	CreateRootSignature();
	InitializeShaders();
    CreateVertexBuffer();
    CreateIndexBuffer();
	CreatePSOs();
}

ID3DBlob* CompileShader(const std::wstring& filename,
	const D3D_SHADER_MACRO* defines,
	const std::string& entrypoint,
	const std::string& target)
{
	UINT compileFlags = 0;

	HRESULT hr = S_OK;

	ID3DBlob* byteCode = nullptr;
	ID3DBlob* errors;
	hr = D3DCompileFromFile(filename.c_str(), defines, D3D_COMPILE_STANDARD_FILE_INCLUDE,
		entrypoint.c_str(), target.c_str(), compileFlags, 0, &byteCode, &errors);

	if (errors != nullptr)
		OutputDebugStringA((char*)errors->GetBufferPointer());

	return byteCode;
}

///////////////////////////////////////////////////////////////////////////////////////


//// include the basic windows header file
//#include <windows.h>
//#include <windowsx.h>
//#include <stdio.h>
//#include <tchar.h>
//#include <stdint.h>
//
//#include <d3d12.h>
//#include "d3dx12.h"
//#include <DXGI1_4.h>		// DXGI (DirectX Graphics Infrastructure) is an API used along with Direct3D. The basic idea of 
//							// DXGI is that some graphics related tasks are common to multipile graphics APIs. Here are some
//							// common functionality, swap chain, display adapters, monitors and supported display modes(
//							// resolution, refresh rate), surface formats
//#include <D3DCompiler.h>
//#include <DirectXMath.h>
//#include <DirectXPackedVector.h>
//#include <DirectXColors.h>
//
//#include <wrl/client.h>
//
//#include <string>
//#include <exception>
//
//namespace Panda {
//	// Helper class for COM exceptions
//	class com_exception : public std::exception {
//	public:
//		com_exception(HRESULT hr) : result(hr) {
//		}
//
//		virtual const char* what() const override {
//			static char s_str[64] = { 0 };
//			sprintf_s(s_str, "Failure with HRESULT of %08X",
//				static_cast<unsigned int>(result));
//			return s_str;
//		}
//
//	private:
//		HRESULT result;
//	};
//
//	// Helper utility converts D3D API failures into exceptions.
//	inline void ThrowIfFailed(HRESULT hr) {
//		if (FAILED(hr)) {
//			throw com_exception(hr);
//		}
//	}
//};
//
//using namespace Panda;
//using namespace DirectX;
//using namespace DirectX::PackedVector;
//using namespace Microsoft::WRL;
//using namespace std;
//
//ComPtr<IDXGIFactory4> 							g_pDXGIFactory = nullptr;			// dxgi factory
//ComPtr<ID3D12Device>							g_pDevice = nullptr;				// the pointer to out Direct3D device interface
//D3D12_FEATURE_DATA_MULTISAMPLE_QUALITY_LEVELS	g_msQualityLevels;	// multi-sample support information
//
//ComPtr<ID3D12CommandQueue>			g_pCommandQueue;				// the pointer to command queue
//ComPtr<ID3D12CommandAllocator>		g_pCommandAllocator;			// the pointer to command buffer allocator
//ComPtr<ID3D12GraphicsCommandList>	g_pCommandList;					// a list to store CPU commands, which will be submmited to GPU to execute
//
//const uint32_t 	g_FrameCount = 2;			// buffer count
//HWND 			g_hWnd;						// window handle
//
//ComPtr<IDXGISwapChain3> 				g_pSwapChain = nullptr;			// the pointer to the swap chain interface 
//uint32_t g_ScreenWidth = 960;
//uint32_t g_ScreenHeight = 480;
//
//ComPtr<ID3D12Fence>	g_pFence;	// the pointer to the fence
//uint32_t g_CurrentFence = 0; 	// current fence value
//uint32_t g_rtvDescriptorSize = 0;
//uint32_t g_dsvDescriptorSize = 0;
//ComPtr<ID3D12DescriptorHeap> g_pRtvHeap;
//ComPtr<ID3D12DescriptorHeap> g_pDsvHeap;
//
//uint32_t g_currentBackBuffer = 0;
//ComPtr<ID3D12Resource> g_SwapChainBuffer[g_FrameCount];
//ComPtr<ID3D12Resource> g_DepthStencilBuffer;
//
//ComPtr<ID3D12Resource> g_pDefaultVertexBuffer;
//D3D12_VERTEX_BUFFER_VIEW	g_DefaultVertexBufferView;
//ComPtr<ID3D12Resource> g_pUploadVertexBuffer;
//D3D12_VERTEX_BUFFER_VIEW g_UploadVertexBufferView;
//
//D3D12_VIEWPORT g_ScreenViewport;
//D3D12_RECT g_ScissorRect;
//
//ComPtr<ID3D12RootSignature>     g_pRootSignature;
//ComPtr<ID3D12PipelineState>     g_pPipelineState;
//
//struct VERTEX
//{
//	XMFLOAT3	Position;
//	XMFLOAT4	Color;
//};
//
//
//LRESULT CALLBACK WindowProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
//void OnResize();
//void FlushCommandQueue();
//
//bool InitMainWindow(HINSTANCE hInstance, int nCmdShow) {
//	WNDCLASSEX wc;
//
//	// clear out the window class for use
//	ZeroMemory(&wc, sizeof(WNDCLASSEX));
//
//	// fill in the struct with the needed information
//	wc.cbSize = sizeof(WNDCLASSEX);
//	wc.style = CS_HREDRAW | CS_VREDRAW;
//	wc.lpfnWndProc = WindowProc;
//	wc.hInstance = hInstance;
//	wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
//	wc.hbrBackground = (HBRUSH)COLOR_WINDOW;
//	wc.lpszClassName = _T("HelloD3D12");
//
//	// register the window class
//	RegisterClassEx(&wc);
//
//	// create the window and use the result as the handle
//	g_hWnd = CreateWindowEx(0,
//		_T("HelloD3D12"),                   // name of the window class
//		_T("Hello Direct 12"),      			// title of the window
//		WS_OVERLAPPEDWINDOW,                  // window style
//		100,                                  // x-position of the window
//		100,                                  // y-position of the window
//		g_ScreenWidth,                         // width of the window
//		g_ScreenHeight,                        // height of the window
//		NULL,                                 // we have no parent window, NULL
//		NULL,                                 // we aren't using menus, NULL
//		hInstance,                            // application handle
//		NULL);                                // used with multiple windows, NULL
//
//	if (!g_hWnd) {
//		MessageBoxW(0, L"Create window failed.", 0, 0);
//		return false;
//	}
//
//	// display the window on the screen
//	ShowWindow(g_hWnd, nCmdShow);
//	UpdateWindow(g_hWnd);
//
//	return true;
//}
//
///***************************************************************************************/
///*
// * DirectX 12
// */
//
//void CreateDevice() {
//	// Try to create hardware device.
//	HRESULT hardwareResult = D3D12CreateDevice(
//		nullptr,	// default adapter
//		D3D_FEATURE_LEVEL_11_0,
//		IID_PPV_ARGS(&g_pDevice));
//
//	// Fallback to WARP(Windows Advanced Rasterization Platform) device
//	ThrowIfFailed(CreateDXGIFactory1(IID_PPV_ARGS(&g_pDXGIFactory)));
//	if (FAILED(hardwareResult)) {
//		ComPtr<IDXGIAdapter> pWarpAdapter;
//		ThrowIfFailed(g_pDXGIFactory->EnumWarpAdapter(IID_PPV_ARGS(&pWarpAdapter)));
//
//		ThrowIfFailed(D3D12CreateDevice(
//			pWarpAdapter.Get(),
//			D3D_FEATURE_LEVEL_11_0,
//			IID_PPV_ARGS(&g_pDevice)));
//	}
//
//	// Check 4X MSAA quality support for our back buffer format.
//	// All Direct3D 11 capable devices support 4X MSAA for all render
//	// target formats, so we only need to check quality support.
//	g_msQualityLevels.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
//	g_msQualityLevels.SampleCount = 4;
//	g_msQualityLevels.Flags = D3D12_MULTISAMPLE_QUALITY_LEVELS_FLAG_NONE;
//	g_msQualityLevels.NumQualityLevels = 0;
//	ThrowIfFailed(g_pDevice->CheckFeatureSupport(
//		D3D12_FEATURE_MULTISAMPLE_QUALITY_LEVELS,
//		&g_msQualityLevels,
//		sizeof(g_msQualityLevels)));
//}
//
//void CreateCommandObjects() {
//	// Command Queue
//	D3D12_COMMAND_QUEUE_DESC queueDesc = {};
//	queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
//	queueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
//	ThrowIfFailed(g_pDevice->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&g_pCommandQueue)));
//
//	// command allocator and command list
//	ThrowIfFailed(g_pDevice->CreateCommandAllocator(
//		D3D12_COMMAND_LIST_TYPE_DIRECT,
//		IID_PPV_ARGS(g_pCommandAllocator.GetAddressOf())));
//
//	ThrowIfFailed(g_pDevice->CreateCommandList(
//		0,
//		D3D12_COMMAND_LIST_TYPE_DIRECT,
//		g_pCommandAllocator.Get(),	// assocated command allocator
//		nullptr, // initial pipelinestateobject
//		IID_PPV_ARGS(g_pCommandList.GetAddressOf())));
//
//	// close command list. 
//	g_pCommandList->Close();
//}
//
//void CreateSwapChain() {
//	// Release the previous swapchain we will be recreating.
//	//if (g_pSwapChain.Get() != nullptr)
//	g_pSwapChain.Reset();
//
//	DXGI_SWAP_CHAIN_DESC sd;
//	sd.BufferDesc.Width = g_ScreenWidth;
//	sd.BufferDesc.Height = g_ScreenHeight;
//	sd.BufferDesc.RefreshRate.Numerator = 60;
//	sd.BufferDesc.RefreshRate.Denominator = 1;
//	sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
//	sd.BufferDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED;
//	sd.BufferDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;
//	sd.SampleDesc.Count = 1;
//	sd.SampleDesc.Quality = 0;
//	sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
//	sd.BufferCount = g_FrameCount;
//	sd.OutputWindow = g_hWnd;
//	sd.Windowed = true;
//	sd.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD; // This flag cannot be used with multisampling and partial presentation
//	sd.Flags = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;
//
//	// Note: Swap chain uses queue to perform flush
//	ComPtr<IDXGISwapChain> swapChain;
//	ThrowIfFailed(g_pDXGIFactory->CreateSwapChain(
//		g_pCommandQueue.Get(),
//		&sd,
//		swapChain.GetAddressOf()));
//	swapChain.As(&g_pSwapChain);
//
//	g_currentBackBuffer = g_pSwapChain->GetCurrentBackBufferIndex();
//}
//
//void CreateDescriptorHeaps() {
//	// fence
//	ThrowIfFailed(g_pDevice->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&g_pFence)));
//
//	// RTVs
//	g_rtvDescriptorSize = g_pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);
//	D3D12_DESCRIPTOR_HEAP_DESC rtvHeapDesc;
//	rtvHeapDesc.NumDescriptors = g_FrameCount;
//	rtvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
//	rtvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
//	rtvHeapDesc.NodeMask = 0;
//	ThrowIfFailed(g_pDevice->CreateDescriptorHeap(
//		&rtvHeapDesc,
//		IID_PPV_ARGS(g_pRtvHeap.GetAddressOf())));
//
//	// DSVs
//	g_dsvDescriptorSize = g_pDevice->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_DSV);
//	D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc;
//	dsvHeapDesc.NumDescriptors = 1;
//	dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
//	dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
//	dsvHeapDesc.NodeMask = 0;
//	ThrowIfFailed(g_pDevice->CreateDescriptorHeap(
//		&dsvHeapDesc,
//		IID_PPV_ARGS(g_pDsvHeap.GetAddressOf())));
//}
//
//void CreateRtvsAndDsv() {
//	// Release the previous resources we will be recreating.
//	for (int i = 0; i < g_FrameCount; ++i)
//		g_SwapChainBuffer[i].Reset();
//	g_DepthStencilBuffer.Reset();
//
//	// Resize the swap chain
//	ThrowIfFailed(g_pSwapChain->ResizeBuffers(
//		g_FrameCount,
//		g_ScreenWidth, g_ScreenHeight,
//		DXGI_FORMAT_R8G8B8A8_UNORM,
//		DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH));
//
//	g_currentBackBuffer = 0;
//
//	// Create render target views.
//	CD3DX12_CPU_DESCRIPTOR_HANDLE rtvHeapHandle(g_pRtvHeap->GetCPUDescriptorHandleForHeapStart());
//	for (int i = 0; i < g_FrameCount; ++i) {
//		ThrowIfFailed(g_pSwapChain->GetBuffer(i, IID_PPV_ARGS(&g_SwapChainBuffer[i])));
//		g_pDevice->CreateRenderTargetView(g_SwapChainBuffer[i].Get(), nullptr, rtvHeapHandle);
//		rtvHeapHandle.Offset(1, g_rtvDescriptorSize);
//	}
//
//	// Create the depth/stencil buffer and view.
//	D3D12_RESOURCE_DESC depthStencilDesc;
//	depthStencilDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
//	depthStencilDesc.Alignment = 0;
//	depthStencilDesc.Width = g_ScreenWidth;
//	depthStencilDesc.Height = g_ScreenHeight;
//	depthStencilDesc.DepthOrArraySize = 1;
//	depthStencilDesc.MipLevels = 1;
//	depthStencilDesc.Format = DXGI_FORMAT_R24G8_TYPELESS;
//	depthStencilDesc.SampleDesc.Count = 1;
//	depthStencilDesc.SampleDesc.Quality = 0;
//	depthStencilDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
//	depthStencilDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;
//
//	// Create clear value.
//	D3D12_CLEAR_VALUE optClear;
//	optClear.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
//	optClear.DepthStencil.Depth = 1.0f;
//	optClear.DepthStencil.Stencil = 0;
//	ThrowIfFailed(g_pDevice->CreateCommittedResource(
//		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
//		D3D12_HEAP_FLAG_NONE,
//		&depthStencilDesc,
//		D3D12_RESOURCE_STATE_COMMON,
//		&optClear,
//		IID_PPV_ARGS(g_DepthStencilBuffer.GetAddressOf())));
//
//	// Create descriptor to mip level 0 of entire resource using the format of the resource.
//	D3D12_DEPTH_STENCIL_VIEW_DESC dsvDesc;
//	dsvDesc.Flags = D3D12_DSV_FLAG_NONE;
//	dsvDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2D;
//	dsvDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
//	dsvDesc.Texture2D.MipSlice = 0;
//	g_pDevice->CreateDepthStencilView(g_DepthStencilBuffer.Get(), &dsvDesc, g_pDsvHeap->GetCPUDescriptorHandleForHeapStart());
//
//	// Transtion the resource from its initial
//	g_pCommandList->ResourceBarrier(
//		1,
//		&CD3DX12_RESOURCE_BARRIER::Transition(
//			g_DepthStencilBuffer.Get(),
//			D3D12_RESOURCE_STATE_COMMON,
//			D3D12_RESOURCE_STATE_DEPTH_WRITE));
//
//}
//
//void CreateVertexBuffer()
//{
//	ThrowIfFailed(g_pCommandList->Reset(g_pCommandAllocator.Get(), nullptr));
//
//	VERTEX OutVertices[] = 
//	{
//		{XMFLOAT3(0.0f, 0.5f, 0.0f), XMFLOAT4(1.0f, 0.0f, 0.0f, 1.0f)},
//		{XMFLOAT3(0.45f, -0.5, 0.0f), XMFLOAT4(0.0f, 1.0f, 0.0f, 1.0f)},
//		{XMFLOAT3(-0.45f, -0.5f, 0.0f), XMFLOAT4(0.0f, 0.0f, 1.0f, 1.0f)}
//	};
//	const UINT vertexBufferSize = sizeof(OutVertices);
//
//	// Create default buffer resource
//	ThrowIfFailed(g_pDevice->CreateCommittedResource(
//		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
//		D3D12_HEAP_FLAG_NONE,
//		&CD3DX12_RESOURCE_DESC::Buffer(vertexBufferSize),
//		D3D12_RESOURCE_STATE_COMMON,
//		nullptr,
//		IID_PPV_ARGS(&g_pDefaultVertexBuffer)));
//
//	ThrowIfFailed(g_pDevice->CreateCommittedResource(
//		&CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
//		D3D12_HEAP_FLAG_NONE,
//		&CD3DX12_RESOURCE_DESC::Buffer(vertexBufferSize),
//		D3D12_RESOURCE_STATE_GENERIC_READ,
//		nullptr,
//		IID_PPV_ARGS(&g_pUploadVertexBuffer)));
//
//	// Describe the data
//	D3D12_SUBRESOURCE_DATA subResourceData = {};
//	subResourceData.pData = OutVertices;
//	subResourceData.RowPitch = vertexBufferSize;
//	subResourceData.SlicePitch = vertexBufferSize;
//	
//	// upload the vertex
//	g_pCommandList->ResourceBarrier(1,
//		&CD3DX12_RESOURCE_BARRIER::Transition(g_pDefaultVertexBuffer.Get(),
//			D3D12_RESOURCE_STATE_COMMON,
//			D3D12_RESOURCE_STATE_COPY_DEST));
//
//	UpdateSubresources<1>(g_pCommandList.Get(), g_pDefaultVertexBuffer.Get(), g_pUploadVertexBuffer.Get(),
//		0, 0, 1, &subResourceData);
//
//	g_pCommandList->ResourceBarrier(1,
//		&CD3DX12_RESOURCE_BARRIER::Transition(g_pDefaultVertexBuffer.Get(),
//			D3D12_RESOURCE_STATE_COPY_DEST,
//			D3D12_RESOURCE_STATE_GENERIC_READ));
//
//	// create the vertex buffer descriptor, or just call it "view"
//	g_DefaultVertexBufferView.BufferLocation = g_pDefaultVertexBuffer->GetGPUVirtualAddress();
//	g_DefaultVertexBufferView.StrideInBytes = sizeof(VERTEX);
//	g_DefaultVertexBufferView.SizeInBytes = vertexBufferSize;
//
//	g_UploadVertexBufferView.BufferLocation = g_pUploadVertexBuffer->GetGPUVirtualAddress();
//	g_UploadVertexBufferView.StrideInBytes = sizeof(VERTEX);
//	g_UploadVertexBufferView.SizeInBytes = vertexBufferSize;
//
//	// Execute the initialization commands.
//	ThrowIfFailed(g_pCommandList->Close());
//	ID3D12CommandList* cmdsLists[] = { g_pCommandList.Get() };
//	g_pCommandQueue->ExecuteCommandLists(_countof(cmdsLists), cmdsLists);
//
//	// Wait until initialization is complete.
//	FlushCommandQueue();
//}
//
//void CreateRootSignature()
//{
//	// create an empty root signature
//	CD3DX12_ROOT_SIGNATURE_DESC rsd;
//	rsd.Init(0, nullptr, 0, nullptr, D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);
//
//	ComPtr<ID3DBlob> signature;
//	ComPtr<ID3DBlob> error;
//	ThrowIfFailed(D3D12SerializeRootSignature(&rsd, D3D_ROOT_SIGNATURE_VERSION_1, &signature, &error));
//	ThrowIfFailed(g_pDevice->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&g_pRootSignature)));
//}
//
//void InitPipeline()
//{
//	ComPtr<ID3DBlob> error;
//
//	ComPtr<ID3DBlob> vertexShader;
//	ComPtr<ID3DBlob> pixelShader;
//
//	D3DCompileFromFile(
//		L"copy.vs",
//		nullptr,
//		D3D_COMPILE_STANDARD_FILE_INCLUDE,
//		"main",
//		"vs_5_0",
//		0,
//		0,
//		&vertexShader,
//		&error);
//	if (error) { OutputDebugString((LPCTSTR)error->GetBufferPointer()); error->Release(); throw std::exception(); }
//
//	D3DCompileFromFile(
//		L"copy.ps",
//		nullptr,
//		D3D_COMPILE_STANDARD_FILE_INCLUDE,
//		"main",
//		"ps_5_0",
//		0,
//		0,
//		&pixelShader,
//		&error);
//	if (error) { OutputDebugString((LPCTSTR)error->GetBufferPointer()); error->Release(); throw std::exception(); }
//
//	// create the input layout object
//	D3D12_INPUT_ELEMENT_DESC ied[] =
//	{
//		{"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
//		{"COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
//	};
//
//	// describe and create the graphics pipeline state object (PSO)
//	D3D12_GRAPHICS_PIPELINE_STATE_DESC psod = {};
//	psod.InputLayout = { ied, _countof(ied) };
//	psod.pRootSignature = g_pRootSignature.Get();
//	psod.VS = { reinterpret_cast<UINT8*>(vertexShader->GetBufferPointer()), vertexShader->GetBufferSize() };
//	psod.PS = { reinterpret_cast<UINT8*>(pixelShader->GetBufferPointer()), pixelShader->GetBufferSize() };
//	psod.RasterizerState = CD3DX12_RASTERIZER_DESC(D3D12_DEFAULT);
//	psod.BlendState = CD3DX12_BLEND_DESC(D3D12_DEFAULT);
//	psod.DepthStencilState.DepthEnable = FALSE;
//	psod.DepthStencilState.StencilEnable = FALSE;
//	psod.SampleMask = UINT_MAX;
//	psod.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
//	psod.NumRenderTargets = 1;
//	psod.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
//	psod.SampleDesc.Count = 1;
//	ThrowIfFailed(g_pDevice->CreateGraphicsPipelineState(&psod, IID_PPV_ARGS(&g_pPipelineState)));
//}
//
//void FlushCommandQueue() {
//	const uint64_t fence = g_CurrentFence;
//	ThrowIfFailed(g_pCommandQueue->Signal(g_pFence.Get(), fence));
//
//	// Advance the fence value to mark commands up to this fence point.
//	g_CurrentFence++;
//
//	// Wait until the GPU has completed commands up to this fence point.
//	if (g_pFence->GetCompletedValue() < fence) {
//		HANDLE eventHandle = CreateEventEx(nullptr, false, false, EVENT_ALL_ACCESS);
//
//		// Set fire event when GPU hits current fence.
//		ThrowIfFailed(g_pFence->SetEventOnCompletion(fence, eventHandle));
//
//		// Wait until the GPU hits current fence event is fired.
//		WaitForSingleObject(eventHandle, INFINITE);
//		CloseHandle(eventHandle);
//	}
//
//	g_currentBackBuffer = g_pSwapChain->GetCurrentBackBufferIndex();
//}
//
//void OnResize() {
//	// We won't resize it until the device,swap chain, command allocator are created.
//	if (g_pDevice == nullptr || g_pSwapChain == nullptr || g_pCommandAllocator == nullptr) {
//		return;
//	}
//
//	// Fluse before changing any resources.
//	FlushCommandQueue();
//
//	ThrowIfFailed(g_pCommandList->Reset(g_pCommandAllocator.Get(), nullptr));
//
//	CreateRtvsAndDsv();
//
//	// Execute the resize commands.
//	ThrowIfFailed(g_pCommandList->Close());
//	ID3D12CommandList* cmdLists[] = { g_pCommandList.Get() };
//	g_pCommandQueue->ExecuteCommandLists(_countof(cmdLists), cmdLists);
//
//	// Wait until resize is complete.
//	FlushCommandQueue();
//
//	// Update the viewport transform to cover the client area.
//	g_ScreenViewport.TopLeftX = 0;
//	g_ScreenViewport.TopLeftY = 0;
//	g_ScreenViewport.Width = static_cast<float>(g_ScreenWidth);
//	g_ScreenViewport.Height = static_cast<float>(g_ScreenHeight);
//	g_ScreenViewport.MinDepth = 0.0f;
//	g_ScreenViewport.MaxDepth = 1.0f;
//
//	g_ScissorRect = { 0, 0, static_cast<LONG>(g_ScreenWidth), static_cast<LONG>(g_ScreenHeight) };
//}
//
//void InitDirect3D12() {
//	CreateDevice();
//
//	CreateCommandObjects();
//
//	CreateSwapChain();
//
//	CreateDescriptorHeaps();
//
//	OnResize();
//
//	CreateRootSignature();
//
//	CreateVertexBuffer();
//
//	InitPipeline();
//
//
//}
//
//void Draw() {
//	// Reuse the memory associated with command recording.
//	// We can only reset when the associated command lists have finished execution on the GPU.
//	ThrowIfFailed(g_pCommandAllocator->Reset());
//
//	// A command list can be reset after it has been added to the command queue via ExecuteCommandList.
//	// Reusing the command list reuses memory.
//	ThrowIfFailed(g_pCommandList->Reset(g_pCommandAllocator.Get(), g_pPipelineState.Get()));
//
//	g_pCommandList->SetGraphicsRootSignature(g_pRootSignature.Get());
//	// Set the viewport and scissor rect. This needs to be reset whenever the command list is reset.
//	g_pCommandList->RSSetViewports(1, &g_ScreenViewport);
//	g_pCommandList->RSSetScissorRects(1, &g_ScissorRect);
//
//	// Indicate a state transition on the resource usage.
//	g_pCommandList->ResourceBarrier(1,
//		&CD3DX12_RESOURCE_BARRIER::Transition(g_SwapChainBuffer[g_currentBackBuffer].Get(),
//			D3D12_RESOURCE_STATE_PRESENT,
//			D3D12_RESOURCE_STATE_RENDER_TARGET));
//
//	// Specify the buffers we are going to render to.
//	g_pCommandList->OMSetRenderTargets(1,
//		&CD3DX12_CPU_DESCRIPTOR_HANDLE(
//			g_pRtvHeap->GetCPUDescriptorHandleForHeapStart(),
//			g_currentBackBuffer,
//			g_rtvDescriptorSize),
//		true,
//		&g_pDsvHeap->GetCPUDescriptorHandleForHeapStart());
//
//	// Clear the back buffer and depth buffer.
//	g_pCommandList->ClearRenderTargetView(
//		CD3DX12_CPU_DESCRIPTOR_HANDLE(
//			g_pRtvHeap->GetCPUDescriptorHandleForHeapStart(),
//			g_currentBackBuffer,
//			g_rtvDescriptorSize),
//		Colors::LightSteelBlue,
//		0,
//		nullptr);
//	g_pCommandList->ClearDepthStencilView(
//		g_pDsvHeap->GetCPUDescriptorHandleForHeapStart(),
//		D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL,
//		1.0f,
//		0,
//		0,
//		nullptr);
//
//	{
//		//g_pCommandList->IASetVertexBuffers(0, 1, &g_UploadVertexBufferView);
//		g_pCommandList->IASetVertexBuffers(0, 1, &g_DefaultVertexBufferView);
//
//		g_pCommandList->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
//
//		g_pCommandList->DrawInstanced(3, 1, 0, 0);
//	}
//
//
//	// Indicate a state transition on the resource usage.
//	g_pCommandList->ResourceBarrier(
//		1,
//		&CD3DX12_RESOURCE_BARRIER::Transition(
//			g_SwapChainBuffer[g_currentBackBuffer].Get(),
//			D3D12_RESOURCE_STATE_RENDER_TARGET,
//			D3D12_RESOURCE_STATE_PRESENT));
//
//	// Done recording commands.
//	ThrowIfFailed(g_pCommandList->Close());
//
//	// Add the command list to the queue for execution.
//	ID3D12CommandList* cmdsLists[] = { g_pCommandList.Get() };
//	g_pCommandQueue->ExecuteCommandLists(_countof(cmdsLists), cmdsLists);
//
//	// swap the back and fornt buffers.
//	ThrowIfFailed(g_pSwapChain->Present(0, 0));
//	
//
//	// Wait until frame commands are complete. This waiting is inefficient and is 
//	// done for simplicity. But it's enough for a show case.
//	FlushCommandQueue();
//}
//
///***************************************************************************************/
//
//int Run() {
//	// this struct holds Windows event messages
//	MSG msg = {};
//
//	// wait for the next message in the queue, store the result in 'msg'
//	while (msg.message != WM_QUIT) {
//		// If there are window message then process them
//		if (PeekMessage(&msg, 0, 0, 0, PM_REMOVE)) {
//			// translate keystroke messages into the right format
//			TranslateMessage(&msg);
//
//			// send the message to the WindowProc function
//			DispatchMessage(&msg);
//		}
//		else {
//			Sleep(100);
//			Draw();
//		}
//	}
//
//	// return this part of the WM_QUIT message to Windows
//	return msg.wParam;
//}
//
//// this is the main message handler for the program
//LRESULT CALLBACK WindowProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) {
//	LRESULT result = 0;
//	bool wasHandled = false;
//
//	// sort through and find what code to run for the message given
//	switch (message)
//	{
//	case WM_CREATE:
//		wasHandled = true;
//		break;
//
//		// case WM_PAINT:
//			// wasHandled = true;
//			// break;
//
//	case WM_SIZE:
//		// Save the new client area dimensions.
//		g_ScreenWidth = LOWORD(lParam);
//		g_ScreenHeight = HIWORD(lParam);
//		if (g_pDevice) {
//			if (wParam == SIZE_MAXIMIZED || wParam == SIZE_RESTORED) {
//				OnResize();
//			}
//		}
//		wasHandled = true;
//		break;
//
//	case WM_DESTROY:
//		FlushCommandQueue();
//		PostQuitMessage(0);
//		wasHandled = true;
//		break;
//	}
//
//	// Handle any messages the switch statement didn't
//	if (!wasHandled) { result = DefWindowProc(hWnd, message, wParam, lParam); }
//	return result;
//}
//
//// the entry point for any Windows program
//int WINAPI WinMain(HINSTANCE hInstance,
//	HINSTANCE hPrevInstance,
//	LPTSTR lpCmdLine,
//	int nCmdShow)
//{
//	// Create window
//	InitMainWindow(hInstance, nCmdShow);
//
//	// Initialize DirectX 12
//	InitDirect3D12();
//
//	OnResize();
//
//	return Run();
//}