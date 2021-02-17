#include <iostream>
#include <fstream>
#include <objbase.h>
#include <d3dcompiler.h>
#include "d3dx12.h"
#include "D3D12GraphicsManager.hpp"
#include "Common/Entry/WindowsApplication.hpp"
#include "SceneManager/SceneManager.hpp"
#include "Core/AssetLoader.hpp"
#include "D3DShaderModule.hpp"

namespace Panda {
#define MSAA_SOLUTION_DEFAULT	1

extern IApplication* g_pApp;

	static void OutputShaderErrorMessage(std::string str)
	{
		std::string filename = "shader-error.txt";
		std::ofstream fout;
		fout.open(filename);

		fout << str;

		fout.close();

		std::cerr << "Error compiling shader. Check shader-error.txt for more infomation." << std::endl;
	}
	
    static void GetHardwareAdapter(IDXGIFactory4* pFactory, IDXGIAdapter1** ppAdapter)
    {
        IDXGIAdapter1* pAdapter = nullptr;
        *ppAdapter = nullptr;

        for (UINT adapterIndex = 0; DXGI_ERROR_NOT_FOUND != pFactory->EnumAdapters1(adapterIndex, &pAdapter); adapterIndex++)
        {
           DXGI_ADAPTER_DESC1 desc;
           pAdapter->GetDesc1(&desc);

           if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE)
           {
               // Don't select the Basic Render Driver adapter.
               continue;
           }

           // Check to see if the adapter supports Direct3D 12, but don't create the
           // actual device yet.
           if (SUCCEEDED(D3D12CreateDevice(pAdapter, D3D_FEATURE_LEVEL_11_0, __uuidof(ID3D12Device), nullptr)))
           {
               break;
           }
        }

        *ppAdapter = pAdapter;
    }

    HRESULT D3D12GraphicsManager::WaitForPreviousFrame()
    {
        // WARNING: Waiting for the frame to complete before continuing is not the best practice.
        // Here is just for simplicity. More advanced samples illustrate how
        // to use fences for efficient resource usage.

        // Signal and increment the fence value.
        HRESULT hr;
        const uint64_t fence = m_FenceValue;
        if (FAILED(hr = m_pCommandQueue->Signal(m_pFence, fence)))
        {
            return hr;
        }

        m_FenceValue++;

        // Wait until the previous frame is finished.
        if (m_pFence->GetCompletedValue() < fence)
        {
            if (FAILED(hr = m_pFence->SetEventOnCompletion(fence, m_FenceEvent)))
            {
                return hr;
            }
            WaitForSingleObject(m_FenceEvent, INFINITE);
        }

        m_FrameIndex = m_pSwapChain->GetCurrentBackBufferIndex();

        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateDescriptorHeaps()
    {
        HRESULT hr;

        /// Describe and create render target view (RTV) descriptor heap
        D3D12_DESCRIPTOR_HEAP_DESC rtvHeapDesc = {};
        rtvHeapDesc.NumDescriptors = k_FrameCount + 1;  // +1 for MSAA Resolver
        rtvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
        rtvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
        if (FAILED(hr = m_pDev->CreateDescriptorHeap(&rtvHeapDesc, IID_PPV_ARGS(&m_pRtvHeap))))
        {
            return hr;
        }

        m_RtvDescriptorSize = m_pDev->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV); // heap size

        /// Describe and create a depth stencil view (DSV) descriptor heap.
        D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
        dsvHeapDesc.NumDescriptors = 1;
        dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
        dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
        if (FAILED(hr = m_pDev->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&m_pDsvHeap))))
        {
            return hr;
        }

        /// Describe and create a Shader Resource View (SRV) and Constant Buffer View(CBV) and Unordered Access View(UAV) descriptor heap
        D3D12_DESCRIPTOR_HEAP_DESC cbvSrvUavHeapDesc = {};
        cbvSrvUavHeapDesc.NumDescriptors = 
            k_FrameCount * (2 * k_MaxSceneObjectCount) +      //  1 perFrame and 1 per DrawBatch
            k_MaxTextureCount;                // k_MaxTextureCount for srv(texture)
        cbvSrvUavHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
        cbvSrvUavHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
        if (FAILED(hr = m_pDev->CreateDescriptorHeap(&cbvSrvUavHeapDesc, IID_PPV_ARGS(&m_pCbvHeap))))
        {
            return hr;
        }

        m_CbvSrvDescriptorSize = m_pDev->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV); // heap size

        /// Describe and create a sampler descriptor heap.
        D3D12_DESCRIPTOR_HEAP_DESC samplerHeapDesc = {};
        samplerHeapDesc.NumDescriptors = k_MaxTextureCount;   // this is the max D3D12 hardware support currently
        samplerHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_SAMPLER;
        samplerHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
        if (FAILED(hr = m_pDev->CreateDescriptorHeap(&samplerHeapDesc, IID_PPV_ARGS(&m_pSamplerHeap))))
        {
            return hr;
        }

        if (FAILED(hr = m_pDev->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&m_pCommandAllocator))))
        {
            return hr;
        }

        return hr;
    }
    HRESULT D3D12GraphicsManager::CreateRenderTarget() 
    {
        HRESULT hr = S_OK;

        D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle = m_pRtvHeap->GetCPUDescriptorHandleForHeapStart();

        // Create a RTV for each frame.
        for (uint32_t i = 0; i < k_FrameCount; i++)
        {
            if (FAILED(hr = m_pSwapChain->GetBuffer(i, IID_PPV_ARGS(&m_pRenderTargets[i])))) {
                return hr;
            }
            m_pDev->CreateRenderTargetView(m_pRenderTargets[i], nullptr, rtvHandle);
            rtvHandle.ptr += m_RtvDescriptorSize;
        }

        // Create intermediate MSAA render target
        D3D12_RENDER_TARGET_VIEW_DESC   renderTargetDesc;
        renderTargetDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
        renderTargetDesc.ViewDimension = D3D12_RTV_DIMENSION_TEXTURE2DMS;

        D3D12_CLEAR_VALUE optimizedClearValue = {};
        optimizedClearValue.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
        optimizedClearValue.Color[0] = 0.0f;
        optimizedClearValue.Color[1] = 0.1f;
        optimizedClearValue.Color[2] = 0.2f;
        optimizedClearValue.Color[3] = 1.0f;

        D3D12_HEAP_PROPERTIES prop = {};
        prop.Type = D3D12_HEAP_TYPE_DEFAULT;
        prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
        prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
        prop.CreationNodeMask = 1;
        prop.VisibleNodeMask = 1;

        D3D12_RESOURCE_DESC textureDesc = {};
        textureDesc.MipLevels = 1;
        textureDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
        textureDesc.Width = g_pApp->GetConfiguration().screenWidth;
        textureDesc.Height = g_pApp->GetConfiguration().screenHeight;
        textureDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET;
        textureDesc.DepthOrArraySize = 1;
        textureDesc.SampleDesc.Count = 4;
        textureDesc.SampleDesc.Quality = DXGI_STANDARD_MULTISAMPLE_QUALITY_PATTERN;
        textureDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;

        if (FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &textureDesc,
            D3D12_RESOURCE_STATE_RENDER_TARGET,
            &optimizedClearValue,
            IID_PPV_ARGS(&m_pMsaaRenderTarget)
        )))
            return hr;

        m_pDev->CreateRenderTargetView(m_pMsaaRenderTarget, &renderTargetDesc, rtvHandle);

#if (!MSAA_SOLUTION_DEFAULT)
		//Describe and create a SRV for the texture
		D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
		srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
		srvDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2DMS;
		srvDesc.Texture2D.MipLevels = 4;
		srvDesc.Texture2D.MostDetailedMip = 0;
		D3D12_CPU_DESCRIPTOR_HANDLE srvHandle;
		size_t textureID = static_cast<uint32_t>(m_TextureIndex.size());
		srvHandle.ptr = m_pCbvHeap->GetCPUDescriptorHandleForHeapStart().ptr + (k_TextureDescStartIndex + textureID) * m_CbvSrvDescriptorSize;
		m_pDev->CreateShaderResourceView(m_pMsaaRenderTarget, &srvDesc, srvHandle);
		m_TextureIndex["MSAA"] = textureID;
#endif 
        

        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateDepthStencil()
    {
        HRESULT hr;

        // Create the depth stencil view
        D3D12_DEPTH_STENCIL_VIEW_DESC depthStencilDesc = {};
        depthStencilDesc.Format = DXGI_FORMAT_D32_FLOAT;
        depthStencilDesc.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2DMS;
        depthStencilDesc.Flags = D3D12_DSV_FLAG_NONE;

        D3D12_CLEAR_VALUE depthOptimizedClearValue = {};
        depthOptimizedClearValue.Format = DXGI_FORMAT_D32_FLOAT;
        depthOptimizedClearValue.DepthStencil.Depth = 1.0f;
        depthOptimizedClearValue.DepthStencil.Stencil = 0;

        D3D12_HEAP_PROPERTIES prop = {};
        prop.Type = D3D12_HEAP_TYPE_DEFAULT;
        prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
        prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
        prop.CreationNodeMask = 1;
        prop.VisibleNodeMask = 1;

        uint32_t width = g_pApp->GetConfiguration().screenWidth;
        uint32_t height = g_pApp->GetConfiguration().screenHeight;
        D3D12_RESOURCE_DESC resourceDesc = {};
        resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;
        resourceDesc.Alignment = 0;
        resourceDesc.Width = width;
        resourceDesc.Height = height;
        resourceDesc.DepthOrArraySize = 1;
        resourceDesc.MipLevels = 1;
        resourceDesc.Format = DXGI_FORMAT_D32_FLOAT;
        resourceDesc.SampleDesc.Count = 4;
        resourceDesc.SampleDesc.Quality = DXGI_STANDARD_MULTISAMPLE_QUALITY_PATTERN;
		//resourceDesc.SampleDesc.Count = 1;
		//resourceDesc.SampleDesc.Quality = 0;
        resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN;
        resourceDesc.Flags = D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;

        if (FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &resourceDesc,
            D3D12_RESOURCE_STATE_DEPTH_WRITE,
            &depthOptimizedClearValue,
            IID_PPV_ARGS(&m_pDepthStencilBuffer)
        )))
        {
            return hr;
        }

        m_pDev->CreateDepthStencilView(m_pDepthStencilBuffer, &depthStencilDesc, m_pDsvHeap->GetCPUDescriptorHandleForHeapStart());

        return hr;
    }

	HRESULT D3D12GraphicsManager::CreateIntervalVertexBuffer()
	{
		HRESULT hr;
		ID3D12Resource* pVertexBufferUploadHeap;
		float fullScreenQuad[] = {
			-1.0f, 1.0f, 0.0f,		0.0f, 0.0f,
			-1.0f, -1.0f, 0.0f,		0.0f, 1.0f,
			1.0f, 1.0f, 0.0f,		1.0f, 0.0f,
			1.0f, -1.0f, 0.0f,		1.0f, 1.0f
		};

		// create vertex GPU heap
		D3D12_HEAP_PROPERTIES prop = {};
		prop.Type = D3D12_HEAP_TYPE_DEFAULT;
		prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
		prop.CreationNodeMask = 1;
		prop.VisibleNodeMask = 1;

		uint32_t size = sizeof(float) * (3 + 2) * 4;
		D3D12_RESOURCE_DESC resourceDesc = {};
		resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
		resourceDesc.Alignment = 0;
		resourceDesc.Width = size;
		resourceDesc.Height = 1;
		resourceDesc.DepthOrArraySize = 1;
		resourceDesc.MipLevels = 1;
		resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
		resourceDesc.SampleDesc.Count = 1;
		resourceDesc.SampleDesc.Quality = 0;
		resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
		resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

		ID3D12Resource* pVertexBuffer;

		if (FAILED(hr = m_pDev->CreateCommittedResource(
			&prop,
			D3D12_HEAP_FLAG_NONE,
			&resourceDesc,
			D3D12_RESOURCE_STATE_COPY_DEST,
			nullptr,
			IID_PPV_ARGS(&pVertexBuffer))))
		{
			return hr;
		}

		prop.Type = D3D12_HEAP_TYPE_UPLOAD;
		prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
		prop.CreationNodeMask = 1;
		prop.VisibleNodeMask = 1;
		
		if (FAILED(hr = m_pDev->CreateCommittedResource(
			&prop,
			D3D12_HEAP_FLAG_NONE,
			&resourceDesc,
			D3D12_RESOURCE_STATE_GENERIC_READ,
			nullptr,
			IID_PPV_ARGS(&pVertexBufferUploadHeap))))
		{
			return hr;
		}

		D3D12_SUBRESOURCE_DATA vertexData = {};
		vertexData.pData = fullScreenQuad;

		UpdateSubresources<1>(m_pCommandList, pVertexBuffer, pVertexBufferUploadHeap, 0, 0, 1, &vertexData);
		D3D12_RESOURCE_BARRIER barrier = {};
		barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		barrier.Transition.pResource = pVertexBuffer;
		barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_DEST;
		barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_VERTEX_AND_CONSTANT_BUFFER;
		barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		m_pCommandList->ResourceBarrier(1, &barrier);

		// initialize the vertex buffer view
		m_VertexBufferViewResolve.BufferLocation = pVertexBuffer->GetGPUVirtualAddress();
		m_VertexBufferViewResolve.StrideInBytes = 20;
		m_VertexBufferViewResolve.SizeInBytes = size;

		m_Buffers.push_back(pVertexBuffer);
		m_Buffers.push_back(pVertexBufferUploadHeap);

		return hr;
	}

    HRESULT D3D12GraphicsManager::CreateVertexBuffer(const SceneObjectVertexArray& vPropertyArray)
    {
        HRESULT hr;

        ID3D12Resource* pVertexBufferUploadHeap;

        // create vertex GPU heap
        D3D12_HEAP_PROPERTIES prop = {};
        prop.Type = D3D12_HEAP_TYPE_DEFAULT;
        prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
        prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
        prop.CreationNodeMask = 1;
        prop.VisibleNodeMask = 1;

        D3D12_RESOURCE_DESC resourceDesc = {};
        resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
        resourceDesc.Alignment = 0;
        resourceDesc.Width = vPropertyArray.GetDataSize();
        resourceDesc.Height = 1;
        resourceDesc.DepthOrArraySize = 1;
        resourceDesc.MipLevels = 1;
        resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
        resourceDesc.SampleDesc.Count = 1;
        resourceDesc.SampleDesc.Quality = 0;
        resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
        resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

        ID3D12Resource* pVertexBuffer;

        if(FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &resourceDesc,
            D3D12_RESOURCE_STATE_COPY_DEST,
            nullptr,
            IID_PPV_ARGS(&pVertexBuffer))))
        {
            return hr;
        }

        prop.Type = D3D12_HEAP_TYPE_UPLOAD;
        prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
        prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
        prop.CreationNodeMask = 1;
        prop.VisibleNodeMask = 1;

        if (FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &resourceDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&pVertexBufferUploadHeap)
        )))
        {
            return hr;
        }

        D3D12_SUBRESOURCE_DATA vertexData = {};
        vertexData.pData = vPropertyArray.GetData();

        UpdateSubresources<1>(m_pCommandList, pVertexBuffer, pVertexBufferUploadHeap, 0, 0, 1, &vertexData);
        D3D12_RESOURCE_BARRIER barrier = {};
        barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        barrier.Transition.pResource = pVertexBuffer;
        barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_DEST;
        barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_VERTEX_AND_CONSTANT_BUFFER;
        barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
        m_pCommandList->ResourceBarrier(1, &barrier);

        // initialize the vertex buffer view
        D3D12_VERTEX_BUFFER_VIEW vertexBufferView;
        vertexBufferView.BufferLocation = pVertexBuffer->GetGPUVirtualAddress();
        vertexBufferView.StrideInBytes = (UINT)(vPropertyArray.GetDataSize() / vPropertyArray.GetVertexCount());
        vertexBufferView.SizeInBytes = (UINT)vPropertyArray.GetDataSize();
        m_VertexBufferView.push_back(vertexBufferView);

        m_Buffers.push_back(pVertexBuffer);
        m_Buffers.push_back(pVertexBufferUploadHeap);

        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateIndexBuffer(const SceneObjectIndexArray& indexArray)
    {
        HRESULT hr;

        ID3D12Resource* pIndexBufferUploadHeap;

        // create index GPU heap
        D3D12_HEAP_PROPERTIES prop = {};
        prop.Type = D3D12_HEAP_TYPE_DEFAULT;
        prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
        prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
        prop.CreationNodeMask = 1;
        prop.VisibleNodeMask = 1;

        D3D12_RESOURCE_DESC resourceDesc = {};
        resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
        resourceDesc.Alignment = 0;
        resourceDesc.Width = indexArray.GetDataSize();
        resourceDesc.Height = 1;
        resourceDesc.DepthOrArraySize = 1;
        resourceDesc.MipLevels = 1;
        resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
        resourceDesc.SampleDesc.Count = 1;
        resourceDesc.SampleDesc.Quality = 0;
        resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
        resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

        ID3D12Resource* pIndexBuffer;

        if (FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &resourceDesc,
            D3D12_RESOURCE_STATE_COPY_DEST,
            nullptr,
            IID_PPV_ARGS(&pIndexBuffer)
        )))
        {
            return hr;
        }

        prop.Type = D3D12_HEAP_TYPE_UPLOAD;

        if (FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &resourceDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&pIndexBufferUploadHeap)
        )))
        {
            return hr;
        }

        D3D12_SUBRESOURCE_DATA indexData = {};
        indexData.pData = indexArray.GetData();

        UpdateSubresources<1>(m_pCommandList, pIndexBuffer, pIndexBufferUploadHeap, 0, 0, 1, &indexData);
        D3D12_RESOURCE_BARRIER barrier = {};
        barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
        barrier.Transition.pResource = pIndexBuffer;
        barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_DEST;
        barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_INDEX_BUFFER;
        barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
        m_pCommandList->ResourceBarrier(1, &barrier);

        // initialize the index buffer view
        D3D12_INDEX_BUFFER_VIEW indexBufferView;
        indexBufferView.BufferLocation = pIndexBuffer->GetGPUVirtualAddress();
        indexBufferView.Format = DXGI_FORMAT_R32_UINT;
        indexBufferView.SizeInBytes = (UINT)indexArray.GetDataSize();
        m_IndexBufferView.push_back(indexBufferView);

        m_Buffers.push_back(pIndexBuffer);
        m_Buffers.push_back(pIndexBufferUploadHeap);

        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateTextureBuffer(SceneObjectTexture& texture)
    {
        HRESULT hr = S_OK;

        auto it = m_TextureIndex.find(texture.GetName());
        if (it == m_TextureIndex.end())
        {
            auto image = texture.GetTextureImage();

            // Describe and create a Texture2D
            D3D12_HEAP_PROPERTIES prop = {};
            prop.Type = D3D12_HEAP_TYPE_DEFAULT;
            prop.CPUPageProperty = D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
            prop.MemoryPoolPreference = D3D12_MEMORY_POOL_UNKNOWN;
            prop.CreationNodeMask = 1;
            prop.VisibleNodeMask = 1;

            D3D12_RESOURCE_DESC textureDesc = {};
            textureDesc.MipLevels = 1;
            textureDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
            textureDesc.Width = image.Width;
            textureDesc.Height = image.Height;
            textureDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
            textureDesc.DepthOrArraySize = 1;
            textureDesc.SampleDesc.Count = 1;
            textureDesc.SampleDesc.Quality = 0;
            textureDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;

            ID3D12Resource* pTextureBuffer;
            ID3D12Resource* pTextureUploadHeap;

            if(FAILED(hr = m_pDev->CreateCommittedResource (
                &prop,
                D3D12_HEAP_FLAG_NONE,
                &textureDesc,
                D3D12_RESOURCE_STATE_COPY_DEST,
                nullptr,
                IID_PPV_ARGS(&pTextureBuffer)
            )))
            {
                return hr;
            }

            const UINT subresourceCount = textureDesc.DepthOrArraySize * textureDesc.MipLevels;
            const UINT64 uploadBufferSize = GetRequiredIntermediateSize(pTextureBuffer, 0, subresourceCount);

            prop.Type = D3D12_HEAP_TYPE_UPLOAD;

            D3D12_RESOURCE_DESC resourceDesc = {};
            resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
            resourceDesc.Alignment = 0;
            resourceDesc.Width = uploadBufferSize;
            resourceDesc.Height = 1;
            resourceDesc.DepthOrArraySize = 1;
            resourceDesc.MipLevels = 1;
            resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
            resourceDesc.SampleDesc.Count = 1;
            resourceDesc.SampleDesc.Quality = 0;
            resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
            resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

            if (FAILED(hr = m_pDev->CreateCommittedResource(
                &prop,
                D3D12_HEAP_FLAG_NONE,
                &resourceDesc,
                D3D12_RESOURCE_STATE_GENERIC_READ,
                nullptr,
                IID_PPV_ARGS(&pTextureUploadHeap)
            )))
            {
                return hr;
            }

            // Copy data to the intermediate upload heap and then schedule a copy
            // from the upload heap to the Texture2D.
            D3D12_SUBRESOURCE_DATA textureData = {};
            if (image.BitCount == 24)
            {
                // DXGI dose not have 24 bit format so we have to extend it to 32bit
                uint32_t newPitch = image.Pitch / 3 * 4;
                size_t dataSize = newPitch * image.Height;
                void* data = g_pMemoryManager->Allocate(dataSize);
                uint8_t* buf = reinterpret_cast<uint8_t*>(data);
                uint8_t* src = reinterpret_cast<uint8_t*>(image.Data);
                for (size_t row = 0; row < image.Height; row++)
                {
                    buf = reinterpret_cast<uint8_t*>(data) + row * newPitch;
                    src = reinterpret_cast<uint8_t*>(image.Data) + row * image.Pitch;
                    for (size_t col = 0; col < image.Width; ++col)
                    {
                        *(uint32_t*)buf = *(uint32_t*)src;
                        buf[3] = 0;
                        buf += 4;
                        src += 3;
                    }
                }

                // we do not need to free the old data because the old data is still referenced by the SceneObject
                // g_pMemoryManager->Free(image.Data, image.DataSize);
                image.Data = data;
                image.DataSize = dataSize;
                image.Pitch = newPitch;
            }
            textureData.pData = image.Data;
            textureData.RowPitch = image.Pitch;
            textureData.SlicePitch = image.Pitch * image.Height;

            UpdateSubresources(m_pCommandList, pTextureBuffer, pTextureUploadHeap, 0, 0, subresourceCount, &textureData);
            D3D12_RESOURCE_BARRIER barrier = {};
            barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
            barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
            barrier.Transition.pResource = pTextureBuffer;
            barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_DEST;
            barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_GENERIC_READ;
            barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
            m_pCommandList->ResourceBarrier(1, &barrier);

            // Describe and create a SRV for the texture
            D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
            srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
            srvDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
            srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
            srvDesc.Texture2D.MipLevels = -1;
            srvDesc.Texture2D.MostDetailedMip = 0;
            D3D12_CPU_DESCRIPTOR_HANDLE srvHandle;
            size_t textureID = static_cast<uint32_t>(m_TextureIndex.size());
            srvHandle.ptr = m_pCbvHeap->GetCPUDescriptorHandleForHeapStart().ptr + (k_TextureDescStartIndex + textureID) * m_CbvSrvDescriptorSize;
            m_pDev->CreateShaderResourceView(pTextureBuffer, &srvDesc, srvHandle);
            m_TextureIndex[texture.GetName()] = textureID;

            m_Buffers.push_back(pTextureUploadHeap);
            m_Textures.push_back(pTextureBuffer);
        }


        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateSamplerBuffer()
    {
        // Describe and create a sampler.
        D3D12_SAMPLER_DESC samplerDesc = {};
        samplerDesc.Filter = D3D12_FILTER_MIN_MAG_MIP_LINEAR;
        samplerDesc.AddressU = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
        samplerDesc.AddressV = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
        samplerDesc.AddressW = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
        samplerDesc.MinLOD = 0;
        samplerDesc.MaxLOD = D3D12_FLOAT32_MAX;
        samplerDesc.MipLODBias = 0.0f;
        samplerDesc.MaxAnisotropy = 1;
        samplerDesc.ComparisonFunc = D3D12_COMPARISON_FUNC_ALWAYS;
        m_pDev->CreateSampler(&samplerDesc, m_pSamplerHeap->GetCPUDescriptorHandleForHeapStart());

        return S_OK;
    }

    HRESULT D3D12GraphicsManager::CreateConstantBuffer()
    {
        HRESULT hr;

		k_SizePerFrameConstantBuffer = ALIGN(sizeof(DrawFrameContext) + sizeof(LightModel) * k_MaxLightCount, 256);
		k_SizePerBatchConstantBuffer = ALIGN(sizeof(DrawBatchContext), 256);
		k_SizeConstantBufferPerFrame = k_SizePerFrameConstantBuffer + k_SizePerBatchConstantBuffer * k_MaxSceneObjectCount;

        D3D12_HEAP_PROPERTIES prop = {
            D3D12_HEAP_TYPE_UPLOAD,
            D3D12_CPU_PAGE_PROPERTY_UNKNOWN,
            D3D12_MEMORY_POOL_UNKNOWN,
            1,
            1
        };

        D3D12_RESOURCE_DESC resourceDesc = {};
        resourceDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
        resourceDesc.Alignment = 0;
        resourceDesc.Width = k_SizeConstantBufferPerFrame * k_FrameCount;
        resourceDesc.Height = 1;
        resourceDesc.DepthOrArraySize = 1;
        resourceDesc.MipLevels = 1;
        resourceDesc.Format = DXGI_FORMAT_UNKNOWN;
        resourceDesc.SampleDesc.Count = 1;
        resourceDesc.SampleDesc.Quality = 0;
        resourceDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
        resourceDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

        ID3D12Resource* pConstantUploadBuffer;
        if (FAILED(hr = m_pDev->CreateCommittedResource(
            &prop,
            D3D12_HEAP_FLAG_NONE,
            &resourceDesc,
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(&pConstantUploadBuffer)
        )))
        {
            return hr;
        }

        // populate descriptor table
        D3D12_CPU_DESCRIPTOR_HANDLE cbvHandle;
        cbvHandle.ptr = m_pCbvHeap->GetCPUDescriptorHandleForHeapStart().ptr;
        for (uint32_t i = 0; i < k_FrameCount; ++i)
        {
            for (uint32_t j = 0; j < k_MaxSceneObjectCount; ++j)
            {
				// Describe and create constant buffer descriptors.
				// 1 per frame and 1 per batch descriptor per object.
				D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc = {};

				// Per frame const buffer descriptors.
				cbvDesc.BufferLocation = pConstantUploadBuffer->GetGPUVirtualAddress() +
					i * k_SizeConstantBufferPerFrame;
				cbvDesc.SizeInBytes = k_SizePerFrameConstantBuffer;
				m_pDev->CreateConstantBufferView(&cbvDesc, cbvHandle);
				cbvHandle.ptr += m_CbvSrvDescriptorSize;

                // Per batch constant buffer descriptors
                cbvDesc.BufferLocation += k_SizePerFrameConstantBuffer + j * k_SizePerBatchConstantBuffer;
                cbvDesc.SizeInBytes = k_SizePerBatchConstantBuffer;
                m_pDev->CreateConstantBufferView(&cbvDesc, cbvHandle);
                cbvHandle.ptr += m_CbvSrvDescriptorSize;
            }
        }

        D3D12_RANGE readRange = {0, 0};
        hr = pConstantUploadBuffer->Map(0, &readRange, reinterpret_cast<void**>(&m_pCbvDataBegin));

        m_Buffers.push_back(pConstantUploadBuffer);

        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateGraphicsResources()
    {
        HRESULT hr;

    #if defined(_DEBUG)
        // Enable the D3D12 debug layer.
        {
            ID3D12Debug* pDebugController;
            if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&pDebugController))))
            {
                pDebugController->EnableDebugLayer();
            }
            SafeRelease(&pDebugController);
        }
    #endif

        IDXGIFactory4* pFactory;
        if (FAILED(hr = CreateDXGIFactory1(IID_PPV_ARGS(&pFactory)))) {
            return hr;
        }

        IDXGIAdapter1* pHardwareAdapter;
        GetHardwareAdapter(pFactory, &pHardwareAdapter);

        if (FAILED(D3D12CreateDevice(pHardwareAdapter,
            D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&m_pDev)))) {

            IDXGIAdapter* pWarpAdapter;
            if (FAILED(hr = pFactory->EnumWarpAdapter(IID_PPV_ARGS(&pWarpAdapter)))) {
                SafeRelease(&pFactory);
                return hr;
            }

            if(FAILED(hr = D3D12CreateDevice(pWarpAdapter, D3D_FEATURE_LEVEL_11_0,
                IID_PPV_ARGS(&m_pDev)))) {
                SafeRelease(&pFactory);
                return hr;
            }
        }


        HWND hWnd = reinterpret_cast<WindowsApplication*>(g_pApp)->GetMainWindow();

        // Describe and create the command queue.
        D3D12_COMMAND_QUEUE_DESC queueDesc = {};
        queueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
        queueDesc.Type  = D3D12_COMMAND_LIST_TYPE_DIRECT;

        if(FAILED(hr = m_pDev->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&m_pCommandQueue)))) {
            SafeRelease(&pFactory);
            return hr;
        }

        // create a struct to hold information about the swap chain
        DXGI_SWAP_CHAIN_DESC1 scd;

        // clear out the struct for use
        ZeroMemory(&scd, sizeof(DXGI_SWAP_CHAIN_DESC1));

        // fill the swap chain description struct
        scd.Width  = g_pApp->GetConfiguration().screenWidth;
        scd.Height = g_pApp->GetConfiguration().screenHeight;
        scd.Format = DXGI_FORMAT_R8G8B8A8_UNORM;     	        // use 32-bit color
        scd.Stereo = FALSE;
        scd.SampleDesc.Count = 1;                               // multi-samples can not be used when in SwapEffect sets to
                                                                // DXGI_SWAP_EFFECT_FLOP_DISCARD
        scd.SampleDesc.Quality = 0;                             // multi-samples can not be used when in SwapEffect sets to
        scd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;      // how swap chain is to be used
        scd.BufferCount = k_FrameCount;                          // back buffer count
        scd.Scaling     = DXGI_SCALING_STRETCH;
        scd.SwapEffect  = DXGI_SWAP_EFFECT_FLIP_DISCARD;        // DXGI_SWAP_EFFECT_FLIP_DISCARD only supported after Win10
                                                                // use DXGI_SWAP_EFFECT_DISCARD on platforms early than Win10
        scd.AlphaMode   = DXGI_ALPHA_MODE_UNSPECIFIED;
        scd.Flags    = DXGI_SWAP_CHAIN_FLAG_ALLOW_MODE_SWITCH;  // allow full-screen transition

        IDXGISwapChain1* pSwapChain;
        if (FAILED(hr = pFactory->CreateSwapChainForHwnd(
                    m_pCommandQueue,                            // Swap chain needs the queue so that it can force a flush on it
                    hWnd,
                    &scd,
                    NULL,
                    NULL,
                    &pSwapChain
                    )))
        {
            SafeRelease(&pFactory);
            return hr;
        }

        SafeRelease(&pFactory);

        m_pSwapChain = reinterpret_cast<IDXGISwapChain3*>(pSwapChain);

        m_FrameIndex = m_pSwapChain->GetCurrentBackBufferIndex();
        
        std::cout << "Creating Descriptor Heaps...";
        if (FAILED(hr = CreateDescriptorHeaps()))
        {
            return hr;
        }
        std::cout << "Done!" << std::endl;

        std::cout << "Creating Render Targets...";
        if (FAILED(hr = CreateRenderTarget()))
        {
            return hr;
        }
        std::cout << "Done!" << std::endl;

		std::cout << "Creating Depth Stencil Buffer ...";
        if (FAILED(hr = CreateDepthStencil())) {
            return hr;
        }
		std::cout << "Done!" << std::endl;

		std::cout << "Creating Root Signatures ...";
        if (FAILED(hr = CreateRootSignature())) {
            return hr;
        }
		std::cout << "Done!" << std::endl;

        std::cout << "Creating PSO ...";
        if (FAILED (hr = InitializePSO()))
            return hr;
        std::cout << "Done!" << std::endl;

        return hr;
    }

    HRESULT D3D12GraphicsManager::CreateRootSignature()
    {
        HRESULT hr = S_OK;

        D3D12_FEATURE_DATA_ROOT_SIGNATURE featureData = {};

        // This is the highest version the sample supports. If CheckFeatureSupport succeeds, the HighestVersion returned will not be greater than this.
        featureData.HighestVersion = D3D_ROOT_SIGNATURE_VERSION_1_1;

        if (FAILED(m_pDev->CheckFeatureSupport(D3D12_FEATURE_ROOT_SIGNATURE, &featureData, sizeof(featureData))))
        {
            featureData.HighestVersion = D3D_ROOT_SIGNATURE_VERSION_1_0;
        }

		// root signature for base pass
		{
			D3D12_DESCRIPTOR_RANGE1 ranges[3] = {
				{D3D12_DESCRIPTOR_RANGE_TYPE_CBV, 2, 0, 0, D3D12_DESCRIPTOR_RANGE_FLAG_DATA_STATIC},
				{D3D12_DESCRIPTOR_RANGE_TYPE_SAMPLER, 1, 0},
				{D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 1, 0, 0, D3D12_DESCRIPTOR_RANGE_FLAG_DATA_STATIC}
			};

			D3D12_ROOT_PARAMETER1 rootParameters[4] = {
				{D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE, {1, &ranges[0]}, D3D12_SHADER_VISIBILITY_ALL},
				{D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE, {1, &ranges[1]}, D3D12_SHADER_VISIBILITY_PIXEL},
				{D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE, {1, &ranges[2]}, D3D12_SHADER_VISIBILITY_PIXEL}
			};
            rootParameters[3].ParameterType = D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS;
            rootParameters[3].Constants = {2, 0, 1};
            rootParameters[3].ShaderVisibility = D3D12_SHADER_VISIBILITY_PIXEL;

			// Allow input layout and deny unecessary access to certain pipeline stages.
			D3D12_ROOT_SIGNATURE_FLAGS rootSignatureFlags =
				D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT |
				D3D12_ROOT_SIGNATURE_FLAG_DENY_HULL_SHADER_ROOT_ACCESS |
				D3D12_ROOT_SIGNATURE_FLAG_DENY_DOMAIN_SHADER_ROOT_ACCESS |
				D3D12_ROOT_SIGNATURE_FLAG_DENY_GEOMETRY_SHADER_ROOT_ACCESS;

			D3D12_ROOT_SIGNATURE_DESC1 rootSignatureDesc = {
				_countof(rootParameters), rootParameters, 0, nullptr, rootSignatureFlags
			};

			D3D12_VERSIONED_ROOT_SIGNATURE_DESC versionedRootSignatureDesc = {
				D3D_ROOT_SIGNATURE_VERSION_1_1
			};

			versionedRootSignatureDesc.Desc_1_1 = rootSignatureDesc;

			ID3DBlob* signature = nullptr;
			ID3DBlob* error = nullptr;
			if (SUCCEEDED(hr = D3D12SerializeVersionedRootSignature(&versionedRootSignatureDesc, &signature, &error)))
			{
				hr = m_pDev->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&m_pRootSignature));
			}

			SafeRelease(&signature);
			SafeRelease(&error);
		}

#if (!MSAA_SOLUTION_DEFAULT)
		// root signature for resolve pass
		{
			D3D12_DESCRIPTOR_RANGE1 ranges[1] = {
				{D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 1, 0, 0, D3D12_DESCRIPTOR_RANGE_FLAG_DATA_STATIC}
			};

			D3D12_ROOT_PARAMETER1 rootParameters[2];
			rootParameters[0] = { D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE, {1, &ranges[0]}, D3D12_SHADER_VISIBILITY_PIXEL };
			rootParameters[1].ParameterType = D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS;
			rootParameters[1].Constants = { 0, 0, 2 };
			rootParameters[1].ShaderVisibility = D3D12_SHADER_VISIBILITY_PIXEL;

			// Allow input layout and deny unnecessary access to certain pipeline stages.
			D3D12_ROOT_SIGNATURE_FLAGS rootSignatureFlags =
				D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT |
				D3D12_ROOT_SIGNATURE_FLAG_DENY_HULL_SHADER_ROOT_ACCESS |
				D3D12_ROOT_SIGNATURE_FLAG_DENY_DOMAIN_SHADER_ROOT_ACCESS |
				D3D12_ROOT_SIGNATURE_FLAG_DENY_GEOMETRY_SHADER_ROOT_ACCESS;

			D3D12_ROOT_SIGNATURE_DESC1 rootSignatureDesc = {
				_countof(rootParameters), rootParameters, 0, nullptr, rootSignatureFlags
			};

			D3D12_VERSIONED_ROOT_SIGNATURE_DESC versionedRootSignatureDesc = {
				D3D_ROOT_SIGNATURE_VERSION_1_1 };

			versionedRootSignatureDesc.Desc_1_1 = rootSignatureDesc;

			ID3DBlob* signature = nullptr;
			ID3DBlob* error = nullptr;
			if (SUCCEEDED(hr = D3D12SerializeVersionedRootSignature(&versionedRootSignatureDesc, &signature, &error)))
			{
				hr = m_pDev->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&m_pRootSignatureResolve));
				if (FAILED(hr))
				{
					printf("failed to create root signature for resolve phase: %x", hr);
				}
			}

			SafeRelease(&signature);
			SafeRelease(&error);
		}
#endif

        return hr;
    }

    // This is the function that loads and prepares the shaders.
    HRESULT D3D12GraphicsManager::InitializePSO()
    {
        HRESULT hr = S_OK;

		// basic pass
		{
			ID3DBlob* pVS = CompileShader(L"../../Asset/Shaders/basic_lighting.hlsl", nullptr, "VSMain", "vs_5_0");
			ID3DBlob* pPS = CompileShader(L"../../Asset/Shaders/basic_lighting.hlsl", nullptr, "PSMain", "ps_5_0");

			D3D12_SHADER_BYTECODE vertexShaderByteCode;
			D3D12_SHADER_BYTECODE pixelShaderByteCode;

			vertexShaderByteCode.pShaderBytecode = reinterpret_cast<BYTE*>(pVS->GetBufferPointer());
			vertexShaderByteCode.BytecodeLength = pVS->GetBufferSize();

			pixelShaderByteCode.pShaderBytecode = reinterpret_cast<BYTE*>(pPS->GetBufferPointer());
			pixelShaderByteCode.BytecodeLength = pPS->GetBufferSize();

			//const char* vsFilename = "Shaders/basic_vs.cso";
			//const char* fsFilename = "Shaders/basic_ps.cso";

			//Buffer vertexShader = g_pAssetLoader->SyncOpenAndReadBinary(vsFilename);
			//Buffer pixelShader = g_pAssetLoader->SyncOpenAndReadBinary(fsFilename);
			
			//vertexShaderByteCode.pShaderBytecode = vertexShader.GetData();
			//vertexShaderByteCode.BytecodeLength = vertexShader.GetDataSize();
			//
			//pixelShaderByteCode.pShaderBytecode = pixelShader.GetData();
			//pixelShaderByteCode.BytecodeLength = pixelShader.GetDataSize();

			// create the input layout object
			D3D12_INPUT_ELEMENT_DESC ied[] =
			{
				{"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
				{"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 1, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
				{"TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 2, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0}
			};

			D3D12_RASTERIZER_DESC rsd = {
				D3D12_FILL_MODE_SOLID, D3D12_CULL_MODE_BACK, TRUE, D3D12_DEFAULT_DEPTH_BIAS, D3D12_DEFAULT_DEPTH_BIAS_CLAMP,
				TRUE, FALSE, FALSE, 0, D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF
			};
			const D3D12_RENDER_TARGET_BLEND_DESC defaultRenderTargetBlend = {
				FALSE, FALSE,
				D3D12_BLEND_ONE, D3D12_BLEND_ZERO, D3D12_BLEND_OP_ADD,
				D3D12_BLEND_ONE, D3D12_BLEND_ZERO, D3D12_BLEND_OP_ADD,
				D3D12_LOGIC_OP_NOOP,
				D3D12_COLOR_WRITE_ENABLE_ALL
			};

			D3D12_BLEND_DESC bld = {
				FALSE, FALSE,
				{
					defaultRenderTargetBlend,
					defaultRenderTargetBlend,
					defaultRenderTargetBlend,
					defaultRenderTargetBlend,
					defaultRenderTargetBlend,
					defaultRenderTargetBlend,
					defaultRenderTargetBlend
				}
			};

			const D3D12_DEPTH_STENCILOP_DESC defaultStencilOp = {
				D3D12_STENCIL_OP_KEEP,
				D3D12_STENCIL_OP_KEEP,
				D3D12_STENCIL_OP_KEEP,
				D3D12_COMPARISON_FUNC_ALWAYS };
			D3D12_DEPTH_STENCIL_DESC dsd = {
				TRUE,
				D3D12_DEPTH_WRITE_MASK_ALL,
				D3D12_COMPARISON_FUNC_LESS,
				FALSE,
				D3D12_DEFAULT_STENCIL_READ_MASK,
				D3D12_DEFAULT_STENCIL_WRITE_MASK,
				defaultStencilOp, defaultStencilOp
			};


			// describe and create the graphics pipeline state object(PSO)
			D3D12_GRAPHICS_PIPELINE_STATE_DESC psod = {};
			psod.pRootSignature = m_pRootSignature;
			psod.VS = vertexShaderByteCode;
			psod.PS = pixelShaderByteCode;
			psod.BlendState = bld;
			psod.SampleMask = UINT_MAX;
			psod.RasterizerState = rsd;
			psod.DepthStencilState = dsd;
			psod.InputLayout = { ied, _countof(ied) };
			psod.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
			psod.NumRenderTargets = 1;
			psod.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
			psod.DSVFormat = DXGI_FORMAT_D32_FLOAT;
			psod.SampleDesc.Count = 4;  // 4X MSAA
			psod.SampleDesc.Quality = DXGI_STANDARD_MULTISAMPLE_QUALITY_PATTERN;
			//psod.SampleDesc.Count = 1;  // 4X MSAA
			//psod.SampleDesc.Quality = 0;

			if (FAILED(hr = m_pDev->CreateGraphicsPipelineState(&psod, IID_PPV_ARGS(&m_pPipelineState))))
			{
				return hr;
			}
		}

#if (!MSAA_SOLUTION_DEFAULT)
		// resolve pass
		{
			const char* vsFilename = "Shaders/msaa_resolver_vs.cso";
			const char* fsFilename = "Shaders/msaa_resolver_ps.cso";

			// load the shaders
			Buffer vertexShader = g_pAssetLoader->SyncOpenAndReadBinary(vsFilename);
			Buffer pixelShader = g_pAssetLoader->SyncOpenAndReadBinary(fsFilename);

			D3D12_SHADER_BYTECODE vertexShaderByteCode;
			vertexShaderByteCode.pShaderBytecode = vertexShader.GetData();
			vertexShaderByteCode.BytecodeLength = vertexShader.GetDataSize();

			D3D12_SHADER_BYTECODE pixelShaderByteCode;
			pixelShaderByteCode.pShaderBytecode = pixelShader.GetData();
			pixelShaderByteCode.BytecodeLength = pixelShader.GetDataSize();

			// create the input layout object
			D3D12_INPUT_ELEMENT_DESC ied[] = {
				{"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0},
				{"TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0}
			};

			D3D12_RASTERIZER_DESC rsd = {
				D3D12_FILL_MODE_SOLID, D3D12_CULL_MODE_BACK, TRUE,
				D3D12_DEFAULT_DEPTH_BIAS, D3D12_DEFAULT_DEPTH_BIAS_CLAMP,
				TRUE, FALSE, FALSE, 0, D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF };

			const D3D12_RENDER_TARGET_BLEND_DESC defaultRenderTargetBlend = { FALSE, FALSE,
				D3D12_BLEND_ONE, D3D12_BLEND_ZERO, D3D12_BLEND_OP_ADD,
				D3D12_BLEND_ONE, D3D12_BLEND_ZERO, D3D12_BLEND_OP_ADD,
				D3D12_LOGIC_OP_NOOP,
				D3D12_COLOR_WRITE_ENABLE_ALL };

			D3D12_BLEND_DESC bld = { FALSE, FALSE,
			{
				defaultRenderTargetBlend,
				defaultRenderTargetBlend,
				defaultRenderTargetBlend,
				defaultRenderTargetBlend,
				defaultRenderTargetBlend,
				defaultRenderTargetBlend,
				defaultRenderTargetBlend
				}
			};

			const D3D12_DEPTH_STENCILOP_DESC defaultStencilOp = {
				D3D12_STENCIL_OP_KEEP,
				D3D12_STENCIL_OP_KEEP,
				D3D12_STENCIL_OP_KEEP,
				D3D12_COMPARISON_FUNC_ALWAYS };

			D3D12_DEPTH_STENCIL_DESC dsd = {
				TRUE,
				D3D12_DEPTH_WRITE_MASK_ALL,
				D3D12_COMPARISON_FUNC_LESS,
				FALSE,
				D3D12_DEFAULT_STENCIL_READ_MASK,
				D3D12_DEFAULT_STENCIL_WRITE_MASK,
				defaultStencilOp, defaultStencilOp };

			D3D12_GRAPHICS_PIPELINE_STATE_DESC psod = {};
			psod.pRootSignature = m_pRootSignatureResolve;
			psod.VS = vertexShaderByteCode;
			psod.PS = pixelShaderByteCode;
			psod.BlendState = bld;
			psod.SampleMask = UINT_MAX;
			psod.RasterizerState = rsd;
			psod.DepthStencilState = dsd;
			psod.InputLayout = { ied, _countof(ied) };
			psod.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
			psod.NumRenderTargets = 1;
			psod.RTVFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;
			psod.DSVFormat = DXGI_FORMAT_UNKNOWN;
			psod.SampleDesc.Count = 1; // no msaa
			psod.SampleDesc.Quality = 0; // no msaa

			if (FAILED(hr = m_pDev->CreateGraphicsPipelineState(&psod, IID_PPV_ARGS(&m_pPipelineStateResolve))))
			{
				return false;
			}
		}
#endif

        if (!m_pCommandList)
        {
            if (FAILED(hr = m_pDev->CreateCommandList(0,
                D3D12_COMMAND_LIST_TYPE_DIRECT,
                m_pCommandAllocator,
                m_pPipelineState,
                IID_PPV_ARGS(&m_pCommandList))))
            {
                return hr;
            }            
        }


        return hr;
    }

    void D3D12GraphicsManager::InitializeBuffers(const Scene& scene)
    {
        HRESULT hr;

		std::cout << "Creating vertex buffer..." << std::endl;
        for (auto _it : scene.GeometryNodes)
        {
			auto pGeometryNode = _it.second;
            if (pGeometryNode && pGeometryNode->Visible())
            {
                auto pGeometry = scene.GetGeometry(pGeometryNode->GetSceneObjectRef());
                assert(pGeometry);
                auto pMesh = pGeometry->GetMesh().lock();
                if (!pMesh) continue;

                // Set the number of vertex properties.
                size_t vertexPropertiesCount = pMesh->GetVertexPropertiesCount();

                // Set the number of vertices in the vertex array.
                auto vertexCount = pMesh->GetVertexCount();

                Buffer buff;

                for (size_t i = 0; i < vertexPropertiesCount; ++i)
                {
                    const SceneObjectVertexArray& vPropertyArray = pMesh->GetVertexPropertyArray(i);
                    CreateVertexBuffer(vPropertyArray);
                }

                // TODO: Implement LOD switching
                // auto indexGroupCount = pMesh->GetIndexGroupCount();
                const SceneObjectIndexArray& indexArray = pMesh->GetIndexArray(0);
                CreateIndexBuffer(indexArray);

                auto materialIndex = indexArray.GetMaterialIndex();
                auto materialKey = pGeometryNode->GetMaterialRef(materialIndex);
                auto material = scene.GetMaterial(materialKey);

                DrawBatchContext dbc;
                dbc.indexCount = (UINT)indexArray.GetIndexCount();
                dbc.propertyCount = vertexPropertiesCount;
                if (material)
                {
                    dbc.material = material;
                }

                dbc.node = pGeometryNode;

                m_DrawBatchContext.push_back(dbc);
            }
        }
		std::cout << "Done!" << std::endl;

#if (!MSAA_SOLUTION_DEFAULT)
		std::cout << "Creating Internal Vertex Buffer..." << std::endl;
		CreateIntervalVertexBuffer();
		std::cout << "Done! " << std::endl;

#endif
		std::cout << "Creating Constant Buffer ...";
		if (FAILED(hr = CreateConstantBuffer())) {
			return;
		}
		std::cout << "Done!" << std::endl;

		std::cout << "Creating Sampler Buffer ...";
		if (FAILED(hr = CreateSamplerBuffer())) {
			return;
		}
		std::cout << "Done!" << std::endl;


		for (auto _it : scene.Materials)
		{
			auto material = _it.second;
			if (material)
			{
				auto color = material->GetBaseColor();
				if (auto texture = color.ValueMap)
				{
					if (FAILED(hr = CreateTextureBuffer(*texture)))
					{
						return;
					}
				}
			}
		}

        if (SUCCEEDED(hr = m_pCommandList->Close()))
        {
            ID3D12CommandList* ppCommandLists[] = {m_pCommandList};
            m_pCommandQueue->ExecuteCommandLists(_countof(ppCommandLists), ppCommandLists);

            if (FAILED(hr = m_pDev->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&m_pFence))))
            {
                return ;
            }

            m_FenceValue = 1;

            m_FenceEvent = CreateEvent(NULL, FALSE, FALSE, NULL);
            if (m_FenceEvent == NULL)
            {
                hr = HRESULT_FROM_WIN32(GetLastError());
                if (FAILED(hr))
                    return ;
            }

            WaitForPreviousFrame();
        }

        hr = PopulateCommandList();

        return ;
    }

    int  D3D12GraphicsManager::Initialize()
    {
        int result = GraphicsManager::Initialize();

        if (!result)
        {
            const GfxConfiguration& config = g_pApp->GetConfiguration();
            m_ViewPort = {0.0f, 0.0f, static_cast<float>(config.screenWidth), static_cast<float>(config.screenHeight), 0.0f, 1.0f};
            m_ScissorRect = {0, 0, static_cast<LONG>(config.screenWidth), static_cast<LONG>(config.screenHeight)};
            result = static_cast<int>(CreateGraphicsResources());
        }

        return result;
    }

    void D3D12GraphicsManager::ClearBuffers()
    {
        SafeRelease(&m_pFence);
        for (auto p : m_Buffers)
        {
            SafeRelease(&p);
        }
		m_Buffers.clear();
        for(auto p : m_Textures)
        {
            SafeRelease(&p);
        }
        m_Textures.clear();
        m_TextureIndex.clear();
        m_VertexBufferView.clear();
        m_IndexBufferView.clear();
        m_DrawBatchContext.clear();
    }

    void D3D12GraphicsManager::Finalize()
    {
        GraphicsManager::Finalize();

        SafeRelease(&m_pRtvHeap);
        SafeRelease(&m_pDsvHeap);
        SafeRelease(&m_pCbvHeap);
        SafeRelease(&m_pSamplerHeap);
        SafeRelease(&m_pRootSignature);
		SafeRelease(&m_pRootSignatureResolve);
        SafeRelease(&m_pCommandQueue);
        SafeRelease(&m_pCommandAllocator);
        SafeRelease(&m_pDepthStencilBuffer);
		SafeRelease(&m_pMsaaRenderTarget);
        for (uint32_t i = 0; i < k_FrameCount; i++) {
            SafeRelease(&m_pRenderTargets[i]);
        }
        SafeRelease(&m_pSwapChain);
        SafeRelease(&m_pDev);
    }

    void D3D12GraphicsManager::Clear()
    {
        GraphicsManager::Clear();
    }

    void D3D12GraphicsManager::Draw()
    {
        PopulateCommandList();

		UpdateConstants();
		RenderBuffers();
        //GraphicsManager::Draw();

        WaitForPreviousFrame();
    }

    HRESULT D3D12GraphicsManager::PopulateCommandList()
    {
        HRESULT hr;

        // command list allocators can only be reset when the associated
        // command lists have finished execution on GPU; apps should use
        // fences to determine GPU execution progress.
        if (FAILED(hr = m_pCommandAllocator->Reset()))
        {
            return hr;
        }


			// however, when ExecuteCommandList() is called on a particular command
			// list, that command list can then be reset at any time and must be before
			// re-recording
			if (FAILED(hr = m_pCommandList->Reset(m_pCommandAllocator, m_pPipelineState)))
			{
				return hr;
			}

			// Indicate that the back buffer will be used as a render target
			D3D12_RESOURCE_BARRIER barrier = {};
			barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
			barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
			barrier.Transition.pResource = m_pMsaaRenderTarget;
#if (!MSAA_SOLUTION_DEFAULT)
			barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE;
#else
			barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
#endif
			barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;
			barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
			m_pCommandList->ResourceBarrier(1, &barrier);

			D3D12_CPU_DESCRIPTOR_HANDLE rtvHandle;
			//rtvHandle.ptr = m_pRtvHeap->GetCPUDescriptorHandleForHeapStart().ptr + m_FrameIndex * m_RtvDescriptorSize;
			rtvHandle.ptr = m_pRtvHeap->GetCPUDescriptorHandleForHeapStart().ptr + k_FrameCount * m_RtvDescriptorSize;
			D3D12_CPU_DESCRIPTOR_HANDLE dsvHandle;
			dsvHandle = m_pDsvHeap->GetCPUDescriptorHandleForHeapStart();
			m_pCommandList->OMSetRenderTargets(1, &rtvHandle, FALSE, &dsvHandle);

			// clear the back buffer to a deep blue
			const FLOAT clearColor[] = { 0.2f, 0.3f, 0.4f, 1.0f };
			m_pCommandList->ClearRenderTargetView(rtvHandle, clearColor, 0, nullptr);
			m_pCommandList->ClearDepthStencilView(dsvHandle, D3D12_CLEAR_FLAG_DEPTH, 1.0f, 0, 0, nullptr);

			// Set necessary state.
			m_pCommandList->SetGraphicsRootSignature(m_pRootSignature);

			ID3D12DescriptorHeap* ppHeaps[] = { m_pCbvHeap, m_pSamplerHeap };
			m_pCommandList->SetDescriptorHeaps(_countof(ppHeaps), ppHeaps);

			// Sampler
			m_pCommandList->SetGraphicsRootDescriptorTable(1, m_pSamplerHeap->GetGPUDescriptorHandleForHeapStart());

			m_pCommandList->RSSetViewports(1, &m_ViewPort);
			m_pCommandList->RSSetScissorRects(1, &m_ScissorRect);
			m_pCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST); // tell D3D our vertices can assemble to triangle list

            // Set num of lights
            auto& DrawFrameContext = m_Frames[0].frameContext;
            m_pCommandList->SetGraphicsRoot32BitConstant(3, static_cast<UINT>(DrawFrameContext.Lights.size()), 0);

			// do 3D rendering on the back buffer here
			int32_t i = 0;
			size_t vertexBufferViewOffset = 0;
			for (auto dbc : m_DrawBatchContext)
			{
				// CBV Per Batch
				D3D12_GPU_DESCRIPTOR_HANDLE cbvSrvHandle;
				uint32_t frameResourceDescriptorOffset = m_FrameIndex * (2 * k_MaxSceneObjectCount); // 2 descriptors for each draw call
				cbvSrvHandle.ptr = m_pCbvHeap->GetGPUDescriptorHandleForHeapStart().ptr +
					(frameResourceDescriptorOffset + i * 2) * m_CbvSrvDescriptorSize;
				m_pCommandList->SetGraphicsRootDescriptorTable(0, cbvSrvHandle);

				// select which vertex buffer(s) to use
				for (uint32_t j = 0; j < dbc.propertyCount; j++)
				{
					m_pCommandList->IASetVertexBuffers(j, 1, &m_VertexBufferView[vertexBufferViewOffset++]);
				}
				// select which index buffer to use
				m_pCommandList->IASetIndexBuffer(&m_IndexBufferView[i]);

				// Texture
				if (dbc.material)
				{
					if (auto texture = dbc.material->GetBaseColor().ValueMap)
					{
						auto textureIndex = m_TextureIndex[texture->GetName()];
						D3D12_GPU_DESCRIPTOR_HANDLE srvHandle;
						srvHandle.ptr = m_pCbvHeap->GetGPUDescriptorHandleForHeapStart().ptr + (k_TextureDescStartIndex + textureIndex) * m_CbvSrvDescriptorSize;
						m_pCommandList->SetGraphicsRootDescriptorTable(2, srvHandle);
					}
				}

				// draw the vertex buffer to the back buffer
				m_pCommandList->DrawIndexedInstanced(dbc.indexCount, 1, 0, 0, 0);
				i++;
			}
#if MSAA_SOLUTION_DEFAULT
		barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		barrier.Transition.pResource = m_pMsaaRenderTarget;
		barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
		barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RESOLVE_SOURCE;
		barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		m_pCommandList->ResourceBarrier(1, &barrier);

		barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		barrier.Transition.pResource = m_pRenderTargets[m_FrameIndex];
		barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
		barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RESOLVE_DEST;
		barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		m_pCommandList->ResourceBarrier(1, &barrier);

		m_pCommandList->ResolveSubresource(m_pRenderTargets[m_FrameIndex], 0, m_pMsaaRenderTarget, 0, DXGI_FORMAT_R8G8B8A8_UNORM);

		barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		barrier.Transition.pResource = m_pMsaaRenderTarget;
		barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RESOLVE_SOURCE;
		barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;
		barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		m_pCommandList->ResourceBarrier(1, &barrier);

		barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		barrier.Transition.pResource = m_pRenderTargets[m_FrameIndex];
		barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RESOLVE_DEST;
		barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;
		barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		m_pCommandList->ResourceBarrier(1, &barrier);
#else
			// msaa resolve pass
			barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
			barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
			barrier.Transition.pResource = m_pMsaaRenderTarget;
			barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
			barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE;
			barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
			m_pCommandList->ResourceBarrier(1, &barrier);

			m_pCommandList->SetPipelineState(m_pPipelineStateResolve);

			// Indicate that the back buffer will be used as a render target.
			barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
			barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
			barrier.Transition.pResource = m_pRenderTargets[m_FrameIndex];
			barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
			barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;
			barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
			m_pCommandList->ResourceBarrier(1, &barrier);

			rtvHandle.ptr = m_pRtvHeap->GetCPUDescriptorHandleForHeapStart().ptr + m_FrameIndex * m_RtvDescriptorSize;
			m_pCommandList->OMSetRenderTargets(1, &rtvHandle, FALSE, nullptr);

			// Set necessary state.
			m_pCommandList->SetGraphicsRootSignature(m_pRootSignatureResolve);

			// Set CBV
			m_pCommandList->SetGraphicsRoot32BitConstants(1, 1, &g_pApp->GetConfiguration().screenWidth, 0);
			m_pCommandList->SetGraphicsRoot32BitConstants(1, 1, &g_pApp->GetConfiguration().screenHeight, 1);

			// Set SRV
			auto textureIndex = m_TextureIndex["MSAA"];
			D3D12_GPU_DESCRIPTOR_HANDLE srvHandle;
			srvHandle.ptr = m_pCbvHeap->GetGPUDescriptorHandleForHeapStart().ptr + (k_TextureDescStartIndex + textureIndex) * m_CbvSrvDescriptorSize;
			m_pCommandList->SetGraphicsRootDescriptorTable(0, srvHandle);

			m_pCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);
			m_pCommandList->IASetVertexBuffers(0, 1, &m_VertexBufferViewResolve);
			m_pCommandList->DrawInstanced(4, 1, 0, 0);

			barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
			barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
			barrier.Transition.pResource = m_pRenderTargets[m_FrameIndex];
			barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
			barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;
			barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
			m_pCommandList->ResourceBarrier(1, &barrier);
		//memset(&barrier, 0x00, sizeof(barrier));
		//barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
		//barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
		//barrier.Transition.pResource = m_pRenderTargets[m_FrameIndex];
		//barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
		//barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;
		//barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
		//m_pCommandList->ResourceBarrier(1, &barrier);
#endif

        hr = m_pCommandList->Close();

        return hr;
    }

    void D3D12GraphicsManager::UpdateConstants()
    {
        GraphicsManager::UpdateConstants();

        // CBV Per Frame
        SetPerFrameShaderParameters();
        int32_t i = 0;
        for (auto dbc : m_DrawBatchContext)
        {
            SetPerBatchShaderParameters(i++);
        }
    }

    VOID D3D12GraphicsManager::RenderBuffers()
    {
        HRESULT hr;

        // execute the command list
        ID3D12CommandList* ppCommandLists[] = {m_pCommandList};
        m_pCommandQueue->ExecuteCommandLists(_countof(ppCommandLists), ppCommandLists);

        // swap the back buffer and the front buffer
        hr = m_pSwapChain->Present(1, 0);

        return;
    }

    bool D3D12GraphicsManager::SetPerFrameShaderParameters()
    {
        //memcpy(m_pCbvDataBegin + m_FrameIndex * k_SizeConstantBufferPerFrame, &m_DrawFrameContext, sizeof(m_DrawFrameContext));
        uint8_t* pHead = m_pCbvDataBegin + m_FrameIndex * k_SizeConstantBufferPerFrame;
        auto& DrawFrameContext = m_Frames[0].frameContext;
        size_t offset = (uint8_t*)&DrawFrameContext.Lights - (uint8_t*)&DrawFrameContext;

        memcpy(pHead,
            &DrawFrameContext,
            offset);
        
		uint32_t t = ALIGN(offset, 16);
        pHead += t;     // 16 bytes alignment

        for (auto light : DrawFrameContext.Lights)
        {
            //size_t size = ALIGN(sizeof(LightContext), 16);
            memcpy(pHead, &light, sizeof(LightModel));
			size_t size = ALIGN(sizeof(LightModel), 16);
            pHead += size;
        }

        return true;
    }

    bool D3D12GraphicsManager::SetPerBatchShaderParameters(int32_t index)
    {
        PerBatchConstants pbc;
        memset(&pbc, 0x00, sizeof(pbc));

        Matrix4f trans;
        trans = *m_DrawBatchContext[index].node->GetCalculatedTransform();

        pbc.objectMatrix = trans;

        if (m_DrawBatchContext[index].material)
        {
            Color color = m_DrawBatchContext[index].material->GetBaseColor();
            if (color.ValueMap)
            {
                pbc.usingDiffuseMap = true;
            }
            else
            {
                pbc.diffuseColor = color.Value;
                pbc.usingDiffuseMap = false;
            }
            
            color = m_DrawBatchContext[index].material->GetSpecularColor();
            if (color.ValueMap)
            {
                pbc.specularColor = Vector4Df(-1.0f);
            }
            else
            {
                pbc.specularColor = color.Value;
            }

            Parameter param = m_DrawBatchContext[index].material->GetSpecularPower();
            pbc.specularPower = param.Value;
        }

        memcpy(m_pCbvDataBegin + m_FrameIndex * k_SizeConstantBufferPerFrame    // offset by frame index
            + k_SizePerFrameConstantBuffer                                      // offset by per frame buffer
            + index * k_SizePerBatchConstantBuffer,                             // offset by object index
            &pbc, sizeof(pbc));

        return true;
    }

	ID3DBlob* D3D12GraphicsManager::CompileShader(const std::wstring& filename,
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
		{
			OutputDebugStringA((char*)errors->GetBufferPointer());
			OutputShaderErrorMessage((char*)errors->GetBufferPointer());
		}

		return byteCode;
	}
}

