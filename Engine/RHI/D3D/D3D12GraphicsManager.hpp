#pragma once
#include <stdint.h>
#include <d3d12.h>
#include <DXGI1_4.h>
#include <vector>
#include <map>
#include <string>
#include "RHI/GraphicsManager.hpp"
#include "Core/Buffer.hpp"
#include "Core/Image.hpp"
#include "SceneManager/SceneObject.hpp"

namespace Panda {
	struct D3dShaderProgram {
		D3D12_SHADER_BYTECODE vertexShaderByteCode;
		D3D12_SHADER_BYTECODE pixelShaderByteCode;
	};

    class D3D12GraphicsManager : public GraphicsManager
    {
    public:
       	virtual int Initialize();
	    virtual void Finalize();

        virtual void Clear();

        virtual void Draw();

    protected:
        bool SetPerFrameShaderParameters();
        bool SetPerBatchShaderParameters(int32_t index);

        // HRESULT InitializeBuffers();
        // HRESULT InitializeShader(const char* vsFilename, const char* psFilename);
        // HRESULT RenderBuffers();

        void UpdateConstants();
        void InitializeBuffers(const Scene& scene);
        void ClearBuffers();
        void RenderBuffers();

    private:
        HRESULT CreateDescriptorHeaps();
        HRESULT CreateRenderTarget();
        HRESULT CreateDepthStencil();
        HRESULT CreateGraphicsResources();
        HRESULT CreateSamplerBuffer();
        HRESULT CreateTextureBuffer(SceneObjectTexture& texture);
        HRESULT CreateConstantBuffer();
        HRESULT CreateIndexBuffer(const SceneObjectIndexArray& indexArray);
        HRESULT CreateVertexBuffer(const SceneObjectVertexArray& vPropertyArray);
        HRESULT CreateRootSignature();
        HRESULT WaitForPreviousFrame();
        HRESULT PopulateCommandList();
		HRESULT InitializePSO();

		HRESULT CreateIntervalVertexBuffer();

		ID3DBlob* CompileShader(const std::wstring& filename,
			const D3D_SHADER_MACRO* defines,
			const std::string& entrypoint,
			const std::string& target);

    private:
        static const uint32_t           k_FrameCount  = 2;
        static const uint32_t           k_MaxSceneObjectCount = 2048;
        static const uint32_t           k_MaxTextureCount = 2048;
		static const uint32_t			k_MaxLightCount = 55;
        static const uint32_t           k_TextureDescStartIndex = k_FrameCount * k_MaxSceneObjectCount * 2;
        ID3D12Device*                   m_pDev       = nullptr;             // the pointer to our Direct3D device interface
        D3D12_VIEWPORT                  m_ViewPort;                         // viewport structure
        D3D12_RECT                      m_ScissorRect;                      // scissor rect structure
        IDXGISwapChain3*                m_pSwapChain = nullptr;             // the pointer to the swap chain interface
        ID3D12Resource*                 m_pRenderTargets[k_FrameCount];      // the pointer to rendering buffer. [descriptor]
        ID3D12Resource*                 m_pDepthStencilBuffer;              // the pointer to the depth stencil buffer
        ID3D12CommandAllocator*         m_pCommandAllocator = nullptr;      // the pointer to command buffer allocator
        ID3D12CommandQueue*             m_pCommandQueue = nullptr;          // the pointer to command queue
        ID3D12RootSignature*            m_pRootSignature = nullptr;         // a graphics root signature defines what resources are bound to the pipeline
        ID3D12DescriptorHeap*           m_pRtvHeap = nullptr;               // an array of descriptors of GPU objects
        ID3D12DescriptorHeap*           m_pDsvHeap = nullptr;               // an array of descriptors of GPU objects
        ID3D12DescriptorHeap*           m_pCbvHeap = nullptr;               // an array of descriptors of GPU objects
        ID3D12DescriptorHeap*           m_pSamplerHeap = nullptr;           // an array of descriptors of GPU objects
        ID3D12PipelineState*            m_pPipelineState = nullptr;         // an object maintains the state of all currently set shaders
                                                                            // and certain fixed function state objects
                                                                            // such as the input assembler, tesselator, rasterizer and output manager
		
        ID3D12GraphicsCommandList*      m_pCommandList = nullptr;           // a list to store GPU commands, which will be submitted to GPU to execute when done

		ID3D12RootSignature*			m_pRootSignatureResolve = nullptr;	// used for resolve msaa render target
		ID3D12PipelineState*			m_pPipelineStateResolve = nullptr;	// used for resolve msaa render target
        ID3D12Resource*                 m_pMsaaRenderTarget;                // render target used for msaa
		D3D12_VERTEX_BUFFER_VIEW		m_VertexBufferViewResolve;

        uint32_t                        m_RtvDescriptorSize;
        uint32_t                        m_CbvSrvDescriptorSize;

        std::vector<ID3D12Resource*>    m_Buffers;                          // the pointer to the vertex buffer
        std::vector<ID3D12Resource*>    m_Textures;                         // the pointer to the vertex buffer
        std::map<std::string, size_t>  m_TextureIndex;
        std::vector<D3D12_VERTEX_BUFFER_VIEW>        m_VertexBufferView;                 // a view of the vertex buffer
        std::vector<D3D12_INDEX_BUFFER_VIEW>         m_IndexBufferView;                  // a view of the index buffer
        ID3D12Resource*                 m_pTextureBuffer = nullptr;         // the pointer to the texture buffer

        struct PerBatchConstants
        {
            Matrix4f    objectMatrix;
            Vector4Df   diffuseColor;
            Vector4Df   specularColor;
            float       specularPower;
            bool        usingDiffuseMap;
            bool        usingNormalMap;
        };

        struct DrawBatchContext
        {
            uint32_t indexCount;
            uint32_t propertyCount;
            std::shared_ptr<SceneGeometryNode> node;
            std::shared_ptr<SceneObjectMaterial> material;
        };

        std::vector<DrawBatchContext> m_DrawBatchContext;

        uint8_t*            m_pCbvDataBegin = nullptr;
		size_t				k_SizePerFrameConstantBuffer;
		size_t				k_SizePerBatchConstantBuffer;
		size_t				k_SizeConstantBufferPerFrame;

        // Synchronization objects
        HANDLE                          m_FenceEvent;
        ID3D12Fence*                    m_pFence = nullptr;
        uint32_t                        m_FenceValue;
    };
}