#include <iostream>
#include <tchar.h>
#include "WindowsApplication.hpp"
#include "D2D/D2DGraphicsManager.hpp"
#include "MemoryManager.hpp"
#include "AssetLoader.hpp"
#include "SceneManager.hpp"
#include "Utility.hpp"
#include "BMP.hpp"
#include "JPEG.hpp"
#include "PNG.hpp"

using namespace Panda;
using namespace std;

namespace Panda {
    class TestGraphicsManager : public D2DGraphicsManager
    {
        public:
            using D2DGraphicsManager::D2DGraphicsManager;
            void DrawImage(const Image image);
        private:
            ID2D1Bitmap* m_pBitmap = nullptr;
    };

    class TestApplication : public WindowsApplication
    {
    public:
        using WindowsApplication::WindowsApplication;

        virtual int Initialize();

        virtual void OnDraw();

    private:
        Image m_Image;
    };
}

namespace Panda {
	GfxConfiguration config(8, 8, 8, 8, 32, 0, 0, 1024, 512, "Texture Load Test (Windows)");
    GameLogic* g_pGameLogic = static_cast<GameLogic*>(new GameLogic);
}

int Panda::TestApplication::Initialize()
{
    int result;

    result = WindowsApplication::Initialize();

    if (result == 0) {
        Buffer buf;

        PngParser pngParser;
        if(m_ArgC > 1)
        {
            buf = g_pAssetLoader->SyncOpenAndReadBinary(m_ppArgV[1]);
        }
        else 
        {
            buf = g_pAssetLoader->SyncOpenAndReadBinary("Textures/eye.png");
        }

        m_Image = pngParser.Parse(buf);
    }

    if (m_Image.BitCount == 24)
    {
        // DXGI does not have 24 bit formats so we have to extend it to 32bit
        uint32_t newPitch = m_Image.Pitch / 3 * 4;
        size_t dataSize = newPitch * m_Image.Height;
        void* data = g_pMemoryManager->Allocate(dataSize);
        uint8_t* buf = reinterpret_cast<uint8_t*>(data);
        uint8_t* src = reinterpret_cast<uint8_t*>(m_Image.Data);
        for (auto row = 0; row < m_Image.Height; ++row)
        {
            buf = reinterpret_cast<uint8_t*>(data) + row * newPitch;
            src = reinterpret_cast<uint8_t*>(m_Image.Data) + row * m_Image.Pitch;
            for (auto col = 0; col < m_Image.Pitch; ++col)
            {
                *(uint32_t*)buf = *(uint32_t*)src;
                buf[3] = 0; // set alpha to 0
                buf += 4;
                src += 3;
            }
        }

        g_pMemoryManager->Free(m_Image.Data, m_Image.DataSize);
        m_Image.Data = data;
        m_Image.DataSize = dataSize;
        m_Image.Pitch = newPitch;
    }

    return result;
}

void Panda::TestApplication::OnDraw()
{
    dynamic_cast<TestGraphicsManager*>(g_pGraphicsManager)->DrawImage(m_Image);
}

void Panda::TestGraphicsManager::DrawImage(const Image image)
{
	HRESULT hr;

    // start build GPU draw command
    m_pRenderTarget->BeginDraw();

    D2D1_BITMAP_PROPERTIES props;
    props.pixelFormat.format = DXGI_FORMAT_R8G8B8A8_UNORM;
    props.pixelFormat.alphaMode = D2D1_ALPHA_MODE_IGNORE;
    props.dpiX = 72.0f;
    props.dpiY = 72.0f;
    SafeRelease(&m_pBitmap);
    hr = m_pRenderTarget->CreateBitmap(D2D1::SizeU(image.Width, image.Height), 
                                                    image.Data, image.Pitch, props, &m_pBitmap);

    D2D1_SIZE_F rtSize = m_pRenderTarget->GetSize();
    D2D1_SIZE_F bmpSize = m_pBitmap->GetSize();

    D2D1_RECT_F source_rect = D2D1::RectF(
                     0,
                     0,
                     bmpSize.width,
                     bmpSize.height
                     );

    float aspect = bmpSize.width / bmpSize.height;
	float dest_height = rtSize.height;
	float dest_width = rtSize.height * aspect;

    D2D1_RECT_F dest_rect = D2D1::RectF(
                     0,
                     0,
                     dest_width,
                     dest_height 
                     );

    m_pRenderTarget->DrawBitmap(m_pBitmap, dest_rect, 1.0f, D2D1_BITMAP_INTERPOLATION_MODE_NEAREST_NEIGHBOR, source_rect);

    // end GPU draw command building
    m_pRenderTarget->EndDraw();
}


