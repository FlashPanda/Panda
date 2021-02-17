#include "TGA.hpp"

namespace Panda
{
    Image TgaParser::Parse(Buffer& buf)
    {
        Image img;
        const uint8_t* pData = buf.GetData();
        const uint8_t* pDataEnd = buf.GetData() + buf.GetDataSize();

        // without any checking?
        std::cout << "Parsing as TGA file:" << std::endl;

        const TGA_FILEHEADER* pFileHeader = reinterpret_cast<const TGA_FILEHEADER*>(pData);
        pData += sizeof(TGA_FILEHEADER);

        #ifdef DEBUG
        std::cout << "ID Length: " << (uint16_t)pFileHeader->IDLength << std::endl;
        std::cout << "Color Map Type: " << (uint16_t)pFileHeader->ColorMapType << std::endl;
        #endif

        if (pFileHeader->ColorMapType)
        {
            std::cout << "Unsupported Color Map. Only Type 0 is supported." << std::endl;
            return img;
        }

        #ifdef DEBUG
        std::cout << "Image Type: " << (uint16_t)pFileHeader->ImageType << std::endl;
        #endif
        if (pFileHeader->ImageType != 2)
        {
            std::cout << "Unsupoorted Image Type. Only Type 2 is suppoerted. " << std::endl;
            return img;
        }

        img.Width = (pFileHeader->ImageSpec[5] << 8) + pFileHeader->ImageSpec[4];
        img.Height = (pFileHeader->ImageSpec[7] << 8) + pFileHeader->ImageSpec[6];
        uint8_t pixelDepth = pFileHeader->ImageSpec[8];
        uint8_t alphaDepth = pFileHeader->ImageSpec[9] & 0x0F;
        #ifdef DEBUG
        std::cout << "Image Width: " << img.Width << std::endl;
        std::cout << "Image Height: " << img.Height << std::endl;
        std::cout << "Image Pixel Depth: " << (uint16_t)pixelDepth << std::endl;
        std::cout << "Image Alpha Depth: " << (uint16_t)alphaDepth << std::endl;
        #endif
        // skin Image 2D
        pData += pFileHeader->IDLength;
        // skip the Color Map. Since we assume the color map type is 0.
        // nothing to skip

        // reading the pixel data
        img.BitCount = 32;
        img.Pitch = (img.Width * (img.BitCount >> 3) + 3) & ~3u; // for GPU address alignment
        img.DataSize = img.Pitch * img.Height;
        img.Data = g_pMemoryManager->Allocate(img.DataSize);

        uint8_t* pOut = (uint8_t*)img.Data;
        for (size_t i = 0; i < img.Height; ++i)
        {
            for (size_t j = 0; j < img.Width; ++j)
            {
                switch(pixelDepth)
                {
                    case 15:
                    {
                        uint16_t color = *(uint16_t*)pData;
                        pData += 2;
                        *(pOut + img.Pitch * i + j * 4) = ((color & 0x7C00) >> 10);     // R
                        *(pOut + img.Pitch * i + j * 4 + 1) = ((color & 0x03E0) >> 5);  // G
                        *(pOut + img.Pitch * i + j * 4 + 2) = (color & 0x001F);         // B
                        *(pOut + img.Pitch * i + j * 4 + 3) = 0xFF;                     // A
                    }
                    break;
                    case 16:
                    {
                        uint16_t color = *(uint16_t*)pData;
                        pData += 2;
                        *(pOut + img.Pitch * i + j * 4) = ((color & 0x7C00) >> 10);     // R
                        *(pOut + img.Pitch * i + j * 4 + 1) = ((color & 0x03E0) >> 5);  // G
                        *(pOut + img.Pitch * i + j * 4 + 2) = (color & 0x001F);         // B
                        *(pOut + img.Pitch * i + j * 4 + 3) = ((color & 0x8000)? 0xFF : 0x00); // A
                    }
                    break;
                    case 24:
                    {
                        *(pOut + img.Pitch * i + j * 4 + 2) = *pData;   // B
                        pData++;
                        *(pOut + img.Pitch * i + j * 4 + 1) = *pData;   // G
                        pData++;
                        *(pOut + img.Pitch * i + j * 4) = *pData;       // R
                        pData++;
                        *(pOut + img.Pitch * i + j * 4 + 3) = 0xFF;     // A
                    }
                    break;
                    case 32:
                    {
                        *(pOut + img.Pitch * i + j * 4 + 2) = *pData;   // B
                        pData++;
                        *(pOut + img.Pitch * i + j * 4 + 1) = *pData;   // G
                        pData++;
                        *(pOut + img.Pitch * i + j * 4) = *pData;       // R
                        pData++;
                        *(pOut + img.Pitch * i + j * 4 + 3) = *pData;   // A
                        pData++;
                    }
                    break;
                    default:
                    break;
                }
            }
        }

        assert(pData <= pDataEnd);

        return img;
    }
}