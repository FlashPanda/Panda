#pragma once
#include <iostream>
#include "../Interface/ImageParser.hpp"

namespace Panda
{
#pragma pack(push, 1)
    typedef struct _BITMAP_FILEHEADER
    {
        uint16_t signature;
        uint32_t size;
        uint32_t reserved;
        uint32_t bitsOffset;
    } BITMAP_FILEHEADER;

#define BITMAP_FILEHEADER_SIZE 14

    typedef struct _BITMAP_HEADER
    {
        uint32_t headerSize;
        int32_t width;
        int32_t height;
        uint16_t planes;
        uint16_t bitCount;
        uint32_t compression;
        uint32_t sizeImage;
        int32_t pelsPerMeterX;
        int32_t pelsPerMeterY;
        uint32_t clrUsed;
        uint32_t clrImportant;
    } BITMAP_HEADER; 
#pragma pack(pop)

    class BmpParser : implements ImageParser
    {
    public:
        virtual Image Parse(Buffer& buf);
    };
}
