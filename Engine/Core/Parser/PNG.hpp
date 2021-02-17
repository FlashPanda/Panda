#pragma once

#include <cstdio>
#include <iostream>
#include <string>
#include <cassert>
#include <queue>
#include <algorithm>
#include "Core/Utility.hpp"
#include "Core/Interface/ImageParser.hpp"
#include "Core/portable.hpp"
#include "zlib/zlib.h"

namespace Panda
{
#pragma pack(push, 1)
    struct PNG_FILEHEADER
    {
        uint64_t Signature;
    };

    ENUM(PNG_CHUNK_TYPE)
    {
        IHDR = "IHDR"_u32,
        PLTE = "PLTE"_u32,
        IDAT = "IDAT"_u32,
        IEND = "IEND"_u32
    };

    static std::ostream& operator<<(std::ostream& out, PNG_CHUNK_TYPE type)
    {
        int32_t n = static_cast<int32_t> (type);
        n = to_endian_native(n);
        char* c = reinterpret_cast<char*>(&n);

        for (size_t i = 0; i < sizeof(int32_t); ++i)
        {
            out << *c++;
        }

        return out;
    }

    struct PNG_CHUNK_HEADER 
    {
        uint32_t    Length;
        PNG_CHUNK_TYPE Type;
    };

    struct PNG_IHDR_HEADER : PNG_CHUNK_HEADER
    {
        uint32_t    Width;
        uint32_t    Height;
        uint8_t     BitDepth;
        uint8_t     ColorType;
        uint8_t     CompressionMethod;
        uint8_t     FilterMethod;
        uint8_t     InterlaceMethod;
    };

    struct PNG_PLTE_HEADER : PNG_CHUNK_HEADER
    {
        Vector<uint8_t, 3>* pEntries;
    };
#pragma pack(pop)

    // report a zlib or i/o error
    static void zerr(int ret);

    class PngParser : implements ImageParser
    {
        protected:
            uint16_t m_Width;
            uint16_t m_Height;
            uint8_t  m_BitDepth;
            uint8_t  m_ColorType;
            uint8_t  m_CompressionMethod;
            uint8_t  m_FilterMethod;
            uint8_t  m_InterlaceMethod;
            size_t   m_ScanLineSize;
            uint8_t  m_BytesPerPixel;

        public:
            virtual Image Parse(Buffer& buf);
    };
}
