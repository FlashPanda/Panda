#pragma once
#include <cstdio>
#include <iostream>
#include <string>
#include <cassert>
#include <queue>
#include <algorithm>
#include "Core/Interface/ImageParser.hpp"
#include "Core/portable.hpp"

namespace Panda
{
#pragma pack(push, 1)
    struct TGA_FILEHEADER
    {
        uint8_t IDLength;
        uint8_t ColorMapType;
        uint8_t ImageType;
        uint8_t ColorMapSpec[5];
        uint8_t ImageSpec[10];
    };
#pragma pack(pop)
    class TgaParser : implements ImageParser
    {
        public:
            virtual Image Parse(Buffer& buf);
    };
}