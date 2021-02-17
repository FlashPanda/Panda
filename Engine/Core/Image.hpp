#pragma once
#include <iostream>
#include "Math/PandaMath.hpp"

namespace Panda
{
    struct Image
    {
        uint32_t Width; // width in pixels
        uint32_t Height; // height in pixels
        void* Data;
        uint32_t BitCount; // bit count per pixel
        uint32_t Pitch; // size of one line, in bytes
        size_t DataSize; // the size of data area, which is pitch * height rather than width * height * bitcount / 8, because of alignments

        Image() : Width(0),
            Height(0),
            Data(nullptr),
            BitCount(0),
            Pitch(0),
            DataSize(0)
            {}
    };

    std::ostream& operator<<(std::ostream& out, const Image& image);
}

