#pragma once
#include "Interface.hpp"
#include "core/Image.hpp"
#include "core/Buffer.hpp"

namespace Panda
{
    Interface ImageParser
    {
    public:
        virtual Image Parse(Buffer& buf) = 0;
    };
}
