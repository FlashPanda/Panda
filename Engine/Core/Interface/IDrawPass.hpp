#pragma once
#include <iostream>
#include "Interface.hpp"
#include "Core/GfxStructures.hpp"

namespace Panda
{
    Interface IDrawPass
    {
        public:
            IDrawPass() = default;
            virtual ~IDrawPass() {}

            virtual void Draw(Frame& frame) = 0;
    };
}