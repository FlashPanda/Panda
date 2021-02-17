#pragma once
#include "Core/Interface/IDrawPass.hpp"

namespace Panda
{
    class ForwardRenderPass : implements IDrawPass
    {
        public:
            ~ForwardRenderPass() = default;
            void Draw(Frame& frame) final;
    };
}