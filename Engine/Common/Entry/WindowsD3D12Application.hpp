#pragma once
#include "WindowsApplication.hpp"

namespace Panda 
{
    class D3D12Application : public WindowsApplication 
    {
        public:
            using WindowsApplication::WindowsApplication;
            void Tick();
    };
}