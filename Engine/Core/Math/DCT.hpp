#pragma once

#include "MathUtility.hpp"

namespace Panda
{
    FORCEINLINE float NormalizeScaleFactor(float a);

    void DCT8x8 (const float in[64], float out[64]);

    void IDCT8x8 (const float in[64], float out[64]);

}