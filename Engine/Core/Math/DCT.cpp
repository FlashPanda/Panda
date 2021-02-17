#include "MathUtility.hpp"
#include "Core/portable.hpp"
#include <math.h>
#include "DCT.hpp"

namespace Panda
{
    FORCEINLINE float NormalizeScaleFactor(float a)
    {
        return (a == 0)? 1.0f / sqrtf(2.0f) : 1.0f;
    }

    void DCT8x8 (const float in[64], float out[64])
    {
        // in case in and out are the same array
        float _in[64];
        memcpy(_in, in, sizeof(float) * 64);

        float oneOverFour = 1.0 / 4.0;
		float PIOverSixteen = PI / 16.0;
		for (int32_t u = 0; u < 8; ++u)
		{
			for (int32_t v = 0; v < 8; ++v)
			{
				float c = 0.0f;
				for (int32_t x = 0; x < 8; ++x)
				{
					for (int32_t y = 0; y < 8; ++y)
					{
						c += _in[x * 8 + y] * cosf((2 * x + 1) * u * PIOverSixteen) * cosf((2 * y + 1) * v * PIOverSixteen);
					}
				}

				c *= NormalizeScaleFactor(u) * NormalizeScaleFactor(v) * oneOverFour;
				out[u * 8 + v] = c;
			}
		}
    }

    void IDCT8x8 (const float in[64], float out[64])
    {
        // in case in and out are the same array
        float _in[64];
        memcpy(_in, in, sizeof(float) * 64);

		float oneOverFour = 1.0f / 4.0f;
		float PIOverSixteen = PI / 16.0f;

        for (int32_t x = 0; x < 8; ++x)
        {
            for (int32_t y = 0; y < 8; ++y)
            {
				float c = 0.0f;
                for (int32_t u = 0; u < 8; ++u)
                {
                    for (int32_t v = 0; v < 8; ++v)
                    {
                        c += _in[u * 8 + v] * NormalizeScaleFactor(u) * NormalizeScaleFactor(v) * cosf((2 * x + 1) * u * PIOverSixteen) * cosf((2 * y + 1) * v * PIOverSixteen);
                    }
                }

                out[x * 8 + y] = c * oneOverFour;
            }
        }
    }
}