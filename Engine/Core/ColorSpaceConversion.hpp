#pragma once
#include <algorithm>
#include "Math/PandaMath.hpp"
#include "portable.hpp"

namespace Panda
{
    typedef Vector<float, 3> RGBf;
    typedef Vector<float, 3> YCbCrf;

    const Matrix4f RGB2YCbCr({
        0.229f, -0.168736f,  0.5f,      0.0f,
        0.587f, -0.331264f, -0.418688f, 0.0f,
        0.114f,  0.5f     , -0.081312f, 0.0f,
        0.0f,    128.0f,     128.0f,    0.0f
    });

    const Matrix4f YCbCr2RGB({
              1.0f,    1.0f     ,    1.0f  ,    0.0f,
              0.0f, -  0.344136f,    1.772f,    0.0f,
            1.402f, -  0.714136f,    0.0f  ,    0.0f,
         -179.456f,  135.458816f, -226.816f,    0.0f,
    });

    inline YCbCrf ConvertRGB2YCbCr(const RGBf& rgb)
    {
		Vector4Df result({ rgb[0], rgb[1], rgb[2], 1.0f });
        TransformCoord(result, RGB2YCbCr);
		return YCbCrf({ std::clamp<float>(result[0] + 0.5f, 0.0f, 255.0f),
					  std::clamp<float>(result[1] + 0.5f, 0.0f, 255.0f),
					  std::clamp<float>(result[2] + 0.5f, 0.0f, 255.0f) });
    }

    inline RGBf ConvertYCbCr2RGB(const YCbCrf& ycbcr)
    {
		Vector4Df result({ ycbcr[0], ycbcr[1], ycbcr[2], 1.0f });
        TransformCoord(result, YCbCr2RGB);
		return RGBf({ std::clamp<float>(result[0] + 0.5f, 0.0f, 255.0f),
					std::clamp<float>(result[1] + 0.5f, 0.0f, 255.0f),
					std::clamp<float>(result[2] + 0.5f, 0.0f, 255.0f) });
    }
}