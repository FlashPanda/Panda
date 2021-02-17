#include "Ray.hpp"

namespace Panda
{
    Vector3Df OffsetRayOrigin(const Vector3Df& p, const Vector3Df& pError,
                              const Vector3Df& n, const Vector3Df& w)
    {
        float d = DotProduct(Abs(n), pError);

        Vector3Df offset = d * n;
        if (DotProduct(w, n) < 0.f) offset = -offset;
        Vector3Df po = p + offset;
        // Round offset point _po_ away from _p_
        for (int32_t i = 0; i < 3; ++i)
        {
            if (offset[i] > 0)
                po[i] = NextFloatUp(po[i]);
            else if (offset[i] < 0.f)
                po[i] = NextFloatDown(po[i]);
        }
        return po;
    }
}