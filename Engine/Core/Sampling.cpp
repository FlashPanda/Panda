#include "Sampling.hpp"

namespace Panda
{
    Point2Df ConcentricSampleDisk(const Point2Df& u)
    {
        // Map uniform random numbers to $[-1,1]^2$
        Point2Df uOffset = 2.f * u - Vector2Df({1.f, 1.f});

        // Handle degeneracy at the origin
        if (uOffset[0] == 0.f && uOffset[1] == 0.f) return Point2Df({0.f, 0.f});

        // Apply concentric mapping to point
        float theta, r;
        if ((std::abs)(uOffset[0]) > (std::abs)(uOffset[1]))
        {
            r = uOffset[0];
            theta = PI_OVER_4 * (uOffset[1] / uOffset[0]);
        }
        else 
        {
            r = uOffset[1];
            theta = PI_OVER_2 - PI_OVER_4 * (uOffset[0] / uOffset[1]);
        }
        return r * Point2Df({std::cosf(theta), std::sinf(theta)});
    }
}