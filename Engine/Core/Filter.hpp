#pragma once

#include "Math/PandaMath.hpp"

namespace Panda
{
    class Filter
    {
        public:
            virtual ~Filter();
            Filter(const Vector2Df& radius)
                : m_Radius(radius), m_InvRadius(Vector2Df({1.f / radius[0], 1.f / radius[1]}))
                {}
            virtual float Evaluate (const Point2Df& p) const = 0;

            const Vector2Df m_Radius, m_InvRadius;
    };
}
