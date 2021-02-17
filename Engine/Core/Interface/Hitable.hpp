#pragma once

#include "Interface.hpp"
#include "Math/Ray.hpp"

namespace Panda
{
    struct HitRecord
    {
        float t;
        Vector3Df p;        // hit position
        Vector3Df normal;   // the normal of hit point 
        std::string matID;  // the ID of a material
    };

    Interface Hitable
    {
        public:
            virtual bool Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const = 0;
    };
}