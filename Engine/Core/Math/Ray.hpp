#pragma once

#include "MathUtility.hpp"
#include "Vector.hpp"
#include "Core/Medium.hpp"

namespace Panda
{
    struct Ray
    {
        Ray() : tMax(FLT_MAX), pMedium(nullptr){}
        Ray(const Vector3Df& origin, const Vector3Df& direction, 
            float _tMax = FLT_MAX, float time = 0.f, const Medium* _pMedium = nullptr) 
            : o(origin), d(direction), tMax(_tMax), time(time), pMedium(_pMedium)
            { }

        Vector3Df operator()(float t) const { return o + t * d;}
        //friend std::ostream& operator<<(std::ostream& out, const Ray& r)
        //{
        //    out << "[o = " << r.o << ", d = " << r.d << ", tMax = " << r.tMax
        //        << ", time = " << r.time << "]";
        //    return out;
        //}

        Vector3Df o; // origin point
        Vector3Df d; // direction
        mutable float tMax; // The param t must be in [0, tMax)
        float time;
        const Medium* pMedium;
    };

    struct RayDifferential : public Ray
    {
        RayDifferential() : Ray(), hasDifferentials(false) {}
        RayDifferential(const Vector3Df& o, const Vector3Df& d, 
            float tMax = FLT_MAX, float time = 0.f, const Medium* pMedium = nullptr)
            : Ray(o, d, tMax, time, pMedium), hasDifferentials(false) {}
        RayDifferential(const Ray& ray) : Ray(ray), hasDifferentials(false) {}
        
        void ScaleDifferentials(float s)
        {
            rxOrigin = o + (rxOrigin - o) * s;
            ryOrigin = o + (ryOrigin - o) * s;
            rxDirection = d + (rxDirection - d) * s;
            ryDirection = d + (ryDirection - d) * s;
        }

        //friend std::ostream& operator<<(std::ostream& out, const RayDifferential& ray)
        //{
        //    out << "[ " << (Ray&)ray << " has differiendtials: "
        //        << (ray.hasDifferentials? "true" : "false") << ", xo = " << ray.rxOrigin
        //        << ", xd = " << ray.rxDirection << ", yo = " << ray.ryOrigin << ", yd = "
        //        << ray.ryDirection;
        //    return out;
        //}

        bool hasDifferentials;
        Vector3Df rxOrigin, ryOrigin;
        Vector3Df rxDirection, ryDirection;
    };

    Vector3Df OffsetRayOrigin(const Vector3Df& p, const Vector3Df& pError,
                              const Vector3Df& n, const Vector3Df& w);
}
