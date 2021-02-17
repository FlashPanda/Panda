#pragma once

#include "MathUtility.hpp"
#include <vector>
#include <cassert>
#include "Vector.hpp"
#include "Ray.hpp"

namespace Panda
{
    template <typename T, int N>
    struct Bounds
    {
        Vector<T, N> pMin, pMax;

        Bounds() 
        {
            T minNum = (std::numeric_limits<T>::lowest)();
            T maxNum = (std::numeric_limits<T>::max)();
            pMin = Vector<T, N>({maxNum, maxNum});
            pMax = Vector<T, N>({minNum, minNum});  // pMin > pMax, which means this is a invalid bound.
        }
        Bounds(const Vector<T, N>& p) : pMin(p), pMax(p) {}
        Bounds(const Vector<T, N>& p1, const Vector<T, N>& p2)
        {
            pMin = Vector<T, N>({std::min(p1.data[0], p2.data[0]), std::min(p1.data[1], p2.data[1])});
            pMin = Vector<T, N>({std::max(p1.data[0], p2.data[0]), std::max(p1.data[1], p2.data[1])});
        }

        Vector<T, N> Diagonal() const {return pMax - pMin;}

        // To 2D, it's an area.
        // To 3D, it's a volume.
        T Space() const 
        {
            Vector<T, N> d = Diagonal();
            T mul = T(1.0);
            for (size_t i = 0; i < N; ++i)
                mul *= d.data[i];
            return mul;
        }

        inline const Vector<T, N>& operator[] (int32_t i ) const
        {
            return (i == 0)? pMin : pMax;
        }
        inline Vector<T, N>& operator[] (int32_t i )
        {
            return (i == 0)? pMin : pMax;
        }

        bool operator==(const Bounds<T, N>& b) const 
        {
            return b.pMin == pMin && b.pMax == pMax;
        }
        bool operator!=(const Bounds<T, N>& b) const
        {
            return b.pMin != pMin || b.pMax != pMax;
        }

        Vector<T, 2> Lerp(const Vector<T, 2>& t) const
        {
            assert(N == 2);
            return Vector<T, 2>({Lerp(t.data[0], pMin.data[0], pMax.data[0]), Lerp(t.data[1], pMin.data[1], pMax.data[1])});
        }

        Vector<T, 2> Offset(const Vector<T, 2>& p) const 
        {
            assert(N == 2);
            Vector<T, 2> o = p - pMin;
            if (pMax.data[0] > pMin.data[0]) o.data[0] /= pMax.data[0] - pMin.data[0];
            if (pMax.data[1] > pMin.data[1]) o.data[1] /= pMax.data[1] - pMin.data[1];
            return o;
        }

        Vector<T, 3> Lerp(const Vector<T, 3>& t) const 
        {
            assert(N == 3);
            return Vector<T, 3> ({
                Lerp(t.data[0], pMin.data[0], pMax.data[0]),
                Lerp(t.data[1], pMin.data[1], pmax.data[1]),
                Lerp(t.data[2], pMin.data[2], pMax.data[2])
            });
        }

        // Return the offset relative to the minimum corner.
        // minimum corner has offset(0, 0, 0), maximum corner has offset(1, 1, 1)
        // Note that the value can larger than the offset of maximum corner.
        Vector<T, 3> Offset(const Vector<T, 3>& p) const 
        {
            Vector<T, 3> o = p - pMin;
            if (pMax.data[0] > pMin.data[0]) o.data[0] /= pMax.data[0] - pMin.data[0];
            if (pMax.data[1] > pMin.data[1]) o.data[1] /= pMax.data[1] - pMin.data[1];
            if (pMax.data[2] > pMin.data[2]) o.data[2] /= pMax.data[2] - pMin.data[2];
            return o;
        }

        // Return the coordinates of one of the eight corners of the bounding box
        Vector<T, 3> Corner(int32_t corner) const 
        {
            assert(N == 3);
            return Vector<T, 3>({
                (*this)[(corner & 0x1)].data[0], 
                (*this)[(corner & 0x10)? 1 : 0].data[1],
                (*this)[(corner & 0x100)? 1 : 0].data[2]
            });
        }

        // Return the surface area of the six faces of the box
        T SurfaceAreas() const 
        {
            assert (N == 3);
            Vector<T, 3> d = Diagonal();
            return 2 * (d.data[0] * d.data[1] + d.data[0] * d.data[2] + d.data[1] * d.data[2]);
        }

        // Return the max index of the three axeses
        int32_t MaximumExtent() const 
        {
            Vector<T, N> d = Diagonal();
            int32_t index = 0;
            T maxValue = d.data[0];
            for (int32_t i = 1; i < N; ++i)
            {
                if (d.data[i] > maxValue)
                {
                    index = i;
                    maxValue = d.data[i];
                }
            }
            
            return index;
        }

        // Return the center and the radius of a sphere which bounds the bounding box.
        void BoundingSphere(Vector<T, 3>& center, float& radius) const 
        {
            assert(N == 3);
            center.data[0] = (pMin.data[0] + pMax.data[0]) / 2;
            center.data[1] = (pMin.data[1] + pMax.data[1]) / 2;
            center.data[2] = (pMin.data[2] + pMax.data[2]) / 2;
            radius = Inside(center, *this)? GetLength(pMax - center) : 0;
        }

        // Calculate the intersections of a ray and a bound.
        // Output parametric ts.
        bool IntersectP(const Ray& ray, float& hitt0, float hitt1) const
        {
            float t0 = 0.f, t1 = ray.tMax;
            for (int32_t i = 0; i < 3; ++i)
            {
                // Update interval for bounding box slab
                float invRayDir = 1.f / ray.d[i];
                float tNear = (pMin[i] - ray.o[i]) * invRayDir;
                float tFar = (pMax[i] - ray.o[i]) * invRayDir;

                // Update parametric interval from slab intersection $t$ values
                if (tNear > tFar) std::swap(tNear, tFar);

                // Update _tFar_ to ensure robust ray-bound intersection
                tFar *= 1 + 2 * Gamma(3); // Make sure the comparison is correct even if _tFar_ error bound overlaps _tNear_'s.
                t0 = tNear > t0? tNear : t0; // if tNear is NaN, > operation will return false which remains t0 unchanged.
                t1 = tFar < t1? tFar : t1; // BUT, t0 > tNear also returns false. So, we can't change the code to 
                                            // t0 > tNear? t0 : tNear;
                if (t0 > t1) return false; // t0 > t1 means that the ray intersect with bound behind its direction.
            }
        }

        // Check for intersection.
        inline bool IntersectP(const Ray& ray, const Vector3Df& invDir, const int dirIsNeg[3]) const
        {
            const Bounds<T, N>& bounds = *this;
            // Check for ray intersection against $x$ and $y$ slabs
            float tMin = (bounds[dirIsNeg[0]].data[0] - ray.o.data[0]) * invDir.data[0];
            float tMax = (bounds[1 - dirIsNeg[0]].data[0] - ray.o.data[0]) * invDir.data[0];
            float tyMin = (bounds[dirIsNeg[1]].data[1] - ray.o.data[1]) * invDir.data[1];
            float tyMax = (bounds[1 - dirIsNeg[1]].data[1] - ray.o.data[1]) * invDir.data[1];

            // Update _tMax_ and _tyMax_ to ensure robust bounds interaction
            tMax *= 1 + 2 * Gamma(3);
            tyMax *= 1 + 2 * Gamma(3);
            if (tMin > tyMax || tyMin > tMax) return false;
            if (tyMin > tMin) tMin = tyMin;
            if (tyMax < tMax) tMax = tyMax;

            // Check for ray intersection against $z$ slab
            float tzMin = (bounds[dirIsNeg[2]].data[2] - ray.o.data[2]) * invDir.data[2];
            float tzMax = (bounds[1 - dirIsNeg[2]].data[2] - ray.o.data[2]) * invDir.data[2];

            // Update _txMax_ to ensure robust bounds intersection
            tzMax *= 1 + 2 * Gamma(3);
            if (tMin > tzMax || tzMin > tMax) return false;
            if (tzMin > tMin) tMin = tzMin;
            if (tzMax < tMax) tMax = tzMax;
            return (tMin < ray.tMax) & (tMax > 0);
        }
    };

    template <typename T, int32_t N>
    bool Overlaps(const Bounds<T, N>& b1, const Bounds<T, N>& b2)
    {
        bool det[N] = {false};
        for (int32_t i = 0; i < N; ++i)
            det[i] = (b1.pMax.data[i] >= b2.pMin.data[i]) && (b1.pMim.data[i] <= b2.pMax.data[i]);
        for (int32_t i = 0; i < N; ++i)
            if (!det[i])
                return false;

        return true;
    }

    // Test a point is or is not inside a bounding box.
    template <typename T, int32_t N>
    bool Inside(const Vector<T, N>& p, const Bounds<T, N>& b)
    {
        return  p.data[0] >= b.pMin.data[0] && p.data[0] <= b.pMax.data[0] &&
                p.data[1] >= b.pMin.data[1] && p.data[1] <= b.pMax.data[1] &&
                p.data[2] >= b.pMin.data[2] && p.data[2] <= b.pMax.data[2];
    }

    // Test a point is or is not exactly in a bounding box(not on the borders)
    template <typename T, int32_t N>
    bool InsideExclusive(const Vector<T, N>& p, const Bounds<T, N>& b)
    {
        return  p.data[0] > b.pMin.data[0] && p.data[0] < b.pMax.data[0] &&
                p.data[1] > b.pMin.data[1] && p.data[1] < b.pMax.data[1] &&
                p.data[2] > b.pMin.data[2] && p.data[2] < b.pMax.data[2];
    }

    // Pad the bounding box with a constant factor in all dimensions.
    template <typename T, int32_t N>
    inline Bounds<T, N> Expand(const Bounds<T, N>& b, T delta)
    {
        return Bounds<T, N>(b.pMin - Vector<T, N>({delta, delta, delta}),
            b.pMax + Vector<T, N>({delta, delta, delta}));
    }

    // Return a new bounding box that encompasses that point as well as the origin box.
    template <typename T>
    Bounds<T, 3> Union (const Bounds<T, 3>& b, const Vector<T, 3>& p)
    {
        return Bounds<T, 3>(
            Vector<T, 3> ({std::min(b.pMin.data[0], p.data[0]),
                            std::min(b.pMin.data[1], p.data[1]),
                            std::min(b.pMin.data[2], p.data[2])}),
            Vector<T, 3> ({std::max(b.pMax.data[0], p.data[0]),
                            std::max(b.pMax.data[1], p.data[1]),
                            std::max(b.pMax.data[2], p.data[2])})
        );
    }

    // Return a new bounding box that encompasses two bounding boxes.
    template <typename T>
    Bounds<T, 3> Union(const Bounds<T, 3>& b1, const Bounds<T, 3>& b2)
    {
        return Bounds<T, 3> (
            Vector<T, 3> ({ std::min(b1.pMin.data[0], b2.pMin.data[0]),
                            std::min(b1.pMin.data[1], b2.pMin.data[1]),
                            std::min(b1.pMin.data[2], b2.pMin.data[2])}),
            Vector<T, 3> ({ std::max(b1.pMax.data[0], b2.pMax.data[0]),
                            std::max(b1.pMax.data[1], b2.pMax.data[1]),
                            std::max(b1.pMax.data[2], b2.pMax.data[2])})
        );
    }

    // Return the intersection of two bounding boxes.
    template <typename T>
    Bounds<T, 3> Intersect(const Bounds<T, 3>& b1, const Bounds<T, 3>& b2)
    {
        return Bounds<T, 3> (
            Vector<T, 3> ({ std::max(b1.pMin.data[0], b2.pMin.data[0]),
                            std::max(b1.pMin.data[1], b2.pMin.data[1]),
                            std::max(b1.pMin.data[2], b2.pMin.data[2])}),
            Vector<T, 3> ({ std::min(b1.pMax.data[0], b2.pMin.data[0]),
                            std::min(b1.pMax.data[1], b2.pMax.data[1]),
                            std::min(b1.pMax.data[2], b2.pMax.data[2])})
        );
    }

    typedef Bounds<float, 2> Bounds2Df;
    typedef Bounds<float, 3> Bounds3Df;
    typedef Bounds<int32_t, 2> Bounds2Di;
    typedef Bounds<int32_t, 3> Bounds3Di;
}