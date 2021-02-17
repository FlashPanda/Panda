#pragma once

#ifndef PI
#define PI 3.14159265358979323846f
#endif

#ifndef TWO_PI
#define TWO_PI 3.14159265358979323846f * 2.0f
#endif

#ifndef HALF_PI
#define HALF_PI 3.14159265358979323846f / 2.0f
#endif

#ifndef INV_PI
#define INV_PI 0.31830988618379067154f
#endif

#ifndef INV_2_PI
#define INV_2_PI 0.15915494309189533577f
#endif

#ifndef INV_4_PI
#define INV_4_PI 0.07957747154594766788f
#endif

#ifndef PI_OVER_2
#define PI_OVER_2 1.57079632679489661923f
#endif

#ifndef PI_OVER_4
#define PI_OVER_4 1.57079632679489661923f
#endif

#ifndef SQRT_2
#define SQRT_2 1.41421356237309504880f
#endif

#ifndef ONE_QUARTER_PI
#define ONE_QUARTER_PI 3.14159265358979323846f / 4.0f
#endif

#include <limits>
#include <cstdint>
#include <cmath>

namespace Panda
{
#ifdef WIN32
// For 32bit app, the value is 1.0 * 2^(-23)
#define MACHINE_EPSILON (std::numeric_limits<float>::epsilon() * 0.5f)
#else
static float MACHINE_EPSILON = (std::numeric_limits<float>::epsilon() * 0.5f);
#endif

#ifdef WIN32
#define MaxFloat std::numeric_limits<float>::max()
#define Infinity std::numeric_limits<float>::infinity()
#else
static float MaxFloat = std::numeric_limits<float>::max();
static float Infinity = std::numeric_limits<float>::infinity();
#endif

#ifndef NOMINMAX
#define NOMINMAX // We want to use std::min and std::max from the standard library
#endif

    static float ShadowEpsilon = 0.0001f;

    enum SortOrder
    {
        Increasing = 0,
        Decreasing = 1
    };

    inline float Lerp(float t, float v1, float v2)
    {
        return (1 - t) * v1 + t * v2;
    }

    inline float Gamma(int32_t n)
    {
        return (n * MACHINE_EPSILON) / (1 - n * MACHINE_EPSILON);
    }

    inline float Radians (float deg) 
    {
        return (PI / 180) * deg;
    }

    inline float Degrees(float rad)
    {
        return (180 / PI) * rad;
    }

	template <typename T, typename U, typename V>
	inline T Clamp(T val, U low, V high) {
		if (val < low)
			return low;
		else if (val > high)
			return high;
		else
			return val;
	}

	template <typename T>
	inline T Abs(const T data)
	{
		return std::abs(data);
	}

    inline uint32_t FloatToBits(float f)
    {
        uint32_t ui = *reinterpret_cast<uint32_t*>(&f);
        return ui;
    }

    inline float BitsToFloat(uint32_t ui)
    {
        float f = *reinterpret_cast<float*>(&ui);
        return f;
    }

    inline uint64_t FloatToBits(double f)
    {
        uint64_t ui = *reinterpret_cast<uint64_t*>(&f);
        return ui;
    }

    inline double BitsToFloat(uint64_t ui)
    {
        double f = *reinterpret_cast<double*>(&ui);
        return f;
    }

    // Get the minnest float which is bigger than input
    inline float NextFloatUp(float v)
    {
        // Handle infinity and negative zero situation
        if (std::isinf(v) && v > 0.f) return v;
        if (v == -0.f) v = 0.f;

        // Advance _v_ to next higher float
        uint32_t ui = FloatToBits(v);
        if (v >= 0)
            ++ui;
        else
            --ui;
        return BitsToFloat(ui);
    }

    // Get the maxest float which is smaller than input
    inline float NextFloatDown(float v)
    {
        // Handle infinity and positive zero situation
        if (std::isinf(v) && v < 0.f) return v;
        if (v == 0.f) v = -0.f;

        uint32_t ui = FloatToBits(v);
        if (v > 0)
            --ui;
        else
            ++ui;
        return BitsToFloat(ui);
    }

    template <typename Predicate>
    int32_t FindInterval(int32_t size, const Predicate& pred)
    {
        int32_t first = 0, len = size;
        while (len > 0)
        {
            int32_t half = len >> 1, middle = first + half;
            // Bisect range based on value of _pred_ at _middle_
            if (pred(middle))
            {
                first = middle + 1;
                len -= half + 1;
            }
            else 
                len = half;
        }
        return Clamp(first - 1, 0, size - 2);
    }
}