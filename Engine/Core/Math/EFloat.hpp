#pragma once

#include <cmath>
#include "MathUtility.hpp"
#include <algorithm>
#include <iostream>

namespace Panda
{
    class EFloat
    {
    public:
        EFloat() {}
        EFloat(float _v, float _err = 0.f)
            : v(_v)
        {
            if (_err == 0.f)
                low = high = v;
            else 
            {
                // Compute conservative bounds by rounding the endpoints away
                // from the middle.
                low = NextFloatDown(v - _err);
                high = NextFloatUp(v + _err);
            }
            #ifndef NDEBUG
            vPrecise = v;
            Check();
            #endif 
        }
        EFloat (const EFloat& ef)
        {   
            ef.Check();
            v = ef.v;
            low = ef.low;
            high = ef.high;
#ifndef NDEBUG
            vPrecise = ef.vPrecise;
#endif
        }
        EFloat& operator=(const EFloat& ef)
        {
            ef.Check();
            if (&ef != this)
            {
                v = ef.v;
                low = ef.low;
                high = ef.high;
#ifndef NDEBUG
                vPrecise = ef.vPrecise;
#endif
            }
            return *this;
        }
        
        #ifndef NDEBUG
        EFloat (float v, long double lD, float err) : EFloat (v, err)
        {
            vPrecise = lD;
            Check();
        }
        #endif

        explicit operator float() const {return v;}
        explicit operator double() const {return v;}
        float GetAbsoluteError() const {return high - low;}
        float UpperBound() const {return high;}
        float LowerBound() const {return low;}
#ifndef NDEBUG
        float GetRelativeError() const 
        {
            return Abs((vPrecise - v) / vPrecise);
        }
        long double PreciseValue() const {return vPrecise;}
#endif 

        EFloat operator+(EFloat ef) const 
        {
            EFloat r;
            r.v = v + ef.v;
#ifndef NDEBUG
            r.vPrecise = vPrecise + ef.vPrecise;
#endif
            // Interval arichemetic addition, with the result rounded away from
            // the value r.v in order to be conservative.
            r.low = NextFloatDown(LowerBound() + ef.LowerBound());
            r.high = NextFloatUp(UpperBound() + ef.UpperBound());
            r.Check();
            return r;
        }

        EFloat operator-(EFloat ef) const
        {
            EFloat r;
            r.v = v - ef.v;
#ifndef NDEBUG
            r.vPrecise = vPrecise - ef.vPrecise;
#endif 
            r.low = NextFloatDown(LowerBound() - ef.UpperBound());
            r.high = NextFloatUp(UpperBound() - ef.LowerBound());
            r.Check();
            return r;
        }

        EFloat operator*(EFloat ef) const
        {
            EFloat r;
            r.v = v * ef.v;
#ifndef NDEBUG
            r.vPrecise = vPrecise * ef.vPrecise;
#endif
            float prod[4] = {
                LowerBound() * ef.LowerBound(), UpperBound() * ef.LowerBound(),
                LowerBound() * ef.UpperBound(), UpperBound() * ef.UpperBound()
            };
            r.low = NextFloatDown(
                (std::min)((std::min)(prod[0], prod[1]), (std::min)(prod[2], prod[3]))
            );
            r.high = NextFloatUp(
                (std::max)((std::max)(prod[0], prod[1]), (std::max)(prod[2], prod[3]))
            );
            r.Check();
            return r;
        }

        EFloat operator/(EFloat ef) const
        {
            EFloat r;
            r.v = v / ef.v;
#ifndef NDEBUG
            r.vPrecise = vPrecise / ef.vPrecise;
#endif
            if (ef.low < 0.f && ef.high > 0.f)
            {
                // Bah. (翻译成中文就是：我擦!) The interval we're dividing by straddles zero,
                // so just return an interval of everything.
                r.low = -Infinity;
                r.high = Infinity;
            }
            else 
            {
                float div[4] = {
                    LowerBound() / ef.LowerBound(), UpperBound() / ef.LowerBound(),
                    LowerBound() / ef.UpperBound(), UpperBound() / ef.UpperBound()
                };
                r.low = NextFloatDown(
                    (std::min)((std::min)(div[0], div[1]), (std::min)(div[2], div[3]))
                );
                r.high = NextFloatUp(
                    (std::max)((std::max)(div[0], div[1]), (std::max)(div[2], div[3]))
                );
            }

            r.Check();
            return r;
        }

        EFloat operator-() const
        {
            EFloat r;
            r.v = -v;
#ifndef NDEBUG
            r.vPrecise = -vPrecise;
#endif
            r.low = -high;
            r.high = -low;
            r.Check();
            return r;
        }

        inline bool operator==(EFloat ef) const {return v == ef.v;}
        inline bool operator!=(EFloat ef) const {return v != ef.v;}

        inline void Check() const {}

        friend std::ostream& operator<<(std::ostream& os, const EFloat& ef)
        {
            os << "v = " << ef.v << " - [" << ef.low << ", " << ef.high << "]\n";
            return os;
        }

    private:
        float v, low, high;
        #ifndef NDEBUG
        long double vPrecise;
        #endif
        friend inline EFloat Sqrt(EFloat fe);
        friend inline EFloat Abs(EFloat fe);
        friend inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat& t0, EFloat& t1);
    };

    // EFloat inline functions
    inline EFloat operator*(float f, EFloat ef)
    {
        return EFloat(f) * ef;
    }

    inline EFloat operator/(float f, EFloat ef)
    {
        return EFloat(f) * ef;
    }

    inline EFloat operator+(float f, EFloat ef)
    {
        return EFloat(f) + ef;
    }

    inline EFloat operator- (float f, EFloat ef)
    {
        return EFloat(f) - ef;
    }

    inline EFloat Sqrt(EFloat ef)
    {
        EFloat r;
        r.v = std::sqrtf(ef.v);
#ifndef NDEBUG
        r.vPrecise = std::sqrt(ef.vPrecise);
#endif
        r.low = NextFloatDown(std::sqrtf(ef.low));
        r.high = NextFloatUp(std::sqrtf(ef.high));
        r.Check();
        return r;
    }

    inline EFloat Abs(EFloat ef)
    {
        if (ef.low >= 0)
            // The entire interval is greater than zero, so we're all set.
            return ef;
        else if (ef.high <= 0)
        {
            // The entire interval is less than zero.
            EFloat r;
            r.v = -ef.v;
#ifndef NDEBUG
            r.vPrecise = -ef.vPrecise;
#endif
            r.low = -ef.high;
            r.high = -ef.low;
            r.Check();
            return r;
        }
        else 
        {
            // The interval straddles zero.
            EFloat r;
            r.v = Abs(ef.v);
#ifndef NDEBUG
            r.vPrecise = Abs(ef.vPrecise);
#endif
            r.low = 0;
            r.high = (std::max)(-ef.low, ef.high);
            r.Check();
            return r;
        }
    }

    inline bool Quadratic (EFloat A, EFloat B, EFloat C, EFloat& t0, EFloat& t1);
    inline bool Quadratic (EFloat A, EFloat B, EFloat C, EFloat& t0, EFloat& t1)
    {
        // Find quadratic discriminant
        double discrim = (double)B.v * (double)B.v - 4. * (double)A.v * (double)C.v;
        if (discrim < 0.) return false;
        double rootDiscrim = std::sqrt(discrim);

        EFloat floatRootDiscrim(rootDiscrim, MACHINE_EPSILON * rootDiscrim);

        // Compute quadratic _t_ values
        EFloat q;
        if ((float)B < 0.f)
            q = -.5f * (B - floatRootDiscrim);
        else 
            q = -.5f * (B + floatRootDiscrim);

        t0 = q / A;
        // In order to avoid cancellation, we use the fomulation of x1 * x2 = c / a to compute t1,
        // which could be found in "Accuracy and Stability of Numerical Algorithms" Section 1.8.
        t1 = C / q;
        if ((float)t0 > (float)t1) std::swap(t0, t1);
        return true;
    }
}