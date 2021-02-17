#pragma once

#include <functional>
#include <atomic>
#include "MathUtility.hpp"

namespace Panda
{
    class AtomicFloat
    {
        public:
            explicit AtomicFloat(float v = 0) { bits = FloatToBits(v);}
            operator float() const {return BitsToFloat(bits);}
            float operator=(float v)
            {
                bits = FloatToBits(v);
                return v;
            }
            void Add(float v)
            {
#ifdef PANDA_FLOAT_AS_DOUBLE
                uint64_t oldBits = bits, newBits;
#else
                uint32_t oldBits = bits, newBits;
#endif
                do{
                    newBits = FloatToBits(BitsToFloat(oldBits) + v);
                }while (!bits.compare_exchange_weak(oldBits, newBits));
            }

        private:
#ifdef PANDA_FLOAT_AS_DOUBLE
            std::atomic<uint64_t> bits;
#else
            std::atomic<uint32_t> bits;
#endif
    };
}
