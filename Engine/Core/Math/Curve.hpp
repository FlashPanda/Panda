#pragma once
#include <vector>
#include "Core/portable.hpp"
#include "MathUtility.hpp"

namespace Panda
{
    ENUM(CurveType)
    {
        kLinear = "LINE"_i32,
        kBezier = "BEZI"_i32
    };

    class CurveBase
    {
        private:
            const CurveType m_kCurveType;

        public:
            CurveBase() = delete;
            CurveBase(CurveType type) : m_kCurveType(type) {}
            virtual ~CurveBase() = default;
            CurveType GetCurveType() const {return m_kCurveType;}
    };

    template <typename VAL, typename PARAM>
    class Curve
    {
        protected:
            std::vector<VAL> m_Knots;

        public:
            virtual PARAM Reverse(VAL t, size_t& index) const = 0;
            virtual VAL Interpolate(PARAM t, const size_t index) const = 0;
            void AddKnot(const VAL knot)
            {
                m_Knots.push_back(knot);
            }
    };
}