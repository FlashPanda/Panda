#pragma once
#include "MathUtility.hpp"
#include "Curve.hpp"
#include "PandaMath.hpp"
#include "Numerical.hpp"

namespace Panda
{
    template<typename VAL, typename PARAM>
    class Linear : public CurveBase, public Curve<VAL, PARAM>
    {
        public:
            Linear() : CurveBase(CurveType::kLinear) {}
            Linear(const std::vector<VAL> knots)
                : Linear()
            {
                Curve<VAL, PARAM>::m_Knots = knots;
            }

            Linear(const VAL* knots, const size_t count)
                : Linear()
            {
                for (size_t i = 0; i < count; ++i)
                {
                    Curve<VAL, PARAM>::m_Knots.push_back(knots[i]);
                }
            }

            PARAM Reverse(VAL t, size_t& index) const final
            {
                VAL t1((float)0), t2((float)0);

                if (Curve<VAL, PARAM>::m_Knots.size() < 2)
                    return PARAM((float)0);

                if (t <= Curve<VAL, PARAM>::m_Knots.front())
                {
                    index = 0;
                    return PARAM((float)0);
                }

                if (t >= Curve<VAL, PARAM>::m_Knots.back())
                {
                    index = Curve<VAL, PARAM>::m_Knots.size();
                    return PARAM((float)1);
                }

                for (size_t i = 1; i < Curve<VAL, PARAM>::m_Knots.size(); ++i)
                {
                    if (t >= Curve<VAL, PARAM>::m_Knots[i])
                        continue;
                    
                    t1 = Curve<VAL, PARAM>::m_Knots[i- 1];
                    t2 = Curve<VAL, PARAM>::m_Knots[i];
                    index = i;
                    break;
                }

                return PARAM((t - t1) / (t2 - t1));
            }

            VAL Interpolate(PARAM s, const size_t index) const final
            {
                if (Curve<VAL, PARAM>::m_Knots.size() == 0)
                    return VAL((float)0);
                else if (Curve<VAL, PARAM>::m_Knots.size() == 1)
                    return VAL((float)1);
                else if (Curve<VAL, PARAM>::m_Knots.size() < index + 1)
                    return Curve<VAL, PARAM>::m_Knots.back();
                else if (index == 0)
                    return Curve<VAL, PARAM>::m_Knots.front();
                else
                {
                    auto t1 = Curve<VAL, PARAM>::m_Knots[index - 1];
                    auto t2 = Curve<VAL, PARAM>::m_Knots[index];

                    return (PARAM(1) - s) * t1 + s * t2;
                }
            }
    };

    template<typename T>
    class Linear<Vector<T, 4>, T> : public CurveBase, public Curve<Vector<T, 4>, T>
    {
        public:
            Linear() : CurveBase(CurveType::kLinear) {}
            Linear(const std::vector<Vector<T, 4>> knots)
                : Linear()
            {
                m_Knots = knots;
            }

            Linear(const Vector<T, 4>* knots, const size_t count)
                : Linear()
            {
                for (size_t i = 0; i < count; ++i)
                    m_Knots.push_back(knots[i]);
            }

            T Reverse(Vector<T, 4> v, size_t& index) const final
            {
                T result = 0;
                assert(0);

                return result;
            }

			Vector<T, 4> Interpolate(T s, const size_t index) const final
            {
				Vector<T, 4> result;
                assert(0);

                return result;
            }
    };

    template<>
    class Linear<Matrix4f, float> : public CurveBase, public Curve<Matrix4f, float>
    {
        public:
            Linear() : CurveBase(CurveType::kLinear) {}
            Linear(const std::vector<Matrix4f> knots)
                : Linear()
            {
                m_Knots = knots;
            }

            Linear(const Matrix4f* knots, const size_t count)
                : Linear()
            {
                for (size_t i = 0; i < count; ++i)
                    m_Knots.push_back(knots[i]);
            }

            float Reverse(Matrix4f v, size_t& index) const final 
            {
                float result = 0.0f;
                assert(0);
                return result;
            }

            Matrix4f Interpolate(float s, const size_t index) const final
            {
                Matrix4f result;
                result.SetIdentity();
                if (m_Knots.size() == 0)
                    return result;
                else if (m_Knots.size() == 1)
                    return m_Knots[0];
                else if (m_Knots.size() < index + 1)
                    return m_Knots.back();
                else if (index == 0)
                    return m_Knots.front();
                else 
                {
                    auto v1 = m_Knots[index - 1];
                    auto v2 = m_Knots[index];

                    Vector3Df translation1, translation2;
                    Vector3Df scalar1, scalar2;
                    Vector3Df rotation1, rotation2;

                    MatrixDecomposition(v1, rotation1, scalar1, translation1);
                    MatrixDecomposition(v2, rotation2, scalar2, translation2);

                    // Interpolate translation
                    Vector3Df translation = (1.0f - s) * translation1 + s * translation2;
                    // Interpolate scalar
                    Vector3Df scalar = (1.0f - s) * scalar1 + s * scalar2;
                    // Interpolate rotation
                    Vector3Df rotation = (1.0f - s) * rotation1 + s * rotation2;

                    // compose the interpolated matrix
                    MatrixComposition(result, rotation, scalar, translation);
                }

                return result;
            }
    };
}