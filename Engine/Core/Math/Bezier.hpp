#pragma once
#include "MathUtility.hpp"
#include "Curve.hpp"
#include "Numerical.hpp"
#include <cassert>
#include <map>

namespace Panda
{
    template<typename VAL, typename PARAM>
    class Bezier : public CurveBase, public Curve<VAL, PARAM>
    {
        private:
            std::map<VAL, VAL> m_IncomingControlPoints;
            std::map<VAL, VAL> m_OutgoingControlPoints;

        public:
            Bezier() :CurveBase(CurveType::kBezier) {}
			Bezier(const std::vector<VAL> knots, const std::vector<VAL> incoming_cp, const std::vector<VAL> outgoint_cp)
				:Bezier()
            {
                assert(knots.size() == incoming_cp.size() && knots.size() == outgoint_cp.size());

				Curve<VAL, PARAM>::m_Knots = knots;

                size_t count = knots.size();

                for (size_t i = 0; i < count; ++i)
                {
					VAL knot = knots[i];
					VAL in_cp = incoming_cp[i];
					VAL out_cp = outgoint_cp[i];
                    AddControlPoints(knot, in_cp, out_cp);
                }
            }

			Bezier(const VAL* knots, const VAL* incoming_cp, const VAL* outgoing_cp, const size_t count)
				:Bezier()
            {
                for (size_t i = 0; i < count; ++i)
                {
					Curve<VAL, PARAM>::m_Knots.push_back(knots[i]);
                    AddControlPoints(knots[i], incoming_cp[i], outgoing_cp[i]);
                }
            }

            void AddControlPoints(const VAL knot, const VAL incoming_cp, const VAL outgoing_cp)
            {
                //assert(incoming_cp <= knot && knot <= outgoing_cp);
                m_IncomingControlPoints.insert({knot, incoming_cp});
                m_OutgoingControlPoints.insert({knot, outgoing_cp});
            }

            PARAM Reverse (const VAL t, size_t& index) const final
            {
                VAL t1, t2;

                if (Curve<VAL, PARAM>::m_Knots.size() < 2)
                {
                    return PARAM((float)0);
                }
                    
                
                if (t <= Curve<VAL, PARAM>::m_Knots.front())
                {
                    index = 0;
                    return PARAM((float)0);
                }

                if (t >= Curve<VAL, PARAM>::m_Knots.back())
                {   
                    index = Curve<VAL, PARAM>::m_Knots.size();
                    return PARAM(1);
                }
                
                for (size_t i = 0; i < Curve<VAL, PARAM>::m_Knots.size(); ++i)
                {
                    if (t >= Curve<VAL, PARAM>::m_Knots[i])
                        continue;

                    t1 = Curve<VAL, PARAM>::m_Knots[i - 1];
                    t2 = Curve<VAL, PARAM>::m_Knots[i];
                    index = i;
                    break;
                }

                VAL c1, c2;
                c1 = m_OutgoingControlPoints.find(t1)->second;
                c2 = m_IncomingControlPoints.find(t2)->second;

                typename NewtonRaphson<VAL, PARAM>::nr_f f = [t2, t1, c2, c1, t](PARAM s){
                    return (t2 - 3.0f * c2 + 3.0f * c1 - t1) * pow(s, 3.0f)
                        + 3.0f * (c2 - 2.0f * c1 + t1) * pow(s, 2.0f)
                        + 3.0f * (c1 - t1) * s
                        + t1 - t;
                };

                typename NewtonRaphson<VAL, PARAM>::nr_fprime fprime = [t2, t1, c2, c1](PARAM s)
                {
                    return 3.0f * (t2 - 3.0f * c2 + 3.0f * c1 - t1) * pow(s, 2.0f)
                        + 6.0f * (c2 - 2.0f * c1 + t1) * s
                        + 3.0f * (c1 - t1);
                };

                return NewtonRaphson<VAL, PARAM>::Solve(0.5f, f, fprime);
            }

            VAL Interpolate(PARAM s, const size_t index) const final
            {
                if (Curve<VAL, PARAM>::m_Knots.size() == 0)
                    return VAL((float)0);
                else if (Curve<VAL, PARAM>::m_Knots.size() == 1)
                    return Curve<VAL, PARAM>::m_Knots.front();
                else if (Curve<VAL, PARAM>::m_Knots.size() < index + 1)
                    return Curve<VAL, PARAM>::m_Knots.back();
                else if (index == 0)
                    return Curve<VAL, PARAM>::m_Knots.front();
                else
                {
                    VAL t1, t2;
                    t1 = Curve<VAL, PARAM>::m_Knots[index - 1];
                    t2 = Curve<VAL, PARAM>::m_Knots[index];

                    VAL c1, c2;
                    c1 = m_OutgoingControlPoints.find(t1)->second;
                    c2 = m_IncomingControlPoints.find(t2)->second;

					VAL result = (t2 - 3.0f * c2 + 3.0f * c1 - t1) * pow(s, 3.0f)
						+ 3.0f * (c2 - 2.0f * c1 + t1) * pow(s, 2.0f)
						+ 3.0f * (c1 - t1) * s
						+ t1;
					return result;
                }
            }
    };

    template<typename T>
    class Bezier<Vector<T, 4>, T> : public CurveBase, public Curve<Vector<T, 4>, T>
    {
        private:
            //std::map<Vector<T, 4>, Vector<T, 4>> m_IncomingControlPoints;
            //std::map<Vector<T, 4>, Vector<T, 4>> m_OutgoingControlPoints;

        public:
            Bezier() : CurveBase(CurveType::kBezier) {}
            Bezier( const std::vector<Vector<T, 4>> knots, const std::vector<Vector<T, 4>> incoming_cp,
                    const std::vector<Vector<T, 4>> outgoing_cp)
                : Bezier()
            {
                assert(knots.size() == incoming_cp.size() && knots.size() == outgoing_cp.size());

                Curve<Vector<T, 4>, T>::m_Knots = knots;

                size_t count = knots.size();

                for (size_t i = 0; i < count; ++i)
                {
                    auto knot = knots[i];
                    auto in_cp = incoming_cp[i];
                    auto out_cp = outgoing_cp[i];
                    AddControlPoints(knot, in_cp, out_cp);
                }
            }

            Bezier(const Vector<T, 4>* knots, const Vector<T, 4>* incoming_cp, const Vector<T, 4>* outgoing_cp, const size_t count)
                : Bezier()
            {
                for (size_t i = 0; i < count; ++i)
                {
                    Curve<Vector<T, 4>, T>::m_Knots.push_back(knots[i]);
                    AddControlPoints(knots[i], incoming_cp[i], outgoing_cp[i]);
                }
            }

            void AddControlPoints(const Vector<T, 4>& knot, const Vector<T, 4>& incoming_cp, const Vector<T, 4>& outgoing_cp)
            {
                //assert(incoming_cp <= knot && knot << outgoing_cp);
                //m_IncomingControlPoints.insert({knot, incoming_cp});
                //m_OutgoingControlPoints.insert({knot, outgoing_cp});
            }

            T Reverse(Vector<T, 4> t, size_t& index) const final
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
    class Bezier<Matrix4f, float> : public CurveBase, public Curve<Matrix4f, float>
    {
        private:
            //std::map<Matrix4f, Matrix4f> m_IncomingControlPoints;
            //std::map<Matrix4f, Matrix4f> m_OutgoingControlPoints;

        public:
            Bezier() : CurveBase(CurveType::kBezier) {}
            Bezier(const std::vector<Matrix4f>& knots, const std::vector<Matrix4f>& incoming_cp, const std::vector<Matrix4f>& outgoing_cp)
                : Bezier()
            {
                assert(knots.size() == incoming_cp.size() && knots.size() == outgoing_cp.size());

                Curve<Matrix4f, float>::m_Knots = knots;

                size_t count = knots.size();

                for (size_t i = 0; i < count; ++i)
                {
                    Matrix4f knot = knots[i];
                    Matrix4f in_cp = incoming_cp[i];
                    Matrix4f out_cp = outgoing_cp[i];
                    AddControlPoints(knot, in_cp, out_cp);
                }
            }

            Bezier(const Matrix4f* knots, const Matrix4f* incoming_cp, const Matrix4f* outgoing_cp, const size_t count)
                : Bezier()
            {
                for (size_t i = 0; i < count; ++i)
                {
                    Curve<Matrix4f, float>::m_Knots.push_back(knots[i]);
                    AddControlPoints(knots[i], incoming_cp[i], outgoing_cp[i]);
                }
            }

            void AddControlPoints(const Matrix4f& knot, const Matrix4f& incoming_cp, const Matrix4f& outgoing_cp)
            {
                //assert(incoming_cp <= knot && knot <= outgoing_up);
                //m_IncomingControlPoints.insert({knot, incoming_cp});
                //m_OutgoingControlPoints.insert({knot, outgoing_cp});
            }

            float Reverse(Matrix4f t, size_t& index) const final
            {
                float result = 0.0f;
                assert(0);

                return result;
            }

            Matrix4f Interpolate(float s, const size_t index) const final
            {
                Matrix4f result;
                result.SetIdentity();

                assert(0);
                return result;
            }
    };
}