#pragma once
#include "MathUtility.hpp"
#include <vector>
#include <cassert>
#include <memory>
#include <unordered_set>

namespace Panda
{
    template <typename T, int N>
    struct Vector
    {
        T data[N] = {0};

        Vector() = default;
        Vector(const T val)
        {
            for (size_t i = 0; i < N; ++i)
                data[i] = val;
        }

        Vector(const std::vector<T>& vec)
        {
            size_t sv = vec.size();
            for (size_t i = 0; i < N && i < sv; ++i)
            {
                data[i] = vec[i];
            }
        }

        // Be careful when using this function.
        Vector(const T* list)
        {
            memcpy_s(data, sizeof(T) * N, list, sizeof(T) * N);
        }

        Vector(const Vector<T, N>& rhs)
        {
            memcpy_s(data, sizeof(T) * N, rhs.data, sizeof(T) * N);
        }

        operator T*() {
            return reinterpret_cast<T*>(this);
        }

        operator const T*() const 
        {
            return reinterpret_cast<const T*>(this);
        }

        T& operator[] (int32_t index)
        {
            return data[index];
        }

        const T operator[](int32_t index) const
        {
            return data[index];
        }

        void Set(const T val)
        {
            for (size_t i = 0; i < N; ++i)
                data[i] = val;
        }

        void Set(const T* list)
        {
            memcpy_s(data, sizeof(T) * N, list, sizeof(T) * N);
        }

        void Set(const std::vector<T>& vec)
        {
            size_t sv = vec.size();
            for (size_t i = 0; i < N && i < sv; ++i)
            {
                data[i] = vec[i];
            }
        }

        Vector& operator=(const T* list)
        {
            Set(list);
            return *this;
        }

        Vector& operator=(const std::vector<T>& vec)
        {
            Set(vec);
            return *this;
        }

        Vector& operator=(const Vector<T, N>& rhs)
        {
            memcpy_s(data, sizeof(T) * N, rhs.data, sizeof(T) * N);
            return *this;
        }

        Vector& operator+=(T scalar)
        {
            for (int32_t i = 0; i < N; ++i)
                data[i] += scalar;
            return *this;
        }

        Vector& operator+=(const Vector<T, N>& rhs)
        {
            for (int32_t i = 0; i < N; ++i)
                data[i] += rhs[i];
            return *this;
        }

        // negative
        Vector<T, N> operator-() const
        {
            Vector<T, N> result;
            for (size_t i = 0; i < N; ++i)
            {
                result.data[i] = -data[i];
            }
            return result;
        }

        // negative
        Vector<T, N> operator-()
        {
            Vector<T, N> result;
            for (size_t i = 0; i < N; ++i)
            {
                result.data[i] = -data[i];
            }
            return result;
        }

        Vector<T, N> operator*= (T scalar)
        {
            for (int32_t i = 0; i < N; ++i)
                data[i] *= scalar;
            return *this;
        }

        Vector<T, N> operator/= (T scalar)
        {
            T inv = (T)1.0 / scalar;
            for (int32_t i = 0; i < N; ++i)
                data[i] *= inv;
            return *this;
        }

		int32_t GetAbsMaxElementIndex()
		{
			int32_t index = 0;
			T value = 0;

			for (int32_t i = 0; i < N; ++i)
			{
				if (std::abs(data[i]) > value) {
					value = std::abs(data[i]);
					index = i;
				}
			}

			return index;
		}
    };

    template <typename T, int N>
    std::ostream& operator<<(std::ostream& out, const Vector<T, N>& vec)
    {
        out << "( ";
        for(size_t i = 0; i < N; ++i)
            out << vec.data[i] << ((i == N -1)? "": ", ");
        out << " )" << std::endl;

        return out;
    }

    template<typename T, int N>
    Vector<T, N> operator+(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec1.data[i] + vec2.data[i];
        }
        return result;
    }

    template <typename T, int N>
    Vector<T, N> operator+(const Vector<T, N>& vec, T scalar)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec.data[i] + scalar;
        }
        return result;
    }

    template <typename T, int N>
    Vector<T, N> operator+(T scalar, const Vector<T, N>& vec)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec.data[i] + scalar;
        }
        return result;
    }

    template <typename T, int N>
    Vector<T, N> operator-(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec1.data[i] - vec2.data[i];
        }
        return result;
    }

    template <typename T, int N>
    Vector<T, N> operator-(const Vector<T, N>& vec, T scalar)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec.data[i] - scalar;
        }
        return result;
    }

    template <typename T, int N>
    Vector<T, N> operator*(const Vector<T, N>& vec, const float scaler)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec.data[i] * scaler;
        }
        return result;
    }

    template <typename T, int N>
    Vector<T, N> operator*(const float scaler, const Vector<T, N>& vec)
    {
        Vector<T, N> result;
        for (size_t i = 0; i < N; ++i)
        {
            result.data[i] = vec.data[i] * scaler;
        }
        return result;
    }

	template<typename T, int N>
	Vector<T, N> operator*(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
	{
		Vector<T, N> result;
		for (size_t i = 0; i < N; ++i)
		{
			result.data[i] = vec1[i] * vec2[i];
		}
		return result;
	}

    template <typename T, int N>
    Vector<T, N> operator/(const Vector<T, N>& vec, const float scaler)
    {
        assert(scaler != 0.0f);
        float t = 1.0f / scaler;
        return vec * t;
    }

	template<typename T, int N>
	Vector<T, N> operator/(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
	{
		Vector<T, N> result;
		for (int32_t i = 0; i < N; ++i)
			result.data[i] = vec1[i] / vec2[i];
		return result;
	}

    template <typename T, int N>
    bool operator>(const Vector<T, N>& vec, const float scalar)
    {
        return GetLength(vec) > scalar;
    }

    template <typename T, int N>
    bool operator>=(const Vector<T, N>& vec, const float scalar)
    {
        return GetLength(vec) >= scalar;
    }

    template <typename T, int N>
    bool operator<(const Vector<T, N>& vec, const float scalar)
    {
        return GetLength(vec) < scalar;
    }

    template <typename T, int N>
    bool operator<=(const Vector<T, N>& vec, const float scalar)
    {
        return GetLength(vec) <= scalar;
    }

	template <typename T>
	T Pow(const T base, const float exponent)
	{
		return std::pow(base, exponent);
	}

	template<typename T, int N>
	Vector<T, N> Pow(const Vector<T, N>& vec, const float exponent)
	{
		Vector<T, N> result;

		for (int32_t i = 0; i < N; ++i)
		{
			result.data[i] = pow(vec[i], exponent);
		}

		return result;
	}

	template <typename T, int N>
	Vector<T, N> Abs(const Vector<T, N>& vec)
	{
		Vector<T, N> result;
		for (int32_t i = 0; i < N; ++i)
			result.data[i] = Abs(vec[i]);
		return result;
	}

    template <typename T, int N>
    bool AlmostZero(const Vector<T, N>& vec)
    {
        bool result = true;
        for (int32_t i = 0; i < N; ++i)
        {
            if (Abs(vec[i]) > std::numeric_limits<T>::epsilon())
            {
                result = false;
                break;
            }
        }

        return result;
    }

    template <typename T, int N>
    bool operator==(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
    {
        return AlmostZero(vec1 - vec2);
    }

    template <typename T, int N>
    bool operator!= (const Vector<T, N>& vec1, const Vector<T, N>& vec2)
    {
        return !(vec1 == vec2);
    }

	template <typename T, int N>
	T DotProduct(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
	{
		T result = 0.0;
		for (size_t i = 0; i < N; ++i)
			result += vec1.data[i] * vec2.data[i];
		return result;
	}

    template <typename T, int N>
    T AbsDotProduct(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
    {
        T result = 0.0;
		for (size_t i = 0; i < N; ++i)
			result += vec1.data[i] * vec2.data[i];
		return std::abs(result);
    }

    template <typename T>
    Vector<T, 3> CrossProduct(const Vector<T, 3>& vec1, const Vector<T, 3>& vec2)
    {
        Vector<T, 3> result;
        
		result.Set({ vec1[1] * vec2[2] - vec1[2] * vec2[1],
			vec1[2] * vec2[0] - vec1[0] * vec2[2],
			vec1[0] * vec2[1] - vec1[1] * vec2[0] });

        return result;
    }
    template <typename T>
    T CrossProduct(const Vector<T, 2>& vec1, const Vector<T, 2>& vec2)
    {
        return vec1.data[0] * vec2.data[1] - vec1.data[1] * vec2.data[0];
    }

    template <typename T, int N>
    Vector<T, N> MulByElement(const Vector<T, N>& vec1, const Vector<T, N>& vec2)
    {
        Vector<T, N> result;
        for (int32_t i = 0; i < N; ++i)
            result.data[i] = vec1[i] * vec2[i];
        return result;
    }

	template <typename T, int N>
	Vector<T, N> Normalize(const Vector<T, N>& vec)
	{
		Vector<T, N> result(vec);
		result = result / GetLength(result);
		return result;
	}

    template <typename T, int N>
    T GetLength(const Vector<T, N>& vec)
    {
        T sum = 0;
        for (int32_t i = 0; i < N; ++i)
        {
            sum += vec[i] * vec[i];
        }
        return std::sqrt(sum);
    }

    template <typename T, int N>
    T GetLengthSquare(const Vector<T, N>& vec)
    {
        T sum = 0;
        for (int32_t i = 0; i < N; ++i)
        {
            sum += vec[i] * vec[i];
        }
        return sum;
    }

    template <typename T, int N>
    inline Vector<T, N> Faceforward(const Vector<T, N>& v1, const Vector<T, N>& v2)
    {
        return (DotProduct(v1, v2) < 0.f)? -v1 : v1;
    }

    template <typename T>
    int32_t MaxDimension(const Vector<T, 3>& v)
    {
        return (v.data[0] > v.data[1])? ((v.data[0] > v.data[2])? 0 : 2) : ((v.data[1] > v.data[2])? 1 : 2);
    }

    template <typename T>
    int32_t MaxAbsDimension(const Vector<T, 3>& v)
    {
        Vector<T, 3> t = Abs(v);
        return (t.data[0] > t.data[1])? ((t.data[0] > t.data[2])? 0 : 2) : ((t.data[1] > t.data[2])? 1 : 2);
    }

    template <typename T>
    T MaxComponent(const Vector<T, 3>& v)
    {
        return (std::max)(v[0], (std::max)(v[1], v[2]));
    }

    template <typename T>
    T MaxAbsComponent(const Vector<T, 3>& v)
    {
        Vector<T, 3> va = Abs(v);
        return (std::max)(va[0], (std::max)(va[1], va[2]));
    }

    template <typename T>
    Vector<T, 3> Permute(const Vector<T, 3>& v, int32_t x, int32_t y, int32_t z)
    {
        return Vector<T, 3>({v[x], v[y], v[z]});
    }

    // TO BE UNDERSTOOD
    template <typename T>
    void CoordinateSystem(const Vector<T, 3>& v1, Vector<T, 3>& v2,
                          Vector<T, 3>& v3)
    {
        if (Abs(v1[0] > Abs(v1[1])))
            v2 = Vector<T, 3>({-v1[2], 0, v1[0]}) / std::sqrtf(v1[0] * v1[0] + v1[2] * v1[2]);
        else
            v2 = Vector<T, 3>({0, v1[2], v1[1]}) / std::sqrtf(v1[1] * v1[1] + v1[2] * v1[2]);
        v3 = CrossProduct(v1, v2);
    }

    template <typename T, int N>
    Vector<T, N> Lerp(float t, const Vector<T, N>& p0, const Vector<T, N>& p1)
    {
        return (1 - t) * p0 + t * p1;
    }
    
    typedef Vector<float, 2> Vector2Df;
    typedef Vector<int32_t, 2> Vector2Di;
    typedef Vector<int16_t, 2> Pixel2D;
    
    typedef Vector<float, 3> Vector3Df;
    typedef Vector<double, 3> Vector3Dd;
    typedef Vector<int32_t, 3> Vector3Di;

    typedef Vector<float, 4> Vector4Df;
    typedef Vector<float, 4> Quaternion;
    typedef Vector<int32_t, 4> Vector4Di;
	typedef Vector<uint8_t, 4> R8G8B8A8Unorm;

    typedef Vector<float, 2> Point2Df;
    typedef Vector<int32_t, 2> Point2Di;
    typedef std::shared_ptr<Point2Df> Point2DPtr;
    typedef std::vector<Point2DPtr> Point2DList;
    typedef Vector<float, 3> Point3Df;
    typedef Vector<int32_t, 3> Point3Di;
    typedef std::shared_ptr<Point3Df> PointPtr;
    typedef std::unordered_set<PointPtr> PointSet;
    typedef std::vector<PointPtr> PointList;
    typedef std::pair<PointPtr, PointPtr> Edge;
}