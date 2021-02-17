#pragma once
#include "MathUtility.hpp"
#include "Vector.hpp"
#include <algorithm>

namespace Panda
{
    template <typename T, int M, int N>
    struct Matrix
    {
		union {
			Vector<T, N> v[M];
			T data[M * N] = { 0 };
			T m[M][N];
		};

        Matrix() {}
        Matrix(const T val)
        {
            size_t mn = M * N;
            for (size_t i = 0; i < mn; ++i)
                data[i] = val;
        }
        Matrix(const std::vector<T>& vec)
        {
            size_t st = vec.size();
            size_t mn = M * N;
            for (size_t i = 0; i < st && i < mn; ++i)
            {
                data[i] = vec[i];
            }
        }

		Matrix(const Matrix& rhs)
		{
			memcpy(data, rhs.data, sizeof(T) * M * N);
		}
        // Be careful when using this function.
        Matrix(const T* list)
        {
            memcpy_s(data, sizeof(T) * M * N, list, sizeof(T) * M * N);
        }

        Vector<T, N>& operator[](int rowIndex)
        {
            return v[rowIndex];
        }

        const Vector<T, N>& operator[] (int rowIndex) const 
        {
            return v[rowIndex];
        }

         operator T*()
        {
            return &data[0];
        }

         operator const T*() const
        {
            return &data[0];
        }

		Matrix& operator=(const Matrix& rhs)
		{
			memcpy(data, rhs.data, sizeof(T) * M * N);
			return *this;
		}

		void Set(const T val)
		{
			size_t mn = M * N;
			for (size_t i = 0; i < mn; ++i)
				data[i] = val;
		}

		void SetRow(const Vector<T, N>& row, int32_t index)
		{
			v[index] = row;
		}

		void SetCol(const Vector<T, M>& col, int32_t index)
		{
			for (int32_t i = 0; i < M; ++i)
			{
				m[i][index] = col.data[i];
			}
		}

		void SetIdentity()
		{
			Set(0);
			int32_t n = (std::min)(M, N);
			for (int32_t i = 0; i < n; ++i)
				m[i][i] = 1;
		}

        Vector<T, N> GetRow(int32_t row) const
        {
            return v[row];
        }

        Vector<T, M> GetCol(int32_t col) const
        {
            Vector<T, M> result;
            for (int32_t i = 0; i < M; ++i)
            {
                result.data[i] = m[i][col];
            }
            return result;
        }
    };

	typedef Matrix<float, 2, 2> Matrix2f;
    typedef Matrix<float, 3, 3> Matrix3f;
    typedef Matrix<float, 4, 4> Matrix4f;

    template <typename T, int M, int N>
    std::ostream& operator<<(std::ostream& out, const Matrix<T, M, N>& matrix)
    {
        out << std::endl;
        out << "{ ";
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                out << matrix.data[i * N + j] << ((j == N - 1)? (i != M - 1? "\n  " : "") : ", ");
            }
        }
        out << " }";
        out << std::endl;

        return out;
    }

    template <typename T, int M, int N>
    Matrix<T, M, N> operator+(const Matrix<T, M, N>& mat1, const Matrix<T, M, N>& mat2)
    {
        Matrix<T, M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i * N + j] = mat1.data[i * N + j] + mat2.data[i * N + j];
        return result;
    }

    template <typename T, int M, int N>
    Matrix<T, M, N> operator-(const Matrix<T, M, N>& mat1, const Matrix<T, M, N>& mat2)
    {
        Matrix<T, M, N> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i * N + j] = mat1.data[i * N + j] - mat2.data[i * N + j];
        return result;
    }

    template <typename T, int M, int N>
    Matrix<T, M, N> operator*(const Matrix<T, M, N>& mat, const float scaler)
    {
        Matrix<T, M, N> result;
        for(size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                result.data[i * N + j] = mat.data[i * N + j] * scaler;
        return result;
    }

    template <typename T, int M, int N>
    Matrix<T, M, N> operator*(const float scaler, const Matrix<T, M, N>& mat)
    {
        return mat * scaler;
    }

    template <typename T, int M, int N>
    Matrix<T, M, N> operator/(const Matrix<T, M, N>& mat, const float scaler)
    {
        float inv = 1.0f / scaler;
        return mat * inv;
    }

    template <typename T, int M, int N, int R>
    Matrix<T, M, R> operator*(const Matrix<T, M, N>& mat1, const Matrix<T, N, R>& mat2)
    {
        Matrix<T, M, R> result;
        for (size_t i = 0; i < M; ++i)
            for (size_t k = 0; k < R; ++k)
                for (size_t j = 0; j < N; ++j)
                    result.data[i * R + k] += mat1.data[i * N + j] * mat2.data[j * R + k];
        return result;
    }

	template <typename T, int M, int N>
	Vector<T, M> operator*(const Matrix<T, M, N>& mat, const Vector<T, N>& vec)
	{
		Vector<T, M> result;
		
		for (int32_t i = 0; i < M; ++i)
		{
			T sum = 0;
			for (int32_t j = 0; j < N; ++j)
			{
				sum += mat.m[i][j] * vec.data[j];
			}
			result.data[i] = sum;
		}

		return result;
	}

	template <typename T, int M, int N>
	Vector<T, N> operator*(const Vector<T, M>& vec, const Matrix<T, M, N>& mat)
	{
		Vector<T, N> result;
		for (int32_t j = 0; j < N; ++j)
		{
			T sum = 0;
			for (int32_t i = 0; i < M; ++i)
			{
				sum += mat[i][j] * vec.data[i];
			}
			result.data[j] = sum;
		}

		return result;
	}

	template<typename T, int M, int N>
	bool AlmostZero(const Matrix<T, M, N>& mat)
	{
		bool result = true;
		for (int32_t i = 0; i < M; ++i)
		{
			if (!result)
				break;

			for (int32_t j = 0; j < N; ++j)
			{
				if (Abs(mat.m[i][j]) > std::numeric_limits<T>::epsilon())
				{
					result = false;
					break;
				}
			}
		}

		return result;
	}

	template <typename T, int M, int N>
	bool operator==(const Matrix<T, M, N>& mat1, const Matrix<T, M, N>& mat2)
	{
		return AlmostZero(mat1 - mat2);
	}

	template <typename T, int M, int N>
	bool operator!=(const Matrix<T, M, N>& mat1, const Matrix<T, M, N>& mat2)
	{
		return !(mat1 == mat2);
	}

	template <typename T, int M, int N>
	void TransposeMatrix(const Matrix<T, M, N>& mat, Matrix<T, N, M>& out)
	{
		Matrix<T, M, N> temp(mat);
		for (int32_t i = 0; i < M; ++i)
		{
			for (int32_t j = 0; j < N; ++j)
			{
				out.m[i][j] = temp.m[j][i];
			}
		}
	}

	template <typename T, int N>
	T GetCofactor(const Matrix<T, N, N>& mat, int32_t i, int32_t j)
	{
		int32_t factor = 1;
		if ((i + j) % 2 == 1)
			factor = -1;
		
		Matrix<T, N - 1, N - 1> newMat;
		int32_t k = 0, l = 0;
		for (int32_t row = 0; row < N; ++row)
		{
			if (row == i)
				continue;
			l = 0;
			for (int32_t col = 0; col < N; ++col)
			{
				if (col == j)
					continue;

				newMat.m[k][l] = mat.m[row][col];
				++l;
			}
			++k;
		}

		return factor * GetDeterminant(newMat);
	}

	template <typename T, int N>
	T GetDeterminant(const Matrix<T, N, N>& mat)
	{
		T sum = 0;
		for (int32_t i = 0; i < N; ++i)
			sum += mat.m[0][i] * GetCofactor(mat, 0, i);
		return sum;
	}

	template <typename T>
	T GetDeterminant(const Matrix<T, 1, 1>& mat)
	{
		return mat.m[0][0];
	}

	template <typename T>
	T GetDeterminant(const Matrix<T, 2, 2>& mat)
	{
		return mat.m[0][0] * mat.m[1][1] - mat.m[0][1] * mat.m[1][0];
	}

	template <typename T, int N>
	bool GetAdjointMatrix(const Matrix<T, N, N>& mat, Matrix<T, N, N>& out)
	{
		Matrix<T, N, N> temp(mat);

		for (int32_t i = 0; i < N; ++i)
		{
			for (int32_t j = 0; j < N; ++j)
			{
				out.m[i][j] = GetCofactor(temp, i, j);
			}
		}

		TransposeMatrix(out, out);
		return true;
	}

	template <typename T, int N>
	bool InverseMatrix(const Matrix<T, N, N>& mat, Matrix<T, N, N>& out)
	{
		T det = GetDeterminant(mat);
		if (det == 0)
		{
			out.SetIdentity();
			return false;
		}

		Matrix<T, N, N> temp;
		GetAdjointMatrix(mat, temp);
		out = temp / det;
		
		return true;
	}
}