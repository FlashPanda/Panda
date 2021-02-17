#pragma once
#include "MathUtility.hpp"
#include "Matrix.hpp"

namespace Panda
{
	/**
	 * Gaussian method: change matrix to row echelon form.
	 */
    template<typename T, int M, int N>
    Matrix<T, M, N> GaussianElimination(const Matrix<T, M, N>& mat)
    {
		Matrix<T, M, N> result(mat);
		// Find the row to start with.
		size_t dealRow = 0;
		size_t dealCol = 0;
		while (dealRow <= M)
		{
			// The rows which row indices bigger than column count will be simply set to 0.
			if (dealCol > N)
			{
				break;
			}

			if (!ChangeRows(result, dealRow, dealCol))
				break;

			// Make the current row pivot to be 1.
			if (result.data[dealRow * N + dealCol] != 1.0f)
			{
				T t = 1.0f / result.data[dealRow * N + dealCol];
				for (size_t col = dealCol; col < N; ++col)
				{
					result.data[dealRow * N + col] *= t;
				}
			}

			if (dealRow == M - 1)
				break; // Finished

			// Make the left rows data[row * N + dealCol] to be zero.
			size_t nextRow = dealRow + 1;
			while (nextRow < M)
			{
				if (result.data[nextRow * N + dealCol] == 0.0f)
				{
					nextRow++;
					continue;
				}
				else
				{
					T t = result.data[nextRow * N + dealCol];
					for (size_t col = dealCol; col < N; ++col)
					{
						result.data[nextRow * N + col] -= t * result.data[dealRow * N + col];
					}
					nextRow++;
				}
			}
			
			dealRow++;
			dealCol++;
		}
			
		return result;
    }

    /**
	 * Gaussian method: change matrix to row echelon form.
	 */
    template<typename T, int M, int N>
    Matrix<T, M, N> GaussianEliminationImprove(const Matrix<T, M, N>& mat)
    {
        Matrix<T, M, N> result(mat);

        std::vector<int32_t> rowOrder(M);
        for(int32_t i = 0; i < M; ++i)
            rowOrder[i] = i;
        
        int32_t pivotIndex = 0;
        int32_t startOrderIndex = 0;
        int32_t rowIndex = 0;
        int32_t col = 0;
        bool cr = ChoosePivot(result, pivotIndex, rowOrder, startOrderIndex, rowIndex, col);
        while(cr)
        {
            if(rowIndex != startOrderIndex)
                std::swap(rowOrder[startOrderIndex], rowOrder[rowIndex]);

			// Make the current row pivot to be 1.
			if (result.m[rowOrder[startOrderIndex]][col] != 1.0f)
			{
				T t = 1.0f / result.m[rowOrder[startOrderIndex]][col];
				for (size_t c = col; c < N; ++c)
				{
					result.m[rowOrder[startOrderIndex]][c] *= t;
				}
			}

            // startOrderIndex and col
            for (int32_t i = startOrderIndex + 1; i < rowOrder.size(); ++i)
            {
                if (result.m[rowOrder[i]][col] == 0.0f)
				{
					continue;
				}
				else
				{
					T t = result.m[rowOrder[i]][col];
					for (size_t c = col; c < N; ++c)
					{
						result.m[rowOrder[i]][c] -= t * result.m[rowOrder[startOrderIndex]][c];
					}
				}
            }

            startOrderIndex++;
            pivotIndex++;
            cr = ChoosePivot(result, pivotIndex, rowOrder, startOrderIndex, rowIndex, col);
        }

		Matrix<T, M, N> r(result);
		for (int32_t i = 0; i < rowOrder.size(); ++i)
		{
			result.v[i] = r.v[rowOrder[i]];
		}
            
        return result;
    }

	/**
	 * Gauss-Jordan Reduction
     * Note: change matrix to reduced row echelon form.
	 */
    template<typename T, int M, int N>
    Matrix<T, M, N> GaussJordanReduction(const Matrix<T, M, N>& mat)
    {
		Matrix<T, M, N> result = GaussianElimination(mat);

		// find the last pivot
		size_t lRow = 0;
		size_t lCol = 0;
		for (size_t i = 0; i < M; ++i)
		{
			size_t col = 0;
			if (!FindPivot(result, i, col))
				break;
			else
			{
				if (i > lRow)
				{
					lRow = i;
					lCol = col;
				}
			}
		}

		// 
		for (size_t i = lRow; i > 0; --i)
		{
			assert(FindPivot(result, i, lCol));
			int prevRow = i - 1;
			while (prevRow >= 0)
			{
				if (result.m[prevRow][lCol] == 0.0f)
				{
					prevRow--;
					continue;
				}
				else
				{
					T t = result.data[prevRow * N + lCol];
					for (size_t col = lCol; col < N; ++col)
					{
						result.data[prevRow * N + col] -= t * result.data[i * N + col];
					}
					prevRow--;
				}
			}
		}

		return result;
    }

	template <typename T, int M, int N>
	bool ChangeRows(Matrix<T, M, N>& mat, size_t row, size_t& col)
	{
		if (row >= N)
			return false;

		for (size_t j = col; j < N; ++j)
		{
			if (mat.data[row * N + j] != 0.0f) // check row element if it is 0
			{
				col = j;
				return true;
			}

			for (size_t i = row; i < M; ++i)
			{
				if (mat.data[i * N + j] != 0.0f)
				{
					// change the two rows
					T temp[N] = { 0.0f };
					memcpy(temp, &mat.data[row * N], sizeof(T) * N);
					memcpy(&mat.data[row * N], &mat.data[i * N], sizeof(T) * N);
					memcpy(&mat.data[i * N], temp, sizeof(T) * N);
					col = j;
					return true;
				}
			}
		}

		return false; // no need to change
	}

	template <typename T, int M, int N>
	bool FindPivot(Matrix<T, M, N>& mat, size_t row, size_t& col)
	{
		if (col > N)
			return false;

		for (size_t i = 0; i < N; ++i)
			if (mat.m[row][i] != 0)
			{
				col = i;
				return true;
			}

		return false;
	}

    /**
     * Choose a profit pivot of the matrix.
     * input parameters:
     *  mat --- the matrix to choose from
     *  pivotIndex --- which pivot to choose
     *  rowOrder --- the row order of the matrix
     *  startIndex --- form which index in rowOrder to choose
     *  rowIndex --- [out param] current pivot row index in rowOrder
     *  col --- [out param] current pivot col
     */
    template <typename T, int M, int N>
    bool ChoosePivot(const Matrix<T, M, N>& mat, int32_t pivotIndex, std::vector<int32_t>& rowOrder, int32_t startIndex, int32_t& rowIndex, int32_t& col)
    {
        // failure case 1: the pivot index is bigger than row count.
        if (pivotIndex >= M)
            return false;
        // failure case 2: the start index is bigger than rowOrder size.
        if (startIndex >= rowOrder.size())
            return false;

        rowIndex = -1;
        col = -1;
        T maxPivot = 0.0f;
        int32_t startCol = pivotIndex;
        while (maxPivot == 0.0f)
        {
            for (int32_t i = startIndex; i < rowOrder.size(); ++i)
            {
                T curPivot = mat.m[rowOrder[i]][startCol];
                if (std::abs(curPivot) > std::abs(maxPivot))
                {
                    maxPivot = curPivot;
                    rowIndex = i;
                    col = startCol;
                }
            }

            if (maxPivot != 0.0f)
                return true;

            // There is no non-zero element in this column. Find in next column.
            ++startCol;
            if (startCol >= N)
                break;
        }

        // failure case 3: there is no non-zero element in the left rows.
        return false;
    }

	/**
	 * LU factorization: make a n by n matrix to a L matrix(which is an unit lower trianglular matrix)
	 *					 multipled by a U matrix(which is a upper triangular).
	 */
	template <typename T, int N>
	bool LUFactorization(const Matrix<T, N, N>& mat, Matrix<T, N, N>& L, Matrix<T, N, N>& U)
	{
		Matrix<T, N, N> result(mat);
		L.SetIdentity();
		U.Set(0);

		int32_t row = 0;
		for (; row < N; ++row)
		{
			// The pivot of the row must NOT be 0.
			if (result.m[row][row] == 0)
				return false;

			for (int32_t i = row + 1; i < N; ++i)
			{
				T l = result.m[i][row] / result.m[row][row];
				L.m[i][row] = l;
				for (int32_t j = 0; j < N; ++j)
				{
					result.m[i][j] -= result.m[row][j] * l;
				}
			}
		}
		U = result;

		return true;
	}

    /***
     * Matrix norms induced by vector norms
     */
    template <typename T, int M, int N>
    T MatrixNormByVectorNorm(const Matrix<T, M, N>& mat, uint32_t factor)
    {
        // Do not support other norms.
        if (factor != 1 && factor != 0xFFFFFFFF)
        {
            return 0;
        }
        
        if (factor == 1)
        {
			T maxSum = 0;
            for (int32_t col = 0; col < N; ++col)
            {
				T sum = 0;
                for (int32_t row = 0; row < M; ++row)
                {
					sum += std::abs(mat.m[row][col]);
                }
                if (sum > maxSum)
                {
                    maxSum = sum;
                }
            }
            return maxSum;
        }

        if (factor == 0xFFFFFFFF)
        {
            T maxSum = 0;
            for (int32_t row = 0; row < M; ++row)
            {
                T sum = 0;
                for (int32_t col = 0; col < N; ++col)
                {
                    sum += std::abs(mat.m[row][col]);
                }
                if (sum > maxSum)
                {
                    maxSum = sum;
                }
            }
            return maxSum;
        }

        return 0;
    }

    /**
     * Matrix norms induced by entrywise
     */
    template <typename T, int M, int N>
    T MatrixNormByEntrywise(const Matrix<T, M, N>& mat, uint32_t factor)
    {
        if(factor == 1)
        {
            T norm = 0;
            for (int32_t i = 0; i < M; ++i)
            {
                for (int32_t j = 0; j < N; ++j)
                {
                    norm += std::abs(mat.m[i][j]);
                }
            }
            return norm;
        }

        if (factor == 0xFFFFFFFF)
        {
            T norm = 0;
            for (int32_t i = 0; i < M; ++i)
            {
                for (int32_t j = 0; j < N; ++j)
                {
                    if (std::abs(mat.m[i][j]) > norm)
                    {
                        norm = std::abs(mat.m[i][j]);
                    }
                }
            }
            return norm;
        }

        T norm = 0;
        for (int32_t i = 0; i < M; ++i)
        {
            for (int32_t j = 0; j < N; ++j)
            {
                norm += std::pow(std::abs(mat.m[i][j]), factor);
            }
        }
        norm = std::pow(norm, 1.0f / factor);

        return norm;
    }

    /**
     * Specially, there is a kind of norm called Frobenius Norm.
     * Which is the special case for enterwise norm with factor is 2.
     */
    template <typename T, int M, int N>
    T FrobeniusNorm(const Matrix<T, M, N>& mat)
    {
        return MatrixNormByEntrywise(mat, 2);
    }

    /// QR Method

    /**
     * The Householder transformation.
     */
    template <typename T, int N>
    bool HouseholderTransformation(const Vector<T, N>& vec, Matrix<T, N, N>& H)
    {
        T vecLength = vec.GetLength();

        T alpha = 0;
        if (vec.data[0] > 0)
            alpha = -vecLength;
        else 
            alpha = vecLength;
        
        Vector<T, N> v(vec);
        v.data[0] = vec.data[0] - alpha;

        Matrix<T, N, N> identity;
        identity.SetIdentity();
        
        Matrix<T, N, N> temp;
        T beta = alpha * (alpha - vec.data[0]);
        T beta_1 = 1.0 / beta;
        for (int32_t i = 0; i < N; ++i)
        {
            for (int32_t j = 0; j < N; ++j)
            {
                temp.m[i][j] = v.data[i] * v.data[j] * beta_1;
            }
        }

        H = identity - temp;

        return true;
    }

    template <typename T, int N>
    bool HouseholderTransformation(const Vector<T, N>& vec, float alpha, Matrix<T, N, N>& H)
    {
        Vector<T, N> v(vec);
        v.data[0] = vec.data[0] - alpha;

        Matrix<T, N, N> identity;
        identity.SetIdentity();
        
        Matrix<T, N, N> temp;
        T beta = alpha * (alpha - vec.data[0]);
        T beta_1 = 1.0 / beta;
        for (int32_t i = 0; i < N; ++i)
        {
            for (int32_t j = 0; j < N; ++j)
            {
                temp.m[i][j] = v.data[i] * v.data[j] * beta_1;
            }
        }

        H = identity - temp;

        return true;
    }

    /**
     * QR method with Householder
     */
    template <typename T, int N>
    bool QRMethodWithHouseholder(const Matrix<T, N, N>& mat, Matrix<T, N, N>& Q, Matrix<T, N, N>& R)
    {
        //Matrix<T, N, N> result(mat);

        //std::vector<Matrix<T, N, N>> Hs;
        //for (int32_t i = 0; i < N; ++i)
        //{
        //    Matrix<T, N - i, N - i> temp;
        //    BuildDownSquareMatrix(result, temp);
        //    Vector<T, N - i> vec(temp.GetCol(0));
        //    Matrix<T, N - i, N - i> H;
        //    HouseholderTransformation(vec, H);

        //    Matrix<T, N, N> wholeH;
        //    BuildWholeMatrixH(H, wholeH);

        //    Hs.push_back(wholeH);
        //}

        //Q.SetIdentity();    
        //for (int32_t i = Hs.size(); i >= 0; ++i)
        //{
        //    Q = Q * Hs[i];
        //}

        //R = result;

        return true;
    }

    template <typename T, int N, int X>
    void BuildDownSquareMatrix(const Matrix<T, N, N>& mat, Matrix<T, X, X>& result)
    {
        int32_t diff = N - X;
        for (int32_t i = 0; i < X; ++i)
        {
            for (int32_t j = 0; j < X; ++j)
            {
                result.m[i][j] = mat.m[i + diff][j + diff];
            }
        }
    }

    template <typename T, int X, int N>
    void BuildWholeMatrixH(const Matrix<T, X, X>& H, Matrix<T, N, N>& wholeH)
    {
        int32_t count = N - X;

        wholeH.SetIdentity();
        for (int32_t i = count; i < N; ++i)
        {
            for (int32_t j = count; j < N; ++j)
            {
                wholeH.m[i][j] = H.m[i - count][j - count];
            }
        }
    }

	/**
	 * mat must be a matrix with rank N.
	 */
	template<typename T, int M, int N>
	bool GramSchmidt(const Matrix<T, M, N>& mat, Matrix<T, M, N>& Q, Matrix<T, N, N>& R)
	{
		Matrix<T, M, N> result(mat);
		Q.Set(0);
		R.Set(0);
		for (int32_t i = 0; i < N; ++i)
		{
			Vector<T, M> col(result.GetCol(i));
			R.m[i][i] = col.GetLength();
			Vector<T, M> q = 1 / R.m[i][i] * col;
			Q.SetCol(q, i);
			for (int32_t j = i + 1; j < N; ++j)
			{
				R.m[i][j] = DotProduct(q,result.GetCol(j));
				Vector<T, M> a(result.GetCol(j));
				a = a - R.m[i][j] * q;
				result.SetCol(a, j);
			}
		}

		return true;
	}

	template <typename T, int M, int N>
	bool LeastSequaresMethod(const Matrix<T, M, N>& Q, const Matrix<T, N, N>& R, const Vector<T, M>& b, Vector<T, N>& x)
	{
		Vector<T, M> b1(b);
		Vector<T, N> c;
		c.Set(0.0f);
		for (int32_t i = 0; i < N; ++i)
		{
			Vector<T, M> col(Q.GetCol(i));
			c.data[i] = DotProduct(col, b1);
			Vector<T, M> col1(c.data[i] * col);
			b1 = b1 - col1;// c.data[i] * col;
		}

		// back substitution process
		for (int32_t i = N - 1; i >= 0; --i)
		{
			T sum = 0;
			for (int32_t j = i + 1; j < N; j++)
			{
				sum += R.m[i][j] * x.data[j];
			}
			x.data[i] = (c.data[i] - sum) / R.m[i][i];
		}

		return true;
	}

    template <typename T, int M, int N>
    void SolveEigenWithPowerMethod(const Matrix<T, M, N>& mat, int32_t iterCount, T& eigenValue, Vector<T, M>& eigenVector)
    {
        //int32_t count = iterCount;
        //const Matrix<T, M, N> tempMat(mat);

        //Vector<T, M> y(tempMat.GetCol(0));
        //Vector<T, M> t;
        //T lambda = 0;
        //
        //while (count > 0)
        //{
        //    y = Normalize(y);
        //    t = tempMat * y;
        //    lambda = DotProduct(t, y);
        //    y = t;

        //    --count;
        //}

        //eigenValue = lambda;
        //eigenVector = Normalize(y);

        /* another algorithm */
		Vector<T, M> U(mat.GetCol(0));
		Vector<T, M> V;
		int32_t count = iterCount;

		while (count > 0)
		{
			V = mat * U;
			int32_t index = V.GetAbsMaxElementIndex();
			U = 1.0f / V.data[index] * V;

			--count;
		}

		eigenVector = Normalize(U);
		eigenValue = V.data[V.GetAbsMaxElementIndex()];
    }

    template <typename T, int N>
    void BuildGivensRotation(int32_t i, int32_t j, T theta, Matrix<T, N, N>& mat)
    {
        T c = cos(theta);
        T s = sin(theta);
        mat.SetIdentity();
        mat.m[i][i] = c;
        mat.m[j][j] = c;
        mat.m[i][j] = -s;
        mat.m[j][i] = s;
    }

	template <typename T, int N>
	bool BuildGivensRotation(const Matrix<T, N, N>& inMat, int32_t i, int32_t j, Matrix<T, N, N>& mat)
	{
		assert(i > j);
		if (inMat.m[i][j] == 0)
		{
			mat.SetIdentity();
			return false;
		}
		T a = inMat.m[j][j];
		T b = inMat.m[i][j];
		T r = std::sqrt(a * a + b * b);
		T c = a / r;
		T s = -b / r;
		mat.SetIdentity();
		mat.m[i][i] = c;
		mat.m[j][j] = c;
		mat.m[i][j] = s;
		mat.m[j][i] = -s;
		return true;
	}

    /**
     * Jacobi Method to calculate eigenvalues and eigenvectors
     * NOTICE: IT IS ONLY APPLIED TO symmetric MATRIX.
     * Algorithm from: https://www.cmi.ac.in/~ksutar/NLA2013/iterativemethods.pdf
	 *				   with problems fixed
     */
    template<typename T, int N>
    Matrix<T, N, N> JacobiMethod(const Matrix<T, N, N>& mat, Vector<T, N>& eigenvalues, Matrix<T, N, N>& eigenvectors)
    {
        Matrix<T, N, N> temp(mat);
        Matrix<T, N, N> vectors;
        vectors.SetIdentity();

        for (int32_t i = 0; i < 19; ++i)
        {
            int32_t li, lj;
            FindLastestOffDiagonalElement(temp, li, lj);

            // calc sine and cosine
            T c, s;
            if (temp.m[li][li] == temp.m[lj][lj])
            {
                c = cos(ONE_QUARTER_PI);
                s = sin(ONE_QUARTER_PI);
            }
            else 
            {
                T theta = 0.5f * atan(2 * temp.m[li][lj] / (temp.m[lj][lj] - temp.m[li][li]));
                c = cos(theta);
                s = sin(theta);
            }

            // Build Givens rotation matrix
            Matrix<T, N, N> G;
            G.SetIdentity();
            G.m[li][li] = G.m[lj][lj] = c;
            G.m[lj][li] = s;
            G.m[li][lj] = -s;

            // Do the Givens transformation.
            Matrix<T, N, N> GT;
            TransposeMatrix(G, GT);
            temp = G * temp * GT;
            vectors = G * vectors;
        }

        for (int32_t i = 0; i < N; ++i)
        {
            eigenvalues.data[i] = temp.m[i][i];
        }
        eigenvectors = vectors;

		return temp;
    }

    // Maybe we can give a tolerance here
    template <typename T, int N>
    bool FindLastestOffDiagonalElement(const Matrix<T, N, N>& mat, int32_t& oi, int32_t& oj)
    {
        T maxValue = 0.0f;
        for (int32_t i = 1; i < N; ++i)
        {
            for (int32_t j = 0; j < i; ++j)
            {
                if (std::abs(mat.m[i][j]) > maxValue)
                {
                    oi = i;
                    oj = j;
                    maxValue = std::abs(mat.m[i][j]);
                }
            }
        }

        return maxValue != 0.0f;
    }
}