#pragma once
#include <functional>
#include <limits>
#include <cmath>
#include <vector>
#include <set>
#include <unordered_set>
#include "Matrix.hpp"
#include "MathUtility.hpp"
#include "core/Utility.hpp"
#include "Ray.hpp"
#include "Bound.hpp"

namespace Panda
{
    inline bool operator==(const Edge& a, const Edge& b)
    {
        return (a.first == b.first && a.second == b.second) || (a.first == b.second && a.second == b.first);
    }
    typedef std::shared_ptr<Edge> EdgePtr;
    inline bool operator==(const EdgePtr& a, const EdgePtr& b)
    {
        return (a->first == b->first && a->second == b->second) || (a->first == b->second && a->second == b->first);
    }
    typedef std::unordered_set<EdgePtr> EdgeSet;
    typedef std::vector<EdgePtr> EdgeList;

    void TransformCoord(Vector4Df& inVec, const Matrix4f& inMat); // Used by the next two methods.
    void TransformCoord(Vector3Df& inVec, const Matrix4f& inMat); // Used to transform point
    void TransformDirection(Vector3Df& inVec, const Matrix4f& inMat); // Used to transform direction
    Ray TransformRay(const Ray& inRay, const Matrix4f& inMat, Vector3Df& oError, Vector3Df& dError);
    Ray TransformRay(const Ray& inRay, const Matrix4f& inMat, const Vector3Df& oErrorIn, const Vector3Df& dErrorIn,
                    Vector3Df& oError, Vector3Df& dError);
    void TransformBounds(Bounds3Df& inBounds, const Matrix4f& inMat);
    
    void TransformPoint(Vector3Df& inVec, const Matrix4f& inMat);
    void TransformPoint(Vector3Df& inVec, const Matrix4f& inMat, Vector3Df& absError);
    void TransformPoint(Vector3Df& inVec, const Matrix4f& inMat, const Vector3Df& inErr, Vector3Df& absError);
    void TransformDirection(Vector3Df& inVec, const Matrix4f& inMat, Vector3Df& absError);
    void TransformDirection(Vector3Df& inVec, const Matrix4f& inMat, const Vector3Df& inErr, Vector3Df& absError);

    void MatrixScale(Matrix4f& outMat, const Vector3Df& vec);
    void MatrixScale(Matrix4f& outMat, const float scalar);
    void MatrixScale(Matrix4f& outMat, const float x, const float y, const float z);

    void MatrixTranslation(Matrix4f& outMat, float x, float y, float z);

    void MatrixTranslation(Matrix4f& outMat, const Vector3Df& inVec);

    void MatrixRotationZ(Matrix4f& outMat, const float angle);

    void MatrixRotationX(Matrix4f& outMat, const float angle);

    void MatrixRotationY(Matrix4f& outMat, const float angle);

    void MatrixRotationYawPitchRoll(Matrix4f& outMat, const float yaw, const float pitch, const float roll);

    void MatrixRotationAxis(Matrix4f& outMat, const Vector3Df& inVec, const float angle);

    void MatrixRotationQuaternion(Matrix4f& outMat, const Quaternion& q);

    void BuildCoordinateSystem(const Vector3Df& v1, Vector3Df& v2, Vector3Df& v3);

    void BuildViewMatrix(Matrix4f& result, const Vector3Df& pos, const Vector3Df& target, const Vector3Df& up, Handness handness = Handness::kHandnessRight);

    void BuildViewMatrixLH(Matrix4f& result, const Vector3Df& pos, const Vector3Df& target, const Vector3Df& up);

    void BuildViewMatrixRH(Matrix4f& result, const Vector3Df& pos, const Vector3Df& target, const Vector3Df& up);

	void BuildOrthographicMatrix(Matrix4f& result, const float near, const float far);

    void BuildPerspectiveFovMatrix(Matrix4f& result, const float fov, const float aspect, const float near, const float far, Handness handness = Handness::kHandnessRight);

    void BuildPerspectiveFovLHMatrix(Matrix4f& result, const float fov, const float aspect, const float near, const float far);

	
    void BuildPerspectiveFovRHMatrix(Matrix4f& result, const float fov, const float aspect, const float near, const float far);
	void BuildPerspectiveFovRHMatrix1(Matrix4f& result, const float fov, const float aspect, const float near, const float far);
	void BuildPerspectiveFovRHMatrix2(Matrix4f& result, const float fov, const float aspect, const float near, const float far);
	void BuildPerspectiveFovRHMatrix(Matrix4f& result, const float l, const float r, float b, float t, const float n, const float f);

    template <typename VAL, typename PARAM>
    struct NewtonRaphson
    {
        typedef std::function<VAL(PARAM)> nr_f;
        typedef std::function<VAL(PARAM)> nr_fprime;

        static inline PARAM Solve(PARAM x0, nr_f f, nr_fprime fprime)
        {
            PARAM x, x1 = x0;

			do
			{
				x = x1;
				VAL fx = f(x);
				VAL fx1 = fprime(x);
				x1 = x - (fx / fx1);
			} while (abs(x1 - x) >= static_cast<PARAM>(10E-6));

            return x1;
        }
    };

    template<typename T, int N>
    void MatrixQRDecomposition(const Matrix<T, N, N>& mat, Matrix<T, N, N>& Q, Matrix<T, N, N>& R)
    {
		//MatrixQRDecompositionWithTypicalGramSchmidt(mat, Q, R);
        MatrixQRDecompositionWithImprovedGramSchmidt(mat, Q, R);
    }

    /**
     * DO NOT FORGET our matrix is row major.
     * SO, it is called QR decompose since A = Q * R. But in our engine, we must 
     * multiple them as A = R * Q
     */
    template <typename T, int N>
    void MatrixQRDecompositionWithImprovedGramSchmidt(const Matrix<T, N, N>& mat, Matrix<T, N, N>& Q, Matrix<T, N, N>& R)
    {
        Q.SetIdentity();
        R.SetIdentity();
        // Improve Gram Schmidt method
        Matrix<T, N, N> temp(mat);
		TransposeMatrix(temp, temp);
        
        for (int32_t i = 0; i < N; ++i)
        {
            Vector<T, N> q(temp.GetCol(i));
			T r = GetLength(q);
            q = Normalize(q);
            temp.SetCol(q, i);
            R.m[i][i] = r;
            for (int32_t j = i + 1; j < N; ++j)
            {
                Vector<T, N> t(temp.GetCol(j));
                r = DotProduct(t, q);
                R.m[i][j] = r;
                t = t - q * r;
                temp.SetCol(t, j);
            }
        }
        
        // Here it is. Our Q matrix is temp variable and R is upper triangle now!
        Q = temp;

        // DO NOT FORGET we need to transpose them!
		TransposeMatrix(Q, Q);
		TransposeMatrix(R, R);
    }

    /**
     * DO NOT FORGET our matrix is row major.
     * SO, it is called QR decompose since A = Q * R. But in our engine, we must 
     * multiple them as A = R * Q
     */
    template <typename T, int N>
    void MatrixQRDecompositionWithTypicalGramSchmidt(const Matrix<T, N, N>& mat, Matrix<T, N, N>& Q, Matrix<T, N, N>& R)
    {
        // Typical Gram Schmidt method
        Matrix<T, N, N> U;
        
        for (int32_t i = 0; i < N; ++i)
        {
			U.v[i] = mat.v[i];
            for (int32_t j = 0; j < i; ++j)
            {
				T numerator, denominator;
				numerator = DotProduct(U.v[i], U.v[j]);
				denominator = GetLength(U.v[j]);
				T coefficient = (denominator) ? numerator / denominator : 0;
				U.v[i] = U.v[i] - coefficient * U[j];
				R.m[i][j] = coefficient;
            }

			R.m[i][i] = GetLength(U.v[i]);
			if (R.m[i][i])
				U.v[i] = U.v[i] / R.m[i][i];
			if (i < N)
				Q.v[i] = U.v[i];
        }
    }

    /**
     * QR Decomposition with Givens Rotation/Refraction
     * DO NOT FORGET our matrix is row major.
     */
    template <typename T, int N>
    void MatrixQRDecompositionWithGivensRotation(const Matrix<T, N, N>& mat, Matrix<T, N, N>& Q, Matrix<T, N, N>& R)
    {
        Matrix<T, N, N> temp(mat);
		TransposeMatrix(temp, temp);

        Q.SetIdentity();
        R.SetIdentity();
        for (int32_t i = 1; i < 3; ++i)
        {
            for (int32_t j = 0; j < i; ++j)
            {
                Matrix<T, N, N> G;
                BuildGivensRotation(temp, i, j, G);
                Q = Q * G;
				TransposeMatrix(G, G);
                temp = G * temp;
            }
        }

        R = temp;

		TransposeMatrix(Q, Q);
		TransposeMatrix(R, R);
    }

    /**
     * Build Givens Rotation matrix.
     * DO NOT FORGET our matrix is row major.
     */
	template <typename T, int N>
	bool BuildGivensRotation(const Matrix<T, N, N>& inMat, int32_t i, int32_t j, Matrix<T, N, N>& mat)
	{
		assert(i > j);
		Matrix<T, N, N> temp(inMat);
		TransposeMatrix(temp, temp);
		if (temp.m[i][j] == 0)
		{
			mat.SetIdentity();
			return false;
		}
		T a = temp.m[j][j];
		T b = temp.m[i][j];
		T r = std::sqrt(a * a + b * b);
		T c = a / r;
		T s = -b / r;
		mat.SetIdentity();
		mat.m[i][i] = c;
		mat.m[j][j] = c;
		mat.m[i][j] = s;
		mat.m[j][i] = -s;

		TransposeMatrix(mat, mat);
		return true;
	}

    /**
     * Calc eigenvalues and eigenvectors
     * The input matrix must be REAL SYMMETRIC MATRIX!
     */
    template <typename T>
    void SolveEigenvaluesAndEigenvectors(const Matrix<T, 3, 3>& realSymmetricMatrix, Vector<T, 3> eigenvalues, Matrix<T, 3, 3>& eigenvectors)
    {
        Matrix<T, 3, 3> A(realSymmetricMatrix);
        T p1 = A.m[0][1] * A.m[0][1] + A.m[0][2] * A.m[0][2] + A.m[1][2] * A.m[1][2];
        if (p1 == 0)
        {
            // A is a diagonal matrix.
            eigenvalues.data[0] = A.m[0][0];
            eigenvalues.data[1] = A.m[1][1];
            eigenvalues.data[2] = A.m[2][2];
            eigenvectors.SetIdentity();
        }
        else
        {
            T q = A.GetTrace() / 3;
            T p2 = pow(A.m[0][0] - q, 2) + pow(A.m[1][1] - q, 2) + pow(A.m[2][2] - q, 2) + 2 * p1;
            T p = sqrt(p2 / 6);
            Matrix3<T> I;
            I.SetIdentity();
            Matrix3<T> B = (1 / p) * (A - q * I);   // I is an identity matrix.
            T r = B.GetDeterminant() / 2;

            T phi;
            if (r <= -1)
                phi = PI / 3;
            else if (r >= 1)
                phi = 0;
            else 
                phi = acos(r) / 3;

            // Determine the eigenvalues
            eigenvalues.data[0] = q + 2 * p * cos(phi);
            eigenvalues.data[2] = q + 2 * p * cos(phi + 2 * PI / 3);
            eigenvalues.data[1] = 3 * q - eigenvalues.data[0] - eigenvalues.data[2]; // The trace of A equals the sum of all eigenvalues.

            // Determine the eigenvectors
            Matrix3<T> M1 = A - eigenvalues.data[0] * I;
            Matrix3<T> M2 = A - eigenvalues.data[1] * I;
            Matrix3<T> M3 = A - eigenvalues.data[2] * I;
            Matrix3<T> tM = M2 * M3;
            eigenvectors.v[0] = tM.v[0];
            tM = M1 * M3;
            eigenvectors.v[1] = tM.v[1];
            tM = M1 * M2;
            eigenvectors.v[2] = tM.v[2];
            eigenvectors.v[0].Normalize();
            eigenvectors.v[1].Normalize();
            eigenvectors.v[2].Normalize();
        }
        
    }

    template <typename T, int N>
    INLINE void PolarDecomposition(const Matrix<T, N, N>& inMat, Matrix<T, N, N>& U, Matrix<T, N, N>& P)
    {
        U = inMat;
        Matrix<T, N, N> UPre;
        Matrix<T, N, N> UInv;

		do
		{
			UPre = U;
			UInv = U;
			if (!InverseMatrix(UInv, UInv))
				assert(0);
			Matrix<T, N, N> UInvTrans;
			TransposeMatrix(UInv, UInvTrans);
			U = (U + UInvTrans) * (T)0.5;
		} while (U != UPre);

        UInv = U;
		if (!InverseMatrix(UInv, UInv)) assert(0);
        P = inMat * UInv;
    }

    INLINE void MatrixComposition(Matrix4f& oMat, const Vector3Df& rotation, const Vector3Df& scalar, const Vector3Df& translation)
    {
        Matrix4f matRotationX, matRotationY, matRotationZ;
        MatrixRotationX(matRotationX, rotation.data[0]);
        MatrixRotationY(matRotationY, rotation.data[1]);
        MatrixRotationZ(matRotationZ, rotation.data[2]);
        Matrix4f matRatation = matRotationX * matRotationY * matRotationZ;
        Matrix4f matScale;
        MatrixScale(matScale, scalar);
        Matrix4f matTranslation;
        MatrixTranslation(matTranslation, translation);
        oMat = matScale * matRatation * matTranslation;
    }

    INLINE void MatrixDecomposition(const Matrix4f& inMat, Vector3Df& rotation, Vector3Df& scalar, Vector3Df& translation)
    {
		translation.Set({ inMat.m[3][0], inMat.m[3][1], inMat.m[3][2] });

        // Use Polar decomposition
		Matrix3f bases({ inMat.m[0][0], inMat.m[0][1], inMat.m[0][2],
			inMat.m[1][0], inMat.m[1][1], inMat.m[1][2],
			inMat.m[2][0], inMat.m[2][1], inMat.m[2][2] }
        );

        Matrix3f U, P;
        PolarDecomposition(bases, U, P);

		scalar.Set({ P.m[0][0], P.m[1][1], P.m[2][2] });

        // decompose the rotation matrix
        float thetaX = atan2(U.m[1][2], U.m[2][2]);
        float thetaY = -asinf(U.m[0][2]);
        float thetaZ = atan2(U.m[0][1], U.m[0][0]);

		rotation.Set({ thetaX, thetaY, thetaZ });
    }

    void BresenhamLineAlgorithm(const Pixel2D& pos1, const Pixel2D& pos2, std::vector<Pixel2D>& result);

    Point2DList BottomFlatTriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3);

	Point2DList TopFlatTriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3);

	Point2DList TriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3);

	bool IsBottomFlat(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3);

	bool IsTopFlat(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3);

    Point2DList BarycentricTriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, Point2Df& pos3);

    inline Vector3Df SphericalDirection(float sinTheta, float cosTheta, float phi)
    {
        return Vector3Df ({sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta});
    }

    inline Vector3Df SphericalDirection(float sinTheta, float cosTheta, float phi,
                                        const Vector3Df& x, const Vector3Df& y, 
                                        const Vector3Df& z)
    {
        return sinTheta * std::cosf(phi) * x + sinTheta * std::sinf(phi) * y +
               cosTheta * z;
    }

    inline float SphericalTheta(const Vector3Df& v)
    {
        return std::acosf(Clamp(v[2], -1.f, 1.f));
    }

    inline float Sphericalphi(const Vector3Df& v)
    {
        float p = std::atan2f(v[1], v[0]);
        return (p < 0.f)? (p + 2 * PI) : p;
    }
}