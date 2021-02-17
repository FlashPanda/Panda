#include "Numerical.hpp"

namespace Panda
{
    void TransformCoord(Vector3Df& inVec, const Matrix4f& inMat)
    {
		Vector4Df temp({ inVec[0], inVec[1], inVec[2], 1.0f });
        TransformCoord(temp, inMat);
        
		inVec.Set({ temp[0], temp[1], temp[2] });

        return;
    }

    void TransformCoord(Vector4Df& inVec, const Matrix4f& inMat)
    {
        Vector4Df temp;
		temp.Set({ inVec[0] * inMat.m[0][0] + inVec[1] * inMat.m[1][0] + inVec[2] * inMat.m[2][0] + inVec[3] * inMat.m[3][0],
			inVec[0] * inMat.m[0][1] + inVec[1] * inMat.m[1][1] + inVec[2] * inMat.m[2][1] + inVec[3] * inMat.m[3][1],
			inVec[0] * inMat.m[0][2] + inVec[1] * inMat.m[1][2] + inVec[2] * inMat.m[2][2] + inVec[3] * inMat.m[3][2],
			inVec[0] * inMat.m[0][3] + inVec[1] * inMat.m[1][3] + inVec[2] * inMat.m[2][3] + inVec[3] * inMat.m[3][3] });
        inVec = temp;
        return;
    }

    void TransformPoint(Vector3Df& inVec, const Matrix4f& inMat)
    {
		Vector4Df temp({ inVec[0], inVec[1], inVec[2], 1.0f });
        TransformCoord(temp, inMat);
        
		inVec.Set({ temp[0], temp[1], temp[2] });

        return;
    }

    void TransformDirection(Vector3Df& inVec, const Matrix4f& inMat)
    {
        Vector4Df temp({ inVec[0], inVec[1], inVec[2], 0.0f });
        TransformCoord(temp, inMat);
        
		inVec.Set({ temp[0], temp[1], temp[2] });
    }

    void TransformPoint(Vector3Df& inVec, const Matrix4f& inMat, Vector3Df& absError)
    {
        // In this function, we treet _inVec_ as a point and transform it to another point.
        // Compute transformed coordinates from _inVec_
        float xp = inVec[0] * inMat.m[0][0] + inVec[1] * inMat.m[1][0] + inVec[2] * inMat.m[2][0] + inMat.m[3][0];
        float yp = inVec[0] * inMat.m[0][1] + inVec[1] * inMat.m[1][1] + inVec[2] * inMat.m[2][1] + inMat.m[3][1];
        float zp = inVec[0] * inMat.m[0][2] + inVec[1] * inMat.m[1][2] + inVec[2] * inMat.m[2][2] + inMat.m[3][2];
        float wp = inVec[0] * inMat.m[0][3] + inVec[1] * inMat.m[1][3] + inVec[2] * inMat.m[2][3] + inMat.m[3][3];

        // Compute absolute error for transformed point
        float xAbsSum = Abs(inVec[0] * inMat.m[0][0]) + Abs(inVec[1] * inMat.m[1][0]) + 
                        Abs(inVec[2] * inMat.m[2][0]) + Abs(inMat.m[3][0]);
        float yAbsSum = Abs(inVec[0] * inMat.m[0][1]) + Abs(inVec[1] * inMat.m[1][1]) +
                        Abs(inVec[2] * inMat.m[2][1]) + Abs(inMat.m[3][1]);
        float zAbsSum = Abs(inVec[0] * inMat.m[0][2]) + Abs(inVec[1] * inMat.m[1][2]) +
                        Abs(inVec[2] * inMat.m[2][2]) + Abs(inMat.m[3][2]);
        absError.data[0] = Gamma(3) * xAbsSum;
        absError.data[1] = Gamma(3) * yAbsSum;
        absError.data[2] = Gamma(3) * zAbsSum;

        assert(wp != 0.f);
        if (wp == 1.f)
            inVec.Set({xp, yp, zp});
        else
            inVec.Set({xp / wp, yp / wp, zp / wp});
    }

    void TransformPoint(Vector3Df& inVec, const Matrix4f& inMat, const Vector3Df& inErr, Vector3Df& absError)
    {
        float xp = inVec[0] * inMat.m[0][0] + inVec[1] * inMat.m[1][0] + inVec[2] * inMat.m[2][0] + inMat.m[3][0];
        float yp = inVec[0] * inMat.m[0][1] + inVec[1] * inMat.m[1][1] + inVec[2] * inMat.m[2][1] + inMat.m[3][1];
        float zp = inVec[0] * inMat.m[0][2] + inVec[1] * inMat.m[1][2] + inVec[2] * inMat.m[2][2] + inMat.m[3][2];
        float wp = inVec[0] * inMat.m[0][3] + inVec[1] * inMat.m[1][3] + inVec[2] * inMat.m[2][3] + inMat.m[3][3];

        absError.data[0] = (Gamma(3) + 1.f) * (Abs(inMat.m[0][0]) * inErr.data[0] + Abs(inMat.m[1][0]) * inErr.data[1] + Abs(inMat.m[2][0]) * inErr.data[2]) +
                            Gamma(3) * (Abs(inVec[0] * inMat.m[0][0]) + Abs(inVec[1] * inMat.m[1][0]) + Abs(inVec[2] * inMat.m[2][0]) + Abs(inMat.m[3][0]));
        absError.data[1] = (Gamma(3) + 1.f) * (Abs(inMat.m[0][1]) * inErr.data[0] + Abs(inMat.m[1][1]) * inErr.data[1] + Abs(inMat.m[2][1]) * inErr.data[2]) +
                            Gamma(3) * (Abs(inVec[0] * inMat.m[0][1]) + Abs(inVec[1] * inMat.m[1][1]) + Abs(inVec[2] * inMat.m[2][1]) + Abs(inMat.m[3][1]));
        absError.data[2] = (Gamma(3) + 1.f) * (Abs(inMat.m[0][2]) * inErr.data[0] + Abs(inMat.m[1][2]) * inErr.data[1] + Abs(inMat.m[2][2]) * inErr.data[2]) + 
                            Gamma(3) * (Abs(inVec[0] * inMat.m[0][2]) + Abs(inVec[1] * inMat.m[1][2]) + Abs(inVec[2] * inMat.m[2][2]) + Abs(inMat.m[3][2]));

        assert(wp != 0.f);
        if (wp == 1.f)
            inVec.Set({xp, yp, zp});
        else
            inVec.Set({xp / wp, yp / wp, zp / wp});
    }

    void TransformDirection(Vector3Df& inVec, const Matrix4f& inMat, Vector3Df& absError)
    {
        float xp = inVec[0] * inMat.m[0][0] + inVec[1] * inMat.m[1][0] + inVec[2] * inMat.m[2][0];
        float yp = inVec[0] * inMat.m[0][1] + inVec[1] * inMat.m[1][1] + inVec[2] * inMat.m[2][1];
        float zp = inVec[0] * inMat.m[0][2] + inVec[1] * inMat.m[1][2] + inVec[2] * inMat.m[2][2];

        absError.data[0] = Abs(inVec[0] * inMat.m[0][0]) + Abs(inVec[1] * inMat.m[1][0]) + Abs(inVec[2] * inMat.m[2][0]);
        absError.data[1] = Abs(inVec[0] * inMat.m[0][1]) + Abs(inVec[1] * inMat.m[1][1]) + Abs(inVec[2] * inMat.m[2][1]);
        absError.data[2] = Abs(inVec[0] * inMat.m[0][2]) + Abs(inVec[1] * inMat.m[1][2]) + Abs(inVec[2] * inMat.m[2][2]);

        inVec.Set({xp, yp, zp});
    }

    void TransformDirection(Vector3Df& inVec, const Matrix4f& inMat, const Vector3Df& inErr, Vector3Df& absError)
    {
        absError.data[0] = (Gamma(3) + 1.f) * (Abs(inMat.m[0][0]) * inErr.data[0] + Abs(inMat.m[1][0]) * inErr.data[1] + Abs(inMat.m[2][0]) * inErr.data[2]) +
                            Gamma(3) * (Abs(inVec[0] * inMat.m[0][0]) + Abs(inVec[1] * inMat.m[1][0]) + Abs(inVec[2] * inMat.m[2][0]));
        absError.data[1] = (Gamma(3) + 1.f) * (Abs(inMat.m[0][1]) * inErr.data[0] + Abs(inMat.m[1][1]) * inErr.data[1] + Abs(inMat.m[2][1]) * inErr.data[2]) +
                            Gamma(3) * (Abs(inVec[0] * inMat.m[0][1]) + Abs(inVec[1] * inMat.m[1][1]) + Abs(inVec[2] * inMat.m[2][1]));
        absError.data[2] = (Gamma(3) + 1.f) * (Abs(inMat.m[0][2]) * inErr.data[0] + Abs(inMat.m[1][2]) * inErr.data[1] + Abs(inMat.m[2][2]) * inErr.data[2]) + 
                            Gamma(3) * (Abs(inVec[0] * inMat.m[0][2]) + Abs(inVec[1] * inMat.m[1][2]) + Abs(inVec[2] * inMat.m[2][2]));

        float xp = inVec[0] * inMat.m[0][0] + inVec[1] * inMat.m[1][0] + inVec[2] * inMat.m[2][0];
        float yp = inVec[0] * inMat.m[0][1] + inVec[1] * inMat.m[1][1] + inVec[2] * inMat.m[2][1];
        float zp = inVec[0] * inMat.m[0][2] + inVec[1] * inMat.m[1][2] + inVec[2] * inMat.m[2][2];

        inVec.Set({xp, yp, zp});
    }

    Ray TransformRay(const Ray& inRay, const Matrix4f& inMat, Vector3Df& oError, Vector3Df& dError)
    {
        Vector3Df o = inRay.o;
        Vector3Df d = inRay.d;

        TransformPoint(o, inMat, oError);
        TransformDirection(d, inMat, dError);

        // offset ray origin to edge of error bounds and compute _tMax_
        float tMax = inRay.tMax;
        float lengthSquared = GetLengthSquare(d);
        if (lengthSquared > 0.f)
        {
            float dt = DotProduct(Abs(d), oError) / lengthSquared;
            o += d * dt;
        }

        return Ray(o, d, tMax, inRay.time, inRay.pMedium);
    }

    Ray TransformRay(const Ray& inRay, const Matrix4f& inMat, const Vector3Df& oErrorIn, const Vector3Df& dErrorIn,
                    Vector3Df& oError, Vector3Df& dError)
    {
        Vector3Df o = inRay.o;
        Vector3Df d = inRay.d;

        TransformPoint(o, inMat, oErrorIn, oError);
        TransformDirection(d, inMat, dErrorIn, dError);

        // offset ray origin to edge of error bounds and compute _tMax_
        float tMax = inRay.tMax;
        float lengthSquared = GetLengthSquare(d);
        if (lengthSquared > 0.f)
        {
            float dt = DotProduct(Abs(d), oError) / lengthSquared;
            o += d * dt;
        }

        return Ray(o, d, tMax, inRay.time, inRay.pMedium);
    }

    // Yeah, we don't know after transforming, which component is the biggest.
    void TransformBounds(Bounds3Df& b, const Matrix4f& inMat)
    {
        Vector3Df v({b.pMin.data[0], b.pMin.data[1], b.pMin.data[2]});
        TransformCoord(v, inMat);
        Bounds3Df ret (v);

        v = {b.pMax.data[0], b.pMin.data[1], b.pMin.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        v = {b.pMin.data[0], b.pMax.data[1], b.pMin.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        v = {b.pMin.data[0], b.pMin.data[1], b.pMax.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        v = {b.pMin.data[0], b.pMax.data[1], b.pMax.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        v = {b.pMax.data[0], b.pMax.data[1], b.pMin.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        v = {b.pMax.data[0], b.pMin.data[1], b.pMax.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        v = {b.pMax.data[0], b.pMax.data[1], b.pMax.data[2]};
        TransformCoord(v, inMat);
        ret = Union(ret, v);

        b.pMin = ret.pMin;
        b.pMax = ret.pMax;
    }

    void MatrixScale(Matrix4f& outMat, const Vector3Df& vec)
    {
        outMat.SetIdentity();
        outMat.m[0][0] = vec.data[0];
        outMat.m[1][1] = vec.data[1];
        outMat.m[2][2] = vec.data[2];
        return;
    }

    void MatrixScale(Matrix4f& outMat, const float x, const float y, const float z)
    {
        outMat.SetIdentity();

        outMat.m[0][0] = x;
        outMat.m[1][1] = y;
        outMat.m[2][2] = z;

        return;
    }

    void MatrixScale(Matrix4f& outMat, const float scalar)
    {
        outMat.SetIdentity();
        outMat.m[0][0] = outMat.m[1][1] = outMat.m[2][2] = scalar;
        return;
    }

    void MatrixTranslation(Matrix4f& outMat, float x, float y, float z)
    {
        outMat.SetIdentity();
        outMat.m[3][0] = x;
        outMat.m[3][1] = y;
        outMat.m[3][2] = z;

        return;
    }

    void MatrixTranslation(Matrix4f& outMat, const Vector3Df& inVec)
    {
        outMat.SetIdentity();
        outMat.m[3][0] = inVec[0];
        outMat.m[3][1] = inVec[1];
        outMat.m[3][2] = inVec[2];

        return;
    }

    void MatrixRotationZ(Matrix4f& outMat, const float angle)
    {
        outMat.SetIdentity();

        float cosValue = cosf(angle);
        float sinValue = sinf(angle);

        outMat.m[0][0] = cosValue;
        outMat.m[0][1] = sinValue;
        outMat.m[1][0] = -sinValue;
        outMat.m[1][1] = cosValue;

        return;
    }

    void MatrixRotationX(Matrix4f& outMat, const float angle)
    {
        outMat.SetIdentity();
        float cosValue = cosf(angle);
        float sinValue = sinf(angle);

        outMat.m[1][1] = cosValue;
        outMat.m[1][2] = sinValue;
        outMat.m[2][1] = -sinValue;
        outMat.m[2][2] = cosValue;

        return;
    }

    void MatrixRotationY(Matrix4f& outMat, const float angle)
    {
        outMat.SetIdentity();
        float cosValue = cosf(angle);
        float sinValue = sinf(angle);

        outMat.m[0][0] = cosValue;
        outMat.m[0][2] = -sinValue;
        outMat.m[2][0] = sinValue;
        outMat.m[2][2] = cosValue;

        return;
    }

	// yaw - y, pitch - x, roll - z
	// order: yaw -> pitch ->roll
    void MatrixRotationYawPitchRoll(Matrix4f& outMat, const float yaw, const float pitch, const float roll)
    {
        outMat.SetIdentity();

        float cYaw, cPitch, cRoll, sYaw, sPitch, sRoll;

        cYaw = cosf(yaw);
        cPitch = cosf(pitch);
        cRoll = cosf(roll);

        sYaw = sinf(yaw);
        sPitch = sinf(pitch);
        sRoll = sinf(roll);

        outMat.m[0][0] = cRoll * cYaw + sRoll * sPitch * sYaw;
        outMat.m[0][1] = sRoll * cPitch;
        outMat.m[0][2] = -cRoll * sYaw + sRoll * sPitch * cYaw;

        outMat.m[1][0] = -sRoll * cYaw + cRoll * sYaw * sPitch;
        outMat.m[1][1] = cPitch * cRoll;
        outMat.m[1][2] = sRoll * sYaw + cRoll * sPitch * cYaw;

        outMat.m[2][0] = cPitch * sYaw;
        outMat.m[2][1] = -sPitch;
        outMat.m[2][2] = cYaw * cPitch;

        return;
    }

    void MatrixRotationAxis(Matrix4f& outMat, const Vector3Df& inVec, const float angle)
    {
        outMat.SetIdentity();

        float c = cosf(angle);
        float s = sinf(angle);
        float one_minus_c = 1.f - c;
        outMat.m[0][0] = c + inVec[0] * inVec[0] * one_minus_c;
        outMat.m[0][1] = s * inVec[2] + inVec[0] * inVec[1] * one_minus_c;
        outMat.m[0][2] = -s * inVec[1] + inVec[0] * inVec[2] * one_minus_c;

        outMat.m[1][0] = -s * inVec[2] + inVec[0] * inVec[1] * one_minus_c;
        outMat.m[1][1] = c + inVec[1] * inVec[1] * one_minus_c;
        outMat.m[1][2] = s * inVec[0] + inVec[1] * inVec[2] * one_minus_c;

        outMat.m[2][0] = s * inVec[1] + inVec[0] * inVec[2] * one_minus_c;
        outMat.m[2][1] = -s * inVec[0] + inVec[1] * inVec[2] * one_minus_c;
        outMat.m[2][2] = c + inVec[2] * inVec[2] * one_minus_c;

        return;
    }

    void MatrixRotationQuaternion(Matrix4f& outMat, const Quaternion& q)
    {
        outMat.SetIdentity();

        outMat.m[0][0] = 1.f - 2.f * q[1] * q[1] - 2.f * q[2] * q[2];
        outMat.m[0][1] = 2.f * q[0] * q[1] + 2.f * q[3] * q[2];
        outMat.m[0][2] = 2.f * q[0] * q[2] - 2.f * q[3] * q[1];

        outMat.m[1][0] = 2.f * q[0] * q[1] - 2.f * q[3] * q[2];
        outMat.m[1][1] = 1.f - 2.f * q[0] * q[0] - 2.f * q[2] * q[2];
        outMat.m[1][2] = 2.f * q[1] * q[2] + 2.f * q[3] * q[0];

        outMat.m[2][0] = 2.f * q[0] * q[2] + 2.f * q[3] * q[1];
        outMat.m[2][1] = 2.f * q[1] * q[2] - 2.f * q[3] * q[0];
        outMat.m[2][2] = 1.f - 2.f * q[0] * q[0] - 2.f * q[1] * q[1];
        
        return;
    }

    // Build coordinate system from only one vector
    // 1. zeroing one of the components of the vector(the origin vector called v1)
    // 2. swap the remaining two
    // 3. negating one of them
    // 4. normalize the vector gotten from step 3(call it v2)
    // 5. cross product with v1 and v2, we will get v3
    // 6. v1, v2, v3 consist a coordinate system
    void BuildCoordinateSystem(const Vector3Df& v1, Vector3Df& v2, Vector3Df& v3)
    {
        Vector3Df tv = Normalize(v1);
        if (std::abs(tv[0]) > std::abs(tv[1]))
        {
            Vector3Df v({-tv[2], 0, tv[0]});
            v = Normalize(v);
            v2[0] = v[0]; v2[1] = v[1]; v2[2] = v[2];
        }
        else
        {
            Vector3Df v({0, tv[2], tv[1]});
            v = Normalize(v);
            v2[0] = v[0]; v2[1] = v[1]; v2[2] = v[2];
        }

        Vector3Df v = CrossProduct(tv, v2);
        v = Normalize(v);
        v3[0] = v[0]; v3[1] = v[1]; v3[2] = v[2];
    }
    
    void BuildViewMatrix(Matrix4f& result, const Vector3Df& pos, const Vector3Df& target, const Vector3Df& up, Handness handness)
    {
        if (handness == Handness::kHandnessRight)
            BuildViewMatrixRH(result, pos, target, up);
        else
            BuildViewMatrixLH(result, pos, target, up);
        
        return;
    }

    void BuildViewMatrixRH(Matrix4f& result, const Vector3Df& pos, const Vector3Df& target, const Vector3Df& up)
    {
        Vector3Df n = pos - target;
		n = Normalize(n);
        Vector3Df u = CrossProduct(up, n);
        u = Normalize(u);
        Vector3Df v = CrossProduct(n, u);

        result.SetIdentity();
        result.m[0][0] = u[0];
        result.m[0][1] = v[0];
        result.m[0][2] = n[0];

        result.m[1][0] = u[1];
        result.m[1][1] = v[1];
        result.m[1][2] = n[1];

        result.m[2][0] = u[2];
        result.m[2][1] = v[2];
        result.m[2][2] = n[2];

        result.m[3][0] = -DotProduct(pos, u);
        result.m[3][1] = -DotProduct(pos, v);
        result.m[3][2] = -DotProduct(pos, n);

        return;
    }

    void BuildViewMatrixLH(Matrix4f& result, const Vector3Df& pos, const Vector3Df& target, const Vector3Df& up)
    {
        Vector3Df n = target - pos;
        n = Normalize(n);
        Vector3Df u = CrossProduct(up, n);
        u = Normalize(u);
        Vector3Df v = CrossProduct(n, u);

        result.SetIdentity();
        result.m[0][0] = u[0];
        result.m[0][1] = v[0];
        result.m[0][2] = n[0];

        result.m[1][0] = u[1];
        result.m[1][1] = v[1];
        result.m[1][2] = n[1];

        result.m[2][0] = u[2];
        result.m[2][1] = v[2];
        result.m[2][2] = n[2];

        result.m[3][0] = -DotProduct(pos, u);
        result.m[3][1] = -DotProduct(pos, v);
        result.m[3][2] = -DotProduct(pos, n);

        return;
    }

	void BuildOrthographicMatrix(Matrix4f& result, const float near, const float far)
	{
		if (g_DepthClipSpace == DepthClipSpace::kDepthClipZeroToOne)
		{
			result.m[0][0] = 1.f; result.m[0][1] = result.m[0][2] = result.m[0][3] = 0.f;
			result.m[1][1] = 1.f; result.m[1][0] = result.m[1][2] = result.m[1][3] = 0.f;
			result.m[2][2] = 1.f / (near - far); result.m[2][0] = result.m[2][1] = result.m[2][3] = 0.f;
			result.m[3][3] = 1.f; result.m[3][2] = near / (near - far); result.m[3][0] = result.m[3][1] = 0.f;
		}
		else // g_DepthClipSpace == DepthClipSpace::kDepthClipNegativeOneToOne
		{
			result.m[0][0] = 1.f; result.m[0][1] = result.m[0][2] = result.m[0][3] = 0.f;
			result.m[1][1] = 1.f; result.m[1][0] = result.m[1][2] = result.m[1][3] = 0.f;
			result.m[2][2] = 1.f / (near - far); result.m[2][0] = result.m[2][1] = result.m[2][3] = 0.f;
			result.m[3][3] = 1.f; result.m[3][2] = (near + far) / (2 * (near - far)); result.m[3][0] = result.m[3][1] = 0.f;
		}
	}

    void BuildPerspectiveFovMatrix(Matrix4f& result, const float fov, const float aspect, const float near, const float far, Handness handness)
    {
        if (handness == Handness::kHandnessRight)
			BuildPerspectiveFovRHMatrix(result, fov, aspect, near, far);
        else
            BuildPerspectiveFovLHMatrix(result, fov, aspect, near, far);

        return;
    }

    void BuildPerspectiveFovRHMatrix(Matrix4f& result, const float xfov, const float aspect, const float near, const float far)
    {
        result.Set(0);

        float cFov = 1.f / tanf(xfov  / 2);
        result.m[0][0] = cFov;
        result.m[1][1] = aspect * cFov;
        result.m[2][3] = -1;

        if (g_DepthClipSpace == DepthClipSpace::kDepthClipZeroToOne)
        {
            result.m[2][2] = far / (near - far);
            result.m[3][2] = near * far / (near - far);
        }
        else /* g_DepthClipSpace == DepthClipSpace::kDepthClipNegativeOneToOne */
        {
            result.m[2][2] = (near + far) / (near - far);
            result.m[3][2] = (2 * near * far) / (near - far);
        }
        
        return;
    }

	void BuildPerspectiveFovRHMatrix1(Matrix4f& result, const float xfov, const float aspect, const float near, const float far)
	{
		result.Set(0);

		float cFov = 1.f / tanf(xfov / 2);
		result.m[0][0] = cFov;
		result.m[1][1] = aspect * cFov;
		result.m[2][3] = -1;

		if (g_DepthClipSpace == DepthClipSpace::kDepthClipZeroToOne)
		{
			result.m[2][2] = far / (far - near);
			result.m[3][2] = near * far / (far - near);
		}
		else /* g_DepthClipSpace == DepthClipSpace::kDepthClipNegativeOneToOne */
		{
			result.m[2][2] = (near + far) / (far - near);
			result.m[3][2] = (2 * near * far) / (far - near);
		}

		return;
	}
	void BuildPerspectiveFovRHMatrix2(Matrix4f& result, const float xfov, const float aspect, const float near, const float far)
	{
		result.Set(0);

		float cFov = 1.f / tanf(xfov / 2);
		result.m[0][0] = -cFov;
		result.m[1][1] = -aspect * cFov;
		result.m[2][3] = 1;

		if (g_DepthClipSpace == DepthClipSpace::kDepthClipZeroToOne)
		{
			result.m[2][2] = far / (far - near);
			result.m[3][2] = near * far / (far - near);
		}
		else /* g_DepthClipSpace == DepthClipSpace::kDepthClipNegativeOneToOne */
		{
			result.m[2][2] = (near + far) / (far - near);
			result.m[3][2] = (2 * near * far) / (far - near);
		}

		return;
	}

    void BuildPerspectiveFovRHMatrix(Matrix4f& result, const float l, const float r, float b, float t, const float n, const float f)
    {
        result.Set(0);

        result.m[0][0] = 2 * n / (r - l);
        result.m[1][1] = 2 * n / (t - b);
        result.m[2][0] = (r + l) / (r - l);
        result.m[2][1] = (t + b) / (t - b);
		result.m[2][3] = -1;
        if (g_DepthClipSpace == DepthClipSpace::kDepthClipZeroToOne)
        {
            result.m[2][2] = f / (n - f);
            result.m[3][2] = n * f / (n - f);
        }
        else /* g_DepthClipSpace == DepthClipSpace::kDepthClipNegativeOneToOne */
        {
            result.m[2][2] = (n + f) / (n - f);
            result.m[3][2] = (2 * n * f) / (n - f);
        }

        return;
    }

    void BuildPerspectiveFovLHMatrix(Matrix4f& result, const float xfov, const float aspect, const float near, const float far)
    {
        result.Set(0);

        float cFov = 1.f / tanf(xfov / 2.f);
        result.m[0][0] = cFov;
        result.m[1][1] = aspect * cFov;
        result.m[2][3] = 1.0f;

        if (g_DepthClipSpace == DepthClipSpace::kDepthClipZeroToOne)
        {
            result.m[2][2] = -far / (near - far);
            result.m[3][2] = near * far / (near - far);
        }
        else
        {
            result.m[2][2] = -(far + near) / (near - far);
            result.m[3][2] = (2 * near * far) / (near - far);
        }

        return;
    }

   void BresenhamLineAlgorithm(const Pixel2D& pos1, const Pixel2D& pos2, std::vector<Pixel2D>& result)
    {
	   result.clear();

        // Makesure that start point is on the left of end point.
        Pixel2D start, end;
        if (pos1.data[0] > pos2.data[0])
        {
            start = pos2;
            end = pos1;
        }
        else 
        {
            start = pos1;
            end = pos2;
        }

        // Some special situations.
        if (start == end)
        {    
            result.push_back(start);
            return;
        }

        if (start.data[0] == end.data[0])
        {
			int16_t diff = 0;
            if (end.data[1] > start.data[1])
                diff = end.data[1] - start.data[1];
            else 
                diff = start.data[1] - end.data[1];
            for (int16_t i = 0; i < diff; ++i)
            {
				int16_t temp = start.data[1] + i;
                Pixel2D p({start.data[0], temp});
                result.push_back(p);
            }
            result.push_back(end);
            return;
        }

        if (start.data[1] == end.data[1])
        {
			int16_t diff = 0;
            if (end.data[0] > start.data[0])
                diff = end.data[0] - start.data[0];
            else 
                diff = start.data[0] - end.data[0];
            for (int16_t i = 0; i < diff; ++i)
            {
				int16_t temp = start.data[0] + i;
                Pixel2D p({temp, start.data[1]});
                result.push_back(p);
            }
            result.push_back(end);
            return;
        }

        // Common situation.
		result.push_back(start);
        float scope = 0.0f;
		scope = std::abs((float)(end.data[1] - start.data[1]) / (end.data[0] - start.data[0]));
        if (std::abs(scope) < 1.0f)
        {
            // Use x increment to determine y increment
			int16_t diff = 0;
            float error = 0.0f;
            diff = end.data[0] - start.data[0];
			int16_t xPos = start.data[0];
			int16_t yPos = start.data[1];
            for(int16_t i = 1; i < diff; ++i)
            {
                xPos++;
                error += scope;
                if (error > 0.5f)
                {
                    if (end.data[1] > start.data[1])
                        yPos++;
                    else
                        yPos--;
                    error -= 1.0f;
                } 
                result.push_back(Pixel2D({xPos, yPos}));
            }
        }
        else 
        {
            // Use y increment to determine x increment
            scope = std::abs((float)(end.data[0] - start.data[0]) / (end.data[1] - start.data[1]));
            // Use x increment to determine y increment
            int16_t diff = 0;
            float error = 0.0f;
			int16_t xPos = start.data[0];
			int16_t yPos = start.data[1];
            diff = end.data[1] - start.data[1];
            if (diff < 0)
                diff = -diff;
            for(int16_t i = 1; i < diff; ++i)
            {
                if (end.data[1] > start.data[1])
                    yPos++;
                else
                    yPos--;
                error += scope;
                if (error > 0.5f)
                {
                    xPos++;
                    error -= 1.0f;
                } 
                result.push_back(Pixel2D({xPos, yPos}));
            }                
        }
        result.push_back(end);
    }

    Point2DList BottomFlatTriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3)
    {
        Point2Df bottomL, bottomR, top;
        if (pos1.data[1] == pos2.data[1])
        {
            bottomL = pos1.data[0] < pos2.data[0]? pos1 : pos2;
            bottomR = pos1.data[0] < pos2.data[0]? pos2 : pos1;
            top = pos3;
        }
        else if (pos1.data[1] == pos3.data[1])
        {
            bottomL = pos1.data[0] < pos3.data[0]? pos1 : pos3;
            bottomR = pos1.data[0] < pos3.data[0]? pos3 : pos1;
            top = pos2;
        }
        else 
        {
            bottomL = pos2.data[0] < pos3.data[0]? pos2 : pos3;
            bottomR = pos2.data[0] < pos3.data[0]? pos3 : pos2;
            top = pos1;
        }

        float invSlope1 = (top.data[0] - bottomL.data[0]) / (top.data[1] - bottomL.data[1]);
        float invSlope2 = (top.data[0] - bottomR.data[0]) / (top.data[1] - bottomR.data[1]);

        float startX = top.data[0];
        float endX = top.data[0];
        Point2DList result;
        for (int32_t y = std::round(top.data[1]); y <= std::round(bottomL.data[1]); ++y)
        {
            for (int32_t x = std::round(startX); x <= std::round(endX); ++x)
            {
                result.push_back(std::make_shared<Point2Df>(Point2Df({(float)x, (float)y})));
            }
            startX += invSlope1;
            endX += invSlope2;
        }

        return result;
    }

	Point2DList TopFlatTriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3)
	{
		Point2Df topL, topR, bottom;
		if (pos1.data[1] == pos2.data[1])
		{
			topL = pos1.data[0] < pos2.data[0] ? pos1 : pos2;
			topR = pos1.data[0] < pos2.data[0] ? pos2 : pos1;
			bottom = pos3;
		}
		else if (pos1.data[1] == pos3.data[1])
		{
			topL = pos1.data[0] < pos3.data[0] ? pos1 : pos3;
			topR = pos1.data[0] < pos3.data[0] ? pos3 : pos1;
			bottom = pos2;
		}
		else
		{
			topL = pos2.data[0] < pos3.data[0] ? pos2 : pos3;
			topR = pos2.data[0] < pos3.data[0] ? pos3 : pos2;
			bottom = pos1;
		}

		float invSlopeL = (bottom.data[0] - topL.data[0]) / (bottom.data[1] - topL.data[1]);
		float invSlopeR = (bottom.data[0] - topR.data[0]) / (bottom.data[1] - topR.data[1]);

		Point2DList result;
		float startX = topL.data[0];
		float endX = topR.data[0];
		for (int32_t y = std::round(topL.data[1]); y <= std::round(bottom.data[1]); ++y)
		{
			for (int32_t x = std::round(startX); x <= std::round(endX); ++x)
			{
				result.push_back(std::make_shared<Point2Df>(Point2Df({ (float)x, float(y) })));
			}
			startX += invSlopeL;
			endX += invSlopeR;
		}
		return result;
	}

	Point2DList TriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3)
	{
		if (IsBottomFlat(pos1, pos2, pos3))
			return BottomFlatTriangleRasterization(pos1, pos2, pos3);

		if (IsTopFlat(pos1, pos2, pos3))
			return TopFlatTriangleRasterization(pos1, pos2, pos3);

		// Sort the points
		Point2Df top, middle, bottom;
		if (pos1.data[1] < pos2.data[1])
		{
			if (pos3.data[1] < pos1.data[1])
			{
				top = pos3;
				middle = pos1;
				bottom = pos2;
			}
			else if (pos3.data[1] < pos2.data[1])
			{
				top = pos1;
				middle = pos3;
				bottom = pos2;
			}
			else
			{
				top = pos1;
				middle = pos2;
				bottom = pos3;
			}
		}
		else
		{
			if (pos3.data[1] < pos2.data[1])
			{
				top = pos3;
				middle = pos2;
				bottom = pos1;
			}
			else if (pos3.data[1] < pos1.data[1])
			{
				top = pos2;
				middle = pos3;
				bottom = pos1;
			}
			else
			{
				top = pos2;
				middle = pos1;
				bottom = pos3;
			}
		}

		// Make the 4th point.
		Point2Df n4(middle);
		n4.data[0] = (bottom.data[0] - top.data[0]) / (bottom.data[1] - top.data[1]) * n4.data[1];
		Point2DList bottomFlat = BottomFlatTriangleRasterization(top, middle, n4);
		Point2DList topFlat = TopFlatTriangleRasterization(middle, n4, bottom);

		// Remove the bottom scanline form "bottomFlat"
		Point2DList::iterator removeStart = bottomFlat.end() - 1;
		for (Point2DList::iterator it = bottomFlat.end() - 1; it >= bottomFlat.begin(); --it)
		{
			if ((*it)->data[1] == (*removeStart)->data[1])
				removeStart = it;
			else
				break;
		}
		bottomFlat.erase(removeStart, bottomFlat.end());

		Point2DList result(bottomFlat.begin(), bottomFlat.end());
		result.insert(result.end(), topFlat.begin(), topFlat.end());
		return result;
	}

	bool IsBottomFlat(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3)
	{
		if (pos1.data[1] == pos2.data[1])
			return pos3.data[1] > pos1.data[1];
		else if (pos1.data[1] == pos3.data[1])
			return pos2.data[1] > pos1.data[1];
		else if (pos2.data[1] == pos2.data[1])
			return pos1.data[1] > pos2.data[1];
		else
			return false;
	}

	bool IsTopFlat(const Point2Df& pos1, const Point2Df& pos2, const Point2Df& pos3)
	{
		if (pos1.data[1] == pos2.data[1])
			return pos3.data[1] < pos1.data[1];
		else if (pos1.data[1] == pos3.data[1])
			return pos2.data[1] < pos1.data[1];
		else if (pos2.data[1] == pos3.data[1])
			return pos1.data[1] < pos2.data[1];
		else
			return false;
	}

    Point2DList BarycentricTriangleRasterization(const Point2Df& pos1, const Point2Df& pos2, Point2Df& pos3)
    {
        Vector2Df AB = {{pos2.data[0] - pos1.data[0], pos2.data[1] - pos1.data[1]}};
        Vector2Df AC = {{pos3.data[0] - pos1.data[0], pos3.data[1] - pos1.data[1]}};
        float aABC = CrossProduct(AB, AC);

        float xmin = pos1.data[0], xmax = pos1.data[0];
        float ymin = pos1.data[1], ymax = pos1.data[1];
        if (pos2.data[0] < xmin) xmin = pos2.data[0];
        if (pos2.data[0] > xmax) xmax = pos2.data[0];
        if (pos2.data[1] < ymin) ymin = pos2.data[1];
        if (pos2.data[1] > ymax) ymax = pos2.data[1];
        if (pos3.data[0] < xmin) xmin = pos3.data[0];
        if (pos3.data[0] > xmax) xmax = pos3.data[0];
        if (pos3.data[1] < ymin) ymin = pos3.data[1];
        if (pos3.data[1] > ymax) ymax = pos3.data[1];

        Point2DList result;
		for (int32_t y = std::round(ymin); y <= std::round(ymax); ++y)
        {
			for (int32_t x = std::round(xmin); x <= std::round(xmax); ++x)
            {
                Vector2Df AP = {{x - pos1.data[0], y - pos1.data[1]}};
                float aAPC = CrossProduct(AP, AC);
                float aABP = CrossProduct(AB, AP);
                float s = aAPC / aABC, t = aABP / aABC;
                if (s >= 0.0f && t >= 0.0f && (s + t) <= 1.0f)
                {
                    result.push_back(std::make_shared<Point2Df>(Point2Df({(float)x, (float)y})));
                }
            }
        }
        return result;
    }
}