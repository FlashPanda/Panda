#include <iostream>
#include <random>
#include "Math/Linear.hpp"
#include "Math/PandaMath.hpp"

using namespace Panda;
using namespace std;

namespace Panda
{
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipZeroToOne;
}

void TestPolar()
{
	Matrix3f A({
		-1.2845, 1.0415, 0.8232,
		25.1124, 34.2865, -4.1928,
		-36.8394, 17.2797, -79.3425
		});
	Matrix3f U, P;

	PolarDecomposition(A, U, P);

	cout.precision(4);
	cout.setf(ios::fixed);

	cout << "Polar Decomposition of matrix A: " << endl;
	cout << A;
	cout << "U: " << endl;
	cout << U;
	cout << "P: " << endl;
	cout << P;
	
}

void TestInterpolate()
{
	cout.precision(4);
	cout.setf(ios::fixed);
	int interpolateCount = 100;

	default_random_engine generator;
	generator.seed(48);
	uniform_real_distribution<float> distribution_r(-1.0f * PI, 1.0f * PI);
	uniform_real_distribution<float> distribution_s(0.1f, 100.0f);
	uniform_real_distribution<float> distribution_t(-1000.0f, 1000.0f);
	auto dice_r = std::bind(distribution_r, generator);
	auto dice_s = std::bind(distribution_s, generator);
	auto dice_t = std::bind(distribution_t, generator);

	// generator start point matrix
	Vector3Df translation_1({ -965.0195, -147.0325, 783.1465 });//({ dice_t(), dice_t(), dice_t() });
	Vector3Df scale_1({ 1.8473,42.7057,89.1682 });// ({ dice_s(), dice_s(), dice_s() });
	Vector3Df rotation_1({ -3.0317,-0.4619,2.4603 });//({ dice_r(), dice_r(), dice_r() });
	Matrix4f matrix_transform_1;
	MatrixComposition(matrix_transform_1, rotation_1, scale_1, translation_1);
	cout << endl;
	cout << "matrix_transform_1 = " << matrix_transform_1 << endl;
	cout << endl;

	Vector3Df v({ -467.0575, -430.2776, 387.5537 });//({ dice_t(), dice_t(), dice_t() });

	// generate end point matrix
	Vector3Df translation_2({ -402.0473,-783.5765,584.0685 });// ({ dice_t(), dice_t(), dice_t() });
	Vector3Df scale_2({ 26.7205,28.5576,69.4083 });// ({ dice_s(), dice_s(), dice_s() });
	Vector3Df rotation_2({ -1.4673,-1.3518,1.2175 });// ({ dice_r(), dice_r(), dice_r() });
	Matrix4f matrix_transform_2;
	MatrixComposition(matrix_transform_2, rotation_2, scale_2, translation_2);
	cout << endl;
	cout << "matrix_transform_2 = " << matrix_transform_2 << endl;
	cout << endl;
	auto v1 = v;
	TransformCoord(v1, matrix_transform_1);
	auto v2 = v;
	TransformCoord(v2, matrix_transform_2);


	cout << "Start Point: " << v1;
	cout << "End Point: " << v2;
	cout << "_________________________" << endl;
	cout << "Start Translation: " << translation_1;
	cout << "End Translation: " << translation_2;
	cout << "____________________" << endl;
	cout << "Start scalar: " << scale_1;
	cout << " End scalar: " << scale_2;
	cout << "_____________________" << endl;
	cout << "Start Rotation:" << rotation_1;
	cout << "End Rotation: " << rotation_2;
	cout << "====================" << endl;

	Linear<Matrix4f, float> linear_interpolator({ matrix_transform_1, matrix_transform_2 });

	cout << endl;
	cout << "Interpolater: " << endl;
	for (int i = 0; i <= interpolateCount; i++)
	{
		auto inter_matrix = linear_interpolator.Interpolate(i * 1.0f / interpolateCount, 1);
		Vector3Df rotation, scalar, translation;
		MatrixDecomposition(inter_matrix, rotation, scalar, translation);
		
		auto v_inter = v;
		cout << endl;
		cout << "v_inter = " << v_inter << endl;
		cout << "inter_matrix = " << inter_matrix << endl;
		TransformCoord(v_inter, inter_matrix);
		cout << "#" << i << endl;
		cout << "___________________" << endl;
		cout << "Interpolated Position: " << v_inter;
		cout << "_____________________" << endl;
		cout << "Rotation: " << rotation;
		cout << "Scalar: " << scalar;
		cout << "Translation " << translation;
		cout << "================" << endl;
	}
}

int main(int argc, char** argv)
{
	//TestPolar();
	TestInterpolate();

    getchar();
    return 0;
}