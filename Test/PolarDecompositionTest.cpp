#include <iomanip>
#include <iostream>
#include <random>
#include "Math/PandaMath.hpp"

using namespace Panda;
using namespace std;

namespace Panda
{
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipZeroToOne;
}

void TestInverse()
{
	Matrix2f A({
		2, 1,
		5, 3
		});
	cout << "A = " << A << endl;
	Matrix2f A_Inv;
	InverseMatrix(A, A_Inv);
	cout << "Inverse A = " << A_Inv << endl;

	Matrix2f B({
		2, 3,
		6, 8
		});
	cout << "B = " << B << endl;
	Matrix2f B_Inv;
	InverseMatrix(B, B_Inv);
	cout << "Inverse B = " << B_Inv << endl;

	Matrix3f C({
		2, 0, 1,
		0, 3, 0,
		5, 0, 3
		});
	cout << "C = " << C << endl;
	Matrix3f C_Inv;
	InverseMatrix(C, C_Inv);
	cout << "Inverse C = " << C_Inv << endl;

	Matrix3f D({
		1, 2, 3,
		-3, 6, 8,
		-3, 1, 3
		});
	cout << "D = " << D << endl;
	Matrix3f D_Inv;
	InverseMatrix(D, D_Inv);
	cout << "Inverse D = " << D_Inv << endl;

	Matrix4f E({
		1, 1, 1, 1,
		1, 2, 3, 4,
		1, 3, 6, 10,
		1, 4, 10, 20
		});
	cout << "E = " << E << endl;
	Matrix4f E_Inv;
	InverseMatrix(E, E_Inv);
	cout << "Inverse E = " << E_Inv << endl;
}

void TestPolar()
{
	default_random_engine generator;
	generator.seed(48);
	uniform_real_distribution<float> distribution_r(-1.0f * PI, 1.0f * PI);
	uniform_real_distribution<float> distribution_s(0.1f, 100.0f);
	uniform_real_distribution<float> distribution_t(-1000.0f, 1000.0f);
	auto dice_r = std::bind(distribution_r, generator);
	auto dice_s = std::bind(distribution_s, generator);
	auto dice_t = std::bind(distribution_t, generator);

	Vector3Df translation({ dice_t(), dice_t(), dice_t() });
	Vector3Df scale({ dice_s(), dice_s(), dice_s() });
	Vector3Df rotation({ dice_r(), dice_r(), dice_r() });
	Matrix4f matrix;
	MatrixComposition(matrix, rotation, scale, translation);

	cerr << "Scale: " << scale;

	//Matrix3f A({
	//    matrix.m[0][0], matrix.m[0][1], matrix.m[0][2],
	//    matrix.m[1][0], matrix.m[1][1], matrix.m[1][2],
	//    matrix.m[2][0], matrix.m[2][1], matrix.m[2][2]
	//});
	Matrix3f A({
		0.1197, -0.1359, 0.1019,
		6.1709, 5.4718, 0.0487,
		-28.4736, 31.4404, 75.3926
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

	Matrix3f A_dash = P * U;
	cout << "U * P: " << A_dash;
	cout << "Error: " << A_dash - A;
}

int main(int , char**)
{
	//TestInverse();
	TestPolar();

    getchar();
    return 0;
}