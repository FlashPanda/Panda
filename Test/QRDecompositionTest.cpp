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


int main(int, char**)
{
    default_random_engine generator;
    generator.seed(48);
    uniform_real_distribution<float> distribution_r(-1.0f * PI, 1.0f * PI);
    uniform_real_distribution<float> distribution_s(0.1f, 1.0f);
    uniform_real_distribution<float> distribution_t(-1000.0f, 1000.0f);
    auto dice_r = std::bind(distribution_r, generator);
    auto dice_s = std::bind(distribution_s, generator);
    auto dice_t = std::bind(distribution_t, generator);

	Vector3Df translation({ dice_t(), dice_t(), dice_t() });
	Vector3Df scale({ dice_s(), dice_s(), dice_s() });
	Vector3Df rotation({ dice_r(), dice_r(), dice_r() });
    Matrix4f matrix;
    MatrixComposition(matrix, rotation, scale, translation);

    //Matrix3f A({
    //    matrix.m[0][0], matrix.m[0][1], matrix.m[0][2],
    //    matrix.m[1][0], matrix.m[1][1], matrix.m[1][2],
    //    matrix.m[2][0], matrix.m[2][1], matrix.m[2][2]
    //});
	Matrix3f A({
		-1.2845, 1.0415, 0.8233,
		25.1114, 34.2873, -4.1928,
		-36.8410, 17.2794, -79.3418
		});
    Matrix3f Q, R;

    MatrixQRDecomposition(A, Q, R);

    cout.precision(4);
    cout.setf(ios::fixed);
    cout << "QR Decomposition of matrix A: " << endl;
    cout << A;
    cout << "Q: " << Q <<endl;
    cout << "R: " << R <<endl;

    Matrix3f A_dash = R * Q;
    cout << "Q * R: " << A_dash;
    cout << "Error: " << A_dash - A; 

    getchar();
    return 0;
}