#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Math/PandaMath.hpp"

using namespace Panda;
using namespace std;

namespace Panda
{
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipZeroToOne;
}

template<class T>
void visualize(T points, const string& name)
{
	cout << name << ':' << endl;

	// dump the result
	for (auto point : points)
	{
		cout << *point;
	}

	// visualize in the console
	// note y is fliped towarding downside for easy print
	int row = 0, col = 0;
	for (auto point : points)
	{
		while (col < point->data[1])
		{
			col++;
			cout << endl;
			row = 0;
		}

		while (row++ < point->data[0])
		{
			cout << ' ';
		}

		cout << '*';
	}

	cout << endl;
}

ostream& operator<<(ostream& out, const std::vector<Pixel2D>& vec)
{
    for (int32_t i = 0; i < vec.size(); ++i)
        out << vec[i] << endl;

	return out;
}

void TestLineAlgorithm()
{
    Pixel2D pos1({0, 0}), pos2({11, 4});
    std::vector<Pixel2D> result;
	BresenhamLineAlgorithm(pos1, pos2, result);
	cout << "start point = " << pos1 << endl;
	cout << "end point = " << pos2 << endl;
	cout << "Line = " << result << endl;
	cout << endl << endl;

	Pixel2D pos3({ 0, 10 }), pos4({ 11, 7 });
	BresenhamLineAlgorithm(pos3, pos4, result);
	cout << "start point = " << pos3 << endl;
	cout << "end point = " << pos4 << endl;
	cout << "Line = " << result << endl;
	cout << endl << endl;

	Pixel2D pos5({ 1, 1 }), pos6({ 3, 6 });
	BresenhamLineAlgorithm(pos5, pos6, result);
	cout << "start point = " << pos5 << endl;
	cout << "end point = " << pos6 << endl;
	cout << "Line = " << result << endl;
	cout << endl << endl;

	Pixel2D pos7({ 1, 7 }), pos8({ 3, 1 });
	BresenhamLineAlgorithm(pos7, pos8, result);
	cout << "start point = " << pos7 << endl;
	cout << "end point = " << pos8 << endl;
	cout << "Line = " << result << endl;
	cout << endl << endl;

}

void TestBottomFlat()
{
	Point2D v1({ 5.0f, 7.0f });
	Point2D v2({ 0.0f, 8.0f });
	Point2D v3({ 21.0f, 8.0f });

	Point2DList points = BottomFlatTriangleRasterization(v1, v2, v3);

	visualize(points, "Bottom Flat Triangle");
}

void TestTopFlat()
{
	Point2D v1 = { { 0.0f, 0.0f } };
	Point2D v2 = { { 21.0f, 0.0f } };
	Point2D v3 = { { 12.0f, 8.0f } };

	Point2DList points = TopFlatTriangleRasterization(v1, v2, v3);

	visualize(points, "General Triangle");
}

void TestNormalTriangle()
{
	Point2D v1 = { { 1.0f, 0.0f } };
	Point2D v2 = { { 16.0f, 9.0f } };
	Point2D v3 = { { 30.0f, 4.0f } };

	Point2DList points = TriangleRasterization(v1, v2, v3);

	visualize(points, "General Triangle");
}

void TestBarycentric()
{
	Point2D v1 = { { 1.0f, 0.0f } };
	Point2D v2 = { { 16.0f, 9.0f } };
	Point2D v3 = { { 30.0f, 4.0f } };

	Point2DList points = BarycentricTriangleRasterization(v1, v2, v3);

	visualize(points, "General Triangle Barycentric");
}

int main()
{
    //TestLineAlgorithm();
	//TestBottomFlat();
	//TestTopFlat();
	//TestNormalTriangle();
	TestBarycentric();

    getchar();
    return 0;
}