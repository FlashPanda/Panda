#include <iostream>
#include "ColorSpaceConversion.hpp"

using namespace std;
using namespace Panda;

namespace Panda
{
	Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipNegativeOneToOne;
}

int main(int argc, const char** argv)
{
    int result = 0;

    RGBf rgb({64, 35, 17});
	cout << "RGB color: " << rgb << endl;;

    YCbCrf ycbcr = ConvertRGB2YCbCr(rgb);
    cout << "When transform to YCbCr: " << ycbcr << endl;

    rgb = ConvertYCbCr2RGB(ycbcr);
    cout << "When transform back to RGB:" << rgb << endl;

	getchar();

    return result;
}