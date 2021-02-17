#include <iostream>
#include <string>
#include "AssetLoader.hpp"
#include "MemoryManager.hpp"
#include "Parser/PNG.hpp"
#include "Parser/JPEG.hpp"

using namespace Panda;
using namespace std;

namespace Panda
{
	Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipNegativeOneToOne;

    MemoryManager* g_pMemoryManager = new MemoryManager();
    AssetLoader*   g_pAssetLoader = new AssetLoader();
}

void TestPng(int argc, const char** argv)
{
	Buffer buf;
	if (argc >= 2)
	{
		buf = g_pAssetLoader->SyncOpenAndReadBinary(argv[1]);
	}
	else
	{
		buf = g_pAssetLoader->SyncOpenAndReadBinary("Textures/eye.png");
	}

	PngParser pngParser;

	Image image = pngParser.Parse(buf);

	cout << image;
}

void TestJpeg(int argc, const char**argv)
{
	Buffer buf1;
	buf1 = g_pAssetLoader->SyncOpenAndReadBinary("Textures/jpeg_decoder_test_8.jpg");
	JfifParser parser1;
	Image image1 = parser1.Parse(buf1);
	cout << image1;

	Buffer buf2;
	buf2 = g_pAssetLoader->SyncOpenAndReadBinary("Textures/jpeg_decoder_test_8.jpg");
	JfifParser parser2;
	Image image2 = parser2.Parser1(buf2);
	cout << image2;

	int x = 0;
}

int main(int argc, const char** argv)
{
    g_pMemoryManager->Initialize();
    g_pAssetLoader->Initialize();

	//TestPng(argc, argv);
	TestJpeg(argc, argv);

    g_pAssetLoader->Finalize();
    g_pMemoryManager->Finalize();

    delete g_pAssetLoader;
    delete g_pMemoryManager;

    return 0;
}