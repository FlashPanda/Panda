#include <iostream>
#include <string>
#include "AssetLoader.hpp"
#include "MemoryManager.hpp"
#include "Parser/JPEG.hpp"

using namespace std;
using namespace Panda;

namespace Panda
{
	Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipNegativeOneToOne;

    MemoryManager* g_pMemoryManager = new MemoryManager();
    AssetLoader* g_pAssetLoader = new AssetLoader();
}

int main(int argc, const char** argv)
{
    g_pMemoryManager->Initialize();
    g_pAssetLoader->Initialize();

    {
        Buffer buf;
        if (argc >= 2)
            buf = g_pAssetLoader->SyncOpenAndReadBinary(argv[1]);
        else 
            buf = g_pAssetLoader->SyncOpenAndReadBinary("Textures/huff_simple0.jpg");
        
        JfifParser jfifParser;

        Image image = jfifParser.Parse(buf);

        cout << image;
    }

    g_pAssetLoader->Finalize();
    g_pMemoryManager->Finalize();

    delete g_pAssetLoader;
    delete g_pMemoryManager;

	getchar();
    return 0;
}