#include <iostream>
#include <string>
#include "AssetLoader.hpp"
#include "MemoryManager.hpp"

using namespace std;
using namespace Panda;

namespace Panda
{
    MemoryManager* g_pMemoryManager = new MemoryManager();
    AssetLoader* g_pAssetLoader = new AssetLoader();
}

int main (int argc, char** argv)
{
    g_pMemoryManager->Initialize();
    g_pAssetLoader->Initialize();

    string shader_pgm = g_pAssetLoader->SyncOpenAndReadFileToString("Shaders/copy.vs");

    cout << shader_pgm;

    g_pAssetLoader->Finalize();
    g_pMemoryManager->Finalize();

    delete g_pAssetLoader;
    delete g_pMemoryManager;

	getchar();
    return 0;
}