#include <iostream>
#include <string>
#include "AssetLoader.hpp"
#include "MemoryManager.hpp"
#include "Parser/OGEX.hpp"

using namespace Panda;
using namespace std;

namespace Panda
{
	Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipZeroToOne;

    MemoryManager* g_pMemoryManager = new MemoryManager();
    AssetLoader* g_pAssetLoader = new AssetLoader();
}

template<typename Key, typename T>
static std::ostream& operator<<(std::ostream& out, std::unordered_map<Key, T> map)
{
    for (auto p : map)
    {
        out << *p.second << endl;
    }

    return out;
}

int main (int , char**)
{
    g_pMemoryManager->Initialize();
    g_pAssetLoader->Initialize();

    string ogex_text = g_pAssetLoader->SyncOpenAndReadFileToString("Scene/cube.ogex");

    OgexParser* ogex_parser = new OgexParser();
    std::unique_ptr<Scene> pScene = ogex_parser->Parse(ogex_text);
    delete ogex_parser;

    std::cout << "Dump of Scene Graph" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << *pScene->SceneGraph << std::endl;

    std::cout << "Dump of Cameras" << std::endl;
    std::cout << "-----------------------" << std::endl;
    std::cout << pScene->Cameras << std::endl;

    std::cout << "Dump of Lights" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << pScene->Lights << std::endl;

    std::cout << "Dump of Geometries" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << pScene->Geometries << std::endl;

    std::cout << "Dump of Materials" << std::endl;
    std::cout << "------------------------" << std::endl;
    std::cout << pScene->Materials << std::endl;

    g_pAssetLoader->Finalize();
    g_pMemoryManager->Finalize();

    delete g_pAssetLoader;
    delete g_pMemoryManager;

	getchar();
    
    return 0;
}