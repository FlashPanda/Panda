#include <iostream>
#include <string>
#include "AssetLoader.hpp"
#include "MemoryManager.hpp"
#include "SceneManager.hpp"

using namespace Panda;
using namespace std;

namespace Panda
{
	Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipNegativeOneToOne;

    MemoryManager* g_pMemoryManager = new MemoryManager();
    AssetLoader* g_pAssetLoader = new AssetLoader();
    SceneManager* g_pSceneManager = new SceneManager();
}

template <typename T>
static ostream& operator<< (ostream& out, unordered_map<string, shared_ptr<T>> map)
{
    for (auto p : map)
        out << *p.second << endl;

    return out;
}

int main(int, char**)
{
    g_pMemoryManager->Initialize();
    g_pSceneManager->Initialize();
    g_pAssetLoader->Initialize();

    g_pSceneManager->LoadScene("Scene/Example.ogex");
    auto& scene = g_pSceneManager->GetScene();

    cout << "Dump of Geometries" << endl;
    cout << "----------------------" << endl;
    for (auto _it : scene.GeometryNodes)
    {
        auto pGeometryNode = _it.second;
        cout << *pGeometryNode << endl;
        if (pGeometryNode)
        {
            cout << *pGeometryNode << endl;
            weak_ptr<SceneObjectGeometry> pGeometry = scene.GetGeometry(pGeometryNode->GetSceneObjectRef());
            auto pObj = pGeometry.lock();
            if (pObj)
                cout << *pObj << endl;
        }
    }

    cout << "Dump of Cameras" << endl;
    cout << "----------------------" << endl;
    for (auto _it : scene.CameraNodes)
    {
        auto pCameraNode = _it.second;
        if (pCameraNode) {
            cout << *pCameraNode << endl;
            weak_ptr<SceneObjectCamera> pCamera = scene.GetCamera(pCameraNode->GetSceneObjectRef());
            auto pObj = pCamera.lock();
            if (pObj)
                cout << *pObj << endl;
        }
    }

    cout << "Dump of Lights" << endl;
    cout << "-----------------------" << endl;
    for (auto _it : scene.LightNodes)
    {
        auto pLightNode = _it.second;
        if (pLightNode) {
            cout << *pLightNode <<endl;
            weak_ptr<SceneObjectLight> pLight = scene.GetLight(pLightNode->GetSceneObjectRef());
            auto pObj = pLight.lock();
            if (pObj)
                cout << *pObj << endl;
        }
    }

    cout << "Dump of metarials" << endl;
    cout << "------------------------" << endl;
    for (auto _it : scene.Materials)
    {
        auto pMaterial = _it.second;
        if (pMaterial)
            cout << *pMaterial << endl;
    }

    cout << "Dump of Bone Nodes" << endl;
    cout << "------------------------" << endl;
    for (auto _it:scene.BoneNodes)
    {
        auto pBone = _it.second.lock();
        if (pBone)
            cout << *pBone << endl;
    }

    g_pSceneManager->Finalize();
    g_pAssetLoader->Finalize();
    g_pMemoryManager->Finalize();

    delete g_pSceneManager;
    delete g_pAssetLoader;
    delete g_pMemoryManager;

    getchar();
    return 0;
}
