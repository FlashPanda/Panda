#include "EditorLogic.hpp"
#include "Interface/IApplication.hpp"
#include "SceneManager.hpp"
#include "GraphicsManager.hpp"

using namespace Panda;
using namespace std;

int EditorLogic::Initialize()
{
    int result;

    std::cout << "[Editor Logic] Editor Logic Initialize" << std::endl;

    // if (g_pApp->GetCommandLineArgumentsCount() > 1)
    // {
    //     auto scene_filename = g_pApp->GetCommandLineArgument(1);
    //     cout << "[EditorLogic] Loading Scene: " << scene_filename << endl;
    //     result = g_pSceneManager->LoadScene(scene_filename);
    // }
    // else
    // {
        cout << "[EditorLogic] Loading Splash Scene" << endl;
		// donuts.ogex
		//result = g_pSceneManager->LoadScene("../../Asset/Scene/apple.dae");
        result = g_pSceneManager->LoadScene("../../Asset/Scene/modern_kitchen.dae");
		//result = g_pSceneManager->LoadScene("Scene/aili.ogex");
    //}

	return 0;
}

void EditorLogic::Finalize()
{
}

void EditorLogic::Tick()
{
}

#ifdef DEBUG
void EditorLogic::DrawDebugInfo()
{
}
#endif

void EditorLogic::OnLeftKeyDown()
{
    auto& scene = g_pSceneManager->GetScene();
    auto pCameraNode = scene.GetFirstCameraNode();
    if (pCameraNode) {
        auto local_axis = pCameraNode->GetLocalAxis(); 
        Vector3Df camera_x_axis;
        memcpy(camera_x_axis.data, local_axis.m[0], sizeof(camera_x_axis));

        // move camera along its local axis x direction
        pCameraNode->MoveBy(camera_x_axis);
    }
}

void EditorLogic::OnRightKeyDown()
{
    auto& scene = g_pSceneManager->GetScene();
    auto pCameraNode = scene.GetFirstCameraNode();
    if (pCameraNode) {
        auto local_axis = pCameraNode->GetLocalAxis(); 
        Vector3Df camera_x_axis;
        memcpy(camera_x_axis.data, local_axis.m[0], sizeof(camera_x_axis));

        // move along camera local axis -x direction
        pCameraNode->MoveBy(camera_x_axis * -1.0f);
    }
}

void EditorLogic::OnUpKeyDown()
{
    auto& scene = g_pSceneManager->GetScene();
    auto pCameraNode = scene.GetFirstCameraNode();
    if (pCameraNode) {
        auto local_axis = pCameraNode->GetLocalAxis(); 
        Vector3Df camera_y_axis;
        memcpy(camera_y_axis.data, local_axis.m[1], sizeof(camera_y_axis));

        // move camera along its local axis y direction
        pCameraNode->MoveBy(camera_y_axis);
    }
}

void EditorLogic::OnDownKeyDown()
{
    auto& scene = g_pSceneManager->GetScene();
    auto pCameraNode = scene.GetFirstCameraNode();
    if (pCameraNode) {
        auto local_axis = pCameraNode->GetLocalAxis(); 
        Vector3Df camera_y_axis;
        memcpy(camera_y_axis.data, local_axis.m[1], sizeof(camera_y_axis));

        // move camera along its local axis -y direction
        pCameraNode->MoveBy(camera_y_axis * -1.0f);
    }
}

void EditorLogic::OnCharKeyDown(uint32_t keyCode)
{
    if (keyCode == KEY_CODE_W)
    {
        auto& scene = g_pSceneManager->GetScene();
        auto pCameraNode = scene.GetFirstCameraNode();
        if (pCameraNode) {
			Vector3Df forward = pCameraNode->GetForwardDirection();
			float scalar = 0.3f;
            
            pCameraNode->MoveBy(forward * scalar);
        }
    }
    else if (keyCode == KEY_CODE_S)
    {
        auto& scene = g_pSceneManager->GetScene();
        auto pCameraNode = scene.GetFirstCameraNode();
		if (pCameraNode) {
			Vector3Df forward = pCameraNode->GetForwardDirection();
			float scalar = -0.3f;

			pCameraNode->MoveBy(forward * scalar);
		}
    }
    else if (keyCode == KEY_CODE_A)
    {
		auto& scene = g_pSceneManager->GetScene();
		auto pCameraNode = scene.GetFirstCameraNode();
		if (pCameraNode) {
			Vector3Df right = pCameraNode->GetRightDirection();
			float scalar = -0.3f;

			pCameraNode->MoveBy(right * scalar);
		}
    }
    else if (keyCode == KEY_CODE_D)
    {
		auto& scene = g_pSceneManager->GetScene();
		auto pCameraNode = scene.GetFirstCameraNode();
		if (pCameraNode) {
			Vector3Df right = pCameraNode->GetRightDirection();
			float scalar = 0.3f;

			pCameraNode->MoveBy(right * scalar);
		}
    }
	else if (keyCode == KEY_CODE_E)
	{
		auto& scene = g_pSceneManager->GetScene();
		auto pCameraNode = scene.GetFirstCameraNode();
		if (pCameraNode) {
			Vector3Df up = pCameraNode->GetUpDirection();
			float scalar = 0.3f;

			pCameraNode->MoveBy(up * scalar);
		}
	}
	else if (keyCode == KEY_CODE_Q)
	{
		auto& scene = g_pSceneManager->GetScene();
		auto pCameraNode = scene.GetFirstCameraNode();
		if (pCameraNode) {
			Vector3Df up = pCameraNode->GetUpDirection();
			float scalar = -0.3f;

			pCameraNode->MoveBy(up * scalar);
		}
	}
	else if (keyCode == KEY_CODE_1)
	{
		if (g_pGraphicsManager->GetCurrentProjectionMethod() == ProjectionMethod::PM_PERSPECTIVE)
			g_pGraphicsManager->UseOrghographicsProjection();
		else
			g_pGraphicsManager->UsePerspectiveProjection();
	}
}

void EditorLogic::OnCharKeyUp(uint32_t keyCode)
{

}
