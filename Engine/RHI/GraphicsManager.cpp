#include <iostream>
#include "GraphicsManager.hpp"
#include "SceneManager/SceneManager.hpp"
#include "Core/Interface/IApplication.hpp"
#include "ForwardRenderPass.hpp"

namespace Panda
{
	int GraphicsManager::Initialize()
	{
		int result = 0;
		m_Frames.resize(k_FrameCount);
		InitConstants();
		m_ProjectionMethod = ProjectionMethod::PM_PERSPECTIVE;

		m_DrawPasses.push_back(std::make_shared<ForwardRenderPass>());
		return result;
	}

	void GraphicsManager::Finalize()
	{
		#ifdef DEBUG
		ClearDebugBuffers();
		#endif
		ClearBuffers();
		//ClearShaders();
	}

	void GraphicsManager::Tick()
	{
		if (g_pSceneManager->IsSceneChanged())
		{
			std::cout << "[GraphicsManager] Detected Scene Change, reinitialize buffers ..." << std::endl;
			ClearBuffers();
			//ClearShaders();

			const Scene& scene = g_pSceneManager->GetScene();
			//InitializeShaders();
			InitializeBuffers(scene);
			g_pSceneManager->NotifySceneIsRenderingQueued();
		}

		UpdateConstants();

		Clear();
		Draw();
		//std::cout << m_DrawFrameContext;
	}

	void GraphicsManager::UpdateConstants()
	{
		// Generate teh view matrix based on the camera's position.
		CalculateCameraMatrix();
		CalculateLights();
	}

	void GraphicsManager::Clear()
	{

	}

	void GraphicsManager::Draw()
	{
		auto& frame = m_Frames[m_FrameIndex];

		for (auto pDrawPass : m_DrawPasses)
		{
			pDrawPass->Draw(frame);
		}
		#ifdef DEBUG
		//RenderDebugBuffers();
		#endif
	}

	void GraphicsManager::InitConstants()
	{
		// Initialize the world/model matrix to the identity matrix.
		m_Frames[m_FrameIndex].frameContext.WorldMatrix.SetIdentity();
	}

	void GraphicsManager::CalculateCameraMatrix()
	{
		auto& scene = g_pSceneManager->GetScene();
		auto pCameraNode = scene.GetFirstCameraNode();
		DrawFrameContext& frameContext = m_Frames[m_FrameIndex].frameContext;
		if (pCameraNode)
		{
			frameContext.ViewMatrix = *pCameraNode->GetCalculatedTransform();
			InverseMatrix(frameContext.ViewMatrix, frameContext.ViewMatrix);
		}
		else 
		{
			// use default camera
			Vector3Df position({0.0f, -5.0f, 0.0f}), lookAt({0.0f, 0.0f, 0.0f}), up({0.0f, 0.0f, 1.0f});
			BuildViewMatrix(frameContext.ViewMatrix, position, lookAt, up);
		}

		float fieldOfView = PI / 3.0f;
		float nearClipDistance = 1.0f;
		float farClipDistance = 100.0f;

		if (pCameraNode)
		{
			auto pCamera = scene.GetCamera(pCameraNode->GetSceneObjectRef());
			// Set the field of view and screen aspect ratio.
			fieldOfView = std::dynamic_pointer_cast<SceneObjectPerspectiveCamera>(pCamera)->GetFov();
			nearClipDistance = pCamera->GetNearClipDistance();
			farClipDistance = pCamera->GetFarClipDistance();
		}

		if (m_ProjectionMethod == ProjectionMethod::PM_PERSPECTIVE)
		{
			const GfxConfiguration& conf = g_pApp->GetConfiguration();
			float screenAspect = (float)conf.screenWidth / conf.screenHeight;

			// Build the perspective projection matrix.
			BuildPerspectiveFovMatrix(frameContext.ProjectionMatrix, fieldOfView, screenAspect, nearClipDistance, farClipDistance, g_ViewHandness);

			// try (l, r, b, t) projection
			//nearClipDistance = 0.69f;
			//BuildPerspectiveFovRHMatrix(m_DrawFrameContext.ProjectionMatrix, -1.0f, 1.0f, -1 / screenAspect, 1 / screenAspect, nearClipDistance, farClipDistance);
		}
		else
		{
			BuildOrthographicMatrix(frameContext.ProjectionMatrix, nearClipDistance, farClipDistance);
		}
	}

	void GraphicsManager::CalculateLights()
	{
		DrawFrameContext& frameContext = m_Frames[m_FrameIndex].frameContext;
		frameContext.AmbientColor = {0.01f, 0.01f, 0.01f};
		frameContext.Lights.clear();
		
		auto& scene = g_pSceneManager->GetScene();
		auto _pLightNode = scene.GetFirstLightNode();
		if (!_pLightNode) 
		{
			LightModel& light = *(new LightModel());
			light.LightPosition = { -1.0f, -5.0f, 0.0f, 1.0f };
			light.LightColor = { 10.0f, 10.0f, 10.0f, 1.0f };
			light.LightDirection = { 0.0f, 0.0f, -1.0f, 0.0f };
			light.LightSize = { 0.0f, 0.0f };
			light.Type = LightType::kLT_Point;
			light.ConstantAtten = 1.0f;
			light.LinearAtten = 0.0f;
			light.QuadraticAtten = 0.0f;
			light.CosOfInnerAngle = 1.0f;
			light.CosOfOuterAngle = 1.0f;

			frameContext.Lights.push_back(light);

			return;
		}

		for (auto pLightNode : scene.LightNodes)
		{
			//LightContext light;
			LightModel& light = *(new LightModel());
			const std::shared_ptr<Matrix4f> transPtr = pLightNode.second->GetCalculatedTransform();
			light.LightPosition = {0.0f, 0.0f, 0.0f, 1.0f};
			TransformCoord(light.LightPosition, *transPtr);
			light.LightDirection = {0.0f, 0.0f, -1.0f, 0.0f};
			TransformCoord(light.LightDirection, *transPtr);

			const std::shared_ptr<SceneObjectLight> pLight = scene.GetLight(pLightNode.second->GetSceneObjectRef());
			if (pLight)
			{
				light.LightColor = pLight->GetColor().Value;
				const AttenCurve& attenCurve = pLight->GetDistanceAttenuation();
				light.ConstantAtten = attenCurve.u.inverseSquareParams1.kc;
				light.LinearAtten = attenCurve.u.inverseSquareParams1.kl;
				light.QuadraticAtten = attenCurve.u.inverseSquareParams1.kq;
				if (pLight->GetType() == SceneObjectType::kSceneObjectTypeLightPoint)
				{
					light.Type = LightType::kLT_Point; // default set to point light
				}
				else if (pLight->GetType() == SceneObjectType::kSceneObjectTypeLightSpot)
				{
					light.Type = LightType::kLT_Spot;
					const std::shared_ptr<SceneObjectSpotLight> _pLight = std::dynamic_pointer_cast<SceneObjectSpotLight>(pLight);
					light.CosOfInnerAngle = _pLight->GetCosOfInnerAngle();
					light.CosOfOuterAngle = _pLight->GetCosOfOuterAngle();
				}
				else if (pLight->GetType() == SceneObjectType::kSceneObjectTypeLightInfinite)
				{
					light.Type = LightType::kLT_Directional;
				}
				else if (pLight->GetType() == SceneObjectType::kSceneObjectTypeLightArea)
				{
					light.Type = LightType::kLT_Area;
					auto plight = std::dynamic_pointer_cast<SceneObjectAreaLight>(pLight);
					light.LightSize = plight->GetDimension();
				}
			}
			else 
			{
				assert(0);
			}

			frameContext.Lights.push_back(std::move(light));
		}
	}

	void GraphicsManager::InitializeBuffers(const Scene& scene)
	{
		std::cout << "[GraphicsManager] GraphicsManager::InitializeBuffers()" << std::endl;
	}
	
	void GraphicsManager::ClearBuffers()
	{
		std::cout << "[GraphicsManager] GraphicsManager::ClearBuffers()" << std::endl;
	}

	void GraphicsManager::UseOrghographicsProjection()
	{
		m_ProjectionMethod = ProjectionMethod::PM_ORTHOGRAPHICS;
	}

	void GraphicsManager::UsePerspectiveProjection()
	{
		m_ProjectionMethod = ProjectionMethod::PM_PERSPECTIVE;
	}

	ProjectionMethod GraphicsManager::GetCurrentProjectionMethod()
	{
		return m_ProjectionMethod;
	}

	void GraphicsManager::UseShaderProgram(const intptr_t shaderProgram)
	{
		std::cout << "[GraphicsManager] UseShaderProgram(" << shaderProgram << ")" << std::endl;
	}

	void GraphicsManager::SetPerFrameConstants(const DrawFrameContext& context)
	{
		std::cout << "[GraphicsManager] SetPerFrameConstants(" << &context << ")" << std::endl;
	}

	void GraphicsManager::DrawBatch(const DrawBatchContext& context)
	{
		std::cout << "[GraphicsManager] DrawBatch(" << &context << ")" << std::endl;
	}

	void GraphicsManager::DrawBatchDepthOnly(const DrawBatchContext& context)
	{
		std::cout << "[GraphicsManager] DrawBatchDepthOnly(" << &context << ")" << std::endl;
	}

	#ifdef DEBUG
	void GraphicsManager::ClearDebugBuffers()
	{
		std::cout << "[GraphicsManager] GraphicsManager::ClearDebugBuffers(void)" << std::endl;
	}

	void GraphicsManager::DrawLine(const Point3Df& from, const Point3Df& to, const Vector3Df& color)
	{
		std::cout << "[GraphicsManager] GraphicsManager::DrawLine(" << from << ", "
			<< to << ", "
			<< color << ")" << std::endl;
	}

	void GraphicsManager::DrawLine(const PointList& vertices, const Vector3Df& color)
	{
		std::cout << "[GraphicsManager] GraphicsManager::DrawLine(" << vertices.size() << ","
			<< color << ")" << std::endl;
	}

	void GraphicsManager::DrawLine(const PointList& vertices, const Matrix4f& trans, const Vector3Df& color)
	{
		std::cout << "[GraphicsManager] GraphicsManager::DrawLine(" << vertices.size() << ","
			<< trans << "," 
			<< color << ")" << std::endl;
	}

	#endif
}
