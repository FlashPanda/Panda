#pragma once
#include <vector>
#include <memory>
#include "Core/GfxStructures.hpp"
#include "Core/Interface/IRuntimeModule.hpp"
#include "Core/Interface/IShaderModule.hpp"
#include "Core/Interface/IDrawPass.hpp"
#include "Core/Math/PandaMath.hpp"
#include "Core/Image.hpp"
#include "SceneManager/Scene.hpp"

namespace Panda {
	ENUM(ProjectionMethod)
	{
		PM_PERSPECTIVE = 0,
		PM_ORTHOGRAPHICS = 1
	};
	class GraphicsManager : implements IRuntimeModule {
		public:
			virtual ~GraphicsManager() {}

			virtual int Initialize();
			virtual void Finalize();

			virtual void Tick();

			virtual void Clear();
			virtual void Draw();

			void UseOrghographicsProjection();
			void UsePerspectiveProjection();
			ProjectionMethod GetCurrentProjectionMethod();

			virtual void UseShaderProgram(const intptr_t shaderProgram);
			virtual void SetPerFrameConstants(const DrawFrameContext& context);
			virtual void DrawBatch(const DrawBatchContext& context);
			virtual void DrawBatchDepthOnly(const DrawBatchContext& context);

			#ifdef DEBUG
			virtual void DrawLine(const Point3Df& from, const Point3Df& to, const Vector3Df& color);
			virtual void DrawLine(const PointList& vertices, const Vector3Df& color);
            virtual void DrawLine(const PointList& vertices, const Matrix4f& trans, const Vector3Df& color);
			virtual void ClearDebugBuffers();
			#endif

		protected:
			// virtual bool InitializeShaders();
			// virtual void ClearShaders();
			virtual void InitializeBuffers(const Scene& scene);
			virtual void ClearBuffers();

			virtual void InitConstants();
			virtual void CalculateCameraMatrix();
			virtual void CalculateLights();
			virtual void UpdateConstants();
			//virtual void RenderBuffers();



		protected:
			ProjectionMethod m_ProjectionMethod;
			static const uint32_t		k_FrameCount = 2;
			static const uint32_t		k_MaxSceneObjectCount = 65535;
			static const uint32_t		k_MaxTextureCount = 2048;

			uint32_t m_FrameIndex = 0;
			std::vector<Frame> m_Frames;
			std::vector<std::shared_ptr<IDrawPass>> m_DrawPasses;
	};

	extern GraphicsManager* g_pGraphicsManager;
}
