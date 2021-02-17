#pragma once
#include "Core/Math/PandaMath.hpp"
#include "Core/Interface/IRuntimeModule.hpp"
#include "Core/Interface/SceneParser.hpp"
#include "Scene.hpp"
#include "BaseSceneNode.hpp"

namespace Panda
{
    class SceneManager : implements IRuntimeModule
    {
        public:
            virtual ~SceneManager();

            virtual int Initialize();
            virtual void Finalize();

            virtual void Tick();

            int LoadScene(const char* sceneFileName);

            bool IsSceneChanged();
            void NotifySceneIsRenderingQueued();

            const Scene& GetScene();
            const Scene& GetSceneForPhysicalSimulation();
            void ResetScene();

            std::weak_ptr<BaseSceneNode> GetRootNode();
            std::weak_ptr<SceneGeometryNode> GetSceneGeometryNode(std::string name);
            std::weak_ptr<SceneObjectGeometry> GetSceneGeometryObject(std::string key);

        protected:
            bool LoadOgexScene(const char* ogexSceneFileName);
			bool LoadDaeScene(const char* daeSceneFileName);

        protected:
            std::shared_ptr<Scene> m_pScene;
            bool m_IsRenderingQueued = false;
            bool m_DirtyFlag = false;
    };

    extern SceneManager* g_pSceneManager;
}
