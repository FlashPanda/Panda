#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include "SceneNode.hpp"
#include "SceneObject.hpp"

namespace Panda
{
    /**
     * Scene class.
     * Contains cameras, lights, materials, geometries etc.
     * Offer the interface to get scene nodes.
     */ 
    class Scene
    {
        public:
            Scene() {
                m_pDefaultMaterial = std::make_shared<SceneObjectMaterial>("default");
            }

            Scene(const char* sceneName) :
                SceneGraph(new BaseSceneNode(sceneName))
            {
            }

            ~Scene() = default;

            const std::shared_ptr<SceneObjectCamera> GetCamera(std::string key) const;
            const std::shared_ptr<SceneCameraNode> GetFirstCameraNode() const;
            //const std::shared_ptr<SceneCameraNode> GetNextCameraNode() const;

            const std::shared_ptr<SceneObjectLight> GetLight(std::string key) const;
            const std::shared_ptr<SceneLightNode> GetFirstLightNode() const;
            //const std::shared_ptr<SceneLightNode> GetNextLightNode() const;

            const std::shared_ptr<SceneObjectMaterial> GetMaterial(std::string key) const;
            const std::shared_ptr<SceneObjectMaterial> GetFirstMaterial() const;
            //const std::shared_ptr<SceneObjectMaterial> GetNextMaterial() const;

            const std::shared_ptr<SceneObjectGeometry> GetGeometry(std::string key) const;
            const std::shared_ptr<SceneGeometryNode> GetFirstGeometryNode() const;
            //const std::shared_ptr<SceneGeometryNode> GetNextGeometryNode() const;

            const std::shared_ptr<SceneObjectShape> GetShape(std::string key) const;

            void LoadResource();
        public:
            std::shared_ptr<BaseSceneNode> SceneGraph;

            std::unordered_map<std::string, std::shared_ptr<SceneCameraNode>>       CameraNodes;
            std::unordered_map<std::string, std::shared_ptr<SceneLightNode>>        LightNodes;
            std::unordered_map<std::string, std::shared_ptr<SceneGeometryNode>>     GeometryNodes;
            std::unordered_map<std::string, std::weak_ptr<SceneBoneNode>>           BoneNodes;
            
            std::unordered_map<std::string, std::shared_ptr<SceneObjectCamera>>     Cameras;
            std::unordered_map<std::string, std::shared_ptr<SceneObjectLight>>      Lights;
            std::unordered_map<std::string, std::shared_ptr<SceneObjectMaterial>>   Materials;
            std::unordered_map<std::string, std::shared_ptr<SceneObjectGeometry>>   Geometries;
            std::unordered_map<std::string, std::shared_ptr<SceneObjectShape>>      Shapes;

            std::vector<std::weak_ptr<BaseSceneNode>>                           AnimatableNodes;
            std::unordered_map<std::string, std::weak_ptr<SceneGeometryNode>>   LUTNameGeometryNode;

        private:
           std::shared_ptr<SceneObjectMaterial> m_pDefaultMaterial;
    };
}