#pragma once
#include "../Interface/SceneParser.hpp"
#include "assimp/scene.h"
#include "assimp/code/ColladaParser.h"

namespace Panda
{
    class DaeParser : implements SceneParser
    {
	private:
		void ConvertDaeStructureToSceneNode(const Assimp::ColladaParser& parser, const Assimp::Collada::Node& inNode, std::shared_ptr<BaseSceneNode>& baseNode, Scene& scene);
		void AddGeometryObject(Scene& scene, const Assimp::Collada::MeshInstance& meshInstance,
			const Assimp::ColladaParser::MeshLibrary& inMeshLib, const Assimp::ColladaParser::MaterialLibrary& inMaterialLib,
			const Assimp::ColladaParser::EffectLibrary& inEffectLib, const Assimp::ColladaParser::ImageLibrary& inImageLib);
		void AppendTransform(std::shared_ptr<BaseSceneNode>& baseNode, const Assimp::Collada::Transform& transform);
		void AddLightObject(Scene& scene, const std::string& inLightName, const Assimp::Collada::Light& inLight);
		void AddCameraObject(Scene& scene, const std::string& inCameraName, const Assimp::Collada::Camera& inCamera);
		void AddMaterial(Scene& scene, const std::string& keyName, const std::string& matName, 
			const Assimp::Collada::Effect& effect, const Assimp::ColladaParser::ImageLibrary& inImageLib);
		void AddDefaultMaterial(Scene& scene);

		const aiLight* FindaiLight(const aiScene& _aiScene, std::string lightName);
		const aiCamera* FindaiCamera(const aiScene& _aiScene, std::string cameraName);
        public:
            DaeParser() = default;
            virtual ~DaeParser() = default;

            virtual std::unique_ptr<Scene> Parse(const std::string& file);
    };
}