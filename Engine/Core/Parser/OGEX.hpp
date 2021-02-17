#pragma once

#include <unordered_map>
#include "OpenGEX/OpenGEX.h"
#include "core/portable.hpp"
#include "SceneManager/SceneNode.hpp"
#include "SceneManager/SceneObject.hpp"
#include "core/Math/Curve.hpp"
#include "SceneManager/Scene.hpp"
#include "core/Interface/SceneParser.hpp"
#include "core/Math/Linear.hpp"

namespace Panda
{
    class OgexParser : implements SceneParser
    {
        private:
            void ConvertOddlStructureToSceneNode(const ODDL::Structure& structure, std::shared_ptr<BaseSceneNode>& baseNode, Scene& scene);

        public:
            OgexParser() = default;
            virtual ~OgexParser() = default;

            virtual std::unique_ptr<Scene> Parse(const std::string& buf);

        private:
            bool m_UpIsYAxis;
    };
}