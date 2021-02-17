#pragma once
#include <vector>
#include <memory>
#include "Core/ParamSet.hpp"
#include "SceneObjectTransform.hpp"
#include "SceneObjectShape.hpp"

namespace Panda
{
    std::vector<std::shared_ptr<SceneObjectShape>> CreateLoopDubDiv(
		const std::shared_ptr<Matrix4f>& objectToWorld, 
		const std::shared_ptr<Matrix4f>& worldToObject,
        const ParamSet& params);
}