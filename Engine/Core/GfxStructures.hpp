#pragma once
#include <vector>
#include "Math/PandaMath.hpp"
#include "SceneManager/Scene.hpp"

namespace Panda
{
    struct LightModel
    {
        Vector4Df LightPosition;
        Vector4Df LightColor;
        Vector4Df LightDirection;
        Vector2Df LightSize;
        LightType Type;
        float ConstantAtten;
        float LinearAtten;
        float QuadraticAtten;
        float CosOfInnerAngle;
        float CosOfOuterAngle;

        LightModel()
        {
            LightPosition = {0.0f, 0.0f, 0.0f, 1.0f};
            LightColor = {1.0f, 1.0f, 1.0f, 1.0f};
            LightDirection = {0.0f, 0.0f, -1.0f, 0.0f};
            LightSize = {0.0f, 0.0f};
            Type = LightType::kLT_Point;
            ConstantAtten = 1.0f;
            LinearAtten = 0.0f;
            QuadraticAtten = 0.0f;
            CosOfInnerAngle = 0.0f;
            CosOfOuterAngle = 0.0f;
        }
    };

    // The context for each FRAME!
    struct DrawFrameContext
    {
        Matrix4f WorldMatrix;
        Matrix4f ViewMatrix;
        Matrix4f ProjectionMatrix;
        Vector3Df AmbientColor;
        //std::vector<LightContext> Lights;
        std::vector<LightModel> Lights;
    };

    // The context for each BATCH!
    // We can draw many batches for each frame.
    // Each RHI can generate its own batch context.
    struct DrawBatchContext
    {
        std::shared_ptr<SceneGeometryNode> node;
        std::shared_ptr<SceneObjectMaterial> material;
        Matrix4f trans;

        virtual ~DrawBatchContext() = default;
    };

    struct Frame
    {
        DrawFrameContext frameContext;
        std::vector<std::shared_ptr<DrawBatchContext>> batchContexts;
    };
}