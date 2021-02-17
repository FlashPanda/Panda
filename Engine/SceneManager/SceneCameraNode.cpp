#include "SceneCameraNode.hpp"

namespace Panda
{
    void SceneCameraNode::GenerateWorldToCameraTransformAndInv()
    {
        Matrix4f viewMatrix = *this->GetCalculatedTransform();
        TransposeMatrix(viewMatrix, viewMatrix);
        m_WorldToCamera = std::make_shared<SceneObjectTransform>(viewMatrix);
        Matrix4f viewMatrixInv;
        InverseMatrix(viewMatrix, viewMatrixInv);
        m_CameraToNDC = std::make_shared<SceneObjectTransform>(viewMatrixInv);
    }

    void SceneCameraNode::GenerateCameraToNDCTransformAndInv()
    {
        float fieldOfView = PI / 3.0f;
		float nearClipDistance = 1.0f;
		float farClipDistance = 100.0f;
        
        float screenAspect = 16.f / 9.f;
        Matrix4f cameraToNDC;
        BuildPerspectiveFovMatrix(cameraToNDC, fieldOfView, screenAspect, nearClipDistance, farClipDistance);

        TransposeMatrix(cameraToNDC, cameraToNDC);
        m_CameraToNDC = std::make_shared<SceneObjectTransform>(cameraToNDC);
        Matrix4f ndcToCamera;
        InverseMatrix(cameraToNDC, ndcToCamera);
        m_NDCToCamera = std::make_shared<SceneObjectTransform>(ndcToCamera);
    }

    // Suppose NDC coordinates are 0 to 1 in z axis.
    void SceneCameraNode::GenerateNDCToRasterTransformAndInv()
    {
        Matrix4f ndcToRaster;
        MatrixTranslation(ndcToRaster, 1.f, -1.f, 0.f);
        Matrix4f scale1;
        MatrixScale(scale1, 0.5f, 0.5f, 1);
        Matrix4f scale2;
        MatrixScale(scale2, 1280.f, 720.f, 1);
        m_NDCToRaster = std::make_shared<SceneObjectTransform>(scale2 * scale1 * ndcToRaster);

        Matrix4f rasterToNDC;
        InverseMatrix(rasterToNDC, ndcToRaster);
        m_RasterToNDC = std::make_shared<SceneObjectTransform>(rasterToNDC);
    }
}