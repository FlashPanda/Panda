#pragma once
#include "BaseSceneNode.hpp"

namespace Panda
{
    class SceneObjectCamera;
    class SceneCameraNode : public BaseSceneNode
    {
        protected:
            Vector3Df m_Target;

        public:
            using BaseSceneNode::BaseSceneNode;

            void SetTarget(Vector3Df& target) {m_Target = target;}
            const Vector3Df& GetTarget() {return m_Target;}
			Vector3Df GetForwardDirection()
			{
				auto pTransform = GetCalculatedTransform();
				Vector4Df forward = Vector4Df({ 0.0f, 0.0f, -1.0f, 0.0f });
				TransformCoord(forward, *pTransform);
				return Normalize(Vector3Df({ forward.data[0], forward.data[1], forward.data[2] }));
			}

			Vector3Df GetRightDirection()
			{
				Vector3Df forward = GetForwardDirection();
				auto pTransform = GetCalculatedTransform();
				Vector4Df up4 = Vector4Df({ 0.0f, 1.0f, 0.0f, 0.0f });
				TransformCoord(up4, *pTransform);
				Vector3Df up3 = Normalize(Vector3Df({ up4.data[0], up4.data[1], up4.data[2] }));

				Vector3Df right = CrossProduct(up3, -forward);
				return Normalize(right);
			}

			Vector3Df GetUpDirection()
			{
				Vector3Df forward = GetForwardDirection();
				Vector3Df right = GetRightDirection();

				return CrossProduct(-forward, right);
			}

            Matrix3f GetLocalAxis()
            {
                Matrix3f result;
                auto pTransform = GetCalculatedTransform();
                Vector3Df target = GetTarget();
                Vector3Df cameraPosition = Vector3Df(0.0f);
                TransformCoord(cameraPosition, *pTransform);
				Vector3Df up({ 0.0f, 0.0f, 1.0f });
                Vector3Df cameraZAxis = cameraPosition - target;
				cameraZAxis = Normalize(cameraZAxis);
                Vector3Df cameraXAxis;
                Vector3Df cameraYAxis;
                cameraXAxis = CrossProduct(cameraZAxis, up);
                cameraYAxis = CrossProduct(cameraXAxis, cameraZAxis);
                result.v[0] = cameraXAxis;
                result.v[1] = cameraYAxis;
                result.v[2] = cameraZAxis;

                return result;
            }

            // For offline-rendering
            void GenerateWorldToCameraTransformAndInv();
            void GenerateCameraToNDCTransformAndInv();
            void GenerateNDCToRasterTransformAndInv();

        public:
            std::shared_ptr<SceneObjectTransform> m_WorldToCamera, m_CameraToWorld;
            std::shared_ptr<SceneObjectTransform> m_CameraToNDC, m_NDCToCamera;
            std::shared_ptr<SceneObjectTransform> m_NDCToRaster, m_RasterToNDC;
    };
}