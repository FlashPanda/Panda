#pragma once

#include "Core/Math/Bound.hpp"
#include "Core/Math/Ray.hpp"
#include "Core/Math/Matrix.hpp"
#include "SceneObjectInteraction.hpp"
#include "BaseSceneObject.hpp"

namespace Panda
{
    class SceneObjectShape : public BaseSceneObject
    {
    public:
		SceneObjectShape(const std::shared_ptr<Matrix4f>& objectToWorld, const std::shared_ptr<Matrix4f>& worldToObject) 
            : BaseSceneObject(SceneObjectType::kSceneObjectTypeMath),
            m_ObjectToWorld(objectToWorld), m_WorldToObject(worldToObject)
            {}
		virtual ~SceneObjectShape() {}
        virtual Bounds3Df ObjectBound() = 0;
        virtual Bounds3Df WorldBound() = 0;
        // Return if there is an intersection and the information of the point
        virtual bool Intersect(const Ray& ray, float& tHit,
			SceneObjectSurfaceInteraction& isect, bool testAlphaTexture = true) = 0;
        // Just return if there is an intersection without any information of the point
        virtual bool IntersectP(const Ray& ray, bool testAlphaTexture = true);
        virtual float Area() = 0;

        /**************************************************************************/
        /* TO BE UNDERSTOOD */
        // Sample a point on the surface of the shape.
        // Return the PDF with respect to area on the surface.
        // virtual Interaction Sample(const Vector2Df& u, float& pdf) = 0;
        // virtual float Pdf(const Interaction& ref) { return 1.f / Area();};

        // Sample a point on the shape given a reference point |ref|.
        // Return the PDF with respect to solid angle from |ref|.
		// virtual Interaction Sample(const Interaction& ref, const Vector2Df& u,
		// 	float& pdf);
		// virtual float Pdf(const Interaction& ref, const Vector3Df& wi);
        /**************************************************************************/

        const std::shared_ptr<Matrix4f> m_ObjectToWorld;
        const std::shared_ptr<Matrix4f> m_WorldToObject;
    };
}