#ifndef SCENE_OBJECT_INTERACTION_HPP
#define SCENE_OBJECT_INTERACTION_HPP
#pragma once

#include "Core/Math/PandaMath.hpp"
#include "Core/Medium.hpp"
#include <memory>
#include "BaseSceneObject.hpp"

namespace Panda
{
	class SceneObjectShape;
    class SceneGeometryNode;
    // The interface of interaction.
    struct SceneObjectInteraction : public BaseSceneObject
    {
        SceneObjectInteraction() : BaseSceneObject(SceneObjectType::kSceneObjectTypeInteraction)
        {}
        SceneObjectInteraction(const Vector3Df& p, const Vector3Df& n, const Vector3Df& _pError,
            const Vector3Df& w, const MediumInterface& medium)
            : BaseSceneObject(SceneObjectType::kSceneObjectTypeInteraction),
              position(p), pError(_pError), wo(w), normal(n), mediumInterface(medium) {}
        SceneObjectInteraction(const Vector3Df& p, const Vector3Df& w,
                    const MediumInterface& medium)
            : BaseSceneObject(SceneObjectType::kSceneObjectTypeInteraction),
              position(p), wo(w), mediumInterface(medium) 
            {}
        SceneObjectInteraction(const Vector3Df& p, const MediumInterface& medium)
            : BaseSceneObject(SceneObjectType::kSceneObjectTypeInteraction),
              position(p), mediumInterface(medium)
            {}

        Ray SpawnRay(const Vector3Df& d)
        {
            Vector3Df o = OffsetRayOrigin(position, pError, normal, d);
            return Ray(o, d, Infinity); // TO BE IMPROVED
        }

        Ray SpawnRayTo(const Vector3Df& p2)
        {
            Vector3Df origin = OffsetRayOrigin(position, pError, normal, p2 - position);
            Vector3Df d = p2 - position;
            return Ray(origin, d, 1 - ShadowEpsilon);
        }

        bool IsSurfaceInteraction() const 
        {
            return normal != Vector3Df();
        }

        Vector3Df position;
        Vector3Df pError;
        Vector3Df wo;
        Vector3Df normal;
        MediumInterface mediumInterface;
    };

    class SceneObjectSurfaceInteraction : public SceneObjectInteraction
    {
        public:
        SceneObjectSurfaceInteraction() : SceneObjectInteraction() {}
        SceneObjectSurfaceInteraction(const Vector3Df& p, const Vector3Df& pError, const Vector2Df& uv,
            const Vector3Df& wo, const Vector3Df& dpdu, const Vector3Df& dpdv,
            const Vector3Df& dndu, const Vector3Df& dndv, 
            SceneObjectShape* sh, int32_t faceIndex = 0);

        void SetShadingGeometry(const Vector3Df& _dpdu, const Vector3Df& _dpdv,
            const Vector3Df& _dndu, const Vector3Df& _dndv);

            
        Vector2Df uv;
        Vector3Df dpdu, dpdv;
        Vector3Df dndu, dndv;
        std::shared_ptr<SceneObjectShape> shape;
        struct 
        {
            Vector3Df normal;
            Vector3Df dpdu, dpdv;
            Vector3Df dndu, dndv;
        } shading;
        std::shared_ptr<SceneGeometryNode> geometry;

        // TODO: add bsdf or bssrdf


        mutable Vector3Df dpdx, dpdy;
        mutable float dudx = 0.f, dvdx = 0.f, dudy = 0.f, dvdy = 0.f;
        
        int32_t faceIndex = 0;
    };

    SceneObjectSurfaceInteraction TransformInteraction(const SceneObjectSurfaceInteraction& si, const Matrix4f& inMat);
}
#endif