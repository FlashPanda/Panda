#include "SceneGeometryNode.hpp"
#include "SceneManager.hpp"
#include "SceneObjectShape.hpp"

namespace Panda
{

    Bounds3Df SceneGeometryNode::WorldBound()
    {
        const std::shared_ptr<SceneObjectShape> pShape = g_pSceneManager->GetScene().GetShape(m_SceneObjectKey);
        return pShape->WorldBound();
    }
    bool SceneGeometryNode::Intersect(const Ray& r, SceneObjectSurfaceInteraction& isect)
    {
        float tHit;
        const std::shared_ptr<SceneObjectShape> pShape = g_pSceneManager->GetScene().GetShape(m_SceneObjectKey);
        if (!pShape->Intersect(r, tHit, isect))
            return false;
        r.tMax = tHit;
        isect.geometry = std::shared_ptr<SceneGeometryNode>(this);

        // Initialize _SurfaceInteraction::mediumInterface_ after _shape_
        // intersection
        // if (m_MediumInterface.IsMediumTransition())
        //     isect.mediumInterface = m_MediumInterface;
        // else 
        //     isect.mediumInterface = MediumInterface(r.pMedium);
        return true;
    }
    bool SceneGeometryNode::IntersectP(const Ray& r)
    {
        const std::shared_ptr<SceneObjectShape> pShape = g_pSceneManager->GetScene().GetShape(m_SceneObjectKey);
        return pShape->IntersectP(r);
    }
}