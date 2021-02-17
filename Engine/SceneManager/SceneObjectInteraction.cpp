#include "SceneObjectInteraction.hpp"
#include "Core/Math/Numerical.hpp"

namespace Panda
{
    SceneObjectSurfaceInteraction TransformInteraction(const SceneObjectSurfaceInteraction& si, const Matrix4f& inMat)
    {
        SceneObjectSurfaceInteraction ret;
        // Transform _p_ and _pError_ in _SurfaceInteraction_
        ret.position = si.position;
        TransformPoint(ret.position, inMat, si.pError, ret.pError);

        // Transform remaining members of _SurfaceInteraction_
        ret.normal = si.normal;
        TransformDirection(ret.normal, inMat);
        Normalize(ret.normal);
        ret.wo = si.wo;
        TransformDirection(ret.wo, inMat);
        Normalize(ret.wo);
        ret.mediumInterface = si.mediumInterface;
        ret.uv = si.uv;
        ret.shape = si.shape;
        ret.dpdu = si.dpdu;
        TransformDirection(ret.dpdu, inMat);
        ret.dpdv = si.dpdv;
        TransformDirection(ret.dpdv, inMat);
        ret.dndu = si.dndu;
        TransformDirection(ret.dndu, inMat);
        ret.dndv = si.dndv;
        TransformDirection(ret.dndv, inMat);
        ret.shading.normal = si.shading.normal;
        TransformDirection(ret.shading.normal, inMat);
        Normalize(ret.shading.normal);

        ret.shading.dpdu = si.shading.dpdu;
        TransformDirection(ret.shading.dpdu, inMat);
        ret.shading.dpdv = si.shading.dpdv;
        TransformDirection(ret.shading.dpdv, inMat);
        ret.shading.dndu = si.shading.dndu;
        TransformDirection(ret.shading.dndu, inMat);
        ret.shading.dndv = si.shading.dndv;
        TransformDirection(ret.shading.dndv, inMat);

        ret.dudx = si.dudx;
        ret.dvdx = si.dvdx;
        ret.dudy = si.dudy;
        ret.dvdy = si.dvdy;

        ret.dpdx = si.dpdx;
        TransformDirection(ret.dpdx, inMat);
        ret.dpdy = si.dpdy;
        TransformDirection(ret.dpdy, inMat);

        // TODO: add bsdf or bssrdf

        ret.shading.normal = Faceforward(ret.shading.normal, ret.normal);
        ret.faceIndex = si.faceIndex;
        return ret;
    }

    SceneObjectSurfaceInteraction::SceneObjectSurfaceInteraction(const Vector3Df& p, const Vector3Df& pError, const Vector2Df& uv,
            const Vector3Df& wo, const Vector3Df& dpdu, const Vector3Df& dpdv,
            const Vector3Df& dndu, const Vector3Df& dndv, 
            SceneObjectShape* sh, int32_t faceIndex)
        : SceneObjectInteraction(p, pError, Normalize(CrossProduct(dpdu, dpdv)), wo, MediumInterface()),
        uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv), shape(sh)
    {
        shading.normal = normal;
        shading.dpdu = dpdu;
        shading.dpdv = dpdv;
        shading.dndu = dndu;
        shading.dndv = dndv;
    }

    void SceneObjectSurfaceInteraction::SetShadingGeometry(const Vector3Df& _dpdu, const Vector3Df& _dpdv,
        const Vector3Df& _dndu, const Vector3Df& _dndv)
    {
        shading.normal = Normalize(CrossProduct(_dpdu, _dpdv));

        shading.dpdu = _dpdu;
        shading.dpdv = _dpdv;
        shading.dndu = _dndu;
        shading.dndv = _dndv;
    }
}