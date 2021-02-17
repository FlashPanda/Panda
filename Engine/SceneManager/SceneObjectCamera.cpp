#include "SceneObjectCamera.hpp"
#include "SceneManager.hpp"
#include "Core/Sampling.hpp"

namespace Panda
{
    float SceneObjectCamera::GenerateRayDifferential(const CameraSample& sample, RayDifferential* rd) const
    {
        float wt = GenerateRay(sample, rd);
        if (wt == 0) return 0;

        // Find camra ray after shifting a fraction of a pixel in the $x$ direction
        float wtx;
        for (float eps : {.05, -.05})
        {
            CameraSample sshift = sample;
            sshift.pFilm[0] += eps;
            Ray rx;
            wtx = GenerateRay(sshift, &rx);
            rd->rxOrigin = rd->o + (rx.o - rd->o) / eps;
            rd->rxDirection = rd->d + (rx.d - rd->d) / eps;
            if (wtx != 0)
                break;
        }

        if (wtx == 0)
            return 0;

        float wty;
        for (float eps : {.05, -.05})
        {
            CameraSample sshift = sample;
            sshift.pFilm[1] += eps;
            Ray ry;
            wty = GenerateRay(sshift, &ry);
            rd->ryOrigin = rd->o + (ry.o - rd->o) / eps;
            rd->ryDirection = rd->d + (ry.d - rd->d) / eps;
            if (wty != 0)
                break;
        }

        if (wty == 0)
            return 0;

        rd->hasDifferentials = true;
        return wt;
    }

    Spectrum SceneObjectCamera::We(const Ray& ray, Point2Df* pRaster2) const
    {
        assert(0);
        return Spectrum(0.f);
    }

    void SceneObjectCamera::Pdf_We(const Ray& ray, float* pdfPos, float* pdfDir) const
    {
        assert(0);
    }

    Spectrum SceneObjectCamera::Sample_Wi(const SceneObjectInteraction& ref, const Point2Df& u,
                                       Vector3Df* wi, float* pdf, Point2Df* pRaster,
                                       VisibilityTester* vis) const
    {
        assert(0);
        return Spectrum(0.f);
    }
     
    float SceneObjectOrthogonalCamera::GenerateRay(const CameraSample& sample, Ray* ray) const
    {
        // Computer raster and camera sample position
        Point3Df pFilm = Point3Df({sample.pFilm[0], sample.pFilm[1], 0});
        const Scene& scene = g_pSceneManager->GetScene();
        auto cameraNode = scene.CameraNodes.find(m_Name);
        if (cameraNode == scene.CameraNodes.end())
        {
            assert(0);
            return 0.f;
        }
        std::shared_ptr<SceneCameraNode> pCameraNode = cameraNode->second;
        Point3Df pCamera = pFilm;
        TransformPoint(pCamera, *pCameraNode->m_RasterToNDC);
        TransformPoint(pCamera, *pCameraNode->m_NDCToCamera);
        *ray = Ray(pCamera, Vector3Df({0.f, 0.f, 1.f}));
        
        // Modify ray for depth of field
        if (m_LensRadius > 0.f)
        {
            // Sample point on lens
            Point2Df pLens = m_LensRadius * ConcentricSampleDisk(sample.pLens);

            // Compute point on plane of focus
            float ft = m_FocalDistance / ray->d[2];
            Point3Df pFocus = (*ray)(ft);
            
            // Update ray for effect of lens
            ray->o = Point3Df({pLens[0], pLens[1], 0.f});
            ray->d = Normalize(pFocus - ray->o);
        }

        ray->time = Lerp(sample.time, m_ShutterOpen, m_ShutterClose);
        ray->pMedium = m_Medium;
        Vector3Df oError, dError;
        *ray = TransformRay(*ray, *pCameraNode->m_CameraToWorld, oError, dError);
        return 1.f;
    }

    float SceneObjectOrthogonalCamera::GenerateRayDifferential(const CameraSample& sample,
                                                               RayDifferential* ray) const 
    {
        // Compute main orthographic viewing ray

        // Compute raster and camera sample positions
        Point3Df pFilm = Point3Df({sample.pFilm[0], sample.pFilm[1], 0});
        Point3Df pCamera = pFilm;
        const Scene& scene = g_pSceneManager->GetScene();
        auto cameraNode = scene.CameraNodes.find(m_Name);
        if (cameraNode == scene.CameraNodes.end())
        {
            assert(0);
            return 0.f;
        }
        std::shared_ptr<SceneCameraNode> pCameraNode = cameraNode->second;
        TransformPoint(pCamera, *pCameraNode->m_RasterToNDC);
        TransformPoint(pCamera, *pCameraNode->m_NDCToCamera);
        *ray = RayDifferential(pCamera, Vector3Df({0.f, 0.f, 1.f}));

        // Modify ray for depth of field 
        if (m_LensRadius > 0)
        {
            // Sample point on lens
            Point2Df pLens = m_LensRadius * ConcentricSampleDisk(sample.pLens);

            // Compute point on plane of focus
            float ft = m_FocalDistance / ray->d[2];
            Point3Df pFocus = (*ray)(ft);

            // Update ray for effect of lens
            ray->o = Point3Df({pLens[0], pLens[1], 0.f});
            ray->d = Normalize(pFocus - ray->o);
        }

        // Compute ray differentials for _OrthogonalCamera_
        if (m_LensRadius > 0)
        {
            // Compute _OrthonogalCamera_ ray differientials accounting for lens

            // Sample point on lens
            Point2Df pLens = m_LensRadius * ConcentricSampleDisk(sample.pLens);
            float ft = m_FocalDistance / ray->d[2];

            Point3Df pFocus = pCamera + m_DxCamera + (ft * Vector3Df({0.f, 0.f, 1.f}));
            ray->rxOrigin = Point3Df({pLens[0], pLens[1], 0});
            ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

            pFocus = pCamera + m_DyCamera + (ft * Vector3Df({0.f, 0.f, 1.f}));
            ray->ryOrigin = Point3Df({pLens[0], pLens[1], 0.f});
            ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
        } else 
        {
            ray->rxOrigin = ray->o + m_DxCamera;
            ray->ryOrigin = ray->o + m_DyCamera;
            ray->rxDirection = ray->ryDirection = ray->d;
        }

        ray->time = Lerp(sample.time, m_ShutterOpen, m_ShutterClose);
        ray->hasDifferentials = true;
        ray->pMedium = m_Medium;
        Vector3Df oError, dError;
        *ray = TransformRay(*ray, *pCameraNode->m_CameraToWorld, oError, dError);
        return 1.f;
    }

    SceneObjectPerspectiveCamera::SceneObjectPerspectiveCamera(const std::string& name, float fov, float shutterOpen, 
                                float shutterClose, float lensr, float focald,
                                Film* film, const Medium* medium
                                )
        : SceneObjectCamera(name, shutterOpen, shutterClose, film, medium)
    {
		m_LensRadius = lensr;
		m_FocalDistance = focald;

        // Compute differential changes in origin for perspective camera rays
        const Scene& scene = g_pSceneManager->GetScene();
        auto cameraNode = scene.CameraNodes.find(m_Name);
        if (cameraNode == scene.CameraNodes.end())
        {
            assert(0);
        }
        std::shared_ptr<SceneCameraNode> pCameraNode = cameraNode->second;
        Point3Df xAxis({1.f, 0.f, 1.f});
        Point3Df origin({0.f, 0.f, 0.f});
        Point3Df yAxis({0.f, 1.f, 0.f});
        TransformPoint(xAxis, *pCameraNode->m_RasterToNDC);
        TransformPoint(xAxis, *pCameraNode->m_NDCToCamera);
        TransformPoint(origin, *pCameraNode->m_RasterToNDC);
        TransformPoint(origin, *pCameraNode->m_NDCToCamera);
        TransformPoint(yAxis, *pCameraNode->m_RasterToNDC);
        TransformPoint(yAxis, *pCameraNode->m_NDCToCamera);
        m_DxCamera = xAxis - origin;
        m_DyCamera = yAxis - origin;
        
        // Compute image plane bounds at $z=1$ for _PerspectiveCamera_
        Point2Di res = film->m_FullResolution;
        Point3Df pMin = origin;
        Point3Df pMax({float(res[0]), float(res[1]), 0});
        TransformPoint(pMax, *pCameraNode->m_RasterToNDC);
        TransformPoint(pMax, *pCameraNode->m_NDCToCamera);
        pMin /= pMin[2];
        pMax /= pMax[2];
        A = std::abs((pMax[0] - pMin[0]) * (pMax[1] - pMin[1]));
    }

    float SceneObjectPerspectiveCamera::GenerateRay(const CameraSample& sample,
        Ray* ray) const 
    {
        // Compute raster and camera sample positions
        Point3Df pFilm = Point3Df({sample.pFilm[0], sample.pFilm[1], 0});
        Point3Df pCamera = pFilm;
        const Scene& scene = g_pSceneManager->GetScene();
        auto cameraNode = scene.CameraNodes.find(m_Name);
        if (cameraNode == scene.CameraNodes.end())
        {
            assert(0);
            return 0.f;
        }
        std::shared_ptr<SceneCameraNode> pCameraNode = cameraNode->second;
        TransformPoint(pCamera, *pCameraNode->m_RasterToNDC);
        TransformPoint(pCamera, *pCameraNode->m_NDCToCamera);
        *ray = Ray(Point3Df({0.f, 0.f, 0.f}), Normalize(pCamera));

        // Modify ray for depth of field
        if (m_LensRadius > 0.f)
        {
            // Sample point on lens
            Point2Df pLens = m_LensRadius * ConcentricSampleDisk(sample.pLens);

            // Compute point on plane of focus
            float ft = m_FocalDistance / ray->d[2];
            Point3Df pFocus = (*ray)(ft);

            // Update ray for effect of lens
            ray->o = Point3Df({pLens[0], pLens[1], 0.f});
            ray->d = Normalize(pFocus - ray->o);
        }

        ray->time = Lerp(sample.time, m_ShutterOpen, m_ShutterClose);
        ray->pMedium = m_Medium;
        Vector3Df oError, dError;
        *ray = TransformRay(*ray, *pCameraNode->m_CameraToWorld, oError, dError);
        return 1.f;
    }

    float SceneObjectPerspectiveCamera::GenerateRayDifferential(const CameraSample& sample, RayDifferential* ray) const
    {
        // Compute raster and camera sample positions
        Point3Df pFilm = Point3Df({sample.pFilm[0], sample.pFilm[1], 0});
        Point3Df pCamera = pFilm;
        const Scene& scene = g_pSceneManager->GetScene();
        auto cameraNode = scene.CameraNodes.find(m_Name);
        if (cameraNode == scene.CameraNodes.end())
        {
            assert(0);
            return 0.f;
        }
        std::shared_ptr<SceneCameraNode> pCameraNode = cameraNode->second;
        TransformPoint(pCamera, *pCameraNode->m_RasterToNDC);
        TransformPoint(pCamera, *pCameraNode->m_NDCToCamera);
        Vector3Df dir = Normalize(pCamera);
        *ray = RayDifferential(Point3Df({0.f, 0.f, 0.f}), dir);

        // Modify ray for depth of field
        // TODO

        // Compute offset rays for SceneObjectPerspectiveCamera_ ray differentials
        // TODO

        ray->time = Lerp(sample.time, m_ShutterOpen, m_ShutterClose);
        ray->pMedium = m_Medium;
        Vector3Df oError, dError;
        *ray = TransformRay(*ray, *pCameraNode->m_CameraToWorld, oError, dError);
        ray->hasDifferentials = true;
        return 1.f;
    }

    Spectrum SceneObjectPerspectiveCamera::We(const Ray& ray, Point2Df* pRaster2 = nullptr) const
    {
        // TODO
    }

    void SceneObjectPerspectiveCamera::Pdf_We(const Ray& ray, float* pdfPos, float* pdfDir) const
    {
        // TODO
    }

    Spectrum SceneObjectPerspectiveCamera::Sample_Wi(const SceneObjectInteraction& ref, const Point2Df& sample,
                        Vector3Df* wi, float* pdf, Point2Df* pRaster,
                        VisibilityTester* vis) const
    {
        // TODO
		return Spectrum();
    }
}