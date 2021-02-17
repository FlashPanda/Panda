#ifndef SCENE_OBJECT_CAMERA_HPP
#define SCENE_OBJECT_CAMERA_HPP
#pragma once
#include "BaseSceneObject.hpp"
#include "Core/Math/PandaMath.hpp"
#include "Core/Medium.hpp"
#include "Core/Film.hpp"
#include "Core/Spectrum.hpp"
#include "SceneObjectInteraction.hpp"
#include "SceneObjectLight.hpp"
#include "SceneObjectTransform.hpp"
#include "Core/ParamSet.hpp"
//#include "SceneManager.hpp"
//#include "Scene.hpp"

namespace Panda
{
    struct CameraSample
    {
        Point2Df pFilm;
        Point2Df pLens;
        float time;
    };

    inline std::ostream& operator<<(std::ostream& out, const CameraSample& cs)
    {
        out << "[pFilm: " << cs.pFilm << ", pLens: " << cs.pLens << ", time: " << cs.time << std::endl;
        return out;
    }

    class SceneCameraNode;
	class Scene;
    // camera
    class SceneObjectCamera : public BaseSceneObject
    {
        protected:
            float m_Aspect;
            float m_NearClipDistance;
            float m_FarClipDistance;
			const float m_ShutterOpen, m_ShutterClose;	// shutter open / close time
			Film* m_Film;
			const Medium* m_Medium;

        public:
            void SetColor(std::string& attrib, Vector4Df& color)
            {
                // TODO: extension
            }

            void SetParam(std::string& attrib, float param)
            {
				if (attrib == "near")
				{
					m_NearClipDistance = param;
				}
				else if (attrib == "far")
				{
					m_FarClipDistance = param;
				}
				else if (attrib == "aspect")
					m_Aspect = param;
            }

            void SetTexture(std::string& attrib, std::string& textureName)
            {
                // TODO: extension
            }

            float GetNearClipDistance() const {return m_NearClipDistance;}
            float GetFarClipDistance() const {return m_FarClipDistance;}
        public:
            SceneObjectCamera() : 
                BaseSceneObject(SceneObjectType::kSceneObjectTypeCamera), m_Aspect(16.0f / 9.0f), m_NearClipDistance(1.0f), m_FarClipDistance(100.0f),
                m_ShutterOpen(0.f), m_ShutterClose(0.f), m_Medium(nullptr),m_Film(nullptr)
                {}
            SceneObjectCamera(const std::string& name) :
                BaseSceneObject(SceneObjectType::kSceneObjectTypeCamera, name), m_Aspect(16.0f / 9.0f), m_NearClipDistance(1.0f), m_FarClipDistance(100.0f),
                m_ShutterOpen(0.f), m_ShutterClose(0.f), m_Medium(nullptr),m_Film(nullptr)
                {}
            SceneObjectCamera(const std::string& name, float shutterOpen, float shutterClose, Film* film, const Medium* medium) :
                BaseSceneObject(SceneObjectType::kSceneObjectTypeCamera, name), m_Aspect(16.0f / 9.0f), m_NearClipDistance(1.0f), m_FarClipDistance(100.0f),
                m_ShutterOpen(shutterOpen), m_ShutterClose(shutterClose), m_Film(film), m_Medium(medium)
                {}
            virtual ~SceneObjectCamera() {}

            virtual float GenerateRay(const CameraSample& sample, Ray* ray) const = 0;
            virtual float GenerateRayDifferential(const CameraSample& sample, RayDifferential* rd) const;
            virtual Spectrum We(const Ray& ray, Point2Df* pRaster2 = nullptr) const;
            virtual void Pdf_We(const Ray& ray, float* pdfPos, float* pdfDir) const;
            virtual Spectrum Sample_Wi(const SceneObjectInteraction& ref, const Point2Df& u,
                                       Vector3Df* wi, float* pdf, Point2Df* pRaster,
                                       VisibilityTester* vis) const;

            friend std::ostream& operator<<(std::ostream& out, const SceneObjectCamera& obj);
    };

    class SceneObjectPorjectiveCamera : public SceneObjectCamera
    {
        public:
            SceneObjectPorjectiveCamera() : SceneObjectCamera(), m_LensRadius(0.f), m_FocalDistance(0.f) {}
            SceneObjectPorjectiveCamera(const std::string& name) : SceneObjectCamera(name) {}
            SceneObjectPorjectiveCamera(const std::string& name, float shutterOpen, float shutterClose, 
                                        float lensr, float focald,
                                        Film* film, const Medium* medium)
                :SceneObjectCamera(name, shutterOpen, shutterClose, film, medium)
            {
                m_LensRadius = lensr;
                m_FocalDistance = focald;
            }

        protected:
            float m_LensRadius;
            float m_FocalDistance;
    };

    class SceneObjectOrthogonalCamera : public SceneObjectCamera
    //class SceneObjectOrthogonalCamera : public SceneObjectCamera
    {
        public:
            //using SceneObjectCamera::SceneObjectCamera;
            SceneObjectOrthogonalCamera() :
				SceneObjectCamera(), m_LensRadius(0.f), m_FocalDistance(0.f) {}
            SceneObjectOrthogonalCamera(const std::string& name) :
				SceneObjectCamera(name) {}
            SceneObjectOrthogonalCamera(const std::string& name, float shutterOpen, float shutterClose, 
                                        float lensr, float focald,
                                        Film* film, const Medium* medium)
                : SceneObjectCamera(name, shutterOpen, shutterClose, film, medium)
            {
				m_LensRadius = lensr;
				m_FocalDistance = focald;
			}

        float GenerateRay(const CameraSample& sample, Ray* ray) const;
        float GenerateRayDifferential(const CameraSample& sample, RayDifferential* rayDiff) const;

        friend std::ostream& operator<<(std::ostream& out, const SceneObjectOrthogonalCamera& obj);

        private:
            Vector3Df m_DxCamera, m_DyCamera;
			float m_LensRadius;
			float m_FocalDistance;
    };

    class SceneObjectPerspectiveCamera : public SceneObjectCamera
    //class SceneObjectPerspectiveCamera : public SceneObjectPorjectiveCamera
    {
        public:
            SceneObjectPerspectiveCamera() : SceneObjectCamera(), m_LensRadius(0.f), m_FocalDistance(0.f) {}
            SceneObjectPerspectiveCamera(const std::string& name, float fov = PI / 2.0f) : SceneObjectCamera(name), m_Fov(fov) {}
            SceneObjectPerspectiveCamera(const std::string& name, float fov, float shutterOpen, 
                                        float shutterClose, float lensr, float focald,
                                        Film* film, const Medium* medium
                                        ) ;

            float GenerateRay(const CameraSample& sample, Ray* ray) const;
            float GenerateRayDifferential(const CameraSample& sample, RayDifferential* rayDiff) const;
            
            Spectrum We(const Ray& ray, Point2Df* pRaster2 = nullptr) const;
            void Pdf_We(const Ray& ray, float* pdfPos, float* pdfDir) const;
            Spectrum Sample_Wi(const SceneObjectInteraction& ref, const Point2Df& sample,
                               Vector3Df* wi, float* pdf, Point2Df* pRaster,
                               VisibilityTester* vis) const;
            void SetParam(std::string& attrib, float param)
            {
                // TODO: handle fovx, fovy
                if (attrib == "fov")
                {
                    m_Fov = param;
                }
                SceneObjectCamera::SetParam(attrib, param);
            }

            
            float GetFov() const {return m_Fov;}
            friend std::ostream& operator<<(std::ostream& out, const SceneObjectPerspectiveCamera& obj);

        protected:
            float m_Fov;
            Vector3Df m_DxCamera, m_DyCamera;
            float A;
			float m_LensRadius;
			float m_FocalDistance;
    };

    class SceneObjectEnvironmentCamera : public SceneObjectCamera
    {
// SceneObjectPorjectiveCamera() : SceneObjectCamera(), m_LensRadius(0.f), m_FocalDistance(0.f) {}
// SceneObjectPorjectiveCamera(const std::string& name) : SceneObjectCamera(name) {}
// SceneObjectPorjectiveCamera(const std::string& name, float shutterOpen, float shutterClose, 
//                         float lensr, float focald,
//                         Film* film, const Medium* medium)
        public:
            SceneObjectEnvironmentCamera() : SceneObjectCamera() {}
            SceneObjectEnvironmentCamera(const std::string& name) : SceneObjectCamera(name) {}
            SceneObjectEnvironmentCamera(const std::string& name, float shutterOpen, float shutterClose, 
                                        Film* film, const Medium* medium)
                : SceneObjectCamera(name, shutterOpen, shutterClose, film, medium)
                {}
            
            float GenerateRay(const CameraSample& sample, Ray* ray) const;
    };
}

#endif
