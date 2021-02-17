#pragma once
#include "BaseSceneObject.hpp"
#include "Core/Math/PandaMath.hpp"

namespace Panda
{
    class SceneObjectTransform : public BaseSceneObject
    {
        protected:
            Matrix4f m_Matrix;
            bool    m_IsSceneObjectOnly;

        public:
            SceneObjectTransform() : BaseSceneObject(SceneObjectType::kSceneObjectTypeTransform) 
			{m_Matrix.SetIdentity(); m_IsSceneObjectOnly = false;}
            SceneObjectTransform(const Matrix4f& matrix, const bool objectOnly = false) : SceneObjectTransform() 
			{m_Matrix = matrix; m_IsSceneObjectOnly = objectOnly;}

			operator Matrix4f() {return m_Matrix;}
			operator const Matrix4f() const {return m_Matrix;}

            void Update(const float amount)
            {
                // should not be used.
                assert(0);
            }

            void Update(const Vector3Df amount)
            {
                // should not be used.
                assert(0);
            }

            void Update(const Vector4Df amount)
            {
                // should not be used.
                assert(0);
            }

            friend std::ostream& operator<<(std::ostream& out, const SceneObjectTransform& obj);
    };

    class SceneObjectTranslation : public SceneObjectTransform
    {
        private:
            char m_Kind = 0;

        public:
			SceneObjectTranslation() { m_Type = SceneObjectType::kSceneObjectTypeTranslate; }
            SceneObjectTranslation(const char axis, const float amount)
				:SceneObjectTranslation()
            {
                m_Kind = axis;

                switch(axis)
                {
                    case 'x':
                        MatrixTranslation(m_Matrix, amount, 0.f, 0.f);
                        break;
                    case 'y':
                        MatrixTranslation(m_Matrix, 0.f, amount, 0.f);
                        break;
                    case 'z':
                        MatrixTranslation(m_Matrix, 0.f, 0.f, amount);
                        break;
                    default:
                        assert(false);
                        break;
                }
            }
            SceneObjectTranslation(const float x, const float y, const float z, const bool object_only = false)
				:SceneObjectTranslation()
            {
                m_Kind = 0;
                MatrixTranslation(m_Matrix, x, y, z);

                m_IsSceneObjectOnly = object_only;
            }
    };

    class SceneObjectRotation : public SceneObjectTransform
    {
        private:
            char m_Kind = 0;

        public:
			SceneObjectRotation() { m_Type = SceneObjectType::kSceneObjectTypeRotate; }
            SceneObjectRotation(const char axis, const float theta, const bool object_only = false)
				:SceneObjectRotation()
            {
                m_Kind = axis;
                m_IsSceneObjectOnly = object_only;

                switch (axis)
                {
                    case 'x':
                        MatrixRotationX(m_Matrix, theta);
                        break;
                    case 'y':
                        MatrixRotationY(m_Matrix, theta);
                        break;
                    case 'z':
                        MatrixRotationZ(m_Matrix, theta);
                        break;
                    default:
                        assert(false);
                        break;
                }
            }

            SceneObjectRotation(Vector3Df& axis, const float theta, const bool object_only = false)
				:SceneObjectRotation()
            {
                m_Kind = 0;
				axis = Normalize(axis);
                MatrixRotationAxis(m_Matrix, axis, theta);
                m_IsSceneObjectOnly = object_only;
            }

            SceneObjectRotation(const Quaternion& q, const bool object_only = false)
				:SceneObjectRotation()
            {
                m_Kind = 0;
                MatrixRotationQuaternion(m_Matrix, q);
                m_IsSceneObjectOnly = object_only;
            }
    };

    class SceneObjectScale : public SceneObjectTransform
    {
    private:
        char m_Kind = 0;

    public:
		SceneObjectScale() { m_Type = SceneObjectType::kSceneObjectTypeScale; }
        SceneObjectScale(const char axis, const float amount, const bool object_only = false)
			:SceneObjectScale()
        {
            m_Kind = axis;
            switch(axis)
            {
                case 'x':
                    MatrixScale(m_Matrix, amount, 1.f, 1.f);
                    break;
                case 'y':
                    MatrixScale(m_Matrix, 1.f, amount, 1.f);
                    break;
                case 'z':
                    MatrixScale(m_Matrix, 1.f, 1.f, amount);
                    break;
                default:
                    assert(0);
            }

			m_IsSceneObjectOnly = object_only;
        }

        SceneObjectScale(const float x, const float y, const float z, const bool object_only = false)
			:SceneObjectScale()
        {
            m_Kind = 0;
            MatrixScale(m_Matrix, x, y, z);
			m_IsSceneObjectOnly = object_only;
        }
    };

    void GenerateViewTransformAndInv(SceneObjectTransform& view, SceneObjectTransform& viewInv);
    void GenerateNDCTransformAndInv(SceneObjectTransform& ndc, SceneObjectTransform& ndcInv);
    void GenerateNDCToRasterTransformAndInv(SceneObjectTransform& ndcToRaster, SceneObjectTransform& rasterToNDC);
}