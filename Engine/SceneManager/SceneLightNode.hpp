#pragma once
#include "BaseSceneNode.hpp"

namespace Panda
{
    class SceneLightNode : public BaseSceneNode
    {
        protected:
            bool    m_IsShadow;

        public:
            using BaseSceneNode::BaseSceneNode;

            void SetIfCastShadow(bool shadow) {m_IsShadow = shadow;}
            const bool CastShadow() {return m_IsShadow;}
    };
}