#pragma once

#include "Interface/Hitable.hpp"

namespace Panda
{
    class HitableList : public Hitable
    {
        public:
            HitableList() {}
            HitableList(Hitable** l, int32_t n) : m_List(l), m_ListSize(n) {}

        public:
            virtual bool Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const;

        private:
            Hitable**   m_List;
            int32_t     m_ListSize;
    };

    bool HitableList::Hit(const Ray& r, float tMin, float tMax, HitRecord& rec) const
    {
        HitRecord tempRec;
        bool hitAnything = false;
        float closestSoFar = tMax;
        for (int32_t i = 0; i < m_ListSize; ++i)
        {
            if (m_List[i]->Hit(r, tMin, closestSoFar, tempRec))
            {
                hitAnything = true;
                closestSoFar = tempRec.t;
                rec = tempRec;
            }
        }

        return hitAnything;
    }
}