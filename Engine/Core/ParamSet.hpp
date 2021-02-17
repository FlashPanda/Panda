#pragma once

#include <map>
#include <vector>
#include <memory>
#include <string>
#include "Math/Vector.hpp"

namespace Panda
{
    template <typename T>
    struct ParamSetItem
    {
        ParamSetItem(const std::string& name, std::unique_ptr<T[]> val, int32_t nValues = 1)
        : name(name), values(std::move(val)), nValues(nValues)
        {}

        const std::string name;
        const std::unique_ptr<T[]> values;
        const int32_t nValues;
        mutable bool lookedUp = false;
    };

    class ParamSet
    {
    public:
        ParamSet() {}
        bool EraseFloat(const std::string& name)
        {
            for (size_t i = 0; i < floats.size(); ++i)
                if(floats[i]->name == name)
                {
                    floats.erase(floats.begin() + i);
                    return true;
                }
            return false;
        }
        void AddFloat(const std::string& name, std::unique_ptr<float[]> v, int32_t nValues = 1)
        {
            EraseFloat(name);
            floats.emplace_back(
                new ParamSetItem<float>(name, std::move(v), nValues)
            );
        }

        bool EraseInt(const std::string& name)
        {
            for (size_t i = 0; i <ints.size(); ++i)
                if (ints[i]->name == name)
                {
                    ints.erase(ints.begin() + i);
                    return true;
                }
            return false;
        }
        void AddInt(const std::string& name, std::unique_ptr<int32_t[]> v, int32_t nValues = 1)
        {
            EraseInt(name);
            ints.emplace_back(
                new ParamSetItem<int32_t>(name, std::move(v), nValues)
            );
        }

        bool EraseBool(const std::string& name)
        {
            for (size_t i = 0; i < bools.size(); ++i)
                if (bools[i]->name == name)
                {
                    bools.erase(bools.begin() + i);
                    return true;
                }
            return false;
        }
        void AddBool(const std::string& name, std::unique_ptr<bool[]> v, int32_t nValues = 1)
        {
            EraseBool(name);
            bools.emplace_back(
                new ParamSetItem<bool>(name, std::move(v), nValues)
            );
        }

        bool ErasePoint2Df(const std::string& name)
        {
            for (size_t i = 0; i < point2Dfs.size(); ++i)
                if (point2Dfs[i]->name == name)
                {
                    point2Dfs.erase(point2Dfs.begin() + i);
                    return true;
                }
            return false;
        }
        void AddPoint2Df(const std::string& name, std::unique_ptr<Vector2Df[]> v, int32_t nValues = 1)
        {
            ErasePoint2Df(name);
            point2Dfs.emplace_back(
                new ParamSetItem<Vector2Df>(name, std::move(v), nValues)
            );
        }

        bool EraseVector2Df(const std::string& name)
        {
            for (size_t i = 0; i < vector2Dfs.size(); ++i)
                if (vector2Dfs[i]->name == name)
                {
                    vector2Dfs.erase(vector2Dfs.begin() + i);
                    return true;
                }
            return false;
        }
        void AddVector2Df(const std::string& name, std::unique_ptr<Vector2Df[]> v, int32_t nValues = 1)
        {
            EraseVector2Df(name);
            vector2Dfs.emplace_back(
                new ParamSetItem<Vector2Df>(name, std::move(v), nValues)
            );
        }

        bool ErasePoint3Df(const std::string& name)
        {
            for (size_t i = 0; i < point3Dfs.size(); ++i)
                if (point3Dfs[i]->name == name)
                {
                    point3Dfs.erase(point3Dfs.begin() + i);
                    return true;
                }
            return false;
        }
        void AddPoint3Df(const std::string& name, std::unique_ptr<Vector3Df[]> v, int32_t nValues = 1)
        {
            ErasePoint3Df(name);
            point3Dfs.emplace_back(
                new ParamSetItem<Vector3Df>(name, std::move(v), nValues)
            );
        }

        bool EraseVector3Df(const std::string& name)
        {
            for (size_t i = 0; i < vector3Dfs.size(); ++i)
                if (vector3Dfs[i]->name == name)
                {
                    vector3Dfs.erase(vector3Dfs.begin() + i);
                    return true;
                }
            return false;
        }
        void AddVector3Df(const std::string& name, std::unique_ptr<Vector3Df[]> v, int32_t nValues = 1)
        {
            EraseVector3Df(name);
            vector3Dfs.emplace_back(
                new ParamSetItem<Vector3Df>(name, std::move(v), nValues)
            );
        }

        bool EraseNormal(const std::string& name)
        {
            for (size_t i = 0; i < normals.size(); ++i)
                if (normals[i]->name == name)
                {
                    normals.erase(normals.begin() + i);
                    return true;
                }
            return false;
        }
        void AddNormal(const std::string& name, std::unique_ptr<Vector3Df[]> v, int32_t nValues = 1)
        {
            EraseNormal(name);
            normals.emplace_back(
                new ParamSetItem<Vector3Df>(name, std::move(v), nValues)
            );
        }

        bool EraseString(const std::string& name)
        {
            for (size_t i = 0; i < strings.size(); ++i)
                if (strings[i]->name == name)
                {
                    strings.erase(strings.begin() + i);
                    return true;
                }
            return false;
        }
        void AddString(const std::string& name, std::unique_ptr<std::string[]> v, int32_t nValues)
        {
            EraseString(name);
            strings.emplace_back(
                new ParamSetItem<std::string>(name, std::move(v), nValues)
            );
        }

        float FindOneFloat(const std::string& name, float d) const 
        {
            for (const auto& t : floats)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        int32_t FindOneInt(const std::string& name, int32_t d) const
        {
            for (const auto& t : ints)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        bool FindOneBool(const std::string& name, bool d) const
        {
            for (const auto& t : bools)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        Vector2Df FindOnePoint2Df(const std::string& name, const Vector2Df& d) const
        {
            for (const auto& t : point2Dfs)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        Vector2Df FindOneVector2Df(const std::string& name, const Vector2Df& d) const
        {
            for (const auto& t : vector2Dfs)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        Vector3Df FindOnePoint3Df(const std::string& name, const Vector3Df& d) const
        {
            for (const auto& t : point3Dfs)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        Vector3Df FindOneVector3Df(const std::string& name, const Vector3Df& d) const
        {
            for (const auto& t : vector3Dfs)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        Vector3Df FindOneNormal(const std::string& name, const Vector3Df& d) const
        {
            for (const auto& t : normals)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        std::string FindOneString(const std::string& name, const std::string& d) const
        {
            for (const auto& t : strings)
                if (t->name == name && t->nValues == 1)
                {
                    t->lookedUp = true;
                    return t->values[0];
                }
            return d;
        }

        const float* FindFloat(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : floats)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const int32_t* FindInt(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : ints)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const bool* FindBool(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : bools)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const Vector2Df* FindPoint2Df(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : point2Dfs)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const Vector2Df* FindVector2Df(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : vector2Dfs)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const Vector3Df* FindPoint3Df(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : point3Dfs)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const Vector3Df* FindVector3Df(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : vector3Dfs)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const Vector3Df* FindNormal(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : normals)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        const std::string* FindString(const std::string& name, int32_t& n) const 
        {
            for (const auto& t : strings)
                if (t->name == name)
                {
                    n = t->nValues;
                    t->lookedUp = true;
                    return t->values.get();
                }
            return nullptr;
        }

        void Clear()
        {
#define DEL_PARAMS(name) (name).erase((name).begin(), (name).end())
            DEL_PARAMS(ints);
            DEL_PARAMS(floats);
            DEL_PARAMS(bools);
            DEL_PARAMS(point2Dfs);
            DEL_PARAMS(vector2Dfs);
            DEL_PARAMS(point3Dfs);
            DEL_PARAMS(vector3Dfs);
            DEL_PARAMS(normals);
            DEL_PARAMS(strings);
#undef DEL_PARAMS
        }
    
    private:
        std::vector<std::shared_ptr<ParamSetItem<bool>>> bools;
        std::vector<std::shared_ptr<ParamSetItem<int32_t>>> ints;
        std::vector<std::shared_ptr<ParamSetItem<float>>> floats;
        std::vector<std::shared_ptr<ParamSetItem<Vector2Df>>> point2Dfs;
        std::vector<std::shared_ptr<ParamSetItem<Vector2Df>>> vector2Dfs;
        std::vector<std::shared_ptr<ParamSetItem<Vector3Df>>> point3Dfs;
        std::vector<std::shared_ptr<ParamSetItem<Vector3Df>>> vector3Dfs;
        std::vector<std::shared_ptr<ParamSetItem<Vector3Df>>> normals;
        std::vector<std::shared_ptr<ParamSetItem<std::string>>> strings;
    };
}