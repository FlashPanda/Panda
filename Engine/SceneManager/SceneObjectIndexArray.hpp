#pragma once
#include "SceneObjectTypeDef.hpp"

namespace Panda
{
    class SceneObjectIndexArray
    {
        protected:
            const uint32_t        m_MaterialIndex;
            const size_t          m_RestartIndex;
            const IndexDataType   m_DataType;

            const void*           m_pData;
            const size_t          m_DataSize;

            public:
                SceneObjectIndexArray(const uint32_t material_index = 0, const uint32_t restart_index = 0, const IndexDataType data_type = IndexDataType::kIndexDataTypeInt16, 
                    const void* data = nullptr, const size_t dataSize = 0) :
                    m_MaterialIndex(material_index), m_RestartIndex(restart_index), m_DataType(data_type),
                    m_pData(data), m_DataSize(dataSize)
                    {}
                SceneObjectIndexArray(SceneObjectIndexArray& arr) = default;
                SceneObjectIndexArray(SceneObjectIndexArray&& arr) = default;

                const uint32_t GetMaterialIndex() const {return m_MaterialIndex;}
                const IndexDataType GetIndexType() const {return m_DataType;}
                const void* GetData() const {return m_pData;}
                size_t GetDataSize() const
                {
                    size_t size = m_DataSize;

                    switch(m_DataType)
                    {
                        case IndexDataType::kIndexDataTypeInt8:
                            size *= sizeof(int8_t);
                            break;
                        case IndexDataType::kIndexDataTypeInt16:
                            size *= sizeof(int16_t);
                            break;
                        case IndexDataType::kIndexDataTypeInt32:
                            size *= sizeof(int32_t);
                            break;
                        case IndexDataType::kIndexDataTypeInt64:
                            size *= sizeof(int64_t);
                            break;
                        default:
                            size = 0;
                            assert(0);
                            break;
                    }

                    return size;
                }

                size_t GetIndexCount() const
                {
                    return m_DataSize;
                }


                friend std::ostream& operator<<(std::ostream& out, const SceneObjectIndexArray& obj);
                    
    };
}