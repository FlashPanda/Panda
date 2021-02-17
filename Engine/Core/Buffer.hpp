#pragma once
#include <memory>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include "MemoryManager.hpp"

namespace Panda
{
    class Buffer
    {
        public:
            Buffer() : m_pData(nullptr), m_Size(0), m_Alignment(alignof(uint32_t)) {}   // alignof:Queries alignment requirements of a type since C++11
            Buffer(size_t size, size_t alignment = 4) : m_Size(size), m_Alignment(alignment) 
            {
                m_pData = reinterpret_cast<uint8_t*> (g_pMemoryManager->Allocate(size, alignment));
            }

            Buffer(const Buffer& rhs)
            {
                m_pData = reinterpret_cast<uint8_t*>(g_pMemoryManager->Allocate(rhs.m_Size, rhs.m_Alignment));
                memcpy(m_pData, rhs.m_pData, rhs.m_Size);
                m_Size = rhs.m_Size;
                m_Alignment = rhs.m_Alignment;
            }

            Buffer(Buffer&& rhs) noexcept
            {
                m_pData = rhs.m_pData;
                m_Size = rhs.m_Size;
                m_Alignment = rhs.m_Alignment;
                rhs.m_pData = nullptr;
                rhs.m_Size = 0;
                rhs.m_Alignment = 4;
            }

            Buffer& operator= (const Buffer& rhs) 
            {
                if (m_Size >= rhs.m_Size && m_Alignment == rhs.m_Alignment)
                {
                    memcpy(m_pData, rhs.m_pData, rhs.m_Size);
                }
                else
                {
                    if (m_pData) g_pMemoryManager->Free(m_pData, m_Size);
                    m_pData = reinterpret_cast<uint8_t*> (g_pMemoryManager->Allocate(rhs.m_Size, m_Alignment));
                    memcpy(m_pData, rhs.m_pData, rhs.m_Size);
                    m_Size = rhs.m_Size;
                    m_Alignment = rhs.m_Alignment;
                }
                
                return *this;
            }

            Buffer& operator=(Buffer&& rhs) noexcept
            {
                if (m_pData) g_pMemoryManager->Free(m_pData, m_Size);
                m_pData = rhs.m_pData;
                m_Size = rhs.m_Size;
                m_Alignment = rhs.m_Alignment;
                rhs.m_pData = nullptr;
                rhs.m_Size = 0;
                rhs.m_Alignment = 4;
                return *this;
            }

            ~Buffer() 
            {
                if (m_pData)
                {
                    g_pMemoryManager->Free(m_pData, m_Size);
                    m_pData = nullptr;
                }
            }

            uint8_t*    GetData() {return m_pData;}
            const uint8_t* GetData() const {return m_pData;}
            size_t GetDataSize() const {return m_Size;}
        public:
            uint8_t* m_pData;
            size_t m_Size;
            size_t m_Alignment;
    };
}