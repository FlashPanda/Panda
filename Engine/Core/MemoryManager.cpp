#include "MemoryManager.hpp"
#include <cstdlib>

#ifndef ALIGN
#define ALIGN(x,a)	(((x) + ((a) - 1)) & ~((a) - 1))
#endif

namespace Panda
{
	Allocator* MemoryManager::m_pAllocators = nullptr;
	uint32_t* MemoryManager::m_pLookUpTable = nullptr;
	bool MemoryManager::m_IsInitialized = false;
	std::vector<void*> MemoryManager::m_LargeBlockPointers = std::vector<void*>();

	int Panda::MemoryManager::Initialize() {
		if (!m_IsInitialized) 
		{
			// 初始化分配器
			m_pAllocators = new Allocator[k_BlockSizeCount];
			for (size_t i = 0; i < k_BlockSizeCount; ++i) 
			{
				m_pAllocators[i].Reset(k_PageSize, k_BlockSizes[i], k_Alignment);
			}

			// 初始化查询表
			m_pLookUpTable = new uint32_t[k_MaxBlockSize + 1];
			size_t j = 0;
			for (size_t i = 0; i <= k_MaxBlockSize; ++i) 
			{
				if (i > k_BlockSizes[j]) 
				{
					++j;
				}
				
				m_pLookUpTable[i] = j;
			}
		
			m_IsInitialized = true;
		}

		return 0;
	}

	void Panda::MemoryManager::Finalize()
	{
		delete[] m_pAllocators;
		delete[] m_pLookUpTable;
		for (int32_t i = 0; i < m_LargeBlockPointers.size(); ++i)
			free(m_LargeBlockPointers[i]);
		m_LargeBlockPointers.clear();
		m_IsInitialized = false;
	}

	void Panda::MemoryManager::Tick()
	{
	}

	void* Panda::MemoryManager::Allocate(size_t inSize)
	{
		Allocator* pAlloc = LookUpAllocator(inSize);
		
		if (pAlloc)
			return pAlloc->Allocate();
		else
		{
			void* p = malloc(inSize);
			m_LargeBlockPointers.push_back(p);
			return p;
		}
			
	}

	void* Panda::MemoryManager::Allocate(size_t inSize, size_t alignment)
	{
		/* Here are "From 0 Series" codes */
		// uint8_t* p;
		// inSize += alignment;
		// Allocator* pAlloc = LookUpAllocator(inSize);
		// if (pAlloc)
		// 	p = reinterpret_cast<uint8_t*> (pAlloc->Allocate());
		// else
		// {
		// 	p = reinterpret_cast<uint8_t*> (malloc(inSize));
		// }

		// p = reinterpret_cast<uint8_t*>(ALIGN(reinterpret_cast<size_t>(p), alignment));

		// return static_cast<void*>(p);
		
		/* Here are my codes */
		uint8_t* p;
		size_t size = ALIGN(inSize, alignment);
		return Allocate(size);
	}

	void Panda::MemoryManager::Free(void* p, size_t inSize) 
	{
		if (m_IsInitialized)
		{
			Allocator* pAlloc = LookUpAllocator(inSize);
			
			if (pAlloc)
				pAlloc->Free(p);
			else
			{
				for (int32_t i = 0; i < m_LargeBlockPointers.size(); ++i)
				{
					if (m_LargeBlockPointers[i] == p)
					{
						m_LargeBlockPointers.erase(m_LargeBlockPointers.begin() + i);
						break;
					}
				}
				free(p);		
			}	
		}
	}

	Panda::Allocator* Panda::MemoryManager::LookUpAllocator(size_t inSize)
	{
		if (inSize <= k_MaxBlockSize)
			return m_pAllocators + m_pLookUpTable[inSize];
		else
			return nullptr;
	}



	
}


