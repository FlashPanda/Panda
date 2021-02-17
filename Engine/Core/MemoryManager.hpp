#pragma once
#include "Interface/IRuntimeModule.hpp"
#include "Allocator.hpp"
#include <new>
#include <vector>
#include <list>
#include <utility>
#include "Utility.hpp"

namespace Panda {
	static const uint32_t k_BlockSizes[] = {
		// 4字节增加
		4,	8,	12,	16,	20,	24,	28,	32,	36,	40,	44,	48,
		52,	56,	60,	64,	68,	72,	76,	80,	84,	88,	92,	96,
		
		// 32字节增加
		128, 160, 192, 224, 256, 288, 320, 352, 384,
		416, 448, 480, 512, 576, 608, 640,
		
		// 64字节增加
		704, 768, 832, 896, 960, 1024
	};
	static const uint32_t k_PageSize = 8192;	// 页尺寸
	static const uint32_t k_Alignment = 4;		// 对齐值
	
	static const uint32_t k_BlockSizeCount = sizeof (k_BlockSizes) / sizeof (k_BlockSizes[0]);	// 预配置的分配器数量
	
	static const uint32_t k_MaxBlockSize = k_BlockSizes[k_BlockSizeCount - 1]; // 预配置的分配器分配的最大内存数
		
	class MemoryManager : implements IRuntimeModule {
	public:
		// C++11的新机制，可变参数模板
		template<typename T, typename... Arguments>
		T* NewOne(Arguments... parameters)
		{
			return new (Allocate(sizeof(T)))T(parameters...);
		}

		template <typename T, typename... Arguments>
		T* New(size_t n, bool runContructor, Arguments... parameters)
		{
			T* ret = (T*)Allocate(n * sizeof(T));
			if (runContructor)
				for (size_t i = 0; i < n; ++i) 
					new (&ret[i]) T(parameters...);
			return ret;
		}
		
		template<typename T>
		void Delete(T* p, size_t n = 1)
		{
			for (size_t i = 0; i < n; ++i)
			{
				p->~T();
				p++;
			}
			
			Free(reinterpret_cast<void*>(p), sizeof(T * n));
		}
	public:
		~MemoryManager() {}
		// virtual ~MemoryManager() {}
		
		virtual int Initialize();
		virtual void Finalize();
		virtual void Tick();
		
		void* Allocate(size_t inSize);
		void* Allocate(size_t inSize, size_t alignment);
		void Free(void* p, size_t inSize);
		
	private:
		static Allocator* LookUpAllocator(size_t inSize);
		
	private:
		static Allocator* m_pAllocators;
		static uint32_t* m_pLookUpTable;
		static bool		m_IsInitialized;
		static std::vector<void*> m_LargeBlockPointers;
	};

	extern MemoryManager* g_pMemoryManager;

	// MemoryLocal is used for locally memory allocation.
	class MemoryLocal
	{
	public:
		MemoryLocal(size_t blockSize = 262144) : m_BlockSize(blockSize) {}
		~MemoryLocal()
		{
			if (m_pCurrentBlock)
				_aligned_free(m_pCurrentBlock);
			for (auto& block : m_UsedBlocks)
				if (block.second)
					_aligned_free(block.second);
			for (auto& block : m_AvailableBlocks)
				if (block.second)
					_aligned_free(block.second);
		}

		void* Alloc(size_t nBytes)
		{
			// Round up _nBytes_ to minimum machine alignment
			const int align = alignof(std::max_align_t);

			static_assert(IsPowerOf2(align), "Minimum alignment not a power of two");

			nBytes = (nBytes + align - 1) & ~(align - 1);
			if (m_CurrentBlockPos + nBytes > m_CurrentAllocSize)
			{
				// Add current block to _m_UsedBlocks_ list
				if (m_pCurrentBlock)
				{
					m_UsedBlocks.push_back(
						std::make_pair(m_CurrentAllocSize, m_pCurrentBlock));
					m_pCurrentBlock = nullptr;
					m_CurrentAllocSize = 0;
				}

				// Get new block of memory for _MemoryLocal_

				// Try to get memory block from _m_AvailableBlocks_
				for (auto iter = m_AvailableBlocks.begin();
					 iter != m_AvailableBlocks.end(); ++iter)
				{
					if (iter ->first >= nBytes)
					{
						m_CurrentAllocSize = iter->first;
						m_pCurrentBlock = iter->second;
						m_AvailableBlocks.erase(iter);
						break;
					}
				}

				if (!m_pCurrentBlock)
				{
					m_CurrentAllocSize = (std::max)(nBytes, m_BlockSize);
					m_pCurrentBlock = (uint8_t*)_aligned_malloc(m_CurrentAllocSize, PANDA_L1_CACHE_LINE_SIZE);
				}
				m_CurrentBlockPos = 0;
			}
			void* ret = m_pCurrentBlock + m_CurrentBlockPos;
			m_CurrentBlockPos += nBytes;
			return ret;
		}

		template <typename T>
		T* Alloc(size_t n = 1, bool runConstructor = true)
		{
			T* ret = (T *)Alloc(n * sizeof(T));
			if (runConstructor)
				for (size_t i = 0; i < n; ++i)
					new (ret + i) T();	// placement new
			return ret;
		}

		void Reset()
		{
			m_CurrentBlockPos = 0;
			m_AvailableBlocks.splice(m_AvailableBlocks.begin(), m_UsedBlocks);
		}

		size_t TotalAllocated() const 
		{
			size_t total = m_CurrentAllocSize;
			for (const auto& alloc : m_UsedBlocks) total += alloc.first;
			for (const auto& alloc : m_AvailableBlocks) total += alloc.first;
			return total;
		}

	private:
		MemoryLocal(const MemoryLocal& ) = delete;
		MemoryLocal& operator=(const MemoryLocal& ) = delete;
		const size_t m_BlockSize;
		size_t m_CurrentBlockPos = 0;
		size_t m_CurrentAllocSize = 0;
		uint8_t* m_pCurrentBlock = nullptr;
		std::list<std::pair<size_t, uint8_t*>> m_UsedBlocks;
		std::list<std::pair<size_t, uint8_t*>> m_AvailableBlocks;
	};
}