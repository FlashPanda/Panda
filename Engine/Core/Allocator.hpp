#pragma once

#include <stddef.h>
#include <stdint.h>

#ifndef ALIGN
#define ALIGN(n, a) (((n) + (a - 1)) & ~((a) - 1))		// 获取对齐的最小值
#endif

namespace Panda{
	struct BlockHeader {
		BlockHeader* pNext;
	};
	
	struct PageHeader {
		PageHeader* pNext;
		BlockHeader* BlockStart() {
			return reinterpret_cast<BlockHeader*>(this + 1);	// PageHeader之后就是块的起始位置
		}
	};
	
	class Allocator{
	public:
		// debug patterns
		static const uint8_t PATTERN_ALIGN = 0xFC;
		static const uint8_t PATTERN_ALLOC = 0xFD;
		static const uint8_t PATTERN_FREE = 0xFE;

		Allocator();
		Allocator(uint32_t inPageSize, uint32_t inBlockSize, uint32_t inAlignment);
		~Allocator();
	
		void Reset(uint32_t inPageSize, uint32_t inBlockSize, uint32_t inAlignment);
		void* Allocate();
		void Free(void* p);
		void FreeAll();
		
	private:
#if defined(_DEBUG)
		// fill a free page with debug patterns
		void FillFreePage(PageHeader* pPage);

		// fill a block with debug patterns
		void FillFreeBlock(BlockHeader* pBlock);

		// fill on allocated block with debug patterns
		void FillAllocatedBlock(BlockHeader* pBlock);
#endif

		BlockHeader* NextBlock(BlockHeader* pCurrentBlock);
		// 禁用拷贝构造函数和赋值操作符
		Allocator(const Allocator& clone);
		Allocator& operator=(const Allocator& rhs);
		
	private:
		size_t		m_BlockSize;			// 块的尺寸
		size_t		m_PageSize;				// 页的尺寸
		size_t		m_BlockCountPerPage;	// 每页的块数
		size_t		m_AlignmentSize;		// 对齐尺寸
		
		
		BlockHeader*	m_pFreeBlockList;		// 空余块列表
		size_t		m_FreeBlockCount;		// 空余的块数
		size_t		m_BlockCount;			// 总块数
		
		PageHeader*		m_pPageList;			// 页列表
		size_t		m_PageCount;			// 页数
	};
}