#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "Allocator.hpp"


#ifndef ALIGN
#define ALIGN(n, a) (((n) + (a - 1)) & ~((a) - 1))		// 获取对齐的最小值
#endif

Panda::Allocator::Allocator()
	:m_BlockSize(0), m_PageSize(0), m_BlockCountPerPage(0), m_AlignmentSize(0), 
	m_pFreeBlockList(nullptr),m_FreeBlockCount(0),m_BlockCount(0),
	 m_pPageList(nullptr), m_PageCount(0)
{
}

Panda::Allocator::Allocator(uint32_t inPageSize, uint32_t inBlockSize, uint32_t inAlignment) 
	:m_pFreeBlockList(nullptr), m_pPageList(nullptr)
{
	Reset(inPageSize, inBlockSize, inAlignment);
}

Panda::Allocator::~Allocator()
{
	FreeAll();
}

void* Panda::Allocator::Allocate() {
	if (m_pFreeBlockList == nullptr) {
		// 分配一页内存
		PageHeader* pNewPage = reinterpret_cast<PageHeader*>(new uint8_t[m_PageSize]);
		m_FreeBlockCount += m_BlockCountPerPage;
		m_BlockCount += m_BlockCountPerPage;
		m_PageCount++;
		
		pNewPage->pNext = m_pPageList;
		m_pPageList = pNewPage;
	 	
		// 将所有内存块串联起来
		BlockHeader* pBlockStart = m_pPageList->BlockStart();
		m_pFreeBlockList = pBlockStart;
		for (size_t i = 0; i < m_BlockCountPerPage - 1; ++i) {
			pBlockStart->pNext = NextBlock(pBlockStart);
			pBlockStart = pBlockStart->pNext;
		}
		pBlockStart->pNext = nullptr;
	}
	
	// 取一块内存返回
	BlockHeader* pBlock = m_pFreeBlockList;
	m_pFreeBlockList = m_pFreeBlockList->pNext;
	m_FreeBlockCount--;
	
	return reinterpret_cast<void*>(pBlock);
}

Panda::BlockHeader* Panda::Allocator::NextBlock(BlockHeader* pCurrentBlock) {
	return reinterpret_cast<BlockHeader*>( reinterpret_cast<uint8_t*> (pCurrentBlock) + m_BlockSize);
}

// 重置分配器
void Panda::Allocator::Reset(uint32_t inPageSize, uint32_t inBlockSize, uint32_t inAlignment) {
	FreeAll();
	
	// 块的大小必须大于块头的大小，并且需要对齐
	size_t minBlockSize = inBlockSize > sizeof(BlockHeader)? inBlockSize : sizeof(BlockHeader);
	m_BlockSize = ALIGN(minBlockSize, inAlignment);
	
	m_AlignmentSize = m_BlockSize - minBlockSize;
	
	// 一页至少一个块
	m_PageSize = inPageSize > (m_BlockSize + sizeof (PageHeader))? inPageSize : (m_BlockSize + sizeof(PageHeader));
	
	m_BlockCountPerPage = (m_PageSize - sizeof(PageHeader)) / m_BlockSize;
	
}

void Panda::Allocator::Free(void* p) 
{
	BlockHeader* pHeader = reinterpret_cast<BlockHeader*>(p);
	pHeader->pNext = m_pFreeBlockList;
	m_pFreeBlockList = pHeader;
	m_FreeBlockCount++;
}

void Panda::Allocator::FreeAll() 
{
	PageHeader* pFree = m_pPageList;
	while (pFree) {
		PageHeader* pTemp = pFree;
		pFree = pFree->pNext;
		
		delete[] reinterpret_cast<uint8_t*>(pTemp);
	}
	
	// 只重置这两个属性因为这两个属性是要用于判断的，其他的都是保存作用
	m_pPageList = nullptr;
	m_pFreeBlockList = nullptr;
}

#if defined(_DEBUG)
void Panda::Allocator::FillFreePage(PageHeader* pPage)
{
	pPage->pNext = nullptr;

	BlockHeader* pBlock = pPage->BlockStart();
	for (uint32_t i = 0; i < m_BlockCountPerPage; ++i)
	{
		FillFreeBlock(pBlock);
		pBlock = NextBlock(pBlock);
	}
}

void Panda::Allocator::FillFreeBlock(BlockHeader* pBlock)
{
	// block header + data
	memset(pBlock, PATTERN_FREE, m_BlockSize - m_AlignmentSize);

	// alignment
	memset(reinterpret_cast<uint8_t*>(pBlock) + m_BlockSize - m_AlignmentSize, PATTERN_ALIGN, m_AlignmentSize);
}

void Panda::Allocator::FillAllocatedBlock(BlockHeader* pBlock)
{
	// block header + data
	memset(pBlock, PATTERN_FREE, m_BlockSize - m_AlignmentSize);

	// alignment
	memset(reinterpret_cast<uint8_t*>(pBlock) + m_BlockSize - m_AlignmentSize, PATTERN_ALIGN, m_AlignmentSize);
}
#endif
