#pragma once
#ifndef _SORT_ALOGS_HPP
#define _SORT_ALOGS_HPP

#include <vector>
#include "MathUtility.hpp"

namespace Panda
{
	int32_t LeftChild(int32_t index);
	int32_t RightChild(int32_t index);

	////////////////////////////////////////////////////////////////////
	template <typename T>
	void MaxHeapify(std::vector<T>& unheaped, int32_t index)
	{
		int32_t maxIndex = unheaped.size() - 1;
		int32_t l = LeftChild(index);
		int32_t r = RightChild(index);
		int32_t largest = index;
		if (l <= maxIndex && unheaped[l] > unheaped[largest])
			largest = l;
		if (r <= maxIndex && unheaped[r] > unheaped[largest])
			largest = r;
		if (largest != index)
		{
			std::swap(unheaped[index], unheaped[largest]);
			MaxHeapify(unheaped, largest);
		}
	}

	template <typename T>
	std::vector<T> BuildMaxHeap(const std::vector<T>& array)
	{
		std::vector<T> temp = array;
		int32_t count = std::floor(array.size() / 2);
		for (int32_t i = count - 1; i >= 0; --i)
			MaxHeapify(temp, i);
		return temp;
	}

	template <typename T>
	T GetMaxHeapValue(const std::vector<T>& heap)
	{
		assert(heap.size() > 0);
		return heap[0];
	}

	template <typename T>
	T ExtractMaxHeapValue(std::vector<T>& heap)
	{
		assert(heap.size() > 0);
		T value = heap[0];
		std::swap(heap[0], heap[heap.size() - 1]);
		heap.pop_back();
		MaxHeapify(heap, 0);
		return value;
	}

	template <typename T>
	std::vector<T> HeapSortDecreasing(const std::vector<T>& array)
	{
		std::vector<T> maxHeap = BuildMaxHeap(array);
		std::vector<T> sortedArray;
		while (maxHeap.size() > 0)
			sortedArray.push_back(ExtractMaxHeapValue(maxHeap));
		return sortedArray;
	}

	///////////////////////////////////////////////////////////////////
	template <typename T>
	void MinHeapify(std::vector<T>& unheaped, int32_t index)
	{
		int32_t maxIndex = unheaped.size() - 1;
		int32_t l = LeftChild(index);
		int32_t r = RightChild(index);
		int32_t smallest = index;
		if (l <= maxIndex && unheaped[l] < unheaped[smallest])
			smallest = l;
		if (r <= maxIndex && unheaped[r] < unheaped[smallest])
			smallest = r;
		if (smallest != index)
		{
			std::swap(unheaped[index], unheaped[smallest]);
			MinHeapify(unheaped, smallest);
		}
	}

	template <typename T>
	std::vector<T> BuildMinHeap(const std::vector<T>& array)
	{
		std::vector<T> temp = array;
		int32_t count = temp.size() / 2;
		for (int32_t i = count - 1; i >= 0; --i)
			MinHeapify(temp, i);
		return temp;
	}

	template <typename T>
	T GetMinHeapValue(const std::vector<T>& heap)
	{
		assert(heap.size() > 0);
		return heap[0];
	}

	template <typename T>
	T ExtractMinHeapValue(std::vector<T>& heap)
	{
		assert(heap.size() > 0);
		T value = heap[0];
		std::swap(heap[0], heap[heap.size() - 1]);
		heap.pop_back();
		MinHeapify(heap, 0);
		return value;
	}

	template <typename T>
	std::vector<T> HeapSortIncreasing(const std::vector<T>& array)
	{
		std::vector<T> minHeap = BuildMinHeap(array);
		std::vector<T> sortedArray;
		while (minHeap.size() > 0)
			sortedArray.push_back(ExtractMinHeapValue(minHeap));
		return sortedArray;
	}

    //////////////////////////////////////////////////////////////
	template <typename T>
	int32_t PositionDecreasing(std::vector<T>& array, int32_t from, int32_t to, int32_t target)
	{
		if (target > to)
			target = to;
		else
			std::swap(array[target], array[to]);
		T key = array[to];
		int32_t i = from - 1;
		for (int32_t j = from; j < to; ++j)
		{
			if (array[j] >= key)
			{
				++i;
				std::swap(array[i], array[j]);
			}
		}
		std::swap(array[i + 1], array[to]);
		return i + 1;
	}

	// In place sorting
	template <typename T>
	void QuickSortDecreasing(std::vector<T>& array, int32_t from, int32_t to)
	{
		if (from < to)
		{
			// To randomize the pivot 
			int32_t pivot = rand() % (to - from + 1) + from;
			int32_t position = PositionDecreasing(array, from, to, pivot);
			QuickSortDecreasing(array, from, position - 1);
			QuickSortDecreasing(array, position + 1, to);
		}
	}

	template <typename T>
	int32_t PositionIncreasing(std::vector<T>& array, int32_t from, int32_t to, int32_t target)
	{
		if (target > to)
			target = to;
		else
			std::swap(array[target], array[to]);
		T key = array[to];
		int32_t i = from - 1;
		for (int32_t j = from; j < to; ++j)
		{
			if (array[j] <= key)
			{
				++i;
				std::swap(array[i], array[j]);
			}
		}
		std::swap(array[i + 1], array[to]);
		return i + 1;
	}

	// In place sorting
	template <typename T>
	void QuickSortIncreasing(std::vector<T>& array, int32_t from, int32_t to)
	{
		if (from < to)
		{
			int32_t pivot = rand() % (to - from + 1) + from;
			int32_t position = PositionIncreasing(array, from, to, pivot);
			QuickSortIncreasing(array, from, position - 1);
			QuickSortIncreasing(array, position + 1, to);
		}
	}
	////////////////////////////////////////////////////////////////////////
	
	void CountingSortIncreasing(const std::vector<int32_t>& array, int32_t maxValue, std::vector<int32_t>& out);

	// Counting sort is used for integers. So, T could be either int32_t or int64_t.
	// There's no way to do with float , double and so on.
	template <typename T>
	void CountingSort(const std::vector<T>& in, std::vector<T>& out, SortOrder order, int32_t maxValue)
	{
		int64_t count = maxValue + 1;
		std::vector<T> countContainter(count, (T)0);

		// Counting the array
		for (int32_t i = 0; i < in.size(); ++i)
			++countContainter[in[i]];

		// Add the count
		if (order == SortOrder::Increasing)
			for (int32_t i = 1; i < countContainter.size(); ++i)
				countContainter[i] = countContainter[i] + countContainter[i - 1];
		else 
			for (int32_t i = countContainter.size() - 2; i >= 0; --i)
				countContainter[i] = countContainter[i] + countContainter[i + 1];

		// place the numbers
		out.assign(in.size(), (T)0);
		if (order == SortOrder::Increasing)
		{
			for (int32_t i = in.size() - 1; i >= 0; --i)
			{
				out[countContainter[in[i]] - 1] = in[i];
				--countContainter[in[i]];				
			}			
		}
		else 
		{
			for (int32_t i = 0; i < in.size(); ++i)
			{
				out[countContainter[in[i]] - 1] = in[i];
				--countContainter[in[i]];
			}
		}
	}

	///////////////////////////////////////////////////////////////////////
	// Note that all the inputs, from, to, target, are all indices.
	template <typename T>
	T RandomizeSelect(std::vector<T>& array, int32_t from, int32_t to, int32_t target)
	{
		if (from == to)
			return array[from];

		int32_t pivot = rand() % (to - from + 1) + from;
		int32_t q = PositionIncreasing(array, from, to, pivot);
		int32_t k = q - from;

		if (target == k)
			return array[q];
		else if (target < k)
			return RandomizeSelect(array, from, q - 1, target);
		else
			return RandomizeSelect(array, q + 1, to, target - (k + 1));
	}

	//////////////////////////////////////////////////////////////////////
	template <typename T>
	void InsertSortingDecreasing(std::vector<T>& array)
	{
		for (int32_t i = 1; i < array.size(); ++i)
		{
			T key = array[i];
			int32_t j = i - 1;
			for (; j >= 0; --j)
			{
				if (key >= array[j])
					std::swap(array[j], array[j + 1]);
				else 
					break;
			}
			array[j + 1] = key;
		}
	}
}

#endif