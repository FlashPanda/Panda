#include "SortAlgos.hpp"

namespace Panda
{
    int32_t LeftChild(int32_t index)
    {
        return 2 * index + 1;
    }
    int32_t RightChild(int32_t index)
    {
        return 2 * index + 2;
    }

	void CountingSortIncreasing(const std::vector<int32_t>& array, int32_t maxValue, std::vector<int32_t>& out)
	{
		int32_t count = maxValue + 1;
		std::vector<int32_t> countContainer(count, 0);

		// Counting the array
		for (int32_t i = 0; i < array.size(); ++i)
			++countContainer[array[i]];

		// Add the count
		for (int32_t i = 1; i < countContainer.size(); ++i)
			countContainer[i] = countContainer[i] + countContainer[i - 1];

		// place the numbers
		out.assign(array.size(), 0);
		for (int32_t i = array.size() - 1; i >= 0; --i)
		{
			out[countContainer[array[i]] - 1] = array[i];
			--countContainer[array[i]];
		}
	}
}