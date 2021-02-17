#ifndef _CORE_UTILITY_HPP_
#define _CORE_UTILITY_HPP_
#pragma once
#include "portable.hpp"

#ifndef PANDA_L1_CACHE_LINE_SIZE
#define PANDA_L1_CACHE_LINE_SIZE 64
#endif

namespace Panda {
//#define DUMP_DETAILS 1
//#define DEBUG

	template<class T>
	inline void SafeRelease(T **ppInterfaceToRelease)
	{
		if (*ppInterfaceToRelease != nullptr)
		{
			(*ppInterfaceToRelease)->Release();

			(*ppInterfaceToRelease) = nullptr;
		}
	}

	ENUM(Handness)
	{
		kHandnessRight,
		kHandnessLeft,
	};

	ENUM(DepthClipSpace)
	{
		kDepthClipZeroToOne,
		kDepthClipNegativeOneToOne,
	};

	extern Handness g_EngineHandness; // DO NOT change this. Engine handness is just a showcase.
	extern Handness g_ViewHandness;
	extern DepthClipSpace g_DepthClipSpace;

	template <typename T>
	inline constexpr bool IsPowerOf2(T v)
	{
		return v && !(v & (v - 1));
	}
}
#endif