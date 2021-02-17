#pragma once
#include <cstdint>
#include <iostream>
#include <stdio.h>
#include <string>

namespace Panda
{
	struct GfxConfiguration
	{
		/// 所有元素的构造函数
		/// \param[in] r 红色分量占多少位
		/// \param[in] g 绿色分量占多少位
		/// \param[in] b 蓝色分量占多少位
		/// \param[in] a alpha分量占多少位
		/// \param[in] d 深度缓存中的深度值占多少位
		/// \param[in] s 模板缓存中的模板值占多少位
		/// \param[in] msaa msaa采样数
		/// \param[in] width 屏幕宽度（以像素为单位）
		/// \param[in] height 屏幕高度（以像素为单位）
		GfxConfiguration(uint32_t r = 8, uint32_t g = 8,
			uint32_t b = 8, uint32_t a = 8,
			uint32_t d = 24, uint32_t s = 0, 
			uint32_t msaa = 0,
			uint32_t width = 1920, uint32_t height = 1080,
			const char* inAppName = "Panda")
			: redBits(r), greenBits(g), blueBits(b),
			  alphaBits(a), depthBits(d), stencilBits(s),
			  msaaSamples(msaa), screenWidth(width), screenHeight(height),
			  appName(inAppName)
			{
			}
			
		uint32_t redBits;		// 红色通道
		uint32_t greenBits;		// 绿色通道
		uint32_t blueBits;		// 蓝色通道
		uint32_t alphaBits;		// alpha通道
		uint32_t depthBits;		// 深度值
		uint32_t stencilBits;	// 模板值
		uint32_t msaaSamples;	// msaa采样数
		uint32_t screenWidth;	// 屏幕宽度（以像素为单位）
		uint32_t screenHeight;	// 屏幕高度（以像素为单位）
		const char* appName;	// 应用名
		
		friend std::ostream& operator<< (std::ostream& out, const GfxConfiguration& conf)
		{
			out << "App Name:" << conf.appName << std::endl;
			out << "GfxConfiguration:" <<
				" R:" << conf.redBits <<
				" G:" << conf.greenBits << 
				" B:" << conf.blueBits << 
				" A:" << conf.alphaBits << 
				" D:" << conf.depthBits << 
				" S:" << conf.stencilBits <<
				" M:" << conf.msaaSamples <<
				" N:" << conf.screenWidth <<
				" H:" << conf.screenHeight <<
				std::endl;
				
			return out;
		}
	};
}
