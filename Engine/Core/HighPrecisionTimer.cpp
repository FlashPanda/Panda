#include <windows.h>
#include "HighPrecisionTimer.hpp"

namespace Panda
{
    HighPrecisionTimer::HighPrecisionTimer()
        :m_MillisecondsPerCount(0.0), m_StartTime(0), m_StopTime(0)
        {
            int64_t countsPerSec;
            QueryPerformanceFrequency((LARGE_INTEGER*)&countsPerSec);
            m_MillisecondsPerCount = 1000.0 / (double)countsPerSec;
        }

	int64_t HighPrecisionTimer::TotalClockCount() const
	{
		return m_StopTime - m_StartTime;
	}

    float HighPrecisionTimer::TotalTime() const
    {
        return (float)((m_StopTime - m_StartTime) * m_MillisecondsPerCount);
    }

    void HighPrecisionTimer::Start()
    {
        QueryPerformanceCounter((LARGE_INTEGER*)&m_StartTime);
    }

    void HighPrecisionTimer::Stop()
    {
        QueryPerformanceCounter((LARGE_INTEGER*)&m_StopTime);
    }
}