#pragma once
#include <cstdint>

namespace Panda
{
    class HighPrecisionTimer
    {
        public:
            HighPrecisionTimer();

			int64_t TotalClockCount() const;
            float TotalTime() const; // in milliseconds
            void Start();
            void Stop();


        private:
            double m_MillisecondsPerCount;

            int64_t m_StartTime;
            int64_t m_StopTime;
    };
}