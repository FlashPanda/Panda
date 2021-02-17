#ifndef _PARALLEL_HPP_
#define _PARALLEL_HPP_

#pragma once
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>


namespace Panda
{
    // Barrier class is much more the same as std::experimental::barrier
    // We defined it because std::experimental::barrier is not in standard 
    // library in C++. However, std::experimental::barrier will become standard
    // in C++ 20.
    class Barrier
    {
    public:
        Barrier(int32_t count) : m_Count (count) {}
        ~Barrier() {}
        void Wait();

    private:
        std::mutex m_Mutex;
        std::condition_variable m_CV;
        int32_t m_Count;
    };

    int32_t MaxThreadCount();
    int32_t NumSystemCores();

    void ParallelInit();
    void ParallelCleanup();

    void ParallelFor(std::function<void(int64_t)> func, int64_t count,
                     int32_t chunkSize = 1);
    extern thread_local int32_t t_ThreadIndex;

}

#endif