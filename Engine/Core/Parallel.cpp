#include "Parallel.hpp"
#include <list>
#include <thread>
#include <algorithm>
#include "Math/Vector.hpp"
#include <iostream>

namespace Panda
{
    static std::vector<std::thread> s_Threads;
    static bool s_ShutdownThreads = false;
    class TaskBlock;
    static TaskBlock* s_WorkList = nullptr;
    static std::mutex s_WorkListMutex;

    class TaskBlock
    {
    public:
        TaskBlock(std::function<void(int64_t)> func1D, int64_t maxIndex,
                  int32_t chunkSize, uint64_t profilerState)
            : m_Func1D(std::move(func1D)),
              m_MaxIndex(maxIndex),
              m_ChunkSize(chunkSize),
              m_ProfilerState(profilerState){}
        TaskBlock(const std::function<void(Vector2Di)>& f, const Vector2Di& count,
                 uint64_t profilerState)
            : m_Func2D(f),
              m_MaxIndex(count[0] * count[1]),
              m_ChunkSize(1),
              m_ProfilerState(profilerState)
        {
            m_nX = count[0];
        }

        bool Finished() const 
        {
            return m_NextIndex >= m_MaxIndex && m_ActiveWorkers == 0;
        }
    public:
        std::function<void(int64_t)> m_Func1D;
        std::function<void(Vector2Di)> m_Func2D;
        const int64_t m_MaxIndex;
        const int32_t m_ChunkSize;
        uint64_t m_ProfilerState;
        int64_t m_NextIndex = 0;
        int32_t m_ActiveWorkers = 0;
        TaskBlock* m_Next = nullptr;
        int32_t m_nX = -1;
    };

    void Barrier::Wait()
    {
        std::unique_lock<std::mutex> lock(m_Mutex);
        if (--m_Count == 0)
            // This is the last thread to reach the barrier; 
            // Wake up all the other ones before exiting.
            m_CV.notify_all();
        else 
            // Otherwise there are still threads that haven't reach it.
            // Give up the lock and wait to be notified.
            m_CV.wait(lock, [this] {return m_Count == 0;});
    }
    static std::condition_variable s_WorkListCondition;

    int32_t MaxThreadCount()
    {
        return NumSystemCores();
    }

    int32_t NumSystemCores()
    {
        return (std::max)(1u, std::thread::hardware_concurrency());
    }

    static void WorkerThreadFunc(int32_t index, std::shared_ptr<Barrier> barrier)
    {
        std::cout << "Started execution in worker thread " << index;
        t_ThreadIndex = index;

        barrier->Wait();

        // Release our reference to the Barrier so that it's freed once all
        // of the threads have cleared it.
        barrier.reset();

        std::unique_lock<std::mutex> lock(s_WorkListMutex);
        while (!s_ShutdownThreads)
        {
            if (!s_WorkList)
            {
                s_WorkListCondition.wait(lock);
            }
            else 
            {
                // Get work from _s_WorkList_ and run task iterations
                TaskBlock& task = *s_WorkList;

                // Run a chunk of task iterations for _task_

                // Find the set of task iterations to run next
                int64_t indexStart = task.m_NextIndex;
                int64_t indexEnd = (std::min)(indexStart + task.m_ChunkSize, task.m_MaxIndex);

                // Update _task_ to reflect iterations this thread will run
                task.m_NextIndex = indexEnd;
                if (task.m_NextIndex == task.m_MaxIndex) s_WorkList = task.m_Next;
                task.m_ActiveWorkers++;

                // Run task indices in _[indexStart, indexEnd)_
                lock.unlock();
                for (int64_t index = indexStart; index < indexEnd; ++index)
                {
                    if (task.m_Func1D)
                        task.m_Func1D(index);
                    else if (task.m_Func2D)
                        task.m_Func2D(Vector2Di({(int32_t)(index % task.m_nX), (int32_t)(index / task.m_nX)}));
                }
                lock.lock();

                // Update _task_ to reflect completion of iterations
                task.m_ActiveWorkers--;
                if (task.Finished()) s_WorkListCondition.notify_all();
            }
        }

        std::cout << "Exiting worker thread " << index;
    }

    thread_local int32_t t_ThreadIndex;
    void ParallelInit()
    {
        int32_t nThreads = MaxThreadCount();
        t_ThreadIndex = 0;

        // Create a barrier so that we can be sure all worker threads are started.
        std::shared_ptr<Barrier> barrier = std::make_shared<Barrier>(nThreads);

        // Launch one fewer worker thread that than the total number we want doing
        // work, since the main thread helps out, too.
        for (int32_t i = 0; i < nThreads - 1; ++i)
            s_Threads.push_back(std::thread(WorkerThreadFunc, i + 1, barrier));
        
        barrier->Wait();
    }

    void ParallelCleanup()
    {
        if (s_Threads.empty()) return;

        {
            std::lock_guard<std::mutex> lock(s_WorkListMutex);
            s_ShutdownThreads = true;
            s_WorkListCondition.notify_all();
        }

        for (std::thread& thread : s_Threads) 
            thread.join();
        s_Threads.erase(s_Threads.begin(), s_Threads.end());
        s_ShutdownThreads = false;
    }

    void ParallelFor(std::function<void(int64_t)> func, int64_t count,
                     int32_t chunkSize)
    {
        assert(s_Threads.size() > 0 || MaxThreadCount() == 1);

        // Run iteration immediately if not using threads or if _count_ is small
        if (s_Threads.empty() || count < chunkSize)
        {
            for(int64_t i = 0; i < count; ++i) func(i);
            return;
        }

        // Create and enqueue _TaskBlck_ for this task
        TaskBlock task(std::move(func), count, chunkSize, 0);

        s_WorkListMutex.lock();
        task.m_Next = s_WorkList;
        s_WorkList = &task;
        s_WorkListMutex.unlock();

        // Notify worker threads of work to be done
        std::unique_lock<std::mutex> lock(s_WorkListMutex);
        s_WorkListCondition.notify_all();

        // Help out with parallel task iterations in the current thread.
        while (!task.Finished())
        {
            // Run a chunk of task iterations for _task_

            // Find the set of task iterations to run next
            int64_t indexStart = task.m_NextIndex;
            int64_t indexEnd = (std::min)(indexStart + task.m_ChunkSize, task.m_MaxIndex);

            // Update _task_ to reflect iterations this thread will run
            task.m_NextIndex = indexEnd;
            if (task.m_NextIndex == task.m_MaxIndex) s_WorkList = task.m_Next;
            task.m_ActiveWorkers++;

            // Run task indices in _[indexStart, indexEnd)_
            lock.unlock();
            for (int64_t index = indexStart; index < indexEnd; ++index)
            {
                if (task.m_Func1D)
                    task.m_Func1D(index);
                else if (task.m_Func2D)
                    task.m_Func2D(Vector2Di({(int32_t)(index % task.m_nX), (int32_t)(index / task.m_nX)}));
            }
            lock.lock();

            // Update _task_ to reflect completion of iterations
            task.m_ActiveWorkers--;
        }
    }
}