// monte-carlo-simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>
#include <cassert>
#include <functional> // For std::bind
#include <random>
#include <thread>
#include <mutex>
#include <condition_variable>


#include "utility.h"

void numCirclesPerThread(unsigned long long N, unsigned int threadId, std::vector<unsigned long long> & numCircles) {
    numCircles[threadId] = numberOfCircles(N);
}


double estimatePiMultiThreaded(unsigned long long OriginalN, unsigned int processor_count) {
    // estimate pi in simple form.
    unsigned long long PerThreadN = OriginalN / processor_count;
    std::vector<unsigned long long> numInCircles(processor_count, 0);

    // preallocate threads.
    std::vector<std::thread> threads;

    // Launch a group of threads
    unsigned long long N = OriginalN, numNSent = 0;
    for (unsigned int i = 0; i < processor_count; ++i) {
        unsigned long long currN = std::min(N, PerThreadN);
        //  numInCircleSimple(unsigned long long N, int threadId, std::vector<unsigned long long> & numCircles)
        threads.emplace_back(std::thread(numCirclesPerThread, currN, i, std::ref(numInCircles)));
        N -= currN;
        numNSent += currN;
    }
    assert(numNSent == OriginalN);

    // Join the threads with the main thread
    unsigned long long totalNumInCircles = 0;
    for (unsigned int i = 0; i < processor_count; ++i) {
        threads[i].join();
        totalNumInCircles += numInCircles[i];
    }
    return estimatePi(totalNumInCircles, OriginalN);
}


std::vector<std::unique_ptr<std::mutex>> mtx;
std::vector<std::unique_ptr<std::condition_variable>> cv;
std::mutex mtx2;
std::condition_variable cv2;
int numberOfResultsReady = 0;
constexpr int RUN_PER_BATCH = 100000;

void numCirclesPerThreadPersistent(unsigned int threadId, std::vector<unsigned long long> & requestQueues, std::vector<unsigned long long> & result) {
    while (true) {
        // wait for signal.
        // locking
        unsigned long long  N;
        {
            std::unique_lock<std::mutex> lock1(*mtx[threadId]);

            // waiting
            cv[threadId]->wait(lock1, [&, threadId] {
                return requestQueues[threadId] >= RUN_PER_BATCH;
            });

            // take request off the queue.
            N = requestQueues[threadId];
            requestQueues[threadId] = 0;
        }
        
        // number of circles.
        unsigned long long  totalNumInCircles = numberOfCircles(N);

        // save result
        result[threadId] = totalNumInCircles;
        
        // notify main.
        std::lock_guard<std::mutex> lock2(mtx2);

        numberOfResultsReady++;

        // signal main.
        cv2.notify_one();
    }
}

void estimatePiContinuously(unsigned int processor_count) {
    std::vector<unsigned long long> requestQueues(processor_count, 0);
    std::vector<unsigned long long> localRequestQueues(processor_count, 0);

    std::vector<unsigned long long> result(processor_count, 0);

    // preallocate threads.
    std::vector<std::thread> threads;

    // initialize synchronization primitive

    // Launch a group of threads
    for (unsigned int i = 0; i < processor_count; ++i) {
        //  numInCircleSimple(unsigned long long N, int threadId, std::vector<unsigned long long> & numCircles)
        threads.emplace_back(std::thread(numCirclesPerThreadPersistent, i, std::ref(requestQueues), std::ref(result)));

        cv.emplace_back(std::make_unique<std::condition_variable>());
        mtx.emplace_back(std::make_unique<std::mutex>());
    }

    unsigned int thread = 0;
    unsigned long long runs = 0, totalNumInCircles = 0;
    while (true) {
        localRequestQueues[thread]++;
        if (localRequestQueues[thread] >= RUN_PER_BATCH) {
            // signal thread to take requests, communicate by one variable.
            //std::cout << "signal thread " << thread << " take request " << localRequestQueues[thread] << std::endl;
            {
                std::lock_guard<std::mutex> lock1(*mtx[thread]);
                requestQueues[thread] = localRequestQueues[thread];
                localRequestQueues[thread] = 0;
                cv[thread]->notify_one();
            }
            
            if (thread == processor_count - 1) {
                // wait for all threads to finish.
                {
                    std::unique_lock<std::mutex> lock2(mtx2);
                    // waiting
                    cv2.wait(lock2, [&] { return numberOfResultsReady == processor_count; });
                }
                
                numberOfResultsReady = 0;
                for (unsigned int threadId = 0; threadId < processor_count; ++threadId) {
                    totalNumInCircles += result[threadId];
                    result[threadId] = 0;
                }
                std::cout << "\rEstimate of Pi = " << estimatePi(totalNumInCircles, runs) << " from " << runs << " runs." << std::endl;
            }
        }

        // increment 
        runs++;
        thread = (thread + 1) % processor_count;
    }
}

int main()
{
    std::cout << std::fixed << std::setprecision(12);

    const auto processor_count = std::thread::hardware_concurrency();
    //unsigned long long OriginalN = 500000000;
    //std::cout << "estimated pi, multi-threaded, N = " << OriginalN << " pi = " << estimatePiMultiThreaded(OriginalN, processor_count) << std::endl;
    
    estimatePiContinuously(processor_count);
}
