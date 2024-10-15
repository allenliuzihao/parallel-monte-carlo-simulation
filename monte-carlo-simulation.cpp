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

#include <algorithm>

#include "utility.h"

void numCirclesPerThread(unsigned long long N, unsigned int threadId, std::vector<unsigned long long> & numCircles) {
    // Create a random number generator
    std::random_device rd;  // Seed
    std::default_random_engine gen(rd()); // Mersenne Twister engine

    // Define the range
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    numCircles[threadId] = numberOfCircles(N, gen, dis);
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
constexpr int RUN_PER_BATCH = 1000000;

void numCirclesPerThreadPersistent(unsigned int threadId, std::vector<unsigned long long> & requestQueues, std::vector<unsigned long long> & result) {
    // Create a random number generator
    std::random_device rd;  // Seed
    std::default_random_engine gen(rd()); // Mersenne Twister engine

    // Define the range
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    std::uniform_real_distribution<double> dis1(0.0, 1.0);

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
        //unsigned long long  totalNumInCircles = numberOfCircles(N, gen, dis);
        //auto  totalNumInCircles = numberOfCirclesStratified(std::sqrt(N), gen, dis1);
        auto unstratified = estimateIntegralSum(N, gen, dis1);

        // save result
        result[threadId] = unstratified;        
        //result[threadId] = totalNumInCircles.stratified;

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
    unsigned long long runs = 0, totalSatisfied = 0;
    while (true) {
        localRequestQueues[thread]++;
        runs++;
        if (localRequestQueues[thread] >= RUN_PER_BATCH) {
            // signal thread to take requests, communicate by one variable.
            {
                std::lock_guard<std::mutex> lock1(*mtx[thread]);
                requestQueues[thread] = localRequestQueues[thread];
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
                unsigned long long currSatisfied = 0, currTotalSatisfied = 0;
                for (unsigned int threadId = 0; threadId < processor_count; ++threadId) {
                    currSatisfied += result[threadId];
                    currTotalSatisfied += localRequestQueues[threadId];
                    result[threadId] = 0;
                    localRequestQueues[threadId] = 0;
                }
                totalSatisfied += currSatisfied;
                // 
                //std::cout << "\rEstimate of Pi = " << estimatePi(currNumInCircles, currTotalNumRuns) << " from current " << currTotalNumRuns << " runs." << std::endl;
                //std::cout << "\rEstimate of Pi = " << estimatePi(totalSatisfied, runs) << " from " << runs << " runs." << std::endl;
                std::cout << "\rEstimate of integral over 0 to 2 = " << estimateIntegral(totalSatisfied, runs) << " from " << runs << " runs." << std::endl;

            }
        }
        // increment 
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
