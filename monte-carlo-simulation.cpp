/*
MIT License

Copyright (c) 2024 Zihao Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <iostream>
#include <format>
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
    std::mt19937 gen(rd()); // Mersenne Twister engine

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
    std::mt19937 gen(rd()); // Mersenne Twister engine

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
        auto integralSum = estimateIntegralSum(N, gen, dis1);

        // save result
        result[threadId] = integralSum;
        //result[threadId] = totalNumInCircles.stratified;

        // notify main.
        {
            std::lock_guard<std::mutex> lock2(mtx2);
            numberOfResultsReady++;
            // signal main.
            cv2.notify_one();
        }
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

    auto start = std::chrono::high_resolution_clock::now();

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
                    // waiting: implements while (!pred()) wait(lock);
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
                //std::cout << "\rEstimate of Pi = " << std::format("{}", estimatePi(totalSatisfied, runs)) << " from " << runs << " runs." << std::endl;
                //std::cout << "\rEstimate of integral = " << estimateIntegral(currSatisfied, currTotalSatisfied) << " from current " << currTotalSatisfied << " runs." << std::endl;
                std::cout << "\rEstimate of integral = " << std::format("{}", estimateIntegral(totalSatisfied, runs)) << " from " << runs << " runs." << std::endl;

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                double seconds = duration.count();
                std::cout << "Duration: " << seconds << " seconds " << " rate (runs/second): " << (currTotalSatisfied / seconds) << std::endl;
                start = end;
            }
        }
        // increment 
        thread = (thread + 1) % processor_count;
    }
}

int main()
{
    //std::cout << std::fixed << std::setprecision(10);

    const auto processor_count = std::thread::hardware_concurrency();
    //unsigned long long OriginalN = 500000000;
    //std::cout << "estimated pi, multi-threaded, N = " << OriginalN << " pi = " << estimatePiMultiThreaded(OriginalN, processor_count) << std::endl;
    
    estimatePiContinuously(processor_count);
}
