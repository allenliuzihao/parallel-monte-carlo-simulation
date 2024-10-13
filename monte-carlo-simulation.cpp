// monte-carlo-simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <random>
#include <thread>

void numInCircleSimple(unsigned long long N, int threadId, std::vector<unsigned long long> & numCircles) {
    // Create a random number generator
    std::random_device rd;  // Seed
    std::mt19937 gen(rd()); // Mersenne Twister engine

    // Define the range
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    
    unsigned long long numInCircle = 0;
    for (int i = 0; i < N; ++i) {
        double rx = dis(gen);
        double ry = dis(gen);

        if (rx * rx + ry * ry <= 1.0) {
            numInCircle++;
        }
    }
    numCircles[threadId] = numInCircle;
}

double estimatePi(unsigned long long numInCircle, unsigned long long N) {
    return 4.0 * (double(numInCircle) / double(N)); 
}

int main()
{
    std::cout << std::fixed << std::setprecision(12);

    const auto processor_count = std::thread::hardware_concurrency();
    std::cout << "processor count: " << processor_count << std::endl;

    // estimate pi in simple form.
    unsigned long long N = 100000000;
    unsigned long long PerThreadN = N / processor_count;
    std::vector<long long int> numInCircles(processor_count, 0);

    // preallocate threads.
    std::vector<std::thread> threads(processor_count);

    // Launch a group of threads
    for (int i = 0; i < processor_count; ++i) {
        unsigned long long currN = std::min(N, PerThreadN);
        //  numInCircleSimple(unsigned long long N, int threadId, std::vector<unsigned long long> & numCircles)
        threads[i] = std::thread(numInCircleSimple, currN, i, numInCircles);
        N -= currN;
    }

    // Join the threads with the main thread
    unsigned long long totalNumInCircles = 0;
    for (int i = 0; i < processor_count; ++i) {
        threads[i].join();
        totalNumInCircles += numInCircles[i];
    }

    double pi = estimatePi(totalNumInCircles, N);
    std::cout << "PI = " << pi << std::endl;
}
