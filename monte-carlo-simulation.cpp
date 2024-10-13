// monte-carlo-simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iomanip>

#include <cstdlib>
#include <random>
#include <thread>

int main()
{
    // Create a random number generator
    std::random_device rd;  // Seed
    std::mt19937 gen(rd()); // Mersenne Twister engine

    // Define the range
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    std::cout << std::fixed << std::setprecision(12);
    const int N = 1000000;
    int numInCircle = 0;
    for (int i = 0; i < N; ++i) {
        float rx = dis(gen);
        float ry = dis(gen);

        if (rx * rx + ry * ry <= 1) {
            numInCircle++;
        }
    }
    double pi = 4.0 * numInCircle / double(N);
    std::cout << "PI = " << pi << std::endl;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
