#pragma once
#include <random>

inline double estimatePi(unsigned long long numInCircle, unsigned long long N) {
    return 4.0 * (double(numInCircle) / double(N));
}

unsigned long long numberOfCircles(unsigned long long N, std::default_random_engine&gen, std::uniform_real_distribution<double> &dis) {
    unsigned long long numInCircle = 0;
    for (int i = 0; i < N; ++i) {
        double rx = dis(gen);
        double ry = dis(gen);

        if (rx * rx + ry * ry <= 1.0) {
            numInCircle++;
        }
    }
    return numInCircle;
}

struct CirclesResult {
    unsigned long long unstratified = 0;
    unsigned long long stratified = 0;
};

CirclesResult numberOfCirclesStratified(unsigned long long sqrtN, std::default_random_engine& gen, std::uniform_real_distribution<double>& dis) {
    unsigned long long unstratified = 0, stratified = 0;

    for (int i = 0; i < sqrtN; ++i) {
        for (int j = 0; j < sqrtN; ++j) {
            double rnd1 = dis(gen);   // 0 - 1
            double rnd2 = dis(gen);   // 0 - 1

            double rx = 2 * rnd1 - 1;
            double ry = 2 * rnd2 - 1;
            if (rx * rx + ry * ry <= 1.0) {
                unstratified++;
            }

            rx = 2 * ((i + rnd1) / sqrtN) - 1;
            ry = 2 * ((j + rnd2) / sqrtN) - 1;
            if (rx * rx + ry * ry <= 1.0) {
                stratified++;
            }
        }
    }

    // return unstratefied and stratified.
    return {
        unstratified, stratified
    };
}