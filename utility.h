#pragma once
#include <random>
#include <algorithm>
#include <cstdlib>

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

struct Result {
    unsigned long long unstratified = 0;
    unsigned long long stratified = 0;
};

Result numberOfCirclesStratified(unsigned long long sqrtN, std::default_random_engine& gen, std::uniform_real_distribution<double>& dis) {
    unsigned long long unstratified = 0, stratified = 0;

    for (int i = 0; i < sqrtN; ++i) {
        for (int j = 0; j < sqrtN; ++j) {
            double rnd1 = dis(gen);   // 0 ~ 1
            double rnd2 = dis(gen);   // 0 ~ 1

            double rx = 2 * rnd1 - 1;       // -1 ~ 1
            double ry = 2 * rnd2 - 1;       // -1 ~ 1
            if (rx * rx + ry * ry <= 1.0) {
                unstratified++;
            }

            // stratified sampling.
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

inline double estimateIntegral(double sum, unsigned long long N) {
    return (sum / double(N));
}

// integral over 0 to 2 is: 2.6666666...
inline double integrandXSquared(double x) {
    return x * x;
}

inline double integrandXPow(double x, double exponent) {
    return std::pow(x, exponent);
}

// integral over 0 to 2 is: 0.903931238481499
inline double integrandSinXPow5(double x) {
    return std::pow(std::sin(x), 5.0);
}

inline double integrandLogSin(double x) {
    return std::log(std::sin(x));
}

double estimateIntegralSum(unsigned long long N, std::default_random_engine& gen, std::uniform_real_distribution<double>& dis) {
    double unstratified = 0;

    constexpr double a = 0, b = 2;

    double deltaX = (b - a) / double(N);
    double x_stratified = 0;

    for (int i = 0; i < N; ++i) {
        double rnd = dis(gen);   // 0 ~ 1
        double x = rnd * (b - a) + a;  // a ~ b
        //unstratified += integrandXSquared(x);
        //unstratified += integrandSinXPow5(x);
        //unstratified += integrandLogSin(x);
        unstratified += integrandXPow(x, 2.5);
    }
    unstratified *= (b - a);

    // return unstratefied and stratified.
    return unstratified;
}
