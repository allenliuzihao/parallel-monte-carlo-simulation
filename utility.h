#pragma once
#include <random>
#include <algorithm>
#include <cstdlib>
#include <numbers>

inline double estimatePi(unsigned long long numInCircle, unsigned long long N) {
    return 4.0 * (double(numInCircle) / double(N));
}

unsigned long long numberOfCircles(unsigned long long N, std::mt19937&gen, std::uniform_real_distribution<double> &dis) {
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

Result numberOfCirclesStratified(unsigned long long sqrtN, std::mt19937& gen, std::uniform_real_distribution<double>& dis) {
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

inline double integrandExponent(double x) {
    return std::exp(x);
}

inline double integrandLogSin(double x) {
    return std::log(std::sin(x));
}

/* uniform sampling strategy */
inline double uniformSample(double a, double b, double rnd) {
    return rnd * (b - a) + a;  // a ~ b
}

inline double uniformPdf(double a, double b) {
    return 1.0 / (b - a);
}

/* linear function sampling strategy */
inline double linearSample(double rnd) {
    return 2.0 * std::sqrt(rnd);  // a ~ b
}

inline double linearPdf(double x) {
    return x / 2.0;
}

/* quadratic sampling strategy */
inline double quadraticSample(double rnd) {
    return 8.0 * std::pow(rnd, 1.0 / 3.0);
}

inline double quadraticPdf(double x) {
    return (3.0 / 8.0) * x * x;
}

struct vec3 {
    double x{}, y{}, z{};
    
    double length() {
        return std::sqrt(x * x + y * y + z * z);
    }

    vec3 normalize() {
        double l = length();
        double inv_l = 1.0 / l;
        if (l == 0) {
            return vec3(0, 0, 0);
        }
        return {
            x * inv_l, y * inv_l, z * inv_l
        };
    }
};

inline vec3 sampleDirectionOnUnitSphereRejection(std::mt19937& gen, std::uniform_real_distribution<double>& dis) {
    while (true) {
        double rnd1 = 2.0 * dis(gen) - 1.0;   // -1 ~ 1
        double rnd2 = 2.0 * dis(gen) - 1.0;   // -1 ~ 1
        double rnd3 = 2.0 * dis(gen) - 1.0;   // -1 ~ 1

        vec3 dir = vec3(rnd1, rnd2, rnd3);
        double l = dir.length();
        if (l <= 1.0 && l > 0.0) {
            return dir.normalize();
        }
    }
    return {};
}

inline vec3 sampleDirectionOnUnitSphere(std::mt19937& gen, std::uniform_real_distribution<double>& dis) {
    double rnd1 = dis(gen);   // 0 ~ 1
    double rnd2 = dis(gen);   // 0 ~ 1

    double phi = 2 * std::numbers::pi * rnd1;
    double sqrtRnd2 = 2 * std::sqrt(rnd2 * (1 - rnd2));
    double x = std::cos(phi) * sqrtRnd2;
    double y = std::sin(phi) * sqrtRnd2;
    double z = 1.0 - 2.0 * rnd2;

    return vec3(x, y, z);
}

inline vec3 sampleDirectionOnHemisphere(std::mt19937& gen, std::uniform_real_distribution<double>& dis) {
    double rnd1 = dis(gen);   // 0 ~ 1
    double rnd2 = dis(gen);   // 0 ~ 1

    double phi = 2 * std::numbers::pi * rnd1;
    double sqrtRnd2 = 2 * std::sqrt(rnd2 * (1 - rnd2));
    double x = std::cos(phi) * sqrtRnd2;
    double y = std::sin(phi) * sqrtRnd2;
    double z = 1.0 - rnd2;

    return vec3(x, y, z);
}

inline vec3 sampleDirectionCosineWeightedOnHemisphere(std::mt19937& gen, std::uniform_real_distribution<double>& dis) {
    double rnd1 = dis(gen);   // 0 ~ 1
    double rnd2 = dis(gen);   // 0 ~ 1

    double phi = 2 * std::numbers::pi * rnd1;
    double sqrtRnd2 = std::sqrt(rnd2);
    double x = std::cos(phi) * sqrtRnd2;
    double y = std::sin(phi) * sqrtRnd2;
    double z = std::sqrt(1.0 - rnd2);

    return vec3(x, y, z);
}

/* uniform direction on unit sphere */
// pi
inline double cosineIntegrand(const vec3& d) {
    return d.z;
}

// 4 * pi / 3
inline double cosineSquaredIntegrand(const vec3 &d) {
    return d.z * d.z;
}

// pi / 2
inline double cosinePow3Integrand(const vec3& d) {
    return d.z * d.z * d.z;
}

inline double unitSphereDirPdf(const vec3 &d) {
    return 1.0 / (4.0 * std::numbers::pi);
}

inline double hemisphereDirPdf(const vec3& d) {
    return 1.0 / (2.0 * std::numbers::pi);
}

inline double cosineWeightedHemisphereDirPdf(const vec3& d) {
    return d.z / std::numbers::pi;
}

double estimateIntegralSum(unsigned long long N, std::mt19937& gen, std::uniform_real_distribution<double>& dis) {
    double unstratified = 0;
    constexpr double a = 0, b = 2;
    for (int i = 0; i < N; ++i) {
        /*
        double rnd = dis(gen);   // 0 ~ 1
        if (rnd == 0) { // x = 0, pdf = 0, x ^ 2 = 0, no contribution.
            continue;
        }
        */
        //double x = uniformSample(a, b, rnd);
        //double pdf = uniformPdf(a, b);
        //double x = linearSample(rnd);
        //double pdf = linearPdf(x);
        //double x = quadraticSample(rnd);
        //double pdf = quadraticPdf(x);
        //vec3 dir = sampleDirectionOnUnitSphere(gen, dis);
        //double pdf = unitSphereDirPdf(dir);
        //vec3 dir = sampleDirectionOnHemisphere(gen, dis);
        //double pdf = hemisphereDirPdf(dir);
        vec3 dir = sampleDirectionCosineWeightedOnHemisphere(gen, dis);
        double pdf = cosineWeightedHemisphereDirPdf(dir);

        //unstratified += integrandXSquared(x) / pdf;
        //unstratified += integrandSinXPow5(x) / pdf;
        //unstratified += integrandLogSin(x);
        //unstratified += integrandXPow(x, 2.5);
        //unstratified += integrandExponent(x) / pdf;
        //unstratified += cosineSquaredIntegrand(dir) / pdf;
        //unstratified += cosinePow3Integrand(dir) / pdf;
        unstratified += cosineIntegrand(dir) / pdf;
    }
    //assert(std::abs(unstratified / N  -  8.0 / 3.0) < 1e-8);

    // return unstratefied and stratified.
    return unstratified;
}
