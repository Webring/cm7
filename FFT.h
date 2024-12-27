#ifndef FFT_H
#define FFT_H
#include "utils.h"

class FFT {
public:
    static ComplexVector computeForward(const ComplexVector &input) {
        size_t N = input.size();
        if (N <= 1) return input;

        if ((N & (N - 1)) != 0) {
            throw invalid_argument("Size of input must be a power of 2.");
        }

        ComplexVector even;
        ComplexVector odd;
        for (size_t i = 0; i < N; ++i) {
            (i % 2 == 0 ? even : odd).push_back(input[i]);
        }

        auto evenFFT = computeForward(even);
        auto oddFFT = computeForward(odd);

        ComplexVector output(N);
        for (size_t k = 0; k < N / 2; ++k) {
            Complex t = exp(Complex(0, -2.0 * M_PI * k / N)) * oddFFT[k];
            output[k] = evenFFT[k] + t;
            output[k + N / 2] = evenFFT[k] - t;
        }
        return output;
    }

    static ComplexVector computeInverse(const ComplexVector &input) {
        size_t N = input.size();
        if (N <= 1) return input;

        ComplexVector conjugatedInput(N);
        for (size_t i = 0; i < N; ++i) {
            conjugatedInput[i] = conj(input[i]);
        }

        auto output = computeForward(conjugatedInput);

        for (size_t i = 0; i < N; ++i) {
            output[i] = conj(output[i]) / static_cast<double>(N);
        }
        return output;
    }
};

ComplexVector computeConvolution(ComplexVector &a, ComplexVector &b) {
    int N = a.size();
    ComplexVector aTransformed = FFT::computeForward(a);
    ComplexVector bTransformed = FFT::computeForward(b);
    ComplexVector result(N);

    for (int i = 0; i < N; i++) {
        result[i] = aTransformed[i] * bTransformed[i];
    }

    return FFT::computeInverse(result);
}

#endif //FFT_H
