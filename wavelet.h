#ifndef WAVELET_H
#define WAVELET_H
#include "FFT.h"


enum class WaveletType {
    haar,
    daubechies6,
    shannon
};

class Wavelet {
private:
    ComplexVector filterU, filterV, inputSignal, reconstructedSignal, highFrequencyCoeffs, lowFrequencyCoeffs;
    ComplexMatrix filtersHigh, filtersLow, psiBasis, phiBasis;

    void computeFilters(int stages) {
        int N = filterU.size();
        ComplexMatrix filtersU(stages), filtersV(stages);
        filtersU[0] = filterU;
        filtersV[0] = filterV;

        for (int i = 1; i < stages; i++) {
            int elementCount = N / static_cast<int>(pow(2, i));
            filtersU[i].resize(elementCount);
            filtersV[i].resize(elementCount);

            for (int n = 0; n < elementCount; n++) {
                int max = static_cast<int>(pow(2, i));
                for (int k = 0; k < max; k++) {
                    filtersU[i][n] += filtersU[0][n + k * N / max];
                    filtersV[i][n] += filtersV[0][n + k * N / max];
                }
            }
        }

        filtersHigh.resize(stages);
        filtersLow.resize(stages);
        filtersHigh[0] = filtersV[0];
        filtersLow[0] = filtersU[0];

        for (int j = 1; j < stages; j++) {
            ComplexVector dilatedU, dilatedV;
            dilateFilter(j, filtersU[j], dilatedU);
            dilateFilter(j, filtersV[j], dilatedV);

            filtersHigh[j] = computeConvolution(filtersLow[j - 1], dilatedV);
            filtersLow[j] = computeConvolution(filtersLow[j - 1], dilatedU);
        }
    }

    void computeBasis(int stage) {
        int N = filterU.size();
        int elementCount = N / static_cast<int>(pow(2, stage));

        if (filtersLow.size() < static_cast<size_t>(stage)) {
            computeFilters(stage + 1);
        }

        psiBasis.resize(elementCount);
        phiBasis.resize(elementCount);

        for (int i = 0; i < elementCount; i++) {
            int index = static_cast<int>(pow(2, stage)) * i;
            shiftSignal(index, filtersHigh[stage - 1], psiBasis[i]);
            shiftSignal(index, filtersLow[stage - 1], phiBasis[i]);
        }
    }

public:
    Wavelet(int N, WaveletType type, const ComplexVector &inputSignal) : inputSignal(inputSignal) {
        filterU.resize(N);
        filterV.resize(N);

        if (type == WaveletType::shannon) {
            // Shannon Wavelet Filters
            filterU[0] = filterV[0] = 1.0 / sqrt(2.0);
            for (int i = 1; i < N; i++) {
                Complex value = Complex(sqrt(2.0) / N * cos(M_PI * i / N) * sin(M_PI * i / 2.0) / sin(M_PI * i / N),
                                        -sqrt(2.0) / N * sin(M_PI * i / N) * sin(M_PI * i / 2.0) / sin(M_PI * i / N));
                filterU[i] = value;
                filterV[i] = pow(-1, i) * value;
            }
        } else if (type == WaveletType::daubechies6) {
            // Daubechies Wavelet Filters
            double a = 1.0 - sqrt(10.0);
            double b = 1.0 + sqrt(10.0);
            double c = sqrt(5.0 + 2 * sqrt(10.0));
            double mult = sqrt(2.0) / 32.0;

            filterU[0] = (b + c) * mult;
            filterU[1] = (2 * a + 3 * b + 3 * c) * mult;
            filterU[2] = (6 * a + 4 * b + 2 * c) * mult;
            filterU[3] = (6 * a + 4 * b - 2 * c) * mult;
            filterU[4] = (2 * a + 3 * b - 3 * c) * mult;
            filterU[5] = (b - c) * mult;

            filterV[0] = -filterU[1];
            filterV[1] = filterU[0];
            filterV[N - 1] = filterU[2];
            filterV[N - 2] = -filterU[3];
            filterV[N - 3] = filterU[4];
            filterV[N - 4] = -filterU[5];
        }
    }

    void runWaveletTransform(int stage) {
        computeBasis(stage);
        int elementCount = psiBasis.size();

        for (int i = 0; i < elementCount; i++) {
            highFrequencyCoeffs.push_back(computeScalarProduct(inputSignal, psiBasis[i]));
            lowFrequencyCoeffs.push_back(computeScalarProduct(inputSignal, phiBasis[i]));
        }

        int N = filterU.size();
        for (int i = 0; i < N; i++) {
            Complex highFreqPart(0, 0), lowFreqPart(0, 0);

            for (int j = 0; j < elementCount; j++) {
                lowFreqPart += lowFrequencyCoeffs[j] * phiBasis[j][i];
                highFreqPart += highFrequencyCoeffs[j] * psiBasis[j][i];
            }

            reconstructedSignal.push_back(lowFreqPart + highFreqPart);
        }
    }

    double calculateError() {
        double error = 0.0;
        for (size_t i = 0; i < inputSignal.size(); ++i) {
            error += pow(inputSignal[i] - reconstructedSignal[i], 2).real();
        }
        return sqrt(error / inputSignal.size());
    }

    void addToFile(string header, ComplexVector data, ofstream &file) {
        file << header << ";";
        for (const auto &value: data) {
            file << value.real() << ";";
        }
        file << endl;
    }

    void saveResults(string &filebegin) {
        auto saveToFile = [](const string &filename, const ComplexVector &data) {
            ofstream file(filename);
            for (const auto &value: data) {
                file << value.real() << "\n";
            }
            file.close();
        };

        saveToFile(filebegin + "_high_frequency.txt", highFrequencyCoeffs);
        saveToFile(filebegin + "_low_frequency.txt", lowFrequencyCoeffs);
        saveToFile(filebegin + "_original_signal.txt", inputSignal);
        saveToFile(filebegin + "_reconstructed_signal.txt", reconstructedSignal);
    }

    ComplexVector getHighFrequencyCoeffs() {
        return highFrequencyCoeffs;
    }

    ComplexVector getLowFrequencyCoeffs() {
        return lowFrequencyCoeffs;
    }
    ComplexVector getInputSignal() {
        return inputSignal;
    }
    ComplexVector getReconstructedSignal() {
        return reconstructedSignal;
    }

private:
    void shiftSignal(int shift, const ComplexVector &data, ComplexVector &result) {
        int N = data.size();
        result.resize(N);

        for (int i = 0; i < N; i++) {
            result[i] = data[(i - shift + N) % N];
        }
    }

    void dilateFilter(int stage, const ComplexVector &data, ComplexVector &result) {
        int factor = static_cast<int>(pow(2, stage));
        int N = data.size() * factor;
        result.resize(N);

        for (int i = 0; i < N; i++) {
            result[i] = (i % factor == 0) ? data[i / factor] : Complex(0, 0);
        }
    }

    Complex computeScalarProduct(const ComplexVector &a, const ComplexVector &b) {
        Complex result(0, 0);
        for (size_t i = 0; i < a.size(); i++) {
            result += a[i] * conj(b[i]);
        }
        return result;
    }
};
#endif //WAVELET_H
