#include <cmath>
#include "iostream"

#include "wavelet.h"

using namespace std;

string getBasisName(WaveletType type) {
    switch (type) {
        case WaveletType::daubechies6:
            return "Добеши (D6)";
        case WaveletType::haar:
            return "Хаара";
        case WaveletType::shannon:
            return "Шенона";
    }
}

int main() {
    setlocale(LC_ALL, "Russian");

    int N = 512;
    int stages = 4;
    double alpha = 1.68;
    double beta = 2.6;

    vector<WaveletType> basises{WaveletType::daubechies6, WaveletType::shannon};


    ComplexVector inputSignal(N);
    for (int n = 0; n < N; n++) {
        if (n >= 128 && n <= 255) {
            inputSignal[n] = Complex(sin(fabs(pow(n - 128, alpha)) / 128.0), 0.0);
        } else if (n >= 384 && n <= 447) {
            inputSignal[n] = Complex(sin(fabs(pow(n - 384, beta)) / 128.0), 0.0);
        } else {
            inputSignal[n] = Complex(0, 0);
        }
    }

    ofstream file_high("high.csv");
    ofstream file_new_z("output.csv");

    for (int stage = 1; stage <= stages; stage++) {
        for (auto basis: basises) {
            string header = to_string(stage) + "-го этапа " + getBasisName(basis);

            Wavelet wavelet(N, basis, inputSignal);
            wavelet.runWaveletTransform(stage);
            wavelet.saveResults(header);

            wavelet.addToFile("Высокие частоты " + header, wavelet.getHighFrequencyCoeffs(), file_high);
            wavelet.addToFile("Востановленный сигнал " + header, wavelet.getReconstructedSignal(), file_new_z);

            cout << header << " ошибка: " << wavelet.calculateError() << endl;
        }
    }
    file_high.close();
    file_new_z.close();


    return 0;
}
