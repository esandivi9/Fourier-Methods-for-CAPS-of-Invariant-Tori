#include <gtest/gtest.h>

#include <fstream>

RealInterval log(RealInterval x) {
    RealInterval log;
    mpfi_log(log.interval, x.interval);
    return log;
}

void plot(vector<RealInterval> data1, vector<RealInterval> data2, vector<RealInterval> data3) {
    if (data1.size() != data2.size() || data1.size() != data3.size()) {
        cerr << "Error: Vectors must have the same size for comparison." << endl;
        return;
    }

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");

    if (gnuplotPipe) {
        // Set the font in the gnuplot script
        fprintf(gnuplotPipe, "set term qt font \"Arial, 10\"\n");
        // fprintf(gnuplotPipe, "set logscale y\n");
        fprintf(gnuplotPipe, "set xrange [-50:50]\n");
        // fprintf(gnuplotPipe, "set yrange [0:0.3]\n");
        fprintf(gnuplotPipe, "plot '-' with lines title 'Data2', '-' with lines title 'Data3'\n");

        // Plot the real parts of both data2 and data3 against data1
        for (size_t i = 0; i < data1.size(); ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", mpfi_get_d(data1[i].interval), mpfi_get_d(data2[i].interval));
        }
        fprintf(gnuplotPipe, "e\n");

        // Separate data sets with a blank line
        fprintf(gnuplotPipe, "\n");

        for (size_t i = 0; i < data1.size(); ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", mpfi_get_d(data1[i].interval), mpfi_get_d(data3[i].interval));
        }
        fprintf(gnuplotPipe, "e\n");

        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        cerr << "Gnuplot not found. Please install Gnuplot on your system." << endl;
    }
}

void plot(vector<double> x, vector<double> y) {
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");

    if (gnuplotPipe) {
        // Set the font in the gnuplot script
        fprintf(gnuplotPipe, "set term qt font \"Arial, 10\"\n");
        fprintf(gnuplotPipe, "plot '-' with dots\n");

        for (size_t i = 0; i < x.size(); ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", x[i], y[i]);
        }

        fprintf(gnuplotPipe, "e\n");
        fflush(gnuplotPipe);
        pclose(gnuplotPipe);
    } else {
        cerr << "Gnuplot not found. Please install Gnuplot on your system." << endl;
    }
}

// Test norms with read torus
TEST(ValidationTests, FourierDecayTest) {
    // Specify the subfolder and file name
    string subfolder = "Inputs/";
    string fileName = subfolder + "K1.231000.txt";

    // Read N from file
    ifstream inputFile(fileName);
    if (!inputFile) {
        cerr << "File Error" << endl;
        // or handle the error accordingly
    }
    int N, lines;
    string line;
    for (lines = 0; getline(inputFile, line); lines++)
        ;
    N = lines;

    // Number of components in our torus
    int m = 2;

    // Allocate the matrices
    ComplexMatrix diff, P1, fP1, P2, fP2, frotP2, Lam, fLam, Id, fId;
    allocAllMatrices(m, m, N, diff, P1, fP1, P2, fP2, frotP2, Lam, fLam, Id, fId);

    // Allocate the vectors
    ComplexVector K, FK, fK, fFK, frotK;
    allocAllVectors(m, N, K, FK, fK, fFK, frotK);

    // Read the data from file
    double lamUnst, lamSt, e, kap;
    vector<vector<double>> Kc(m, vector<double>(N));
    vector<vector<vector<double>>> P0(m, vector<vector<double>>(m, vector<double>(N)));
    readData(fileName, m, N, e, kap, Kc, P0, K, P1, Lam, Id, lamUnst, lamSt);
    RealInterval eps = e, kappa = kap;

    vectorFFT(fK, K);
    fK[0][N / 2] = 0;
    fK[1][N / 2] = 0;

    vector<RealInterval> indices;
    vector<RealInterval> coeffMods;
    vector<RealInterval> fit;
    RealInterval c = "8.4e-8";
    RealInterval rho = "6.3e-5";

    for (int k = 0; k < N; k++) {
        if (mod(fK[0][k]) < 1.e-12)
            continue;
        coeffMods.push_back(log(mod(fK[0][k])));
        fit.push_back(c * exp(DPI * (-1) * abs(k) * rho));
        // print(coeffMods[k]);
        /*if (k > N / 2)
            indices.push_back(k - N);
        else*/
        indices.push_back(k);
    }

    // for (int k = 0; k < coeffMods.size(); k++)
    // print(coeffMods[k]);

    // plot(indices, coeffMods, fit);
    //performLinearRegressionPlot(indices, coeffMods);
}

// Test norms with read torus
TEST(ValidationTests, Torus) {
    ifstream inputFile("Inputs/K0.100000.txt");
    int N = 1024;

    ComplexVector torus;
    alloc(torus, 2, N);

    vector<double> secondColumn;
    vector<double> thirdColumn;
    double value2, value3, value4, value5;

    string line;
    while (getline(inputFile, line)) {
        istringstream iss(line);
        // Variables to store values from the second and third columns
        double discardValue;  // We'll discard the values from other columns
        if (iss >> discardValue >> discardValue >> discardValue >> discardValue >> value2 >> value3 >> value4 >> value5) {
            secondColumn.push_back(value2);  // Store the value from the second column
            thirdColumn.push_back(value3);   // Store the value from the third column
        }
    }

    inputFile.close();

    for (int i = 0; i < N; i++) {
        torus[0][i] = ComplexInterval(secondColumn[i], 0);
        torus[1][i] = ComplexInterval(thirdColumn[i], 0);
    }

    // F(torus, torus, 1, 0.1);
    ComplexVector coeffs = torus;
    vectorFFT(coeffs, torus);
    double rho = 1.e-5;
    RealInterval rhoInterval = rho;
    auto fNorm = fourierNorm(coeffs, rhoInterval);
    auto rNorm = rhoNorm(torus, rho);
    auto sNorm = supNorm(torus);

    mpfr_t midPoint1, midPoint2;
    mpfr_init2(midPoint1, prec);
    mpfr_init2(midPoint2, prec);
    mpfi_mid(midPoint1, rNorm.interval);
    mpfi_mid(midPoint2, fNorm.interval);
    EXPECT_TRUE(mpfr_cmp(midPoint1, midPoint2) < 0);
}

// Test norms with flat torus
TEST(ValidationTests, FlatTorus) {
    int N = 1024;

    ComplexVector torus;
    alloc(torus, 2, N);

    for (int k = 0; k < N; k++) {
        torus[0][k] = ComplexInterval((cos(2 * M_PI * k) / N) + 1. / 2, 0);
        torus[1][k] = ComplexInterval((sin(4 * M_PI * k) / N), 0);
    }

    // F(torus, torus, 1, 0.1);
    ComplexVector coeffs = torus;
    vectorFFT(coeffs, torus);
    double rho = 1.e-5;
    RealInterval rhoHat = 1.e-2;
    RealInterval rhoInterval = rho;
    auto fNorm = fourierNorm(coeffs, rhoInterval);
    auto rNorm = rhoNorm(torus, rho);
    auto sNorm = supNorm(torus);

    EXPECT_TRUE(rNorm <= fNorm);
}

// Test norms with toy torus
TEST(ValidationTests, TestTorus) {
    int N = 4096;

    ComplexVector torus;
    alloc(torus, 1, N);

    for (int k = 0; k < N; k++) {
        torus[0][k] = ComplexInterval(cos((2 * M_PI * k) / N) + 1. / 2, 0);
    }

    ComplexVector coeffs = torus;
    vectorFFT(coeffs, torus);
    double rho = 1.e-5;
    RealInterval rhoHat = 1.e-2;
    RealInterval rhoInterval = rho;
    auto fNorm = fourierNorm(coeffs, rhoInterval);
    auto rNorm = rhoNorm(torus, rho);
    auto sNorm = supNorm(torus);
    mpfr_t midPoint1, midPoint2;
    mpfr_init2(midPoint1, prec);
    mpfr_init2(midPoint2, prec);
    mpfi_mid(midPoint1, rNorm.interval);
    mpfi_mid(midPoint2, fNorm.interval);
    EXPECT_TRUE(mpfr_cmp(midPoint1, midPoint2) < 0);
}

// Test the Fourier norm increases as the complex width increases
TEST(ValidationTests, fourierNormDecrease) {
    int N = 4096;

    ComplexVector torus;
    alloc(torus, 1, N);
    // Generate a torus
    for (int k = 0; k < N; k++) {
        torus[0][k] = ComplexInterval(cos((2 * M_PI * k) / N) + 1. / 2, 0);
    }

    ComplexVector coeffs = torus;
    vectorFFT(coeffs, torus);

    // Check the Fourier norm increases as rho increases
    RealInterval fNorm, normRho, fNormaux = 0;
    for (double rhoTest = 0; rhoTest <= 1.01e-3; rhoTest += 1.e-4) {
        normRho = rhoTest;
        fNorm = fourierNorm(coeffs, normRho);
        EXPECT_TRUE(fNormaux <= fNorm);
        fNormaux = fNorm;
    }
}

// Test that Fourier coefficients get less precise the more times the transform is performed and the bigger the N
TEST(ValidationTests, fourierTransformError) {
    for (int N = 8; N <= 8192; N *= 2) {
        vector<ComplexInterval> torus(N);
        // Generate a torus
        for (int k = 0; k < N; k++) {
            torus[k] = ComplexInterval(cos((2 * M_PI * k) / N) + 1. / 2, 0);
        }

        vector<ComplexInterval> coeffs = torus, testCoeffs(N);
        fft(coeffs, torus);

        ifft(torus, coeffs);
        fft(testCoeffs, torus);
        // The coefficients are less refined the more times the fft is performed
        for (int k = 0; k < N; k++)
            EXPECT_TRUE(mpfi_cmp(coeffs[k].real, testCoeffs[k].real) <= 0 && mpfi_cmp(coeffs[k].imag, testCoeffs[k].imag) <= 0);

        // The bigger the order of N, the bigger the difference between coefficients
        ComplexInterval difference = 0;
        if (N == 8) {
            difference = testCoeffs[N / 4] - coeffs[N / 4];
        } else if (N == 8192) {
            ComplexInterval bigDifference = testCoeffs[N / 4] - coeffs[N / 4];
            EXPECT_TRUE(mpfi_cmp(difference.real, bigDifference.real) <= 0 && mpfi_cmp(difference.imag, bigDifference.imag) <= 0);
        }
    }
}