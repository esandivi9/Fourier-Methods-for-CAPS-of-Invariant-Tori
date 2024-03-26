#include "IntervalClasses.h"

// Read the initial data from file and store in intervals
void readData(const string &fileName, unsigned m, unsigned N, double &eps, double &kappa,
              vector<vector<double>> &Kc,
              vector<vector<vector<double>>> &P0,
              ComplexVector &K,
              ComplexMatrix &P1,
              ComplexMatrix &Lam,
              ComplexMatrix &Id,
              double &lam0, double &lam1) {
    ifstream inputFile(fileName);
    if (!inputFile.is_open()) {
        cerr << "Error opening file: " << fileName << endl;
        // Handle the error appropriately
        return;
    }

    for (unsigned k = 0; k < N; ++k) {
        string line;
        getline(inputFile, line);

        istringstream iss(line);

        // Ignore the first two columns
        string ignoreValue;
        iss >> ignoreValue;
        iss >> ignoreValue;

        iss >> eps >> kappa >> Kc[0][k] >> Kc[1][k];

        // Read the desired columns
        for (unsigned i = 0; i < m; ++i) {
            for (unsigned j = 0; j < m; ++j) {
                iss >> P0[i][j][k];
            }
        }
        iss >> lam0 >> lam1;

        // Ignore the last two columns
        for (int i = 0; i < 2; i++) {
            iss >> ignoreValue;
        }

        // Additional error checking if needed
        if (iss.fail()) {
            cerr << "Error reading values from line: " << line << endl;
            // Handle the error appropriately
            return;
        }

        for (unsigned i = 0; i < m; ++i) {
            K[i][k] = ComplexInterval(Kc[i][k], 0);
            for (unsigned j = 0; j < m; ++j) {
                P1[i][j][k] = ComplexInterval(P0[i][j][k], 0);
                if (i == j) {
                    if (i == 0) {
                        Lam[i][j][k] = ComplexInterval(lam0, 0);
                    } else {
                        Lam[i][j][k] = ComplexInterval(lam1, 0);
                    }
                    Id[i][j][k] = ComplexInterval(1, 0);
                } else {
                    Id[i][j][k] = 0;
                    Lam[i][j][k] = 0;
                }
            }
        }
    }

    inputFile.close();
}

// Print Fourier coefficients of a ComplexVector
void printFourier(string fileName, ComplexVector &fK) {
    int N = fK[0].size();
    // Extracting the number from the fileName
    size_t pos = fileName.find_last_of("K") + 1;
    double number = stod(fileName.substr(pos));

    int precision = 4;
    if (fmod(number * 100, 10) < 1.e-14)
        precision = 1;

    // Creating the output file name
    stringstream outputFile;
    outputFile << "Outputs/Fourier" << fixed << setprecision(precision) << number << ".txt";

    // Opening the output file
    FILE *output = fopen(outputFile.str().c_str(), "w");

    for (int k = N / 2; k < N; k++) {
        fprintf(output, "%d\t% .16le\t% .16le\t% .16le\t% .16le\n", k - N, mpfi_get_d(fK[0][k].real), mpfi_get_d(fK[0][k].imag),
                mpfi_get_d(fK[1][k].real), mpfi_get_d(fK[1][k].imag));
    }
    for (int k = 0; k < N / 2; k++) {
        fprintf(output, "%d\t% .16le\t% .16le\t% .16le\t% .16le\n", k, mpfi_get_d(fK[0][k].real), mpfi_get_d(fK[0][k].imag),
                mpfi_get_d(fK[1][k].real), mpfi_get_d(fK[1][k].imag));
    }
}

// Print Fourier coefficients of a ComplexVector
void printFourier(string fileName, ComplexMatrix &M) {
    int N = M[0][0].size();
    // Extracting the number from the fileName
    size_t pos = fileName.find_last_of("K") + 1;
    double number = stod(fileName.substr(pos));

    int precision = 4;
    if (fmod(number * 100, 10) < 1.e-14)
        precision = 1;

    // Creating the output file name
    stringstream outputFileSt;
    outputFileSt << "Outputs/FourierPSt" << fixed << setprecision(precision) << number << ".txt";

    // Opening the output file for the stable part
    FILE *outputSt = fopen(outputFileSt.str().c_str(), "w");

    for (int k = N / 2; k < N; k++) {
        fprintf(outputSt, "%d\t% .16le\t% .16le\t% .16le\t% .16le\n", k - N, mpfi_get_d(M[0][0][k].real), mpfi_get_d(M[0][0][k].imag),
                mpfi_get_d(M[1][0][k].real), mpfi_get_d(M[1][0][k].imag));
    }
    for (int k = 0; k < N / 2; k++) {
        fprintf(outputSt, "%d\t% .16le\t% .16le\t% .16le\t% .16le\n", k, mpfi_get_d(M[0][0][k].real), mpfi_get_d(M[0][0][k].imag),
                mpfi_get_d(M[1][0][k].real), mpfi_get_d(M[1][0][k].imag));
    }

    // Creating the output file name
    stringstream outputFileUnst;
    outputFileUnst << "Outputs/FourierPUnst" << fixed << setprecision(precision) << number << ".txt";

    // Opening the output file
    FILE *outputUnst = fopen(outputFileUnst.str().c_str(), "w");

    for (int k = N / 2; k < N; k++) {
        fprintf(outputUnst, "%d\t% .16le\t% .16le\t% .16le\t% .16le\n", k - N, mpfi_get_d(M[0][1][k].real), mpfi_get_d(M[0][1][k].imag),
                mpfi_get_d(M[1][1][k].real), mpfi_get_d(M[1][1][k].imag));
    }
    for (int k = 0; k < N / 2; k++) {
        fprintf(outputUnst, "%d\t% .16le\t% .16le\t% .16le\t% .16le\n", k, mpfi_get_d(M[0][1][k].real), mpfi_get_d(M[0][1][k].imag),
                mpfi_get_d(M[1][1][k].real), mpfi_get_d(M[1][1][k].imag));
    }
}

// Increase grid size for functions by adding zeroes
int increaseGridSize(int factor, ComplexVector &fK, ComplexVector &K, ComplexVector &FK, ComplexVector &fFK,
                     ComplexVector &frotK, ComplexMatrix &fP1, ComplexMatrix &fP2, ComplexMatrix &fLam,
                     ComplexMatrix &fId, ComplexMatrix &diff, ComplexMatrix &P1, ComplexMatrix &P2, ComplexMatrix &frotP2,
                     ComplexMatrix &Lam, ComplexMatrix &Id) {
    int N = fK[0].size();
    int m = P1.size();
    int newSize = N * factor;
    for (int i = 0; i < m; i++) {
        fK[i].insert(fK[i].begin() + N / 2, newSize - N, ComplexInterval(0));
        K[i].resize(newSize);
        FK[i].resize(newSize);
        fFK[i].resize(newSize);
        frotK[i].resize(newSize);
        for (int j = 0; j < m; j++) {
            fP1[i][j].insert(fP1[i][j].begin() + N / 2, newSize - N, ComplexInterval(0));
            fP2[i][j].insert(fP2[i][j].begin() + N / 2, newSize - N, ComplexInterval(0));
            fLam[i][j].insert(fLam[i][j].begin() + N / 2, newSize - N, ComplexInterval(0));
            fId[i][j].insert(fId[i][j].begin() + N / 2, newSize - N, ComplexInterval(0));
            diff[i][j].resize(newSize);
            P1[i][j].resize(newSize);
            P2[i][j].resize(newSize);
            frotP2[i][j].resize(newSize);
            Lam[i][j].resize(newSize);
            Id[i][j].resize(newSize);
        }
    }

    return newSize;
}

// Reduce grid size for functions by removing noise in the middle
int reduceGridSize(double factor, ComplexVector &fK, ComplexVector &K, ComplexVector &FK, ComplexVector &fFK,
                   ComplexVector &frotK, ComplexMatrix &fP1, ComplexMatrix &fP2, ComplexMatrix &fLam,
                   ComplexMatrix &fId, ComplexMatrix &diff, ComplexMatrix &P1, ComplexMatrix &P2, ComplexMatrix &frotP2,
                   ComplexMatrix &Lam, ComplexMatrix &Id) {
    int N = fK[0].size();
    int m = P1.size();
    int originalSize = N;
    int numElementsToRemove = N * factor;
    int newSize = originalSize - numElementsToRemove;

    // Calculate the position to start erasing from
    size_t eraseStart = (originalSize - numElementsToRemove) / 2;

    // Erase the middle N/2 elements
    for (int i = 0; i < m; i++) {
        fK[i].erase(fK[i].begin() + eraseStart, fK[i].begin() + eraseStart + numElementsToRemove);
        K[i].resize(newSize);
        FK[i].resize(newSize);
        fFK[i].resize(newSize);
        frotK[i].resize(newSize);
        for (int j = 0; j < m; j++) {
            fP1[i][j].erase(fP1[i][j].begin() + eraseStart, fP1[i][j].begin() + eraseStart + numElementsToRemove);
            fP2[i][j].erase(fP2[i][j].begin() + eraseStart, fP2[i][j].begin() + eraseStart + numElementsToRemove);
            fLam[i][j].erase(fLam[i][j].begin() + eraseStart, fLam[i][j].begin() + eraseStart + numElementsToRemove);
            fId[i][j].erase(fId[i][j].begin() + eraseStart, fId[i][j].begin() + eraseStart + numElementsToRemove);
            diff[i][j].resize(newSize);
            P1[i][j].resize(newSize);
            P2[i][j].resize(newSize);
            frotP2[i][j].resize(newSize);
            Lam[i][j].resize(newSize);
            Id[i][j].resize(newSize);
        }
    }

    return newSize;
}

// Discrete Fourier Transform
void dft(vector<ComplexInterval> &coef, vector<ComplexInterval> &grid) {
    int N = grid.size();
    ComplexInterval sum, exp;
    for (int k = 0; k < N; k++) {
        sum = 0;
        for (int j = 0; j < N; j++) {
            exp = cExp(-2.0 * k * (double)j / N * PI);  // Note: -2*k*j/N is perfectly represented by the computer, but the product by PI is not
            sum += (grid[j] * exp);
        }
        coef[k] = (1. / N) * sum;  // Note: 1/N is perfectly represented by the computer
    }
}

// Inverse Discrete Fourier Transform
void idft(vector<ComplexInterval> &grid, vector<ComplexInterval> &coef) {
    int N = coef.size();
    ComplexInterval sum, exp;
    for (int k = 0; k < N; k++) {
        sum = 0;
        for (int j = 0; j < N; j++) {
            exp = cExp(2.0 * k * (double)j / N * PI);
            sum += coef[j] * exp;
        }
        grid[k] = sum;
    }
}

// Fast Fourier Transform using a grid sized as a power of 2
void fft_rec(vector<ComplexInterval> &x, int N) {
    // Check if it is split enough
    if (N <= 1) {
        return;
    }

    // Split even and odd
    vector<ComplexInterval> odd(N / 2);
    vector<ComplexInterval> even(N / 2);
    for (int i = 0; i < N / 2; i++) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // Split into tasks
    fft_rec(even, N / 2);
    fft_rec(odd, N / 2);

    // Calculate DFT
    for (int k = 0; k < N / 2; k++) {
        ComplexInterval t = (1. / 2 * odd[k]) * cExp(RealInterval(-2 * PI * (double)k / N));
        x[k] = 1. / 2 * even[k] + t;
        x[N / 2 + k] = 1. / 2 * even[k] - t;
    }
}

// Copy vectors and perform transform
void fft(vector<ComplexInterval> &coef, vector<ComplexInterval> &grid) {
    int N = grid.size();
    coef = grid;
    // Start recursion
    fft_rec(coef, N);
}

// Inverse Fast Fourier Transform
void ifft_rec(vector<ComplexInterval> &x, int N) {
    // Check if it is split enough
    if (N <= 1) {
        return;
    }

    // Split even and odd
    vector<ComplexInterval> odd(N / 2);
    vector<ComplexInterval> even(N / 2);
    for (int i = 0; i < N / 2; i++) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // Split into tasks
    ifft_rec(even, N / 2);
    ifft_rec(odd, N / 2);

    // Calculate DFT
    for (int k = 0; k < N / 2; k++) {
        ComplexInterval t = odd[k] * cExp(RealInterval(DPI * (double)k / N));
        x[k] = even[k] + t;
        x[N / 2 + k] = even[k] - t;
    }
}

// Copy vectors and perform inverse transform
void ifft(vector<ComplexInterval> &grid, vector<ComplexInterval> &coef) {
    int N = coef.size();
    grid = coef;
    // Start recursion
    ifft_rec(grid, N);
}

// Transform vectors evaluated over a grid into vectors of Fourier coefficients
void vectorFFT(ComplexVector &coef, ComplexVector &grid) {
    int n = grid.size();
    for (int i = 0; i < n; i++)
        fft(coef[i], grid[i]);
}

// Transform vectors of Fourier coefficients into vectors evaluated over a grid
void vectorIFFT(ComplexVector &grid, ComplexVector &coef) {
    int n = coef.size();
    for (int i = 0; i < n; i++)
        ifft(grid[i], coef[i]);
}

// Transform matrices evaluated over a grid into matrices of Fourier coefficients
void matrixFFT(ComplexMatrix &coef, ComplexMatrix &grid) {
    int n = grid.size();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fft(coef[i][j], grid[i][j]);
}

// Transform matrices of Fourier coefficients into matrices evaluated over a grid
void matrixIFFT(ComplexMatrix &grid, ComplexMatrix &coef) {
    int n = coef.size();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ifft(grid[i][j], coef[i][j]);
}

// Apply Standard Map to each entry of a vector of several components
void F(ComplexVector &FX, ComplexVector &X, RealInterval kappa, RealInterval eps) {
    int N = X[0].size();
    for (int j = 0; j < N; j++) {
        FX[1][j] = X[1][j] - (kappa / DPI) * sin(DPI * X[0][j]) - eps * sin((1. * j) / N * DPI);
        FX[0][j] = X[0][j] + FX[1][j];
    }
}

// Apply Standard Map over boxes of the grid
void Fbox(ComplexVector &FX, ComplexVector &X, RealInterval kappa, RealInterval eps, double complexWidth) {
    int N = X[0].size();
    vector<ComplexInterval> theta(N);
    for (int j = 0; j < N; j++) {
        // Although inflationRho is double, we are creating an interval around 0, so 0 will be for sure included
        setBounds(theta[j], (1. * j) / N - 1. / (2 * N), (1. * j) / N + 1. / (2 * N), -complexWidth, complexWidth);
        FX[1][j] = X[1][j] - (kappa / DPI) * sin(DPI * X[0][j]) - eps * sin(DPI * theta[j]);
        FX[0][j] = X[0][j] + FX[1][j];
    }
}

// Evaluate the Jacobian of the Standard Map over each entry of a vector with several components
void DF(ComplexMatrix &dif, ComplexVector &X, RealInterval kappa) {
    int N = X[0].size();
    for (int k = 0; k < N; k++) {
        dif[0][0][k] = ComplexInterval(1, 0) - (kappa * cos(DPI * X[0][k]));
        dif[0][1][k] = ComplexInterval(1, 0);
        dif[1][0][k] = kappa * (-1) * cos(DPI * X[0][k]);
        dif[1][1][k] = ComplexInterval(1, 0);
    }
}

// Evaluate the Hessian of the Standard Map over each entry of a vector with several components
// Since each component of F has the same Hessian, no need to compute the maximum
RealInterval DF2_supNorm(ComplexVector &X, RealInterval kappa) {
    int N = X[0].size();
    RealInterval norm, supr = 0;
    for (int k = 0; k < N; k++) {
        norm = mod((DPI * kappa) * sin(DPI * X[0][k]));
        supr = sup(supr, norm);
    }
    return supr;
}

// Compute C_N Fourier error bound
RealInterval CN(RealInterval rho, RealInterval rhoHat, int N) {
    RealInterval S1, S2, T, a, b, c;

    a = exp((-2 * N) * PI * rhoHat) / (one - exp(-2 * N * PI * rhoHat));
    b = (exp(-2 * PI * (rhoHat + rho)) + 1) / (exp(-2 * PI * (rhoHat + rho)) - 1);
    c = one - exp(N * PI * (rhoHat + rho));
    S1 = a * b;
    S1 = S1 * c;

    b = (exp(DPI * (rhoHat - rho)) + 1) / (exp(DPI * (rhoHat - rho)) - 1);
    c = one - exp(-N * PI * (rhoHat - rho));
    S2 = a * b;
    S2 = S2 * c;

    c = exp(-N * PI * (rhoHat - rho));
    T = b * c;

    return S1 + S2 + T;
}

// Find the rho and rho hat that minimize C_N for a given N
pair<double, double> findRhos(int N) {
    const double step = 1.e-2;
    const double maxRho = 1.0;
    double goodRho = 0.0, goodRhoHat = 0.0;
    RealInterval min = numeric_limits<double>::infinity();

    for (double rho = 0.0; rho <= maxRho; rho += step) {
        for (double rhoHat = rho; rhoHat <= maxRho; rhoHat += step) {
            RealInterval rhoInt = rho;
            RealInterval rhoHatInt = rhoHat;

            RealInterval C_N = CN(rhoInt, rhoHatInt, N);

            if (C_N < min) {
                min = C_N;
                goodRho = rho;
                goodRhoHat = rhoHat;
                // print(min);
            }
        }
    }

    return make_pair(goodRho, goodRhoHat);
}

// Rotate Fourier coefficients
void fourierRot(vector<ComplexInterval> &frotx, vector<ComplexInterval> &fx, RealInterval &omega) {
    int N = fx.size();
    ComplexInterval compExp;
    for (int k = 0; k < N; k++) {
        int idx = (k < (float)N / 2) ? k : (k - N);

        compExp = cExp(idx * DPI * omega);
        frotx[k] = fx[k] * compExp;
    }
}

// Rotate a vector by a real factor
void fourierRot(ComplexVector &frotx, ComplexVector &fx, RealInterval &omega) {
    int n = fx.size();

    for (int i = 0; i < n; i++)
        fourierRot(frotx[i], fx[i], omega);
}

// Rotate a matrix by a real factor
void fourierRot(ComplexMatrix &frotx, ComplexMatrix &fx, RealInterval &omega) {
    int n = fx.size();
    int m = fx[0].size();

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            fourierRot(frotx[i][j], fx[i][j], omega);
}

// Rotate a function by a complex factor, used to calculate K over a complex box
void fourierRot(vector<ComplexInterval> &frotx, vector<ComplexInterval> &fx, ComplexInterval &phi) {
    int N = fx.size();
    RealInterval auxImag, auxReal;

    mpfi_mul(auxImag.interval, DPI.interval, phi.imag);
    mpfi_mul(auxReal.interval, DPI.interval, phi.real);

    for (int k = 0; k < N; k++) {
        int idx = (k < (float)N / 2) ? k : (k - N);
        frotx[k] = (one / exp(idx * auxImag)) * fx[k] * cExp(idx * auxReal);
    }
}

// Rotate a vector by a complex factor, used to calculate K over a complex box
void fourierRot(ComplexVector &frotx, ComplexVector &fx, ComplexInterval &phi) {
    int n = fx.size();

    for (int i = 0; i < n; i++)
        fourierRot(frotx[i], fx[i], phi);
}

// Rotate a matrix by a complex factor, used to evaluate matrices on a complex box
void fourierRot(ComplexMatrix &frotx, ComplexMatrix &fx, ComplexInterval &phi) {
    int n = fx.size();
    int m = fx[0].size();

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            fourierRot(frotx[i][j], fx[i][j], phi);
}

// Thickens a function on a complex strip by applying the boxes method on an input on grid space
vector<ComplexInterval> thickenGridInput(vector<ComplexInterval> &gridInput, double realWidth, double complexWidth) {
    int N = gridInput.size();
    vector<ComplexInterval> fourierInput(N);
    // Bring our function to Fourier space
    fft(fourierInput, gridInput);

    // Create the desired box around 0
    ComplexInterval phi;
    setBounds(phi, -realWidth, realWidth, -complexWidth, complexWidth);

    // Evaluating a function over a box is the same as rotating it in Fourier space by the box around 0 (phi)
    vector<ComplexInterval> thickFourier(N);
    fourierRot(thickFourier, fourierInput, phi);

    // fourierBoxes contains now the function evaluated on the boxes but in Fourier space. We must bring it back to grid
    vector<ComplexInterval> thickGrid(N);
    ifft(thickGrid, thickFourier);

    return thickGrid;
}

ComplexVector thickenGridInput(ComplexVector &gridInput, double realWidth, double complexWidth) {
    ComplexVector thickInput;
    alloc(thickInput, gridInput.size(), gridInput[0].size());
    for (int i = 0; i < gridInput.size(); i++)
        thickInput[i] = thickenGridInput(gridInput[i], realWidth, complexWidth);
    return thickInput;
}

// Thickens a function on a complex strip by applying the boxes method on an input in Fourier space
vector<ComplexInterval> thickenFourierInput(vector<ComplexInterval> &fourierInput, double realWidth, double complexWidth) {
    int N = fourierInput.size();

    // Create the desired box around 0
    ComplexInterval phi;
    setBounds(phi, -realWidth, realWidth, -complexWidth, complexWidth);

    // Evaluating a function over a box is the same as rotating it in Fourier space by the box around 0 (phi)
    vector<ComplexInterval> thickFourier(N);
    fourierRot(thickFourier, fourierInput, phi);

    // fourierBoxes contains now the function evaluated on the boxes but in Fourier space. We must bring it back to grid
    vector<ComplexInterval> thickGrid(N);
    ifft(thickGrid, thickFourier);

    return thickGrid;
}

ComplexVector thickenFourierInput(ComplexVector &fourierInput, double realWidth, double complexWidth) {
    ComplexVector thickInput;
    alloc(thickInput, fourierInput.size(), fourierInput[0].size());
    for (int i = 0; i < fourierInput.size(); i++)
        thickInput[i] = thickenFourierInput(fourierInput[i], realWidth, complexWidth);
    return thickInput;
}

// Compute Fourier Norm
RealInterval fourierNorm(vector<ComplexInterval> &coef, RealInterval &rho) {
    RealInterval sum = 0;
    int N = coef.size();
    for (int k = 0; k < N; k++) {
        int idx = (k < (float)N / 2) ? k : (k - N);
        sum += mod(coef[k]) * exp(abs(idx) * (DPI * rho));
    }
    return sum;
}

// Compute Fourier Norm of a vector
RealInterval fourierNorm(ComplexVector &x, RealInterval &rho) {
    RealInterval supr = fourierNorm(x[0], rho);
    for (int i = 1; i < x.size(); i++) {
        RealInterval rowSup = fourierNorm(x[i], rho);
        supr = sup(rowSup, supr);
    }
    return supr;
}

// Compute Fourier Norm of a matrix
RealInterval fourierNorm(ComplexMatrix &x, RealInterval &rho) {
    RealInterval maxSum = 0;
    for (int i = 0; i < x.size(); i++) {
        RealInterval rowSum = 0;
        for (int j = 0; j < x[i].size(); j++) {
            rowSum += fourierNorm(x[i][j], rho);
        }
        maxSum = sup(maxSum, rowSum);
    }

    return maxSum;
}

// Compute Supremum Norm of an array of complex intervals
RealInterval supNorm(const vector<ComplexInterval> &x) {
    RealInterval supr;
    supr = 0;
    for (const ComplexInterval &interval : x) {
        RealInterval modulus = mod(interval);
        supr = sup(supr, modulus);
    }
    return supr;
}

// Compute Supremum Norm of a vector
RealInterval supNorm(const ComplexVector &x) {
    RealInterval supr = supNorm(x[0]);

    for (int i = 1; i < x.size(); i++) {
        RealInterval rowSup = supNorm(x[i]);
        supr = sup(rowSup, supr);
    }

    return supr;
}

// Compute Supremum Norm of a matrix
RealInterval supNorm(ComplexMatrix &x) {
    RealInterval maxSum = 0;
    for (int i = 0; i < x.size(); i++) {
        RealInterval rowSum = 0;
        for (int j = 0; j < x[i].size(); j++) {
            rowSum += supNorm(x[i][j]);
        }
        maxSum = sup(maxSum, rowSum);
    }

    return maxSum;
}

// Compute the rho norm on a complex strip by applying the boxes method
RealInterval rhoNorm(vector<ComplexInterval> &function, double complexWidth) {
    int N = function.size();
    // Calculates the supremum of all the entries of the thickened function
    return supNorm(thickenGridInput(function, 1. / (2 * N), complexWidth));
}

// Compute the rho norm of a vector of functions on a complex strip by applying the boxes method
RealInterval rhoNorm(ComplexVector &gridInput, double complexWidth) {
    RealInterval supr;
    supr = rhoNorm(gridInput[0], complexWidth);

    for (int i = 1; i < gridInput.size(); i++) {
        RealInterval rowSup = rhoNorm(gridInput[i], complexWidth);
        supr = sup(rowSup, supr);
    }

    return supr;
}

// Compute the rho norm of a matrix of functions on a complex strip by applying the boxes method
RealInterval rhoNorm(ComplexMatrix &gridInput, double complexWidth) {
    RealInterval maxSum = 0;
    for (int i = 0; i < gridInput.size(); i++) {
        RealInterval rowSum = 0;
        for (int j = 0; j < gridInput[i].size(); j++) {
            rowSum += rhoNorm(gridInput[i][j], complexWidth);
        }
        maxSum = sup(maxSum, rowSum);
    }

    return maxSum;
}
