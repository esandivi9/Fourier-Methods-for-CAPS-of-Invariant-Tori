#include <gtest/gtest.h>

#include <fstream>

TEST(FourierTests, DFT) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 1), ComplexInterval(0, 0), ComplexInterval(1, 1), ComplexInterval(0, 0)};
    vector<ComplexInterval> fourierCoeffs(4), grid(4), expectedCoeffs = {ComplexInterval(0.5, 0.5), ComplexInterval(0, 0), ComplexInterval(0.5, 0.5), ComplexInterval(0, 0)};

    dft(fourierCoeffs, vec);
    idft(grid, fourierCoeffs);

    EXPECT_TRUE(grid.size() == vec.size());

    // Check if the Fourier coefficients are correct and if applying the inverse transform gives us the original data
    for (int i = 0; i < vec.size(); i++) {
        EXPECT_TRUE(expectedCoeffs[i].isContained(fourierCoeffs[i]));
        EXPECT_TRUE(vec[i].isContained(grid[i]));
    }

    vec = {ComplexInterval(2, 0), ComplexInterval(1, 0), ComplexInterval(1, 0), ComplexInterval(2, 0)};
    expectedCoeffs = {ComplexInterval("1.5", "0"), ComplexInterval(1. / 4, 1. / 4), ComplexInterval(0, 0), ComplexInterval(1. / 4, -1. / 4)};

    fft(fourierCoeffs, vec);
    ifft(grid, fourierCoeffs);

    EXPECT_TRUE(grid.size() == vec.size());

    for (int i = 0; i < vec.size(); i++) {
        EXPECT_TRUE(expectedCoeffs[i].isContained(fourierCoeffs[i]));
        EXPECT_TRUE(vec[i].isContained(grid[i]));
    }
}

TEST(FourierTests, FFT) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 1), ComplexInterval(0, 0), ComplexInterval(1, 1), ComplexInterval(0, 0)};
    vector<ComplexInterval> fourierCoeffs(4), grid(4), expectedCoeffs = {ComplexInterval(0.5, 0.5), ComplexInterval(0, 0), ComplexInterval(0.5, 0.5), ComplexInterval(0, 0)};

    fft(fourierCoeffs, vec);
    ifft(grid, fourierCoeffs);

    EXPECT_TRUE(grid.size() == vec.size());

    // Check if the Fourier coefficients are correct and if applying the inverse transform gives us the original data
    for (int i = 0; i < vec.size(); i++) {
        EXPECT_TRUE(expectedCoeffs[i].isContained(fourierCoeffs[i]));
        EXPECT_TRUE(vec[i].isContained(grid[i]));
    }

    vec = {ComplexInterval(2, 0), ComplexInterval(1, 0), ComplexInterval(1, 0), ComplexInterval(2, 0)};
    expectedCoeffs = {ComplexInterval("1.5", "0"), ComplexInterval(1. / 4, 1. / 4), ComplexInterval(0, 0), ComplexInterval(1. / 4, -1. / 4)};

    fft(fourierCoeffs, vec);
    ifft(grid, fourierCoeffs);

    EXPECT_TRUE(grid.size() == vec.size());

    for (int i = 0; i < vec.size(); i++) {
        EXPECT_TRUE(expectedCoeffs[i].isContained(fourierCoeffs[i]));
        EXPECT_TRUE(vec[i].isContained(grid[i]));
    }
}

TEST(FourierTests, VectorFourierTests) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 1), ComplexInterval(0, 0), ComplexInterval(1, 1), ComplexInterval(0, 0)};
    vector<ComplexInterval> expectedCoeffs = {ComplexInterval(0.5, 0.5), ComplexInterval(0, 0), ComplexInterval(0.5, 0.5), ComplexInterval(0, 0)};

    ComplexVector vector, coefVector, expectedCoefVector, gridVector;
    allocAllVectors(2, 4, vector, coefVector, expectedCoefVector, gridVector);

    // Assign data to the vector
    for (int i = 0; i < 2; i++) {
        vector[i] = vec;
        expectedCoefVector[i] = expectedCoeffs;
    }

    vectorFFT(coefVector, vector);
    vectorIFFT(gridVector, coefVector);

    // Check if the Fourier coefficients are correct and if applying the inverse transform gives us the original data
    // for every entry in our vector
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < vec.size(); k++) {
            EXPECT_TRUE(expectedCoefVector[i][k].isContained(coefVector[i][k]));
            EXPECT_TRUE(vector[i][k].isContained(gridVector[i][k]));
        }
}

TEST(FourierTests, MatrixFourierTests) {
    vector<ComplexInterval> vec = {ComplexInterval(1, 1), ComplexInterval(0, 0), ComplexInterval(1, 1), ComplexInterval(0, 0)};
    vector<ComplexInterval> expectedCoeffs = {ComplexInterval(0.5, 0.5), ComplexInterval(0, 0), ComplexInterval(0.5, 0.5), ComplexInterval(0, 0)};

    ComplexMatrix matrix, coefMatrix, expectedCoefMatrix, gridMatrix;
    allocAllMatrices(2, 2, 4, matrix, coefMatrix, expectedCoefMatrix, gridMatrix);

    // Assign data to the matrix
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            matrix[i][j] = vec;
            expectedCoefMatrix[i][j] = expectedCoeffs;
        }

    matrixFFT(coefMatrix, matrix);
    matrixIFFT(gridMatrix, coefMatrix);

    // Check if the Fourier coefficients are correct and if applying the inverse transform gives us the original data
    // for every entry in our matrix
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < vec.size(); k++) {
                EXPECT_TRUE(expectedCoefMatrix[i][j][k].isContained(coefMatrix[i][j][k]));
                EXPECT_TRUE(matrix[i][j][k].isContained(gridMatrix[i][j][k]));
            }
}

TEST(StandardMapTests, F) {
    ComplexVector K, FK;
    RealInterval kappa = 1, eps = 0;
    alloc(K, 2, 2);
    alloc(FK, 2, 2);

    // Check whether (1/2, 0) is a fixed point
    K[0][0] = ComplexInterval(0.5, 0);
    K[0][1] = ComplexInterval(0.5, 0);
    K[1][0] = 0;
    K[1][1] = 0;
    F(FK, K, kappa, eps);
    ComplexVector expectedResult = K;
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(FK[i][k]));

    // Check whether (3/2, 0) is a fixed point
    K[0][0] = ComplexInterval("1.5", "0");
    K[0][1] = ComplexInterval("1.5", "0");
    F(FK, K, kappa, eps);
    expectedResult = K;
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(FK[i][k]));

    // Check whether (1, 1) maps to (2, 1)
    K[0][0] = ComplexInterval(1, 0);
    K[0][1] = ComplexInterval(1, 0);
    K[1][0] = ComplexInterval(1, 0);
    K[1][1] = ComplexInterval(1, 0);
    F(FK, K, kappa, eps);
    expectedResult = K;
    expectedResult[0][0] = ComplexInterval(2, 0);
    expectedResult[0][1] = ComplexInterval(2, 0);
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(FK[i][k]));

    // Check whether (1+i, 1-i) maps to (2, 1-i) for kappa = 0
    kappa = 0;
    K[0][0] = ComplexInterval(1, 1);
    K[0][1] = ComplexInterval(1, 1);
    K[1][0] = ComplexInterval(1, -1);
    K[1][1] = ComplexInterval(1, -1);
    F(FK, K, kappa, eps);
    expectedResult[0][0] = ComplexInterval(2, 0);
    expectedResult[0][1] = ComplexInterval(2, 0);
    expectedResult[1][0] = ComplexInterval(1, -1);
    expectedResult[1][1] = ComplexInterval(1, -1);
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(FK[i][k]));

    // Test the perturbed Standard Map with eps = 1 and preimage (1,1)
    K.clear();
    FK.clear();
    alloc(K, 2, 4);
    alloc(FK, 2, 4);
    eps = 1;
    K[0][0] = ComplexInterval(1, 0);
    K[0][1] = ComplexInterval(1, 0);
    K[0][2] = ComplexInterval(1, 0);
    K[0][3] = ComplexInterval(1, 0);

    K[1][0] = ComplexInterval(1, 0);
    K[1][1] = ComplexInterval(1, 0);
    K[1][2] = ComplexInterval(1, 0);
    K[1][3] = ComplexInterval(1, 0);

    F(FK, K, kappa, eps);
    expectedResult.clear();
    alloc(expectedResult, 2, 4);
    expectedResult[0][0] = ComplexInterval(2, 0);
    expectedResult[0][1] = ComplexInterval(1, 0);
    expectedResult[0][2] = ComplexInterval(2, 0);
    expectedResult[0][3] = ComplexInterval(3, 0);

    expectedResult[1][0] = ComplexInterval(1, 0);
    expectedResult[1][1] = ComplexInterval(0, 0);
    expectedResult[1][2] = ComplexInterval(1, 0);
    expectedResult[1][3] = ComplexInterval(2, 0);
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 4; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(FK[i][k]));
}

TEST(StandardMapTests, FBoxes) {
    ComplexVector K, FK;
    RealInterval kappa = 1, eps = 1;
    double complexWidth = 0.01;
    alloc(K, 2, 2);
    alloc(FK, 2, 2);

    // The real part of the boxes (theta) will be [-1/4, 1/4] and [1/4, 3/4] (enclosures for 0 and 1/2),
    // which when combined with the fixed point (1/2, 0) will give the image of the box around it, [-1/2, 3/2] and [-1, 1]
    K[0][0] = ComplexInterval(0.5, 0);
    K[0][1] = ComplexInterval(0.5, 0);
    K[1][0] = 0;
    K[1][1] = 0;
    Fbox(FK, K, kappa, eps, complexWidth);
    ComplexVector expectedResult;
    alloc(expectedResult, 2, 2);
    for (int k = 0; k < 2; k++) {
        setBounds(expectedResult[0][k], -1. / 2, 3 * 1. / 2, 0, 0);
        setBounds(expectedResult[1][k], -1, 1, 0, 0);
    }

    for (int i = 0; i < 2; i++)
        for (int k = 0; k < 2; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(FK[i][k]));
}

TEST(StandardMapTests, DF) {
    ComplexMatrix dif, expectedMat;
    allocAllMatrices(2, 2, 2, dif, expectedMat);

    ComplexVector K;
    alloc(K, 2, 2);
    K[0][0] = ComplexInterval(0.5, 0);
    K[0][1] = ComplexInterval(0.5, 0);
    K[1][0] = 0;
    K[1][1] = 0;

    RealInterval kappa = 1;
    DF(dif, K, kappa);

    for (int k = 0; k < 2; k++) {
        expectedMat[0][0][k] = ComplexInterval(2, 0);
        expectedMat[0][1][k] = ComplexInterval(1, 0);
        expectedMat[1][0][k] = ComplexInterval(1, 0);
        expectedMat[1][1][k] = ComplexInterval(1, 0);
    }

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                EXPECT_TRUE(expectedMat[i][j][k].isContained(dif[i][j][k]));
}

TEST(StandardMapTests, DF2_norm) {
    ComplexVector K;
    RealInterval kappa = 1;
    alloc(K, 2, 2);
    K[0][0] = ComplexInterval(0.25, 0);
    K[0][1] = ComplexInterval(1, 0);
    K[1][0] = 0;
    K[1][1] = 0;
    RealInterval norm = DF2_supNorm(K, kappa);
    RealInterval expectedValue = DPI;
    EXPECT_TRUE(expectedValue.isContained(norm));
}

// The following test is to play around with rho, rhoHat and N and see how CN changes.
TEST(FourierBoundsTests, CNBound) {
    RealInterval rho = 1e-3, rhoHat = 1e-2;
    int N = 4096;
    auto CNBound = CN(rho, rhoHat, N);
    // cout << "CN: ";
    // print(CNBound);
}

TEST(FourierBoundsTests, CNBoundPlay) {
    //int N = 1024;
    //auto pair = findRhos(N);
    //cout << pair.first << "\t" << pair.second << endl;
}

TEST(FourierBoundsTests, FourierRotation) {
    int N = 4;
    vector<ComplexInterval> vec = {ComplexInterval(1, 1), ComplexInterval(1, -1), ComplexInterval(2, 0), ComplexInterval(-3, 4)};
    vector<ComplexInterval> rotated(N);
    vector<ComplexInterval> expectedResult = vec;

    // First test identity rotation
    RealInterval omega = 0;
    fourierRot(rotated, vec, omega);
    for (int i = 0; i < N; i++)
        EXPECT_TRUE(expectedResult[i].isContained(rotated[i]));

    omega = "1.5";
    fourierRot(rotated, vec, omega);
    expectedResult = {ComplexInterval(1, 1), ComplexInterval(-1, 1), ComplexInterval(2, 0), ComplexInterval(3, -4)};
    for (int i = 0; i < N; i++)
        EXPECT_TRUE(expectedResult[i].isContained(rotated[i]));

    omega = "3.7";
    fourierRot(rotated, vec, omega);
    expectedResult[0] = 1;
    mpfi_interv_d(expectedResult[1].real, -1.27, -1.26);
    mpfi_interv_d(expectedResult[1].imag, -0.65, -0.64);
    mpfi_interv_d(expectedResult[2].real, -1.62, -1.61);
    mpfi_interv_d(expectedResult[2].imag, -1.18, -1.17);
    mpfi_interv_d(expectedResult[3].real, -2.88, -2.87);
    mpfi_interv_d(expectedResult[3].imag, -4.09, -4.08);

    for (int i = 0; i < N; i++)
        EXPECT_TRUE(rotated[i].isContained(expectedResult[i]));
}

TEST(FourierBoundsTests, FourierComplexRotationOnVector) {
    int N = 4;
    ComplexVector vec;
    alloc(vec, 2, N);
    vec[0] = {ComplexInterval(1, 1), ComplexInterval(1, -1), ComplexInterval(2, 0), ComplexInterval(-3, 4)};
    vec[1] = {ComplexInterval(1, 1), ComplexInterval(1, -1), ComplexInterval(2, 0), ComplexInterval(-3, 4)};
    ComplexVector rotated;
    alloc(rotated, 2, N);
    ComplexVector expectedResult = vec;

    // First test identity rotation
    ComplexInterval phi = 0;
    fourierRot(rotated, vec, phi);
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < N; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(rotated[i][k]));

    phi = ComplexInterval("1.5", "0");
    fourierRot(rotated, vec, phi);
    expectedResult[0] = {ComplexInterval(1, 1), ComplexInterval(-1, 1), ComplexInterval(2, 0), ComplexInterval(3, -4)};
    expectedResult[1] = {ComplexInterval(1, 1), ComplexInterval(-1, 1), ComplexInterval(2, 0), ComplexInterval(3, -4)};
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < N; k++)
            EXPECT_TRUE(expectedResult[i][k].isContained(rotated[i][k]));

    phi = ComplexInterval("3.7", "-0.3");
    fourierRot(rotated, vec, phi);
    for (int i = 0; i < 2; i++) {
        expectedResult[i][0] = 1;
        mpfi_interv_d(expectedResult[i][1].real, -8.36, -8.29);
        mpfi_interv_d(expectedResult[i][1].imag, -4.28, -4.21);
        mpfi_interv_d(expectedResult[i][2].real, -0.04, -0.03);
        mpfi_interv_d(expectedResult[i][2].imag, -0.03, -0.02);
        mpfi_interv_d(expectedResult[i][3].real, -0.44, -0.43);
        mpfi_interv_d(expectedResult[i][3].imag, -0.63, -0.62);
    }

    for (int i = 0; i < 2; i++)
        for (int k = 0; k < N; k++)
            EXPECT_TRUE(rotated[i][k].isContained(expectedResult[i][k]));
}

TEST(FourierBoundsTests, FourierComplexRotationOnMatrix) {
    int N = 4;
    ComplexMatrix vec;
    alloc(vec, 2, 2, N);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            vec[i][j] = {ComplexInterval(1, 1), ComplexInterval(1, -1), ComplexInterval(2, 0), ComplexInterval(-3, 4)};
    ComplexMatrix rotated;
    alloc(rotated, 2, 2, N);
    ComplexMatrix expectedResult = vec;

    // First test identity rotation
    ComplexInterval phi = 0;
    fourierRot(rotated, vec, phi);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < N; k++)
                EXPECT_TRUE(expectedResult[i][j][k].isContained(rotated[i][j][k]));

    phi = ComplexInterval("1.5", "0");
    fourierRot(rotated, vec, phi);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            expectedResult[i][j] = {ComplexInterval(1, 1), ComplexInterval(-1, 1), ComplexInterval(2, 0), ComplexInterval(3, -4)};

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < N; k++)
                EXPECT_TRUE(expectedResult[i][j][k].isContained(rotated[i][j][k]));

    phi = ComplexInterval("3.7", "-0.3");
    fourierRot(rotated, vec, phi);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            expectedResult[i][j][0] = 1;
            mpfi_interv_d(expectedResult[i][j][1].real, -8.36, -8.29);
            mpfi_interv_d(expectedResult[i][j][1].imag, -4.28, -4.21);
            mpfi_interv_d(expectedResult[i][j][2].real, -0.04, -0.03);
            mpfi_interv_d(expectedResult[i][j][2].imag, -0.03, -0.02);
            mpfi_interv_d(expectedResult[i][j][3].real, -0.44, -0.43);
            mpfi_interv_d(expectedResult[i][j][3].imag, -0.63, -0.62);
        }

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < N; k++)
                EXPECT_TRUE(rotated[i][j][k].isContained(expectedResult[i][j][k]));
}

TEST(FourierNormTests, FourierNormOfSeries) {
    vector<ComplexInterval> coeffs = {ComplexInterval("1.5", "0"), ComplexInterval(1. / 4, 1. / 4), ComplexInterval(0, -1), ComplexInterval(1. / 4, -1. / 4)};
    RealInterval rho = "1.e-3";
    auto norm = fourierNorm(coeffs, rho);
    RealInterval expectedResult;
    mpfi_interv_d(expectedResult.interval, 3.22, 3.23);
    EXPECT_TRUE(norm.isContained(expectedResult));

    coeffs.resize(3);
    coeffs = {ComplexInterval(-2, 0), ComplexInterval(-1, 4), ComplexInterval(-5, 0)};
    rho = 1;
    norm = fourierNorm(coeffs, rho);
    mpfi_interv_d(expectedResult.interval, 4887.34, 4887.35);
    EXPECT_TRUE(norm.isContained(expectedResult));
}

TEST(FourierNormTests, FourierNormOfVectorOfSeries) {
    ComplexVector vec;
    alloc(vec, 2, 4);
    vec[0] = {ComplexInterval("1.5", "0"), ComplexInterval(1. / 4, 1. / 4), ComplexInterval(0, -1), ComplexInterval(1. / 4, -1. / 4)};
    vec[1] = {ComplexInterval("0.5", "0"), ComplexInterval(1. / 8, 1. / 4), ComplexInterval(0, -1), ComplexInterval(1. / 4, -1. / 8)};
    RealInterval rho = "1.e-3";
    auto norm = fourierNorm(vec, rho);
    RealInterval expectedResult;
    mpfi_interv_d(expectedResult.interval, 3.22, 3.23);
    EXPECT_TRUE(norm.isContained(expectedResult));

    vec.clear();
    alloc(vec, 4, 3);
    rho = 1;
    vec[0] = {ComplexInterval(-2, 0), ComplexInterval(-1, 4), ComplexInterval(-5, 0)};
    vec[1] = {ComplexInterval(-1, 0), ComplexInterval(-1, 2), ComplexInterval(-4, 0)};
    vec[2] = {ComplexInterval(-1, 1), ComplexInterval(-1. / 2, 3), ComplexInterval(2, 2)};
    vec[3] = {ComplexInterval(0, 0), ComplexInterval(0, -1), ComplexInterval(3, -1)};
    norm = fourierNorm(vec, rho);
    mpfi_interv_d(expectedResult.interval, 4887.34, 4887.35);
    EXPECT_TRUE(norm.isContained(expectedResult));
}

TEST(FourierNormTests, FourierNormOfMatrixOfSeries) {
    ComplexMatrix mat;
    alloc(mat, 2, 2, 4);
    mat[0][0] = {ComplexInterval(0, 0), ComplexInterval(1, -1), ComplexInterval(2, -3), ComplexInterval(-20, 0)};
    mat[0][1] = {ComplexInterval(2, 3), ComplexInterval(5, 5), ComplexInterval(1, -4), ComplexInterval(0, 10)};
    mat[1][0] = {ComplexInterval(-1, 3), ComplexInterval(1, 2), ComplexInterval(3, -4), ComplexInterval(5, -5)};
    mat[1][1] = {ComplexInterval(1, 0), ComplexInterval(-2, -3), ComplexInterval(0, -1), ComplexInterval(3, 7)};
    RealInterval rho = "1.e-3";
    auto norm = fourierNorm(mat, rho);
    RealInterval expectedResult;
    mpfi_interv_d(expectedResult.interval, 50.15, 50.16);
    EXPECT_TRUE(norm.isContained(expectedResult));

    mat.clear();
    alloc(mat, 4, 2, 3);
    rho = 1;
    mat[0][0] = {0, ComplexInterval("1.2", "3.2"), -1};
    mat[0][1] = {1, ComplexInterval("0", "4.5"), 2};
    mat[2][0] = {-4, ComplexInterval(-2, 3), -1};
    mat[2][1] = {2, ComplexInterval(-1, 4), -3};
    mat[3][0] = {3, ComplexInterval(3, 2), -1};
    mat[3][1] = {1, ComplexInterval(1, 3), -5};
    norm = fourierNorm(mat, rho);
    mpfi_interv_d(expectedResult.interval, 8173.56, 8173.58);
    EXPECT_TRUE(norm.isContained(expectedResult));
}

TEST(SupremumNormTests, SupNormOfFunction) {
    vector<ComplexInterval> vec = {ComplexInterval(0, 0), ComplexInterval(1, -1), ComplexInterval(2, -3), ComplexInterval(-20, 0)};
    auto norm = supNorm(vec);
    RealInterval expectedResult = 20;
    EXPECT_TRUE(expectedResult.isContained(norm));

    vec.resize(3);
    vec = {-2, ComplexInterval(0, 4), -5};
    norm = supNorm(vec);
    mpfi_interv_d(expectedResult.interval, 7.07, 7.08);
    EXPECT_TRUE(norm.isContained(expectedResult));
}

TEST(SupremumNormTests, SupNormOfVectorOfFunctions) {
    ComplexVector vec;
    alloc(vec, 2, 4);
    vec[0] = {ComplexInterval(0, 0), ComplexInterval(1, -1), ComplexInterval(2, -3), ComplexInterval(-20, 0)};
    vec[1] = {ComplexInterval(2, 3), ComplexInterval(5, 5), ComplexInterval(1, -4), ComplexInterval(0, 10)};
    auto norm = supNorm(vec);
    RealInterval expectedResult = 20;
    EXPECT_TRUE(expectedResult.isContained(norm));

    vec.clear();
    alloc(vec, 4, 3);
    vec[0] = {0, ComplexInterval("1.2", "3.2"), -1};
    vec[1] = {1, ComplexInterval("0", "4.5"), 2};
    vec[2] = {-4, ComplexInterval(-2, 3), -1};
    vec[3] = {4, ComplexInterval(0, 4), -5};

    norm = supNorm(vec);
    mpfi_interv_d(expectedResult.interval, 7.07, 7.08);
    EXPECT_TRUE(norm.isContained(expectedResult));
}

TEST(SupremumNormTests, SupNormOfMatrixOfFunctions) {
    ComplexMatrix mat;
    alloc(mat, 2, 2, 4);
    mat[0][0] = {ComplexInterval(0, 0), ComplexInterval(1, -1), ComplexInterval(2, -3), ComplexInterval(-20, 0)};
    mat[0][1] = {ComplexInterval(2, 3), ComplexInterval(5, 5), ComplexInterval(1, -4), ComplexInterval(0, 10)};
    mat[1][0] = {ComplexInterval(-1, 3), ComplexInterval(1, 2), ComplexInterval(3, -4), ComplexInterval(5, -5)};
    mat[1][1] = {ComplexInterval(1, 0), ComplexInterval(-2, -3), ComplexInterval(0, -1), ComplexInterval(3, 7)};

    auto norm = supNorm(mat);
    RealInterval expectedResult = 30;
    EXPECT_TRUE(expectedResult.isContained(norm));

    mat.clear();
    alloc(mat, 4, 2, 3);
    mat[0][0] = {0, ComplexInterval("1.2", "3.2"), -1};
    mat[0][1] = {1, ComplexInterval("0", "4.5"), 2};
    mat[2][0] = {-4, ComplexInterval(-2, 3), -1};
    mat[2][1] = {2, ComplexInterval(-1, 4), -3};
    mat[3][0] = {3, ComplexInterval(3, 2), -1};
    mat[3][1] = {1, ComplexInterval(1, 3), -5};

    norm = supNorm(mat);
    mpfi_interv_d(expectedResult.interval, 11.31, 11.33);
    EXPECT_TRUE(norm.isContained(expectedResult));
}

TEST(RhoNorms, RhoNorm) {
    vector<ComplexInterval> vec = {ComplexInterval(0, 0), ComplexInterval(1, -1), ComplexInterval(2, -3), ComplexInterval(-20, 0)};
    double rho = 1.e-3;
    RealInterval rNorm = rhoNorm(vec, rho);
    vector<ComplexInterval> coeffs(vec.size());
    fft(coeffs, vec);
    RealInterval rhoInterval = rho;
    auto fNorm = fourierNorm(coeffs, rhoInterval);

    // All the processes involved in calculating the rho norm make the intervals wider to accommodate for errors. Although this value actually is smaller
    // than the Fourier norm, the right extreme of its interval can be bigger than that of the Fourier norm interval, giving thus a wrong comparison.
    // That is why it can be more meaningful to check for the center of the intervals for comparison.
    mpfr_t midPoint1, midPoint2;
    mpfr_init2(midPoint1, prec);
    mpfr_init2(midPoint2, prec);
    mpfi_mid(midPoint1, rNorm.interval);
    mpfi_mid(midPoint2, fNorm.interval);
    EXPECT_TRUE(mpfr_cmp(midPoint1, midPoint2) < 0);
    mpfr_clear(midPoint1);
    mpfr_clear(midPoint2);

    coeffs = {ComplexInterval("1.5", "0"), ComplexInterval(1. / 4, 1. / 4), ComplexInterval(0, -1), ComplexInterval(1. / 4, -1. / 4)};
    fNorm = fourierNorm(coeffs, rhoInterval);
    ifft(vec, coeffs);
    rNorm = rhoNorm(vec, rho);

    mpfr_init2(midPoint1, prec);
    mpfr_init2(midPoint2, prec);
    mpfi_mid(midPoint1, rNorm.interval);
    mpfi_mid(midPoint2, fNorm.interval);
    EXPECT_TRUE(mpfr_cmp(midPoint1, midPoint2) < 0);
}