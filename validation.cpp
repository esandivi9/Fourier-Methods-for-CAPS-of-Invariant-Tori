#include <fstream>

#include "Functions.h"

int main() {
    // Specify the subfolder and file name
    string subfolder = "Inputs/";
    string fileName = subfolder + "K1.234200.txt";

    // Read N from file
    ifstream inputFile(fileName);
    if (!inputFile) {
        cerr << "File Error" << endl;
        return 0;  // or handle the error accordingly
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
    double LamUnst, LamSt, e, kap;
    vector<vector<double>> Kc(m, vector<double>(N));
    vector<vector<vector<double>>> P0(m, vector<vector<double>>(m, vector<double>(N)));
    readData(fileName, m, N, e, kap, Kc, P0, K, P1, Lam, Id, LamUnst, LamSt);
    RealInterval eps = e, kappa = kap;
    double lamSt = abs(LamSt);
    double lamUnst = abs(1. / LamUnst);

    // Generate P2 by inverting P1
    P2 = inverse(P1);

    // Realify by setting Nyquist term to 0, this will be the object of the validation
    // First turn to Fourier space
    vectorFFT(fK, K);
    matrixFFT(fP1, P1);
    matrixFFT(fP2, P2);
    // Since the Lambda matrix will be our input, we can adapt it as needed
    matrixFFT(fLam, Lam);
    // Turn also the Identity matrix to Fourier space as will be needed later
    matrixFFT(fId, Id);

    // Now set Nyquist term to 0
    for (int i = 0; i < m; i++) {
        fK[i][N / 2] = 0;
        for (int j = 0; j < m; j++) {
            fP1[i][j][N / 2] = 0;
            fP2[i][j][N / 2] = 0;
            fLam[i][j][N / 2] = 0;
            fId[i][j][N / 2] = 0;
        }
    }

    RealInterval zero = 0;
    // Check everything is alright by checking that the coefficients are complex conjugates
    for (int k = 1; k < N / 2; k++) {
        for (int i = 0; i < m; i++) {
            // If the imaginary parts add up to 0 they are conjugate
            if (mpfi_cmp((fK[i][k] + fK[i][N - k]).imag, zero.interval) == 0)
                continue;
            else
                cout << "Fourier coeffs of K not complex conjugate" << endl;
            for (int j = 0; j < m; j++) {
                if (mpfi_cmp((fP1[i][j][k] + fP1[i][j][N - k]).imag, zero.interval) == 0 && mpfi_cmp((fP2[i][j][k] + fP2[i][j][N - k]).imag, zero.interval) == 0)
                    continue;
                else
                    cout << "Fourier coeffs of P1 or P2 not complex conjugate" << endl;
            }
        }
    }

    //====================================================================================
    // Duplicate size of vectors by adding 0s in the middle if needed
    // N = increaseGridSize(2, fK, K, FK, fFK, frotK, fP1, fP2, fLam, fId, diff, P1, P2, frotP2, Lam, Id);

    //====================================================================================
    // Remove elements from the middle of the vectors to remove noise

    // For epsilon = 1.2342 reduce 3./4. For 0.5 or 1 reduce 15./16
    N = reduceGridSize(3. / 4, fK, K, FK, fFK, frotK, fP1, fP2, fLam, fId, diff, P1, P2, frotP2, Lam, Id);

    //====================================================================================

    // Go back to grid space
    vectorIFFT(K, fK);
    matrixIFFT(P1, fP1);
    matrixIFFT(P2, fP2);

    //====================================================================================
    // Print Fourier coefficients of K and P

    printFourier(fileName, fK);
    printFourier(fileName, fP1);

    //====================================================================================
    // Invariance Error

    // Start by rotating the torus
    RealInterval omega;
    omega = 5;
    omega = (sqrt(omega) - 1) / 2;

    fourierRot(frotK, fK, omega);
    F(FK, K, kappa, eps);
    vectorFFT(fFK, FK);
    // Calculate the difference
    ComplexVector fE;
    alloc(fE, m, N);
    fE = fFK - frotK;

    // Compute the norm of the difference
    double rho = 7.e-4;
    RealInterval rhoInterv = rho;
    RealInterval fENorm = fourierNorm(fE, rhoInterv);

    // Create thicker K0 by boxes method and apply F over it
    double rhoHat = 4.e-3;
    ComplexVector thickK = thickenFourierInput(fK, 1. / (2 * N), rhoHat);

    // Note: we avoid computing the fft as much as possible, as it worsens the coefficients, especially for bigger N,
    // hence we use the Fourier coefficients we have from the beginning

    ComplexVector thickFK;
    alloc(thickFK, m, N);

    RealInterval rhoHatInterv = rhoHat;
    Fbox(thickFK, thickK, kappa, eps, rhoHat);
    auto supr = supNorm(thickFK);
    RealInterval C_N = CN(rhoInterv, rhoHatInterv, N);
    cout << "Epsilon: " << e << endl;
    cout << "N: " << N << endl;
    cout << "Rho: " << rho << endl;
    cout << "RhoHat: " << rhoHat << endl;
    cout << "lambda_s: " << lamSt << endl;
    cout << "lambda_u: " << LamUnst << endl;
    double lam = 1 - 1. / max(1. / (1 - lamSt), lamUnst / (1 - lamUnst));
    cout << "lambda: " << lam << endl;
    cout << "C_N: ";
    print(C_N);
    RealInterval invarianceError;
    invarianceError = C_N * supr + fENorm;
    /*cout << "||F(K_0)||_rhoHat: ";
    print(supr);
    cout << "||F^tilde(K_0) - K_0||_Fourier, rho: ";
    print(fENorm);
    cout << "Invariance error: ";
    print(invarianceError);*/

    //====================================================================================

    // Reducibility Error

    // First term
    // Start by computing P2(theta + omega) by applying the rotation by omega in Fourier space
    fourierRot(frotP2, fP2, omega);
    // Compute Fourier norm since object is in Fourier space
    RealInterval rotP2Norm = fourierNorm(frotP2, rhoHatInterv);

    // Now compute M_0(theta), which is the Differential of F but evaluated on a complex strip, that is, a thickened torus thickK
    ComplexMatrix thickDFK;
    alloc(thickDFK, m, m, N);
    DF(thickDFK, thickK, kappa);

    // Compute rhoHat norm since the object is in grid space
    RealInterval M0Norm = supNorm(thickDFK);

    // Compute Fourier norm of P1 since we have it in Fourier space
    RealInterval P1Norm = fourierNorm(fP1, rhoHatInterv);

    // Second term
    // Bring the rotated P2 to grid space. We already have P1
    ComplexMatrix rotP2;
    alloc(rotP2, m, m, N);
    matrixIFFT(rotP2, frotP2);

    // Calculate DF without the thickened torus
    DF(diff, K, kappa);

    // Perform the multiplication
    ComplexMatrix matMult;
    alloc(matMult, m, m, N);
    matMult = (rotP2 * diff) * P1;
    // Turn the product to Fourier space
    ComplexMatrix fmatMult;
    alloc(fmatMult, m, m, N);
    matrixFFT(fmatMult, matMult);

    // Substract lambda from the product of matrices, if performed on grid, use the constant lambda, otherwise transform lambda
    for (int k = 0; k < N; k++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                fmatMult[i][j][k] -= fLam[i][j][k];

    // Compute the Fourier norm of the result, since it's in Fourier space
    RealInterval matMultNorm = fourierNorm(fmatMult, rhoInterv);

    RealInterval redError;
    redError = (C_N * rotP2Norm * M0Norm * P1Norm) + matMultNorm;
    cout << "Reducibility error: ";
    print(redError);

    //====================================================================================

    // Invertibility Error

    // We already have the Fourier norm of P1, let's compute that of P2
    RealInterval P2Norm = fourierNorm(fP2, rhoHatInterv);

    // For the second summand we need to multiply P2 and P1 in grid space
    ComplexMatrix prodInverses;
    alloc(prodInverses, m, m, N);
    prodInverses = P2 * P1;
    // Take it to Fourier space
    ComplexMatrix fprodInverses;
    alloc(fprodInverses, m, m, N);
    matrixFFT(fprodInverses, prodInverses);
    // Substract the identity matrix in Fourier space
    for (int k = 0; k < N; k++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                fprodInverses[i][j][k] -= fId[i][j][k];

    // Compute Fourier norm
    RealInterval prodInvNorm = fourierNorm(fprodInverses, rhoInterv);

    // Add the two terms together
    RealInterval invertibilityError;
    invertibilityError = (C_N * P2Norm * P1Norm) + prodInvNorm;
    cout << "Invertibility error: ";
    print(invertibilityError);

    //====================================================================================
    // Now we can compute sigma

    // First we check that (redError + invertibilityError)/(1-lambda)<1
    cout << "Is the condition for hyperbolicity satisfied? ";
    bool hyperbolicity = ((redError + invertibilityError) / (1 - lamSt)) < 1;
    cout << (hyperbolicity ? "Yes" : "No") << endl;
    if (!hyperbolicity) {
        cout << "The initial torus is not hyperbolic, and thus cannot be validated." << endl;
        return 0;
    }

    RealInterval sigma;
    // Let's try compute norms of P1 and P2 using boxes method (norm in grid rather than Fourier)
    P1Norm = rhoNorm(P1, rho);
    P2Norm = rhoNorm(P2, rho);
    // Print the norms of P1 and P2
    cout << "Norm of P1: ";
    print(P1Norm);
    cout << "Norm of P2: ";
    print(P2Norm);
    sigma = P1Norm * one / (one - (RealInterval(lam) + redError + invertibilityError)) * P2Norm;
    cout << "Hyperbolicity bound sigma: ";
    print(sigma);

    //====================================================================================

    // In order to compute the norm of the second differential, since such supremum is taken for y and theta in the domain,
    // y is R close to the torus and theta is in the torus of thickness rho. This means we have to evaluate the second
    // differential over a thickened torus with radius R and a thickened grid of radius rho.

    // Let's thicken the grid and evaluate the torus over it. For that we use the original torus in Fourier space
    // and thicken it in the usual way with rho as complex width
    ComplexInterval phi;
    setBounds(phi, -1. / (2 * N), 1. / (2 * N), -rho, rho);
    ComplexVector thickFourier;
    alloc(thickFourier, m, N);
    fourierRot(thickFourier, fK, phi);

    // Now simply expand the torus on each side by R

    double R = 1.5e-2;
    ComplexInterval Rinterv;
    setBounds(Rinterv, -R, R, -R, R);

    thickFourier[0][0] += Rinterv;
    thickFourier[1][0] += Rinterv;

    // Notice that we could also have brought the thick torus back to grid and add the R box to each entry of the torus,
    // adding it to just the first value of the Fourier coefficients is the same

    // Bring it back to grid
    ComplexVector thickestK;
    alloc(thickestK, m, N);

    vectorIFFT(thickestK, thickFourier);
    // Since theta does not appear in the second differential, we do not need to thicken it
    RealInterval b;
    b = DF2_supNorm(thickestK, kappa);
    cout << "b bound: ";
    print(b);

    //====================================================================================
    // Put it all together to check the conditions fro validation

    auto rMinus = (one - sqrt(one - 2 * sigma * sigma * b * invarianceError)) / (sigma * b);
    auto rPlus = (one + sqrt(one - 2 * sigma * sigma * b * invarianceError)) / (sigma * b);
    /*cout << "Candidate r- = ";
    print(rMinus);
    cout << "Candidate r+ = ";
    print(rPlus);
    cout << "Condition 2 in r-: " << endl;
    print(sigma * b * rMinus);
    cout << "Condition 2 in r+: " << endl;
    print(sigma * b * rPlus);*/

    // Since in rMinus we have the first r to satisfy the first condition and in rPlus the last to satisfy it,
    // we can run the loop from rMinus to rPlus to find the r's that satisfy both conditions. For that we
    // can also base the step on the width of the interval rPlus - rMinus. 
    double rMin = -1, rMax = -1;
    mpfr_t rMi;
    mpfr_init2(rMi, prec);
    mpfi_get_left(rMi, rMinus.interval);
    mpfr_t rPl;
    mpfr_init2(rPl, prec);
    mpfi_get_right(rPl, rPlus.interval);
    double rM = mpfr_get_d(rMi, MPFR_RNDN), rP = mpfr_get_d(rPl, MPFR_RNDN);
    for (double r = rM; r <= rP && r <= R; r += abs(rP - rM) / 1.e6) {
        auto Condition1 = 1. / 2 * sigma * b * r * r - r + sigma * invarianceError <= 0;
        auto Condition2 = sigma * b * r < 1;

        if (Condition1 && Condition2) {
            if (rMin == -1)
                rMin = r;
            rMax = r;
        }
    }

    if (rMin != -1) {
        cout << "Torus is validated with" << endl;
        cout << "R: " << R << endl
             << "r-: " << rMin << endl
             << "r+: " << rMax << endl;
    } else
        cout << "The torus could not be validated." << endl;

    return 0;
}