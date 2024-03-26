# Fourier Methods for CAPS of Invariant Tori

This project contains the code to execute the validation of a fiberwise hyperbolic invariant torus in a skew-product dynamical system using Fourier methods. There is a validation.cpp which executes the validation and allows for a choice of input values. There is also a test.cpp that uses googletests to make sure all the implementarions are correct. The code uses the interval arithmetics package MPFI for the accurate results. In such a way two classes have been created, ComplexInterval for the easy manipulation of complex intervals (that one can think as boxes) and RealInterval, that allows to work with single intervals. Classic algorithms for many procedures have been adapted to these classes, such as the all important FFT using a Radix-2 DIT, needing then inputs over powers-of-two-sized grids.

- The IntervalClasses.h header file contains the implementations and methods for such classes, as well as other useful functions that require both, like elementary operations.

- The Functions.h header file includes all the necessary functions for the validation to work, such as computations of iterations of maps, jacobians, FFTs and IFFTs, norms, etc.

The validation also creates output files with the Fourier coefficients of certain objects in case one would like to study them by plotting or fitting. That's what can be found in the Output folder.
