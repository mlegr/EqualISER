# Equal[iseR]

## Introduction
Equal[iseR] is a frequency-domain solver for non-smooth systems of ODE/DAE targetting the vibration analysis of structural mechanical systems subject to unilateral and frictional contact occurrences. It is an adaptation of the well-known Harmonic Balance Method. It relies on [Equal]ity-based versions of the governing equations solved in a weighted [Resi]dual sense.

Unilateral contact and friction commonly take the form of differential inclusions which can be recast into more conventional equality-based identities through the use of non-smooth projections. Such equations are then approximately transformed using classical weighted residual techniques. The resulting system of nonsmooth/nonlinear equations is solved using the classical Hybrid Powell solver.

## Supporting paper
The paper [A compact, equality-based weighted residual formulation for periodic solutions of systems undergoing frictional occurrences](https://doi.org/10.25518/2684-6500.190) (see also the version available on HAL: [hal-04189699](https://hal.science/hal-04189699)) introduces the methodology and attentant equations.

## History
2024.03.27 first public release (version 1.0)

### Version 1.0
The project comes with two versions of the solution procedure

* `Matlab` Three versions are available (Tested on Matlab R2021b):
    * __ImplicitEuler__ Time-marching solution method based on the first-order Implicit Euler scheme for comparison purposes with the frequency-domain solution strategies.
    * __WeightedResiduals__ Frequency-domain solution strategy based on the proposed weighted-residual formulation. Integrals are computed using classical quadrature scheme and the nonlinear equations are solved using the class Hybrid Powell (or dog-leg) nonlinear solver. Used in the initial version of the above paper.
    * __FFT__ Same as __WeightedResiduals__ but the integrals are computed using the FFT algorithm. This version of the solver is the fastest.
* `Python` One version is available (Tested on Python 3.8 and 3.11):
   * __FFT__ Essentially the same as the Matlab __FFT__ even though the iFFT operation involved in the algorithm is computed differently. 



## Authors
They are listed in the licence.
