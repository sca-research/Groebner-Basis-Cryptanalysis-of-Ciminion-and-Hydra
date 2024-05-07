# Gröbner Basis Cryptanalysis Ciminion and Hydra
This repository contains implementations of polynomials models for [Ciminion](https://doi.org/10.1007/978-3-030-77886-6_1) and [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9) analyzed in [this paper](dummy).

- For [Ciminion](https://doi.org/10.1007/978-3-030-77886-6_1):
    - The [SageMath](https://www.sagemath.org/) implementation contains a function that computes a DRL Gröbner basis over arbitrary fields and for arbitrary number of rounds.
- For [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9):
    - The [SageMath](https://www.sagemath.org/) implementation contains a function that verifies that a [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9) instance is in generic coordinates for arbitrary number of rounds and arbitrary number of samples.
    - With the [OSCAR](https://www.oscar-system.org/) implementation one can compute the step degrees in F4 for Gröbner basis computations of [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9).
    - With both implementations one can extract a zero-dimensional DRL Gröbner basis after an affine change of coordinates.

## Requirements
- [SageMath](https://www.sagemath.org/) 10.3.
- [OSCAR](https://www.oscar-system.org/) 1.0.2.
