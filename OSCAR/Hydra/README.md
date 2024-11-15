# Hydra

- [Hydra.jl](./Hydra.jl) contains the [OSCAR](https://www.oscar-system.org/) implementation for [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9).
- [Hydra_polynomial_model.jl](./Hydra_polynomial_model.jl) contains the [OSCAR](https://www.oscar-system.org/) polynomial model for [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9).
- [Hydra_step_degrees_with_Gröbner_basis.ipynb](./Hydra_step_degrees_with_Gröbner_basis.ipynb), a notebook that first extracts a zero-dimensional DRL Gröbner basis from the Hydra polynomial via an affine change of coordinates, and then records the step degrees in F4 for low round numbers.
- [Hydra_step_degrees.jl](./Hydra_step_degrees.jl) is a script to redcord the step degrees in F4 for low round numbers with or without the Gröbner basis.
- [Hydra_characteristic_polynomial.jl](./Hydra_characteristic_polynomial.jl) is a script to record the running time of the computation of the characteristic polynomial for low round numbers.
- [Hydra_fglm.jl](./Hydra_fglm.jl) is a script to record the running time of FGLM for low round numbers.

## Requirements
- [OSCAR](https://www.oscar-system.org/) 1.2.0.

## Usage
For usage examples see [Hydra.ipynb](./Hydra.ipynb).
