# ConvexHullPricing

## Overview
The `ConvexHullPricing.jl` package implements various primal and dual optimization methods to compute Convex Hull Prices, relying on a structure model to load thermal generators and unit commitment instances. It is primarily designed to benchmark optimization methods for this problem but one can also compute various metrics that are useful to practitioners.

## Current version solvers
The current version is still in the making. The proposed solvers are:

1. Dual methods
   - Bundle Level Method
    > N. Stevens and A. Papavasiliou, “Application of the level method for computing locational convex hull prices,” IEEE Transactions on Power Systems, (2022)
    > C. Lemaréchal, A. Nemirovskii, and Y. Nesterov, New variants of bundle methods, Math.Program., (1995)
   - D-Adaptation
    > A. Defazio, K. Mishchenko, Learning-rate-free learning by D-adaptation. arXiv preprint arXiv:2301.07733, (2023)
   - DowG
    > A. Khaled, K. Mishchenko, C. Jin, DoWG Unleashed: An Efficient Universal Parameter-Free Gradient Descent Method, arXiv preprint arXiv:2305.16284, (2023)
   - Polyak Subgradient Method
    > B. Polyak, Gradient methods for the minimisation of functionals. USSR Computational Mathematics and Mathematical Physics, (1963)
   - Subgradient Method
    > C. Wang, T. Peng, P. B. Luh, P. Gribik, and L. Zhang, The subgradient simplex cutting plane method for extended locational marginal prices, IEEE Transactions on Power Systems, (2013)
    > S. Boyd, L. Xiao, and A. Mutapcic, “Subgradient methods,” lecture notes of EE392o, Stanford University, Autumn Quarter, (2003)

2. Primal methods
   - Column Generation (Dantzig-Wolfe decomposition)
    > P. Andrianesis, D. Bertsimas, M. C. Caramanis, and W. W. Hogan, Computation of convex hull prices in electricity markets with non- convexities using Dantzig-Wolfe decomposition, IEEE Transactions on Power Systems, (2021)
   - Column-and-Row Generation (Dantzig-Wolfe decomposition on the extended formulation of the primal)
    > M. Paquet, Comparison between row generation, column generation and column-and-row generation for computing convex hull prices in day-ahead electricity markets. Ecole polytechnique de Louvain, Université catholique de Louvain, (2021)
   - Row Generation (Benders decomposition)
    > B. Knueven, J. Ostrowski, A. Castillo, and J.-P. Watson, A computationally efficient algorithm for computing convex hull prices, Computers & Industrial Engineering, (2022)

## Installation

Download Julia, and follow the instructions described [here](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-unregistered-packages).

## Citing this package

If you use `ConvexHullPricing.jl` for published work, we encourage you to cite the software using the following Bibtex citation:
```bibtex

```
## License
`ConvexHullPricing.jl` is released under a MIT License. `ConvexHullPricing.jl` has been developed as part of the ITN-ETN project TraDE-OPT funded by the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 861137.

## Get in touch
Comments and suggestions are more than welcome, get in touch via [mail](mailto:sofiane.tanji@uclouvain.be) or an issue!
