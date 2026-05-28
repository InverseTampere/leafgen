# LeafGen: Foliage generation for QSMs and point clouds

## Description

LeafGen is a MATLAB package containing methods for generating virtual foliage on 3D tree models based on parametric leaf distributions. The tree models can be either quantitative structure models (QSM) or point clouds of trees. The user can choose parametric distribution functions on variables of tree structure for the leaf area density distribution (LADD), leaf orientation distribution (LOD), and leaf size distribution (LSD).

![Example illustration on QSM](readme-figure-qsm.png)

![Example illustration on canopy hull](readme-figure-ch.png)

## Documentation

The general LeafGen framework is presented in the following article

- Mönkkönen P, Van den Broeck W, Ali-Löytty S, Calders K, and Raumonen P (2025). *LeafGen*: Foliage Generation in 3D Tree Models. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.70080.

For technical details and guidance on usage, see the [LeafGen documentation].

## Quick start guide

To quickly test running the method:

1. Download or copy the LeafGen repository to your computer.
2. Open MATLAB and set the main path to the `src` directory of the LeafGen repository
3. Open and run one of the three main files in MATLAB (preferably start with either `main_qsm_direct.m` or `main_canopy_hull.m`, as the `main_leaf_cylinder_library.m` requires a long computation.)

To view a detailed tutorial on each of these files and more, see the  [LeafGen documentation].

[LeafGen documentation]: https://leafgen-docs.github.io/

## Licensing

### Code

The code in this repository is licensed under GPLv3. See `LICENSE` for more infromation.

### Example data

The files in `src/example-data/` include material derived from:

**Dataset:** Calders, K., Verbeeck, H., Burt, A., Origo, N., Nightingale, J., Malhi, Y., Wilkes, P., Raumonen, P., Bunce, R., & Disney, M. (2022). *Terrestrial laser scanning data Wytham Woods: individual trees and quantitative structure models (QSMs)* (Version 1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7307956

**Related article:** Calders, K., Verbeeck, V., Burt, A., et al. *Laser scanning reveals potential underestimation of biomass carbon in temperate forest.* Ecological Solutions and Evidence. https://doi.org/10.1002/2688-8319.12197

**License:** Creative Commons Attribution 4.0 International (CC BY 4.0). https://creativecommons.org/licenses/by/4.0/

These files remain licensed under CC BY 4.0.

