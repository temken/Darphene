[![Build Status](https://github.com/temken/Darphene/workflows/Build%20Status/badge.svg)](https://github.com/temken/Darphene/actions)
[![codecov](https://codecov.io/gh/temken/Darphene/branch/main/graph/badge.svg)](https://codecov.io/gh/temken/template_cpp_cmake_obscura)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# Darphene - Modelling general dark matter - electron interactions in graphene targets

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7774374.svg)](https://doi.org/10.5281/zenodo.7774374)
[![arXiv](https://img.shields.io/badge/arXiv-2303.15497-B31B1B.svg)](https://arxiv.org/abs/2303.15497)


<img width="500" src="https://user-images.githubusercontent.com/29034913/227547647-9d2825dd-9cfd-4fde-9573-04d6c885a99f.png">

*Darphene* calculates predicted signal rates for dark matter detection experiments using graphene as target material using the [tight-binding approximation](https://en.wikipedia.org/wiki/Tight_binding). Details on the physics can be found in the corresponding [publication](https://arxiv.org/abs/2303.15497).

<img width="250" src="https://user-images.githubusercontent.com/29034913/226607786-762d3a59-88d8-4927-bd16-eb2a577391d2.png">

## General notes

- The material properties of the graphene target are modeled using the tight-binding approxmation, which uses atomic orbitals of carbon as a basis for the initial wavefunction of the electrons in graphene.
- The final electron states are approximated as plane-waves.
- *Darphene* is written in C++ and build using CMake.
- Parts of the calculations, in particular tabulation of event rates, are parallelized using MPI.

<details><summary>Repository content</summary>
<p>

The included folders are:

- *bin/*: This folder contains the executable after successful installation together with the configuration files.
- *data/*: Contains additional data necessary for the simulations, e.g. the solar model tables.
- *external/*: This folder will only be created and filled during the build with CMake and will contain the [obscura](https://github.com/temken/obscura) library necessary for all direct detection computations.
- *include/*: All header files of Darphene can be found here.
- *results/*: Each run of Darphene generates result files in a dedicated sub-folder named after the run's simulation ID string, which is specified in the configuration file.
- *src/*: Here you find the source code of Darphene.
- *tests/*: All code and executable files of the unit tests are stored here.

</p>
</details>


## Getting started

<details><summary>1. Dependencies</summary>
<p>

Before we can install *Darphene*, we need to make sure that a few dependencies are taken care of.

- [arb](https://arblib.org/): For the evaluation of hypergeometric functions.
- [boost](https://www.boost.org/): For numerical integration (used by *libphysica*).
- [CMake](https://cmake.org/): *Darphene* as well as the libraries *libphysica* and *obscura* are built with CMake.
- [eigen](https://eigen.tuxfamily.org): For the numerical procedure to find eigenvalues and eigenvectors of the generalized eigen problem.
- [libconfig](https://github.com/hyperrealm/libconfig): For the configuration files, *Darphene* uses the *libconfig* library (required version at least 1.6). This will be installed by *libphysica*, if it is not already installed.
- [libphysica](https://github.com/temken/libphysica): Automatically downloaded to */external/obscura/external/*, compiled, and linked by *CMake*.
- [MPI](https://www.mpi-forum.org/): The tabulation of DM observables is accelerated via parallelization using MPI.
- [obscura](https://github.com/temken/obscura): Automatically downloaded to */external/*, compiled, and linked by *CMake*.

On macOS, the dependencies can be installed using [homebrew](https://brew.sh/).

```
>brew install arb boost cmake eigen libconfig open-mpi
```

On Linux machines, we can use APT:

```
>sudo apt-get update -y
>sudo apt-get install -y arb libboost-all-dev libeigen3-dev libconfig-dev openmpi-bin openmpi-doc libopenmpi-dev
```

<details><summary>Notes on libconfig</summary>
<p>

The installation of `libconfig` is optional, since `libphysica` will install it automatically, if it is not available. Alternatively, it can be built from the source files via

```
>wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
>tar -xvzf libconfig-1.7.2.tar.gz
>pushd libconfig-1.7.2
>./configure
>make
>sudo make install
>popd
```
</p>
</details>

</p>
</details>

<details><summary>2. Download & Installation</summary>
<p>
The *Darphene* source code can be downloaded by cloning this git repository:

```
>git clone https://github.com/temken/Darphene.git 
>cd Darphene
```

The code is compiled and the executable is created using CMake.

```
>cmake -E make_directory build
>cd build
>cmake -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF ..
>cmake --build . --config Release
>cmake --install .
```

If everything worked well, there should be the executable *Darphene* in the */bin/* folder.

</p>
</details>

<details><summary>3. Usage</summary>
<p>
Once *Darphene* is installed, it can run by running the following command from the */bin/* folder:

```
>./Darphene Darphene.cfg
```

Alternative, one can use MPI to speed up calculations.

```
>mpirun -n N Darphene Darphene.cfg
```

where *N* is the number of desired MPI processes.

</p>
</details>

## Version History

- 27.03.2023: Release of version 0.1.0

## Everything else

<details><summary>Citing Darphene</summary>
<p>

If you decide to use this code, please cite the latest archived version,

> Emken, T., 2023, Darphene [Code, v0.1.0], [[DOI:10.5281/zenodo.7774374]](https://doi.org/10.5281/zenodo.7774374)

<details><summary>BibTeX entry</summary>
<p>

```
@software{Darphene,
  author = {Emken, Timon},
  title = {{Darphene [Code, v0.1.0]}},
  year         = {2023},
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {DOI:10.5281/zenodo.7774374},
  url          = {https://doi.org/10.5281/zenodo.7774374},
  howpublished={The code can be found under \url{https://github.com/temken/Darphene}. Version 0.1.0 is archived as \href{https://doi.org/10.5281/zenodo.7774375}{DOI:10.5281/zenodo.7774375}}
}
```
</p>
</details>

As well as the original publication,

> R. Catena, T. Emken, M. Matas, N.A. Spaldin, E. Urdshals , 2023,  **Direct searches for general dark matter-electron interactions with graphene detectors: Part I. Electronic structure calculations**, [[arXiv:2303.15497]](https://arxiv.org/abs/2303.15497).

</p>
</details>

<details><summary>Author & Contact</summary>
<p>

The author of *Darphene* is Timon Emken.

For questions, bug reports or other suggestions please open an [issue](https://github.com/temken/Darphene/issues).
</p>
</details>

<details><summary>License</summary>
<p>

This project is licensed under the MIT License - see the LICENSE file.

</p>
</details>