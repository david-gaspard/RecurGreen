# RecurGreen

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
	- [Physical model](#physical-model)
	- [Purposes](#purposes)
* [USAGE AND OPTIONS](#usage-and-options)
	- [`recurgreen` executable](#recurgreen-executable)
	- [`waveguide` namelist](#waveguide-namelist)
	- [Tasks](#tasks)
	- [`plothisto` executable](#plothisto-executable)
* [REFERENCES](#references)

## PRESENTATION

RecurGreen is a Fortran 2008 program to compute the Green function $G^+(x,y;x',y')$ corresponding to the propagation of a scalar wave in a 2D disordered waveguide starting from some position $(x',y')$.
The algorithm used is a variant of the _recursive Green method_ [[1]](#1), hence the name of the program, expressed in the basis of the transverse eigenmodes of the waveguide.
The focus of this program is the computation of the transmission matrix and the distribution of transmission eigenvalues when averaging over the realizations of the disorder.

### Physical model

The Green function equation has the form

<p>$$ \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + k^2 + \mathrm{i}\varepsilon - U(x,y) \right) \psi(x,y) = \delta(x-x')\delta(y-y') $$</p>

where $\psi(x,y)$ is the sought wavefunction, $k$ is the central wavenumber, $\varepsilon$ a complex shift that may be used to include absorption, and $U(x,y)$ is a random Gaussian potential whose probability distribution is given by

<p>$$ \mathcal{P}[U(x,y)] \propto \exp\left( -\frac{1}{2\alpha} \int_0^L dx \int_0^W dy\: U(x,y)^2 \right) $$</p>

$L$ being the length of the disordered region, $W$ its width, and $\alpha$ is the disorder strength, which is related to the scattering mean free path by $\ell_{\rm s}=k/(\pi\nu\alpha)$ under the independent scattering approximation.
This potential has zero correlation length, so that the independent scatterings are isotropic and the transport and scattering mean free path are equal.
The wave equation is discretized along the longitudinal coordinate $x$ and expressed on the basis of the transverse modes of the waveguide instead of a spatial discretization along the $y$ coordinate.

### Purposes

The main focus of the program is the computation of the _transmission matrix_ $\mathsf{t}$ according to the Fisher-Lee relation [[2]](#2)

<p>$$ t_{ij} = 2\mathrm{i} \sqrt{k_{x,i} k_{x,j}} G^+_{ij}(L;0) $$</p>

where $G^+_{ij}(x;x')$ is the projection of the Green function on the basis of waveguide modes $\chi_i(y)$,

<p>$$ G^+_{ij}(x;x') = \int_0^W dy \int_0^W dy'\: \chi^*_i(y) G^+(x,y:x',y') \chi_j(y') $$</p>

and $k_{x,i}=\sqrt{k^2 - p_{y,i}^2}$ is the longitudinal component of the wavevector, and $p_{y,i}$ is the transverse component.

The program computes in particular the _distribution of transmission eigenvalues_, $\rho(T)$, which are the eigenvalues of $\mathsf{t}^\dagger \mathsf{t}$.
This distribution is defined by

<p>$$ \rho(T) = \frac{1}{N} \sum_{n=1}^N \langle\delta(T - T_n)\rangle = \frac{1}{N} \mathrm{Tr}\langle\delta(T - \mathsf{t}^\dagger \mathsf{t})\rangle $$</p>

where $\langle\cdot\rangle$ denotes the average over the disorder, and $N$ is the number of transmission eigenvalues.

## USAGE AND OPTIONS

### `recurgreen` executable

The syntax of the executable `recurgreen` reads:
```
recurgreen settings.nml
```
This program reads and executes the instructions given by the configuration file `settings.nml`.

This file must contain at least three [Fortran namelists](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/namelist.html):

* `settings`: The only entry of this namelist is `task` which determines the quantity to compute. See [Tasks](#tasks) for the list of available tasks.
* `waveguide`: This namelist contains the parameters of the disordered waveguide. See [waveguide](#waveguide-namelist) for the list of parameters.
* A namelist which depends on the task prescribed in the `settings` namelist. See [Tasks](#tasks) for the list of available tasks.

### `waveguide` namelist

This namelist contains several parameters which characterize the disordered waveguide:

* `boundtype`: Type of boundary conditions in the transverse direction. Either `periodic` or `dirichlet`.
* `wol`: Width-to-wavelength ratio $W/\lambda$ (periodic boundary conditions). Integrer values are not valid because the highest mode has zero momentum ($k_x=0$), causing divisions by zeros. A typical value is `wol=50.5`. 
* `exfac`: Multiplication factor for the number of modes to improve the disorder quality in the transverse direction (>1 includes evanescent modes).
* `aratio`: Aspect ratio of the disordered region, $L/W$. A typical value is `aratio=1`.
* `dscat`: Scattering thickness of the waveguide, i.e., length over the mean free path, $L/\ell_{\rm s}$.
           Note that it is only the _target_ value of the scattering thickness assuming the independent-scattering approximation.
           The real value of the scattering thickness may be slightly larger for $L/\ell_{\rm s}\gg 1$.
* `dabso`: Absorption thickness of the waveguide, i.e., length over the absorption mean free path, $L/\ell_{\rm a}$.
           Absorption is realized by putting $\varepsilon=k/\ell_{\rm a}$ in the Green function equation above. Absorption is thus uniformly distributed in the medium.
* `nxdiso`: Number of "x" points in the disordered region used in the discretization. Recommended value is `nxdiso = 2*wol*exfac*aratio` or larger.
* `nxfree`: Number of "x" points in the free region on one side (same on other side). This parameter should only be used when plotting the wavefunction and is useless for other tasks. Recommended value is `nxfree=0`.
* `napera`: Fraction of excited modes at the input lead (port "a"). Also numerical aperature or filtering parameter. It must be between 0 and 1.
* `naperb`: Fraction of observed modes at output lead (port "b"). Also numerical aperature or filtering parameter. It must be between 0 and 1.

### Tasks

The `task` entry of the namelist `settings` may contain any of the following instructions and must be accompanied with a namelist of the same title.

#### `tdistrib`

Computes the distribution of transmission eigenvalues $\rho(T)$ defined by 

$$ \rho(T) = \frac{1}{N} \sum_{n=1}^N \langle\delta(T - T_n)\rangle $$

The corresponding namelist must contain the following entries:
* `nseed`: Total number of realizations of the disorder (seeds are used as input of the random number generator). Recommended value is `nseed=50000` to reduce the noise but it of course depends on the number of bins.
* `nbins`: Number of bins of the histogram. Recommended value is `nbins=256` to resolve the structures of the distribution.
* `nthreads`: Number of threads used to compute the transmission-eigenvalue distribution in parallel. Recommended value is the number of CPU cores.
* `save_eigvals`: If true, saves the raw transmission eigenvalues in a file. Note that this may generate big files (tens of MB).
    This option can be used to replot the distribution with other binning parameters (see also [plothisto](#plothisto-executable)).

#### `wavefunction`

Computes the wavefunction for a given input plane wave in a given realization of the disorder.
The parameters of the corresponding namelist are:
* `iseed`: Seed index of the realization of the disorder.
* `ipos`: Index of the initial/incident position. It must be in 1..nx. `ipos=1` corresponds to the leftmost edge.
* `imode`: Index of the initial/incident mode. It must be in 1..nmu. `imode=1` always corresponds to the frontal mode for both `boundtype=periodic` and `boundtype=dirichlet`.

#### `tprofile`

Computes the disorder-averaged profile of transmission eigenstates $\langle|\psi_T(x)|^2\rangle$ where 

<p>$$ \psi_T(x,y) = 2\mathrm{i}\sqrt{k} \sum_{m,m'} \chi_m(y) G^+_{mm'}(x;x_0) \sqrt{k_{x,m'}} u_{T,m'} $$</p>

and $u_T$ is the transmission eigenvector (eigenvector of $\mathsf{t}^\dagger\mathsf{t}$) normalized by $|u_T|^2=1$.
$\psi_T(x,y)$ is the transmission eigenfunction associated to transmission eigenvalue $T$.
This quantity is not only averaged in the transverse direction ("y" direction), but also over all the eigenvalues $T$ located in the search interval $[T_{\min}, T_{\max}]$.
The parameters of the corresponding namelist are:
* `tmin`: Minimum transmission eigenvalue of the search interval $[T_{\min}, T_{\max}]$.
* `tmax`: Maximum transmission eigenvalue of the search interval $[T_{\min}, T_{\max}]$.
* `nseed`: Number of realizations of the disorder used to compute $\langle|\psi_T(x)|^2\rangle$. Recommended value is `nseed=10000`.
* `nthreads`: Number of threads used to compute the transmission-eigenvalue distribution in parallel. Recommended value is number of CPU cores.

#### `tprofile_group`

Same as `tprofile` but computes 4 transmission eigenstate profiles at the same time, which saves a factor 4 of computation time.
The value 4 is chosen for convenience, it can be easily tweaked by modifying the line `nprofile=4` in the source code.
The parameters of the corresponding namelist are:
* `tm`: List of central values for the transmission eigenvalue. Recommended argument is `0.998, 0.50, 0.100, 0.0010` as it produces best results.
* `dt`: Corresponding list of half-width for the transmission eigenvalue. The intervals are thus: $[T_{\min}, T_{\max}] = [T_{\rm m}-\delta t, T_{\rm m}+\delta t]$. Recommended argument is `0.001, 0.01, 0.002, 0.0001` as it produces best results.
* `nseed`: Number of realizations of the disorder used to compute $\langle|\psi_T(x)|^2\rangle$. Recommended value is `nseed=10000`.
* `nthreads`: Number of threads used to compute the transmission-eigenvalue distribution in parallel. Recommended value is number of CPU cores.

#### `tpattern`

Computes the radiance within the disordered medium, i.e., the intensity in each direction, given by the square modulus of the deposition matrix, $|Z_{ij}|^2$, where $Z_{ij}=2\mathrm{i}\sqrt{k_{x,i}k_{x,j}}G^+_{ij}(x;x')$.
This task can be used to test the validity of the _isotropy assumption_ (diffusion approximation) used in the DMPK theory, i.e., the relative significance of the radiance directionality.
The parameters of the corresponding namelist are:
* `nseed`: Number of realizations of the disorder used to average the radiance.
* `imode`: Index of the initial/incident mode. It must be in 1..nmu. `imode=1` always corresponds to the frontal mode for both `boundtype=periodic` and `boundtype=dirichlet`.
* `nthreads`: Number of threads used to compute the radiation pattern in parallel. Recommended is number of CPU cores.
* `probatype`: Specifies the observable more precisely.
    - `incident`: Compute the transmission probability $|t_{ij}|^2$ as a function of $i$ for a fixed incident channel $j$.
    - `total`: Compute the total transmission probability $\sum_{i} |t_{ij}|^2$ of channel $j$. The index `imode` of the incident wave is then ignored.
    - `deposition`: Compute the deposition probability at the center of the medium, i.e., $|Z_{ij}|^2$ as a function of $i$ for a fixed incident channel $j$.

### `plothisto` executable

The syntax of the executable `plothisto` reads:
```
plothisto plothisto_settings.nml
```
This program reads and executes the instructions given by the configuration file `plothisto_settings.nml`.
Its purpose is to replot the transmission eigenvalue distribution computed by the task `tdistrib` with different binning parameters.
The settings file `plothisto_settings.nml` contains a single Fortran namelist with the following entries:
* `bintype`: Type of lattice used for the bins. Either `linear` or `chebyshev`. 
* `nbins`: Number of desired bins for the distribution.
* `xmin`: Minimum abscissa of the distribution.
* `xmax`: Maximum abscissa of the distribution.

## REFERENCES

<a id="1">[1]</a>
Baranger, Harold U. and DiVincenzo, David P. and Jalabert, Rodolfo A. and Stone, A. Douglas,
_Classical and quantum ballistic-transport anomalies in microjunctions_,
[Phys. Rev. B **44**, 10637-10675 (1991)](https://doi.org/10.1103/PhysRevB.44.10637).

<a id="2">[2]</a>
Fisher, Daniel S. and Lee, Patrick A.,
_Relation between conductivity and transmission matrix_,
[Phys. Rev. B **23**, 6851-6854 (1981)](https://doi.org/10.1103/PhysRevB.23.6851).
