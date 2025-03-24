# RecurGreen

* [PRESENTATION](#presentation)
* [INSTALLATION](#installation)
	- [Dependencies](#dependencies)
	- [Installation guide](#installation-guide)
* [USAGE AND OPTIONS](#usage-and-options)
	- [recurgreen](#recurgreen)
	- [plothisto](#plothisto)
    - [integrate_rho.py](#integrate_rho.py)
* [CONFIGURATION FILE](#configuration-file)
	
* [EXAMPLE](#example)

## PRESENTATION

RecurGreen is a Fortran program to solve the wave equation in a 2D disordered waveguide using the recursive Green method.
The wave equation has the form 

<p>$$ \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + k^2 - U(x,y) \right) \psi(x,y) = 0 $$</p>

where $\psi(x,y)$ is the sought wavefunction, $k$ is the central wavenumber, and $U(x,y)$ is a random Gaussian potential whose probability distribution is given by

<p>$$ \mathcal{P}[U(x,y)] \propto \exp( -\frac{1}{2\alpha} \int_0^L \mathrm{d}{x} \int_0^W \mathrm{d}{y}\: U(x,y)^2 ) $$</p>

$L$ being the length of the disordered region, and $W$ its width.
The wave equation is discretized along the longitudinal coordinate $x$ and expressed on the basis of the transverse modes of the waveguide instead of a spatial discretization along the $y$ coordinate.


## INSTALLATION

### Dependencies

### Installation guide

## USAGE AND OPTIONS

### recurgreen

### plothisto

### integrate_rho.py

## CONFIGURATION FILE

## EXAMPLE
