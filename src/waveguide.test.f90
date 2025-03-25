!!**************************************************************************************************
!! Created on 2024-05-17 at 17:41:04 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran program to test the waveguide type and in particular the recursive Green method.
!!**************************************************************************************************
program waveguide_test
    use waveguide
	implicit none
	
    call test_init_waveguide()
    call test_unitarity()
    call test_green_v1_v2()
    
	contains
    !!**********************************************************************************************
    !! Test the waveguide initialization with the given parameters.
    !!**********************************************************************************************
    subroutine test_init_waveguide()
        type(waveguide_t) :: wguide
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb
        integer :: nxdiso, nxfree, nseed, nbins, nthreads
        
        print '(a)', "====== TEST INIT_WAVEGUIDE ======"
        
        !! Physical settings:
        boundtype = "dirichlet"  !! Types of boundary conditions in the transverse direction. Either "periodic" or "dirichlet".
        wol = 4.30_wp    !! Width-to-wavelength ratio W/lambda. Typically: wol={25.70, 25.98, 26.02}.
        aratio = 1.0     !! Aspect ratio of the disordered region, L/W.
        exfac = 1.2_wp   !! Multiplication factor for the number of modes taken to improve the representation of the disorder in the transverse direction (this includes evanescent modes).
        dscat = 0.2_wp   !! Optical depth of scattering of the waveguide (length over mean free path, L/lscat).
        dabso = 0.0      !! Optical depth of absorption of the waveguide (length over mean free path, L/labso).
        napera = 1.0     !! Numerical aperture at the input lead.
        naperb = 1.0     !! Numerical aperture at the output lead.
        
        !! Numerical settings:
        nxdiso = 10     !! Number of slices in the "x" direction in the disordered region.
        nxfree = 2     !! Number of slices in the "x" direction in the free region.
        nseed = 5000   !! Total number of seeds which is also the total number of random configurations (Computation speed for wol=25.50 is about 250 seeds/min).
        nbins = 256    !! Number of bins of the histogram (typically nbins=100, recommended is nbins=256 to match the settings of the Eilenberger solver).
        nthreads = 8   !! Number of threads used to compute the transmission-eigenvalue distribution in parallel.
        
        call wguide%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide%print_parameters(stdout, "|")
        call wguide%del()
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the unitarity of the reflection and transmission matrices using the recursive Green method
    !! for several random configurations of the disorder.
    !! This uses version v1 of the recursive Green algorithm.
    !!**********************************************************************************************
    subroutine test_unitarity()
        type(waveguide_t) :: wguide
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb
        integer :: nxdiso, nxfree, iseed, nseed
        
        print '(a)', "====== TEST UNITARITY V1 ======"
        
        boundtype = "dirichlet"  !! Types of boundary conditions in the transverse direction. Either "periodic" or "dirichlet".
        wol = 4.30_wp    !! Width-to-wavelength ratio W/lambda. Typically: wol={25.70, 25.98, 26.02}.
        aratio = 1.0     !! Aspect ratio of the disordered region, L/W.
        exfac = 1.2_wp   !! Multiplication factor for the number of modes taken to improve the representation of the disorder in the transverse direction (this includes evanescent modes).
        dscat = 0.2_wp   !! Optical depth of scattering of the waveguide (length over mean free path, L/lscat).
        dabso = 0.0      !! Optical depth of absorption of the waveguide (length over mean free path, L/labso).
        nxdiso = 100     !! Number of slices in the "x" direction in the disordered region.
        nxfree = 20      !! Number of slices in the "x" direction in the free region.
        nseed = 5        !! Total number of seeds which is also the total number of random configurations (Computation speed for wol=25.50 is about 250 seeds/min).
        napera = 1.0     !! Numerical aperture at the input lead.
        naperb = 1.0     !! Numerical aperture at the output lead.
        
        call wguide%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide%print_parameters(stdout, "|")
        
        do iseed = 1, nseed
            call wguide%test_unitarity(iseed, 1)
        end do
        
        print '(a)', "====== TEST UNITARITY V2 ======"
        
        do iseed = 1, nseed
            call wguide%test_unitarity(iseed, 2)
        end do
        
        call wguide%del()
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compare the reflection and transmission matrices predicted by the two version of the
    !! recursive Green algorithm.
    !!**********************************************************************************************
    subroutine test_green_v1_v2()
        type(waveguide_t) :: wguide
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb
        integer :: nxdiso, nxfree, iseed, nseed
        
        print '(a)', "====== TEST GREEN V1/V2 ======"
        
        boundtype = "dirichlet"  !! Types of boundary conditions in the transverse direction. Either "periodic" or "dirichlet".
        wol = 4.30_wp    !! Width-to-wavelength ratio W/lambda. Typically: wol={25.70, 25.98, 26.02}.
        aratio = 1.0     !! Aspect ratio of the disordered region, L/W.
        exfac = 1.2_wp   !! Multiplication factor for the number of modes taken to improve the representation of the disorder in the transverse direction (this includes evanescent modes).
        dscat = 0.2_wp   !! Optical depth of scattering of the waveguide (length over mean free path, L/lscat).
        dabso = 0.0      !! Optical depth of absorption of the waveguide (length over mean free path, L/labso).
        nxdiso = 100     !! Number of slices in the "x" direction in the disordered region.
        nxfree = 20      !! Number of slices in the "x" direction in the free region.
        nseed = 5        !! Total number of seeds which is also the total number of random configurations (Computation speed for wol=25.50 is about 250 seeds/min).
        napera = 1.0     !! Numerical aperture at the input lead.
        naperb = 1.0     !! Numerical aperture at the output lead.
        
        call wguide%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide%print_parameters(stdout, "|")
        do iseed = 1, nseed
            call wguide%test_green_v1_v2(iseed)
        end do
        call wguide%del()
        
    end subroutine
    
end program
