!!**************************************************************************************************
!! Created on 2024-05-15 at 09:47:43 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran program to test the implementation of the random normal generator. 
!!**************************************************************************************************
program random_normal_test
    use random_normal
	implicit none
	
    real(wp), allocatable :: array(:)
    integer :: nn, seed
    
    nn = 1000000  !! Number of random numbers used in the tests.
    seed = 1      !! Seed value to ensure reproducibility.
    
    allocate(array(nn))  !! Allocate the required space for the sequence of random numbers.
	call rand_normal(array, seed)  !! Generate a random sequence (normally distributed).
    
    !! Run the test over the array:
    !!call test_cshift()
	call test_normal_moments()
    call test_normal_correlation()
	call test_randi_ansic()
    
	contains
	
	!!**********************************************************************************************
	!! Subroutine to test the first moments of the normal random generator.
	!!**********************************************************************************************
	subroutine test_normal_moments()
		real(wp) :: mean, sigma, asymm, kurto, sampledev   !! Moments of the Gaussian distribution.
        
        print '(a)', "====== TEST NORMAL MOMENTS ======"
        
        sampledev = 1./sqrt(real(nn))  !! Sample deviation.
        mean = sum(array)/size(array)
        sigma = sqrt(sum(array*array)/size(array) - mean*mean)
        asymm = sum(array**3)/size(array)
        kurto = sum(array**4)/size(array)
        
        print '(a,g0.6,a)', "Sample dev = ", sampledev, " (fluctuations allowed)"
        print '(a,g0.6,a)', "Mean  = ", mean , " (should be 0)"
        print '(a,g0.6,a)', "Sigma = ", sigma, " (should be 1)"
        print '(a,g0.6,a)', "Asymm = ", asymm, " (should be 0)"
        print '(a,g0.6,a)', "Kurto = ", kurto, " (should be 3)"
        
	end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test Fortran's cshift() subroutine.
	!!**********************************************************************************************
    subroutine test_cshift()
        integer :: a(5)
        
        print '(a)', "====== TEST CSHIFT ======"
        
        a = [ 1, 2, 3, 4, 5 ]
        print '(a,"[",5(1x,i0,1x),"]")', "a  = ", a
        print '(a,"[",5(1x,i0,1x),"]")', "a' = ", cshift(a, 1)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the correlations of the normal random generator.
	!!**********************************************************************************************
	subroutine test_normal_correlation()
		real(wp) :: corr   !! Correlations in the normal rnadom sequence.
        integer, parameter :: ncorr = 5  !! Number of computed correlations.
        integer :: i
        
        print '(a)', "====== TEST NORMAL CORRELATION ======"
        
        do i = 1, ncorr
            corr = sum(cshift(array, i)*array)/size(array)
            print '(a,i0,a,g0.6,a)', "Corr(", i, ") = ", corr , " (should be 0)"
        end do
        
	end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the uniformity of the basic seed generator iurand_lcg().
	!!**********************************************************************************************
	subroutine test_randi_ansic()
        integer, allocatable :: intarray(:)   !! Seed array.
        integer, allocatable :: histo(:)
        integer, parameter :: ns = 100000   !! Number of samples.
        integer, parameter :: nbin = 10    !! Number of bins for the histogram.
        integer :: i, b, seed
        real(wp) :: x
        character(len=30) :: style
        
        print '(a)', "====== TEST RANDI_ANSIC ======"
        
        !! 1. Prepare the sequence:
        allocate(intarray(ns))
        seed = 1
        call randi_ansic(intarray, seed)
        
        !! 2. Prepare the histogram:
        allocate(histo(nbin))
        histo = 0
        do i = 1, ns
            x = real(intarray(i))/real(huge(nn))  !! This is a number in the interval [-1, 1].
            b = ceiling(nbin*(x + 1)/2)
            !!print '(i0)', b
            histo(b) = histo(b) + 1
        end do
        
        write (style, '(a,i0,a)') "(a,", nbin, "(1x,i0,1x),a)"
        print style, "Histo = [", histo, "]"
        
	end subroutine
	
end program
