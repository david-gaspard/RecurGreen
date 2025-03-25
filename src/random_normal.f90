!!**************************************************************************************************
!! Created on 2024-05-14 at 17:34:31 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module to generate normally distributed random numbers.
!!**************************************************************************************************
module random_normal
    use base_utils
	implicit none
    
	contains
	!!**********************************************************************************************
    !! Subroutine to fill the given array with a pseudorandom sequence of 32-bits integers of
    !! uniform distribution. This subroutine is based on the simple linear congruential generator (LCG)
    !! used by ANSI C. This quick-and-dirty generator is NOT RECOMMENDED for use in scientific computations
    !! except for generating the seed of Fortran's random_seed(put=seedarr) subroutine because it is not
    !! related to the sequence itself.
    !! See also: https://en.wikipedia.org/wiki/Linear_congruential_generator
    !! See also: https://citeseer.ist.psu.edu/viewdoc/download?doi=10.1.1.53.3686&rep=rep1&type=pdf
	!!**********************************************************************************************
    subroutine randi_ansic(intarray, seed)
        integer, intent(out) :: intarray(:)
        integer, intent(in) :: seed
        integer :: i, sample
        
        sample = seed
        do i = 1, size(intarray)
            sample = 1103515245*sample + 12345   !! Uniform deviate sequence for signed 32-bit integers (modulo 2^31).
            intarray(i) = sample
        end do
        
    end subroutine
	
	!!**********************************************************************************************
	!! Subroutine to fill the given array with pseudorandom normal numbers with zero mean and unit variance (mu=0, sigma=1).
    !! This subroutine uses the Box-Muller method to generate the normal variables.
    !! array = Output array with normal random numbers (real array).
    !! seed  = Seed index (32-bit integer). The number of different random sequences produced by
    !!         this subroutine is thus limited to 2^32.
	!!**********************************************************************************************
	subroutine rand_normal(array, seed)
        real(wp), intent(out) :: array(:)
		integer, intent(in) :: seed
        integer :: i, n  !! n = Desired number of normal random numbers.
        integer :: nseedarr  !! Number of actual seeds used by the random_number() subroutine.
        integer, allocatable :: seedarr(:)  !! Actual seeds used by the random_number() subroutine.
        real(wp), allocatable :: uniforms(:)   !! Uniform numbers used to compute the normal numbers.
        !!character(len=30) :: style
        
        !! 1. Prepare the seeds:
        call random_seed(size=nseedarr)  !! Get the expected length of the seed (8 * 32-bit integers = 256 bits for xoshiro256).
        allocate(seedarr(nseedarr))
        call randi_ansic(seedarr, seed) !! Initialize the seed array using the LCG of ANSI C to get maximum entropy.
        !!write (style, '(a,i0,a)') "(a,", nseedarr, "(i0,1x))"
        !!print style, "seedarr = ", seedarr
        
        !! 2. Generate the uniform numbers:
        n = size(array)   !! Now 'n' is the desired number of normal random numbers.
        allocate(uniforms(2*n))  !! Each normal random number needs two uniform random numbers with the Box-Muller method.
        call random_seed(put=seedarr)  !! Initialize the seed of the random generator.
        call random_number(uniforms)   !! Generate the sequence of uniform random numbers (xoshiro256 for GNU Fortran).
        
        !! 3. Use the Box-Muller transformation:
        do i = 1, n
            array(i) = sqrt(-2.0_wp*log(uniforms(2*i-1)))*cos(2.0_wp*pi*uniforms(2*i))
        end do
        
	end subroutine
	
end module
