!!**************************************************************************************************
!! Created on 2024-05-15 at 16:16:34 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran module to define the disordered waveguide type. We assume that the waveguide is two-dimensional
!! with periodic boundary conditions in the transverse direction.
!! Regarding the implementation of the recursive Green method, see Appendix G of David Gaspard's notes
!! on the Nazarov technique.
!!**************************************************************************************************
module waveguide
    use random_normal
	implicit none
    
    type :: waveguide_t
		private
            !! Physical settings:
            integer :: nxdiso   !! Number of points of the disordered region in the longitudinal "x" direction. Recommended is nxdiso = 2*wol*exfac*aratio (or larger).
            integer :: nxfree   !! Number of points of the free region on one side (the other side also has "nxfree" too). Recommended is 0 for intensive simulations.
            integer :: nx       !! Total number of points in the "x" direction used in the discretization (nx = nxdiso + 2*nxfree).
            real(wp) :: wol     !! Width-to-wavelength ratio W/lambda. Typically: wol={25.50, 25.98, 26.02}.
            real(wp) :: exfac   !! Multiplication factor for the number of modes taken to improve the representation of the disorder in the transverse direction (this includes evanescent modes).
            real(wp) :: aratio  !! Aspect ratio of the disordered region, L/W. This parameter should not affect the results significantly.
            real(wp) :: dscat   !! Optical depth of scattering of the waveguide (length over mean free path, L/lscat).
            real(wp) :: dabso   !! Optical depth of absorption of the waveguide (length over mean free path, L/labso).
            real(wp) :: napera  !! Fraction of controlled modes at port "a" (left-hand input port).
            real(wp) :: naperb  !! Fraction of observed modes at port "b" (right-hand output port).
            
            !! Internal parameters:
            integer :: seed     !! Current value of the seed. The seed uniquely determines the realization of the disorder (reproducibility).
            real(wp) :: kdx     !! Value of k*dx = 2*pi*dx/lambda, with dx = L/nxdiso.
            real(wp) :: dos     !! Density of states in the waveguide, defined as: nu = (1/V) sum_p delta(k^2 - p^2), V being the volume of the system, and p the momentum (wavevector).
            integer :: nmu         !! Total number of modes, also length of mumesh.
            integer :: imin, imax  !! Index range of observed modes in the transmission matrix: t(i,j). Affected by filtering parameters.
            integer :: jmin, jmax  !! Index range of controlled modes in the transmission matrix: t(i,j). Affected by filtering parameters.
            integer :: nprop       !! Number of propagating modes (WARN: Is it used finally ????).
            integer :: neigv       !! Number of transmission eigenvalue corresponding to a single realization of the disorder.
            integer :: ny          !! Number of points on the y mesh (mesh in the transverse direction). Typically, ny >= nmu.
            character(len=50) :: boundtype           !! Type of boundary oncidtion in the transverse direction. Either "periodic" or "dirichlet".
            real(wp), allocatable :: xmesh(:)        !! Normalized positions of points in the longitudinal direction, x_i/L. Length=nx.
            real(wp), allocatable :: ymesh(:)        !! Normalized positions of points in the transverse direction, y_i/W. Length=ny.
            real(wp), allocatable :: kpmesh(:)       !! List of direction sines of the modes: kp_n = k_perp,n/k = sin(theta_n). Length=nmu.
            complex(wp), allocatable :: mumesh(:)    !! List of direction cosines of the modes: mu_n = k_x,n/k = sqrt(1 - (k_perp,n/k)^2) = cos(theta_n). Length=nmu.
            complex(wp), allocatable :: modal(:,:)   !! Modal matrix containing the modes in columns, modal(iy, imu). Dimensions=(ny, nmu).
            complex(wp), allocatable :: upot(:,:)    !! Dimensionless potential U(x,y)*dx/k. upot(ix,iy). Dimensions=(nxdiso, ny).
            complex(wp), allocatable :: bmat(:)      !! Diagonal of the B matrix (length = nmu).
            complex(wp), allocatable :: amat(:,:,:)  !! List of A matrices. A(ix, m, n): ix = Position index (x), (m,n) = Mode indices.
            
		contains
            !! Initialization:
            procedure :: init => init_waveguide                                 !! Initializes the waveguide type.
            procedure :: init_copy => init_copy_waveguide                       !! Initializes a copy of the given waveguide type.
            procedure :: del => del_waveguide                                   !! Deletes the content of the waveguide type.
            !! Computation:
            procedure :: compute_refl_trans_v1 => compute_refl_trans_v1_waveguide  !! Compute the reflection and transmission matrices.
            procedure :: compute_tdistrib => compute_tdistrib_waveguide         !! Compute and returns the transmission eigenvalues.
            procedure :: compute_wavefunction => compute_wavefunction_waveguide !! Compute and returns the wavefunction in the waveguide.
            procedure :: compute_deposition => compute_deposition_waveguide     !! Compute the deposition matrix (transmission matrix in the bulk).
            procedure :: compute_tprofile => compute_tprofile_waveguide         !! Compute the intensity profile of transmission-eigenstate, |psi_T(x)|^2.
            procedure :: compute_tprofile_group => compute_tprofile_group_waveguide !! Compute a group of transmission profiles for given values of T.
            !! Tests:
            procedure :: test_unitarity => test_unitarity_waveguide             !! Test the recursive Green algo v1 using unitarity relation t'*t + r'*r = 1.
            procedure :: test_green_v1_v2 => test_green_v1_v2_waveguide         !! Compare the predictions of the two versions of the recursive Green algo.
            !! Getters:
            procedure :: get_nx => get_nx_waveguide                             !! Returns the total number of point in the "x" direction.
            procedure :: get_xmesh => get_xmesh_waveguide                       !! Returns the "xmesh", the array of points "x".
            procedure :: get_nmu => get_nmu_waveguide                           !! Returns the number of simulated modes (dimensions of the transmission matrix). 
            procedure :: get_neigv => get_neigv_waveguide                       !! Returns the number of transmission eigenvalues taking into account the filtering.
            procedure :: get_parameters => get_parameters_waveguide             !! Returns all the base parameters.
            !! Output:
            procedure :: print_parameters => print_parameters_waveguide         !! Print the summary of parameters to the given output.
            procedure :: latex_parameters => latex_parameters_waveguide         !! Print the summary of parameters in LaTeX format.
            procedure :: generate_directory => generate_directory_waveguide     !! Generate an appropriate file path where to store the results.
            
    end type
    
	contains
	
	!!**********************************************************************************************
	!! Subroutine to initialize the waveguide type and allocate the arrays.
    !! After use, the destructor self%del() should be used to free the memory.
	!!**********************************************************************************************
	subroutine init_waveguide(self, boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        class(waveguide_t), intent(inout) :: self
        character(len=*), intent(in) :: boundtype
        real(wp), intent(in) :: wol, exfac, aratio, dscat, dabso, napera, naperb
        integer, intent(in) :: nxdiso, nxfree
        integer :: nxmin, i, n
        real(wp) :: kl
        real(wp), parameter :: exfacmax = 1e4   !! Maximum exfac used in simulations.
        
        !! Check for possible errors in the transverse direction (modal space):
        if (wol < epsilon(wol)) then
            write (stderr, '(a)') tag_error // "Invalid width-to-wavelength ratio."
            stop errmsg_invalid_arg
        else if (aratio < epsilon(aratio)) then
            write (stderr, '(a)') tag_error // "Invalid aspect ratio."
            stop errmsg_invalid_arg
        else if (exfac < 1. .or. exfac > exfacmax) then
            write (stderr, '(a,g0.6,a)') tag_error // "Invalid modal basis extension factor. It must be in [1, ", exfacmax, "]."
            stop errmsg_invalid_arg
        end if
        self%wol = wol
        self%aratio = aratio
        self%exfac = exfac
        kl = 2.0_wp*pi*wol*aratio    !! Product k*L = 2*pi*L/lambda (L=Length of the disordered region).
        
        !! Check for possible errors in the number of x points:
        nxmin = ceiling(kl/pi)  !! Minimum number of x points (Shannon limit).
        if (nxdiso < nxmin) then
            write (stdout, '(a,i0,a)') tag_warn // "Few x points compared to Shannon limit (choose at least nxdiso=", nxmin, ")."
        end if
        self%nxdiso = nxdiso
        if (nxfree < 0) then
            write (stderr, '(a)') tag_error // "Invalid number of x points, cannot be negative."
            stop errmsg_invalid_arg
        end if
        self%nxfree = nxfree
        self%nx = nxdiso + 2*nxfree   !! Total number of "x" points is "nxdiso" plus "nxfree" at the two edges.
        self%kdx = kl/nxdiso   !! Discretization step (assuming the conductor length is the unit length).
        
        !! Check the optical depths:
        if (dscat < 0. .or. dabso < 0.) then
            write (stderr, '(a)') tag_error // "Invalid optical depth, cannot be negative."
            stop errmsg_invalid_arg
        else if (dscat > kl .or. dabso > kl) then
            write (stdout, '(a,g0.3,a)') tag_warn // "The mean free path is smaller than the wavelength &
                &in contradiction with the weak disorder assumption (max is dscat=", kl, ")."
        end if
        self%dscat = dscat
        self%dabso = dabso
        
        !! Check the fraction of modes (filtering parameters):
        if (napera <= 0. .or. napera > 1.) then
            write (stderr, '(a,g0.3,a)') tag_error // &
                "Invalid filtering parameter napera=", napera, " (expected in [0, 1]), aborting..."
            stop errmsg_invalid_arg
        else if (naperb <= 0. .or. naperb > 1.) then
            write (stderr, '(a,g0.3,a)') tag_error // &
                "Invalid filtering parameter naperb=", naperb, " (expected in [0, 1]), aborting..."
            stop errmsg_invalid_arg
        end if
        self%napera = napera
        self%naperb = naperb
        
        !! Prepare the positional "x" mesh [xmesh(i) = x_i/L]:
        allocate(self%xmesh(self%nx))
        self%xmesh = [( self%aratio*(i - 0.5_wp - self%nxfree)/self%nxdiso , i = 1, self%nx )]
        !!call print_real_array("xmesh", "f0.3", self%xmesh)
        
        !! Prepare the modal quantities:
        select case (boundtype)
            case ("periodic")
                
                !! Check for possible errors:
                if (mod(wol, 1.0_wp) == 0.0_wp) then
                    write (stderr, '(a,f0.3,a)') tag_error // "Invalid width wol=", wol, &
                        ", cannot be an integer for '" // trim(boundtype) // "' boundary conditions..."
                    stop errmsg_invalid_arg
                end if
                
                !! Compute the number of modes (nmu, ny, nprop):
                self%nmu   = 2*floor(exfac*wol) + 1       !! Total number of modes to simulate (for periodic boundary conditions).
                self%ny    = self%nmu                     !! Number of point on the y mesh. Typically equal to nmu, but can be larger.
                self%imin  = self%nmu/2 + 1 - floor(naperb*wol)    !! Min index of the observed modes (output, port "b", right-hand port).
                self%imax  = self%nmu/2 + 1 + floor(naperb*wol)    !! Max index of the observed modes (output, port "b", right-hand port).
                self%jmin  = self%nmu/2 + 1 - floor(napera*wol)    !! Min index of the controlled modes (input, port "a", left-hand port).
                self%jmax  = self%nmu/2 + 1 + floor(napera*wol)    !! Max index of the controlled modes (input, port "a", left-hand port).
                
                self%neigv = 2*floor(min(napera, naperb)*wol) + 1  !! Number of transmission eigenvalues (or singular values) per disorder realization.
                self%nprop = 2*floor(wol) + 1                      !! Number of propagating modes.
                
                !! Prepare the positional "y" mesh [ymesh(i) = y_i/W]:
                allocate(self%ymesh(self%ny))
                self%ymesh = [( (i - 0.5_wp)/self%ny , i = 1, self%ny )]
                !!call print_real_array("ymesh", "f0.3", self%ymesh)
                
                !! Compute the transverse wavevector kpmesh(imu) = k_perp,n/k:
                allocate(self%kpmesh(self%nmu))
                self%kpmesh = [( (n - self%nmu/2 - 1)/wol, n = 1, self%nmu )]
                !!call print_real_array("kpmesh", "f0.3", self%kpmesh)
                
                !! Compute the modal matrix modal(iy, imu) = chi_imu(y_iy):
                allocate(self%modal(self%ny, self%nmu))
                do n = 1, self%nmu
                    do i = 1, self%ny
                        self%modal(i, n) = exp(iu * 2 * pi * wol * self%kpmesh(n) * self%ymesh(i))
                    end do
                end do
                self%modal = self%modal/sqrt(real(self%ny, kind=wp))   !! Normalize the matrix to get unitarity.
                !!call print_complex_matrix("modal", "f0.3", self%modal)
                
            case ("dirichlet")
                
                !! Check for possible errors:
                if (mod(wol, 0.5_wp) == 0.0_wp) then
                    write (stderr, '(a,f0.3,a)') tag_error // "Invalid width wol=", wol, &
                        ", cannot be a half integer for '" // trim(boundtype) // "' boundary conditions..."
                    stop errmsg_invalid_arg
                end if
                
                !! Compute the number of modes (nmu, ny, nprop):
                self%nmu   = floor(2*exfac*wol)   !! Number of modes to simulate (Dirichlet boundary conditions).
                self%ny    = self%nmu             !! Number of point on the y mesh. Typically equal to nmu, but can be larger.
                self%imin  = 1                     !! Min index of the observed modes (output, port "b", right-hand port).
                self%imax  = floor(2*naperb*wol)    !! Max index of the observed modes (output, port "b", right-hand port).
                self%jmin  = 1                     !! Min index of the controlled modes (input, port "a", left-hand port).
                self%jmax  = floor(2*napera*wol)    !! Max index of the controlled modes (input, port "a", left-hand port).
                
                self%neigv = floor(2*min(napera, naperb)*wol)   !! Number of transmission eigenvalues (or singular values) per disorder realization.
                self%nprop = floor(2*wol)                     !! Number of propagating modes.
                
                !! Prepare the positional "y" mesh [ymesh(i) = y_i/W]:
                allocate(self%ymesh(self%ny))
                self%ymesh = [( (i - 0.5_wp)/self%ny , i = 1, self%ny )]
                !!call print_real_array("ymesh", "f0.3", self%ymesh)
                
                !! Compute the transverse wavevector kpmesh(imu) = k_perp,n/k:
                allocate(self%kpmesh(self%nmu))
                self%kpmesh = [( 0.5_wp*n/wol, n = 1, self%nmu )]
                !!call print_real_array("kpmesh", "f0.3", self%kpmesh)
                
                !! Compute the modal matrix modal(iy, imu) = chi_imu(y_iy):
                allocate(self%modal(self%ny, self%nmu))
                do n = 1, self%nmu
                    do i = 1, self%ny
                        self%modal(i, n) = sin(2 * pi * wol * self%kpmesh(n) * self%ymesh(i))
                    end do
                end do
                self%modal = sqrt(2.0_wp/self%ny) * self%modal  !! Normalize the matrix to get quasi-unitarity.
                !!call print_complex_matrix("modal", "f0.3", self%modal)
                
            case default
                write (stderr, '(a)') tag_error // "Unknown boundary condition type '" // &
                    trim(boundtype) // "', aborting..."
                stop errmsg_invalid_arg
                
        end select
        self%boundtype = boundtype   !! Save the boundary type.
        
        !! Compute the mumesh and the DOS:
        allocate(self%mumesh(self%nmu))   !! Allocate space for the array of directions cosines of modes.
        self%mumesh = sqrt(1.0_wp + iu*dabso/kl - self%kpmesh**2)
        !!call print_complex_array("mumesh", "f0.3", self%mumesh)
        self%dos = sum( real(1.0_wp/sqrt((1.0_wp, 0.0_wp) - self%kpmesh**2)) )/(wol*(2.0_wp*pi)**2)
        !!print '(a,f0.6)', tag_info // "DOS = ", self%dos
        
        !! Prepare the B matrix since it is constant:
        allocate(self%bmat(self%nmu))
        self%bmat = self%mumesh/sin(self%mumesh*self%kdx)
        
        !! Allocate space for the other arrays:
        allocate(self%amat(self%nx, self%nmu, self%nmu))
        allocate(self%upot(self%nxdiso, self%ny))
        self%upot = 0.0_wp
        
	end subroutine
    
    !!**********************************************************************************************
    !! Initialize a copy of the given waveguide type "wguide_src".
    !! Initialize the current waveguide using the settings of the given waveguide "wguide_src".
    !! After use, the destructor self%del() should be used to free the memory.
    !!**********************************************************************************************
    subroutine init_copy_waveguide(self, wguide_src)
        class(waveguide_t), intent(inout) :: self
        class(waveguide_t), intent(in) :: wguide_src
        character(len=50) :: boundtype
        integer :: nxdiso, nxfree
        real(wp) :: wol, exfac, aratio, dscat, dabso, napera, naperb
        call wguide_src%get_parameters(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)  !! Extract the base parameters.
        call self%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)  !! Initialize "self" with the base parameters.
    end subroutine
    
	!!**********************************************************************************************
	!! Frees the memory allocated by the arrays of the muspace object "self".
	!!**********************************************************************************************
    subroutine del_waveguide(self)
        class(waveguide_t), intent(inout) :: self
        !!print '(a)', tag_info // "Deleting the waveguide content."
        if (allocated(self%xmesh)) deallocate(self%xmesh)
        if (allocated(self%ymesh)) deallocate(self%ymesh)
        if (allocated(self%kpmesh)) deallocate(self%kpmesh)
        if (allocated(self%mumesh)) deallocate(self%mumesh)
        if (allocated(self%modal)) deallocate(self%modal)
        if (allocated(self%upot)) deallocate(self%upot)
        if (allocated(self%bmat)) deallocate(self%bmat)
        if (allocated(self%amat)) deallocate(self%amat)
    end subroutine
    
    !!**********************************************************************************************
    !! Returns the number of simulated modes, also the dimensions of the transmission matrix (ignoring the filtering).
    !!**********************************************************************************************
    function get_nmu_waveguide(self) result(nmu)
        class(waveguide_t), intent(in) :: self
        integer :: nmu
        nmu = self%nmu
    end function
    
    !!**********************************************************************************************
    !! Returns the number of propagating modes, also the number of transmission eigenvalues contained
    !! in the array self%teigvals.
    !!**********************************************************************************************
    function get_neigv_waveguide(self) result(neigv)
        class(waveguide_t), intent(in) :: self
        integer :: neigv
        neigv = self%neigv
    end function
    
    !!**********************************************************************************************
    !! Returns the number of points in the mesh along the longitudinal direction of the waveguide ("x" direction).
    !!**********************************************************************************************
    function get_nx_waveguide(self) result(nx)
        class(waveguide_t), intent(in) :: self
        integer :: nx
        nx = self%nx
    end function
    
    !!**********************************************************************************************
    !! Returns the "xmesh", the position of points in the "x" direction.
    !!**********************************************************************************************
    subroutine get_xmesh_waveguide(self, xmesh) 
        class(waveguide_t), intent(in) :: self
        real(wp), intent(out) :: xmesh(:)
        xmesh = self%xmesh
    end subroutine
    
    !!**********************************************************************************************
    !! Returns all the base parameters of the current waveguide type.
    !! This subroutine can be used to create a deep copy of the current waveguide type.
    !!**********************************************************************************************
    subroutine get_parameters_waveguide(self, boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        class(waveguide_t), intent(in) :: self
        character(len=*), intent(out) :: boundtype
        integer, intent(out) :: nxdiso, nxfree
        real(wp), intent(out) :: wol, exfac, aratio, dscat, dabso, napera, naperb
        boundtype = self%boundtype
        nxdiso = self%nxdiso
        nxfree = self%nxfree
        wol = self%wol
        exfac = self%exfac
        aratio = self%aratio
        dscat = self%dscat
        dabso = self%dabso
        napera = self%napera
        naperb = self%naperb
    end subroutine
	
    !!**********************************************************************************************
    !! Set or reset the disordered Hamiltonian for the given seed.
    !!**********************************************************************************************
    subroutine set_hamiltonian_waveguide(self, iseed)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        real(wp), allocatable :: randseq(:)
        integer :: ix, nnr_tot
        real(wp) :: stdu
        
        !! 1. Prepare the dynamical matrix elements of the Hamiltonian (A_1, A_2, ..., A_nx) which are not related to the disorder:
        self%amat = 0.0_wp   !! Reset all the elements to zero. Necessary because the same object can be used for several disorder configurations.
        call add_diag(iu*self%mumesh - self%mumesh/tan(self%mumesh*self%kdx), self%amat(1, :, :))  !! Matrices with boundary conditions.
        self%amat(self%nx, :, :) = self%amat(1, :, :)    !! Copy the boundary condition matrix at the rightmost edge.
        call add_diag(-2.0_wp*self%mumesh/tan(self%mumesh*self%kdx), self%amat(2, :, :))  !! Bulk matrices.
        do ix = 3, self%nx-1
            self%amat(ix, :, :) = self%amat(2, :, :)  !! Copy the matrix within the bulk since they are the same.
        end do
        
        !! 2. Prepare a sequence of normal random numbers:
        self%seed = iseed   !! Stores the seed in the current waveguide type "self".
        nnr_tot = self%ny * self%nxdiso   !! Total number of normal random numbers required for one realization of the disorder. 
        allocate(randseq(nnr_tot))        !! Allocate space for the random sequence.
        call rand_normal(randseq, iseed)  !! Generate all the required normal random numbers for one realization of the disorder.
        stdu = sqrt(self%dscat*self%ny/(2.0_wp*pi*pi*self%dos*self%wol*self%nxdiso))  !! Standard deviation of the samples of U(x_s)*Deltax/k.
        randseq = stdu * randseq   !! Rescale the random elements by the standard deviation.
        self%upot = reshape(randseq, [self%nxdiso, self%ny])   !! Saves the disorder in the potential.
        
        !!call print_complex_matrix("Upot", "f0.3", self%upot)
        
        !! 3. Project the disordered potential onto the modal basis and add to the matrices A_1, A_2, ... A_nx:
        do ix = 1, self%nxdiso
            self%amat(ix+self%nxfree, :, :) = self%amat(ix+self%nxfree, :, :) &
                - matmul(transpose(conjg(self%modal)), diag_mul_left(self%upot(ix, :), self%modal))
        end do
        
        deallocate(randseq)  !! Free the array of normal random numbers since they are now stored in self%amat.
        
    end subroutine
    
    !!**********************************************************************************************
    !! Computes the transmission eigenvalues corresponding to the given seed "iseed" (must be 1..nchunk),
    !! and save them to the array "tallvals" in the chunk of index "iseed".
    !! This subroutine uses the recursive Green method to solve the wave equation within the disordered waveguide.
    !! This is the main computational routine of the waveguide type (no other routine call is needed).
    !!**********************************************************************************************
    subroutine compute_tdistrib_waveguide(self, iseed, tallvals)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        real(wp), intent(inout) :: tallvals(:)
        complex(wp), allocatable :: Gss(:, :), Gs1(:, :), tmat(:, :), sqrtmumesh(:)
        real(wp), allocatable :: teigvals(:)
        integer :: nchunk, cmin, cmax
        integer :: ix
        
        allocate(Gss(self%nmu, self%nmu))  !! Allocate the matrices.
        allocate(Gs1(self%nmu, self%nmu))
        allocate(sqrtmumesh(self%nmu))
        allocate(tmat(self%imax-self%imin+1, self%jmax-self%jmin+1))
        allocate(teigvals(self%neigv))
        
        !! 1. First check for possible errors:
        nchunk = size(tallvals)/self%neigv   !! Expected number of chunks.
        if (nchunk == 0) then
            write (stdout, '(a,i0,a)') tag_error // "Invalid eigenvalue array length, &
                &received ", size(tallvals), ", expected at least ", self%neigv, "..."
            stop errmsg_invalid_arg
        else if (nchunk*self%neigv /= size(tallvals)) then
            write (stdout, '(a,i0,a,i0,a)') tag_error // "Invalid eigenvalue array length, &
                &received ", size(tallvals), ", expected ", nchunk*self%neigv, "..."
            stop errmsg_invalid_arg
        else if (iseed < 1 .or. nchunk < iseed) then
            write (stdout, '(a,i0,a,i0,a)') tag_error // "Invalid chunk index, &
                &received ", iseed, ", expected in [1, ", nchunk, "]..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Prepare the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 3. Compute the transmission eigenvalues: Initialize the first slice, Gss = Gs1 = inv(A_1):
        Gss = self%amat(1, :, :)
        call inverse_complex(Gss)
        Gs1 = Gss
        
        do ix = 2, self%nx  !! Recursive Green method - Loop over the slices of the disordered region:
            
            !! Gss = inverse(Cs) = inv(A_ix - B * Gss * B)
            Gss = self%amat(ix, :, :) - diag_mul_right(diag_mul_left(self%bmat, Gss), self%bmat)
            call inverse_complex(Gss)
            
            !! Gs1 = -Gss * B * Gs1
            Gs1 = -matmul(Gss, diag_mul_left(self%bmat, Gs1))
            
        end do
        
        !! 4. Use the Fisher-Lee relation to get the transmission matrix:
        sqrtmumesh = sqrt(self%mumesh)
        tmat = 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh(self%imin:self%imax), &
            Gs1(self%imin:self%imax, self%jmin:self%jmax)), sqrtmumesh(self%jmin:self%jmax))
        
        !! 5. Compute the transmission eigenvalues using LAPACK's SVD:
        call singular_values(tmat, teigvals)  !! First compute the transmission coefficients using singular value decomposition.
        teigvals = teigvals * teigvals  !! The transmission eigenvalues are the SQUARES of the transmission coefficients.
        
        !! 6. Saves the transmission eigenvalues at the right place in "tallvals":
        cmin = (iseed - 1) * self%neigv + 1
        cmax = cmin + self%neigv - 1
        tallvals(cmin:cmax) = teigvals    !! Copy the eigenvalues at the right place.
        
        deallocate(Gss)  !! Deallocate memory.
        deallocate(Gs1)
        deallocate(tmat)
        deallocate(sqrtmumesh)
        deallocate(teigvals)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the reflection and transmission matrices for a given realization "iseed" of the disorder.
    !! The reflection and transmission matrices are written to "rmat" and "tmat", respectively.
    !! This subroutine uses version v1 of the recursive Green algorithm.
    !!**********************************************************************************************
    subroutine compute_refl_trans_v1_waveguide(self, iseed, rmat, tmat)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        complex(wp), intent(out) :: rmat(:, :), tmat(:, :)
        complex(wp), allocatable :: Gss(:, :), Gs1(:, :), G1s(:, :), G11(:, :)
        complex(wp) :: sqrtmumesh(self%nmu)
        integer :: ix
        
        allocate(Gss(self%nmu, self%nmu))  !! Allocate the matrices.
        allocate(Gs1(self%nmu, self%nmu))
        allocate(G1s(self%nmu, self%nmu))
        allocate(G11(self%nmu, self%nmu))
        
        !! Setup the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! Initialize the first slice, Gss = Gs1 = inv(A_1):
        Gss = self%amat(1, :, :)
        call inverse_complex(Gss)
        Gs1 = Gss
        G1s = Gss
        G11 = Gss
        
        do ix = 2, self%nx  !! Recursive Green method - Loop over the slices of the disordered region:
            
            !! Gss = inverse(Cs) = inv(A_ix - B * Gss * B)
            Gss = self%amat(ix, :, :) - diag_mul_right(diag_mul_left(self%bmat, Gss), self%bmat)
            call inverse_complex(Gss)
            
            !! Gs1 = -Gss * B * Gs1
            Gs1 = -matmul(Gss, diag_mul_left(self%bmat, Gs1))
            
            !! G11 = G11 - G1s * B * Gs1
            G11 = G11 - matmul(G1s, diag_mul_left(self%bmat, Gs1))
            
            !! G1s = -G1s * B * Gss
            G1s = -matmul(G1s, diag_mul_left(self%bmat, Gss))
            
        end do
        
        !! Use the Fisher-Lee relation to get the transmission and reflection matrices:
        sqrtmumesh = cmplx(sqrt(real(self%mumesh)), kind=wp)  !! The real part removes the evanescent modes.
        tmat = 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh, Gs1), sqrtmumesh)
        rmat = -identity_complex(self%nmu) + 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh, G11), sqrtmumesh)
        
        deallocate(Gss)  !! Deallocate memory.
        deallocate(Gs1)
        deallocate(G1s)
        deallocate(G11)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the reflection and transmission matrices for a given realization "iseed" of the disorder.
    !! The reflection and transmission matrices are written to "rmat" and "tmat", respectively.
    !! This subroutine uses version v2 of the recursive Green algorithm.
    !!**********************************************************************************************
    subroutine compute_refl_trans_v2_waveguide(self, iseed, rmat, tmat)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        complex(wp), intent(out) :: rmat(:, :), tmat(:, :)
        complex(wp), allocatable :: ind(:, :, :), green(:, :, :)
        complex(wp) :: sqrtmumesh(self%nmu)
        
        !! 1. Setup the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 2. Prepare the independent term:
        allocate(ind(self%nx, self%nmu, self%nmu))
        allocate(green(self%nx, self%nmu, self%nmu))
        ind = (0.0_wp, 0.0_wp)                      !! All the elements of the independent term are zero,
        ind(1, :, :) = identity_complex(self%nmu)   !! except the first block, which is the identity matrix.
        
        !! 3. Solves for the Green function "green":
        call solve_block_tridiag(self%amat, self%bmat, ind, green)
        
        !! 4. Use the Fisher-Lee relations to get the reflection and transmission matrices:
        sqrtmumesh = cmplx(sqrt(real(self%mumesh)), kind=wp)  !! The real part removes the evanescent modes.
        tmat = 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh, green(self%nx, :, :)), sqrtmumesh)
        rmat = -identity_complex(self%nmu) + 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh, green(1, :, :)), sqrtmumesh)
        
        deallocate(ind)
        deallocate(green)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the version v1 of the recursive Green algorithm using the unitarity of the reflection and transmission matrices.
    !! This subroutine is intended for testing and performs more computations than for the transmission
    !! eigenvalues only. Therefore, the use of this subroutine for computational purposes is not recommended.
    !! iseed = Seed index of the random realization of the disorder.
    !!**********************************************************************************************
    subroutine test_unitarity_waveguide(self, iseed, version)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed, version
        complex(wp), allocatable :: tmat(:, :), rmat(:, :), diff(:, :)
        
        allocate(tmat(self%nmu, self%nmu))
        allocate(rmat(self%nmu, self%nmu))
        allocate(diff(self%nmu, self%nmu))
        
        !! Test the unitarity:
        select case (version)
            case (1)
                call compute_refl_trans_v1_waveguide(self, iseed, rmat, tmat)
            case (2)
                call compute_refl_trans_v2_waveguide(self, iseed, rmat, tmat)
            case default
                write (stderr, '(a)') tag_error // "Unknown recursive Green version (must be 1 or 2), aborting..."
        end select
        
        diff = matmul(conjg(transpose(tmat)), tmat) + matmul(conjg(transpose(rmat)), rmat) - identity_complex(self%nmu)
        
        !!call print_complex_matrix("Diff", "f11.6", diff)
        
        write (stdout, '(a,i5,a,g0.6,a)') tag_test // "#", self%seed, &
            " | norm(t'*t + r'*r - 1) = ", norm2_complex(diff), " (expected 0)."
        
        deallocate(tmat)
        deallocate(rmat)
        deallocate(diff)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compare the Green matrices computed using the first and second version of the recursive Green
    !! algorithm. This subroutine is intended for testing purposes only. Therefore, its use for computation
    !! is not appropriate.
    !!**********************************************************************************************
    subroutine test_green_v1_v2_waveguide(self, iseed)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        complex(wp), allocatable :: rmatv1(:, :), tmatv1(:, :), rmatv2(:, :), tmatv2(:, :)
        
        allocate(rmatv1(self%nmu, self%nmu))
        allocate(tmatv1(self%nmu, self%nmu))
        allocate(rmatv2(self%nmu, self%nmu))
        allocate(tmatv2(self%nmu, self%nmu))
        
        !! Reflection and transmission matrices according to the two versions of the recursive Green algo:
        call compute_refl_trans_v1_waveguide(self, iseed, rmatv1, tmatv1)
        call compute_refl_trans_v2_waveguide(self, iseed, rmatv2, tmatv2)
        
        write (stdout, '(a,i5,a,g0.6,a)') tag_test // "#", self%seed, &
            " | norm(r_v1 - r_v2) = ", norm2_complex(rmatv1 - rmatv2), " (expected 0)."
        
        write (stdout, '(a,i5,a,g0.6,a)') tag_test // "#", self%seed, &
            " | norm(t_v1 - t_v2) = ", norm2_complex(tmatv1 - tmatv2), " (expected 0)."
        
        deallocate(rmatv1)
        deallocate(tmatv1)
        deallocate(rmatv2)
        deallocate(tmatv2)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the wavefunction in the waveguide for a given realization "iseed" of the disorder.
    !! The data are written to the real 2D array "wfundata" in columns: [x, y, Re(psi), Im(psi), Upot].
    !! Coordinates x=[0, aratio] and y=[-0.5, 0.5] correspond to the disordered region.
    !! 
    !! Arguments:
    !! iseed = (IN) Seed index of the realization of the disorder.
    !! ipos  = (IN) Index of the initial/incident position. It must be in [1, nx].
    !! imode = (IN) Index of the initial/incident mode. It must be in [1, nmu].
    !! wfundata = (OUT) 2D real array containing the wavefunction in columns: [x, y, Re(psi), Im(psi), Upot].
    !!**********************************************************************************************
    subroutine compute_wavefunction_waveguide(self, iseed, ipos, imode, wfundata)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed, ipos, imode
        real(wp), intent(out), allocatable :: wfundata(:, :)
        complex(wp), allocatable :: ind(:, :, :), green(:, :, :), psi(:)
        integer :: ix, iy, imu_mode
        
        allocate(wfundata(self%nx * self%ny, 5))  !! Allocation of "wfundata" is here. Deallocation must take place out of this subroutine.
        
        !! 1. Check for possible errors:
        if (imode < 1 .or. self%nmu < imode) then
            write (stderr, '(3(a,i0),a)') tag_error // "Invalid initial mode index, received imode=", imode, &
                ", expected in 1..", self%nmu, ", aborting..."
            stop errmsg_invalid_arg
        else if (ipos < 1 .or. self%nx < ipos) then
            write (stderr, '(2(a,i0),a)') tag_error // "Invalid initial position index, received ipos=", ipos, &
                ", expected in 1..", self%nx, ", aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Shift the mode index in the periodic case to get imode=1 correspond to the ground state.
        if (self%boundtype == "periodic") then
            imu_mode = mod(imode + self%nmu/2 - 1, self%nmu) + 1
        else
            imu_mode = imode
        end if
        
        !! 3. Setup the realization of the disorder:
        print '(a)', tag_info // "Setting up the disordered potential..."
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 4. Solves the linear system for the Green function "green":
        print '(a)', tag_info // "Solving the linear system for the Green function..."
        allocate(ind(self%nx, self%nmu, 1))
        allocate(green(self%nx, self%nmu, 1))
        allocate(psi(self%ny))   !! Temporary variable "psi" contains the values of psi(x,y) at the current value of "x".
        ind = (0.0_wp, 0.0_wp)                     !! All the elements of the independent term are zero,
        ind(ipos, imu_mode, 1) = (1.0_wp, 0.0_wp)  !! except the block at "ipos", which is the identity matrix.
        call solve_block_tridiag(self%amat, self%bmat, ind, green)
        
        !! 5. Compute and save the wave function at each pixel:
        print '(a)', tag_info // "Computing the wavefunction by Fourier transform..."
        do ix = 1, self%nx
            
            psi = matmul(self%modal, green(ix, :, 1))  !! Apply the DFT to the modal components of the Green function to switch to position basis.
            
            do iy = 1, self%ny
                wfundata((ix-1)*self%ny + iy, :) = [self%xmesh(ix), self%ymesh(iy), psi(iy)%re, psi(iy)%im, 0.0_wp]
            end do
        end do
        
        !! 6. Save also the potential Upot(x,y) = U(x,y)*dx/k:
        do ix = 1, self%nxdiso
            do iy = 1, self%ny
                wfundata((ix + self%nxfree - 1)*self%ny + iy, 5) = self%upot(ix, iy)%re
            end do
        end do
        
        deallocate(ind)
        deallocate(green)
        deallocate(psi)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the deposition matrix at some final position "fpos" for a given initial position "ipos".
    !! The deposition matrix is defined the same way as the transmission matrix but when the position of observation is taken withing the disordered medium.
    !! 
    !! Arguments:
    !! iseed  = (IN) Seed index of the realization of the disorder.
    !! fpos   = (IN) Index of the final position of observation along the "x" direction. Must be in [1, nx].
    !! ipos   = (IN) Index of the position along the "x" direction. Must be in [1, nx].
    !! depmat = (OUT) Deposition matrix. Dimensions: (nmu, nmu).
    !!**********************************************************************************************
    subroutine compute_deposition_waveguide(self, iseed, fpos, ipos, depmat)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed, ipos, fpos
        complex(wp), intent(out) :: depmat(:, :)
        complex(wp), allocatable :: ind(:, :, :), green(:, :, :)
        complex(wp) :: sqrtmumesh(self%nmu)
        
        !! 1. Check for possible errors:
        if (fpos < 1 .or. self%nx < fpos) then
            write (stderr, '(2(a,i0),a)') tag_error // "Invalid final position index, received fpos=", fpos, &
                ", expected in 1..", self%nx, ", aborting..."
            stop errmsg_invalid_arg
        else if (ipos < 1 .or. self%nx < ipos) then
            write (stderr, '(2(a,i0),a)') tag_error // "Invalid initial position index, received ipos=", ipos, &
                ", expected in 1..", self%nx, ", aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Setup the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 3. Solves the linear system for the Green function "green":
        allocate(ind(self%nx, self%nmu, self%nmu))
        allocate(green(self%nx, self%nmu, self%nmu))
        
        ind = (0.0_wp, 0.0_wp)                       !! All the elements of the independent term are zero,
        ind(ipos, :, :) = identity_complex(self%nmu) !! except the first block, which is the identity matrix.
        call solve_block_tridiag(self%amat, self%bmat, ind, green)
        
        sqrtmumesh = cmplx(sqrt(real(self%mumesh)), kind=wp)  !! The real part removes the evanescent modes.
        depmat(:, :) = 2.0 * iu * diag_mul_right(diag_mul_left(sqrtmumesh, green(fpos, :, :)), sqrtmumesh)
        
        deallocate(ind)
        deallocate(green)
        
    end subroutine
    
    !!**********************************************************************************************
    !! This subroutine computes the transmission-eigenstate intensity profile formally defined by
    !! 
    !!    I_T(x) = sum_n |psi_n(x,y)|^2 delta(T - T_n),
    !! 
    !! where psi_n(x,y) is the wavefunction of the n-th transmission eigenstate (whose eigenvalue is "T_n"),
    !! and "T" is the selected transmission eigenvalue for which the profile is desired.
    !! This quantity is averaged in the transverse direction ("y" direction).
    !! It should be noted that, in order to approximate the Dirac delta, a rectangular window [Tmin, Tmax] must be chosen.
    !! WARNING: This subroutine is now considered obsolete because it uses another normalization convention for the transmission eigenstate.
    !! 
    !! Arguments:
    !! iseed = (IN) Seed index of the realization of the disorder.
    !! tmin  = (IN) Minimum transmission eigenvalue of the search range for transmission eigenstates. Must be in [0, 1].
    !! tmax  = (IN) Maximum transmission eigenvalue of the search range for transmission eigenstates. Must be in [0, 1].
    !! itx   = (OUT) Real array containing the intensity profile function, I_T(x).
    !!**********************************************************************************************
    subroutine compute_itx_profile_waveguide(self, iseed, tmin, tmax, itx)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        real(wp), intent(in) :: tmin, tmax
        real(wp), intent(out) :: itx(:)
        complex(wp) :: ind(self%nx, self%nmu, self%nmu), green(self%nx, self%nmu, self%nmu), tmat(self%nmu, self%nmu)
        complex(wp) :: sqrtmumesh(self%nmu), U(self%nmu, self%nmu), vec(self%nmu)
        real(wp) :: tvals(self%nmu)
        integer :: iv, ix, nt
        
        !! 1. Check for possible errors:
        if (tmin < 0.0_wp .or. tmin > 1.0_wp) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid Tmin, received Tmin=", tmin, ", expected in [0, 1], aborting..."
            stop errmsg_invalid_arg
        else if (tmax < 0.0_wp .or. tmax > 1.0_wp) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid Tmax, received Tmax=", tmax, ", expected in [0, 1], aborting..."
            stop errmsg_invalid_arg
        else if (tmin >= tmax) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid range (Tmin >= Tmax), aborting..."
            stop errmsg_invalid_arg
        else if (size(itx) /= self%nx) then
            write (stderr, '(2(a,i0),a)') tag_error // "Invalid size for itx, received ", size(itx), &
                ", expected ", self%nx, ", aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Setup the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 3. Prepare the independent term:
        ind = (0.0_wp, 0.0_wp)                      !! All the elements of the independent term are zero,
        ind(1, :, :) = identity_complex(self%nmu)   !! except the first block, which is the identity matrix.
        
        !! 4. Solves for the Green function "green":
        call solve_block_tridiag(self%amat, self%bmat, ind, green)
        
        !! 5. Use the Fisher-Lee relations to get the transmission matrix:
        sqrtmumesh = cmplx(sqrt(real(self%mumesh)), kind=wp)  !! The real part removes the evanescent modes.
        tmat = 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh, green(self%nx, :, :)), sqrtmumesh)
        
        !! 6. Compute the singular values and singular vectors for the transmission matrix using LAPACK:
        call diagonalize_ahca(tmat, tvals, U)
        
        !! 7. Compute the intensity profile I_T(x) by summing over the transmission states in the interval [Tmin, Tmax]:
        itx(:) = 0.0_wp    !! Initialize to zero because it is incremented.
        nt = 0  !! Number of transmission eigenvalues found in the interval [Tmin, Tmax].
        
        do iv = 1, self%nmu  !! Loop on the transmission eigenvalues.
            if (tmin <= tvals(iv) .and. tvals(iv) <= tmax) then
                do ix = 1, self%nx  !! Loop on the slice.
                    vec = matmul(green(ix, :, :), sqrtmumesh*U(:, iv))
                    itx(ix) = itx(ix) + sum(real(vec*conjg(vec)))
                end do
                nt = nt + 1
            end if
        end do
        
        itx = itx * 2.0_wp/(pi*self%wol*abs(tmax-tmin))   !! Normalize by the correct factor.
        
        !!print '(2(a,i0),a)', tag_info // "Seed=", iseed, ", found ", nt, " transmission states..."
        
    end subroutine
    
    !!**********************************************************************************************
    !! This subroutine computes the transmission-eigenfunction profile defined by
    !! 
    !!    |psi_T(x)|^2 , where psi_T(x,y) = 2i*sqrt(k) sum_m,m' chi_m(y) G^+_mm'(x|x_a) sqrt(k_x,m') (u_T)_m'
    !! 
    !! psi_T(x,y) is the transmission eigenfunction associated to transmission eigenvalue T.
    !! This quantity is not only averaged in the transverse direction ("y" direction), but also over all
    !! the eigenvalues T located in the search interval [Tmin, Tmax].
    !! 
    !! Arguments:
    !! iseed    = (IN) Seed index of the realization of the disorder.
    !! tmin     = (IN) Minimum transmission eigenvalue of the search range for transmission eigenstates. Must be in [0, 1].
    !! tmax     = (IN) Maximum transmission eigenvalue of the search range for transmission eigenstates. Must be in [0, 1].
    !! tprofile = (OUT) Real array containing the intensity profile function, |psi_T(x)|^2.
    !! nt       = (OUT) Number of transmission states found in the given interval [Tmin, Tmax]. When nt=0, profile=0 and should be ignored for further averaging.
    !!**********************************************************************************************
    subroutine compute_tprofile_waveguide(self, iseed, tmin, tmax, tprofile, nt)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        real(wp), intent(in) :: tmin, tmax
        real(wp), intent(out) :: tprofile(:)
        integer, intent(out) :: nt
        complex(wp) :: ind(self%nx, self%nmu, self%nmu), green(self%nx, self%nmu, self%nmu), tmat(self%nmu, self%nmu)
        complex(wp) :: sqrtmumesh(self%nmu), U(self%nmu, self%nmu), vec(self%nmu)
        real(wp) :: tvals(self%nmu)
        integer :: iv, ix
        
        !! 1. Check for possible errors:
        if (tmin < 0.0_wp .or. tmin > 1.0_wp) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid Tmin, received Tmin=", tmin, ", expected in [0, 1], aborting..."
            stop errmsg_invalid_arg
        else if (tmax < 0.0_wp .or. tmax > 1.0_wp) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid Tmax, received Tmax=", tmax, ", expected in [0, 1], aborting..."
            stop errmsg_invalid_arg
        else if (tmin >= tmax) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid range (Tmin >= Tmax), aborting..."
            stop errmsg_invalid_arg
        else if (size(tprofile) /= self%nx) then
            write (stderr, '(2(a,i0),a)') tag_error // "Invalid size for tprofile, received ", size(tprofile), &
                ", expected ", self%nx, ", aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Setup the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 3. Prepare the independent term:
        ind = (0.0_wp, 0.0_wp)                      !! All the elements of the independent term are zero,
        ind(1, :, :) = identity_complex(self%nmu)   !! except the first block, which is the identity matrix.
        
        !! 4. Solves for the Green function "green":
        call solve_block_tridiag(self%amat, self%bmat, ind, green)
        
        !! 5. Use the Fisher-Lee relations to get the transmission matrix:
        sqrtmumesh = cmplx(sqrt(real(self%mumesh)), kind=wp)  !! The real part removes the evanescent modes.
        tmat = 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh, green(self%nx, :, :)), sqrtmumesh)
        
        !! 6. Compute the singular values and singular vectors for the transmission matrix using LAPACK:
        call diagonalize_ahca(tmat, tvals, U)
        
        !! 7. Compute the intensity profile |psi_T(x)|^2 by summing over the transmission states in the interval [Tmin, Tmax]:
        tprofile(:) = 0.0_wp    !! Initialize to zero because it is incremented.
        nt = 0  !! Number of transmission eigenvalues found in the interval [Tmin, Tmax].
        
        do iv = 1, self%nmu  !! Loop on the transmission eigenvalues.
            if (tmin <= tvals(iv) .and. tvals(iv) <= tmax) then
                do ix = 1, self%nx  !! Loop on the slice.
                    vec = matmul(green(ix, :, :), sqrtmumesh*U(:, iv))
                    tprofile(ix) = tprofile(ix) + sum(real(vec*conjg(vec)))
                end do
                nt = nt + 1
            end if
        end do
        
        if (nt /= 0) then
            tprofile = tprofile * (4.0_wp/nt)   !! Normalize by the correct factor.
        end if
        
        !!print '(2(a,i0),a)', tag_info // "Seed=", iseed, ", found ", nt, " transmission states..."
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the same transmission-eigenfunction profile as in the subroutine 
    !! compute_tprofile_quartet_waveguide() but for 4 transmission eigenvalues at the same time.
    !! This subroutine thus saves time considerably compared to the previous one.
    !! This subroutine also computes the average transmission probability Tavg in order to estimate the
    !! effective scattering thickness of the medium. This is necessary for L/lscat >> 1.
    !! 
    !! Arguments:
    !! iseed    = (IN) Seed index of the realization of the disorder.
    !! tm       = (IN) Real array of central transmission eigenvalues. Dimensions=nprofile.
    !! dt       = (IN) Real array of half-widths of the earch range of transmission eigenvalues. The intervals are thus [Tmin, Tmax] = [Tm-dt, Tm+dt]. Same length as "tm".
    !! tprofile = (OUT) Real array containing the intensity profile functions, |psi_T(x)|^2, for each value of Tm in the array "tm". Dimensions=(nx, nprofile).
    !! nt       = (OUT) Array of number of transmission states found in the given interval [Tmin, Tmax]. When nt=0, profile=0 and should be ignored for further averaging. Dimensions=nprofile.
    !! tavg     = (OUT) Average transmission probability computed as the average of transmission eigenvalues.
    !!**********************************************************************************************
    subroutine compute_tprofile_group_waveguide(self, iseed, tm, dt, tprofile, nt, tavg)
        class(waveguide_t), intent(inout) :: self
        integer, intent(in) :: iseed
        real(wp), intent(in) :: tm(:), dt(:)  !! Dimensions=nprofile.
        real(wp), intent(out) :: tprofile(:, :)  !! Dimensions=(nx, nprofile).
        integer, intent(out) :: nt(:)  !! Dimensions=nprofile.
        real(wp), intent(out) :: tavg
        complex(wp) :: ind(self%nx, self%nmu, self%nmu), green(self%nx, self%nmu, self%nmu)
        complex(wp) :: sqrtmumesh(self%nmu), vec(self%nmu)
        complex(wp), allocatable :: tmat(:, :), U(:, :)
        real(wp), allocatable :: tvals(:)
        integer :: ip, iv, ix, nprofile
        
        !! 1. Check for possible errors:
        nprofile = size(tm)   !! Number of transmission profiles to compute.
        if (any(tm < 0.0_wp) .or. any(tm > 1.0_wp)) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid transmission value, expected in [0, 1], aborting..."
            stop errmsg_invalid_arg
        else if (any(dt <= 0.0_wp) .or. any(dt > 0.5_wp)) then
            write (stderr, '(a,g0.3,a)') tag_error // "Invalid transmission half-width, expected 0 < dt <= 1/2, aborting..."
            stop errmsg_invalid_arg
        else if (size(tprofile, 1) /= self%nx .or. size(tprofile, 2) /= nprofile) then
            write (stderr, '(4(a,i0),a)') tag_error // "Invalid size for tprofile, received (", size(tprofile, 1), &
                ", ", size(tprofile, 2), "), expected (", self%nx, ", ", nprofile, "), aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Setup the realization of the disorder:
        call set_hamiltonian_waveguide(self, iseed)
        
        !! 3. Prepare the independent term:
        ind = (0.0_wp, 0.0_wp)                      !! All the elements of the independent term are zero,
        ind(1, :, :) = identity_complex(self%nmu)   !! except the first block, which is the identity matrix.
        
        !! 4. Solves for the Green function "green":
        call solve_block_tridiag(self%amat, self%bmat, ind, green)
        
        !! 5. Use the Fisher-Lee relations to get the transmission matrix:
        allocate(tmat(self%imax-self%imin+1, self%jmax-self%jmin+1))
        sqrtmumesh = cmplx(sqrt(real(self%mumesh)), kind=wp)  !! The real part removes the evanescent modes.
        tmat = 2.0_wp * iu * diag_mul_right(diag_mul_left(sqrtmumesh(self%imin:self%imax), &
            green(self%nx, self%imin:self%imax, self%jmin:self%jmax)), sqrtmumesh(self%jmin:self%jmax))
        
        !! 6. Compute the singular values and singular vectors for the transmission matrix using LAPACK:
        allocate(tvals(self%neigv))
        allocate(U(self%neigv, self%neigv))
        call diagonalize_ahca(tmat, tvals, U)
        tavg = sum(tvals)/self%neigv  !! Compute the average transmission probability as the average of transmission eigenvalues.
        
        !! 7. Compute the intensity profile |psi_T(x)|^2 by summing over the transmission states in the interval [Tmin, Tmax]:
        do ip = 1, nprofile  !! Loop on the transmission profiles to compute.
            tprofile(:, ip) = 0.0_wp  !! Initialize to zero because it is incremented.
            nt(ip) = 0  !! Number of transmission eigenvalues found in the interval [Tmin, Tmax].
            do iv = 1, self%neigv  !! Loop on the transmission eigenvalues.
                if (abs(tvals(iv) - tm(ip)) <= dt(ip)) then   !! If the eigenvalue is located in the search interval, then add it.
                    do ix = 1, self%nx  !! Loop on the slice.
                        vec = matmul(green(ix, :, self%jmin:self%jmax), &
                            sqrtmumesh(self%jmin:self%jmax)*U(:, iv))
                        tprofile(ix, ip) = tprofile(ix, ip) + sum(real(vec*conjg(vec)))
                    end do
                    nt(ip) = nt(ip) + 1
                end if
            end do
            if (nt(ip) /= 0) then  !! Avoid dividing by zero.
                tprofile(:, ip) = tprofile(:, ip) * (4.0_wp/nt(ip))  !! Normalize by the correct factor.
            end if
        end do
        
        !! 8. Free memory:
        deallocate(tmat)
        deallocate(tvals)
        deallocate(U)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Print a summary of the parameters of the current waveguide type in the output "fp" (typically
    !! the standard output or a file), with the given "prefix" (typically a comment character).
    !!**********************************************************************************************
    subroutine print_parameters_waveguide(self, fp, prefix)
        class(waveguide_t), intent(in) :: self
        integer, intent(in) :: fp
        character(len=*), intent(in) :: prefix
        
        write (fp, '(2(a,f0.6),2(a,i0))') &
            trim(prefix) // " Waveguide 2D with boundtype='" // trim(self%boundtype) // "', W/lambda=", self%wol, &
            ", exfac=", self%exfac, ", Nmode=", self%nmu, ", Nprop=", self%nprop
        write (fp, '(2(a,f0.6,a,i0),a,i0)') &
            trim(prefix) // " Filtering: NAper_a=", self%napera, " (Nmode_a=", (self%jmax-self%jmin+1), &
            "), NAper_b=", self%naperb, " (Nmode_b=", (self%imax-self%imin+1), "), Neigv=", self%neigv
        write (fp, '(2(a,f0.6),2(a,i0),a,f0.6)') &
            trim(prefix) // " Medium: L/lscat=", self%dscat, ", L/labso=", self%dabso, &
            ", Nx_diso=", self%nxdiso, ", Nx_free=", self%nxfree, ", aratio=L/W=", self%aratio
        
    end subroutine
    
    !!**********************************************************************************************
    !! Print a summary of the parameters of the current waveguide type in the string "latexstr"
    !! in the LaTeX format. This will be used for figure title.
    !!**********************************************************************************************
    subroutine latex_parameters_waveguide(self, latexstr)
        class(waveguide_t), intent(in) :: self
        character(len=*), intent(out) :: latexstr
        
        write (latexstr, '(2(a,f0.6),2(a,i0),a)') &
            "Waveguide 2D with boundtype='" // trim(self%boundtype) // "', $W/\lambda=", self%wol, &
            "$, $e_{\rm f}=", self%exfac, "$, $N_{\rm mode}=", self%nmu, "$, $N_{\rm prop}=", self%nprop, "$, "
        write (latexstr, '(2(a,f0.6,a,i0),a,i0,a)') trim(latexstr) // &
            "$m_{\rm a}=", self%napera, "$ ($N_{\rm mode,a}=", (self%jmax-self%jmin+1), &
            "$), $m_{\rm b}=", self%naperb, "$ ($N_{\rm mode,b}=", (self%imax-self%imin+1), &
            "$), $N_{\rm eigv}=", self%neigv, "$, "
        write (latexstr, '(2(a,f0.6),2(a,i0),a,f0.6,a)') trim(latexstr) // &
            "$L/\ell_{\rm s}=", self%dscat, "$, $L/\ell_{\rm a}=", self%dabso, "$, $N_{x,{\rm diso}}=", self%nxdiso, &
            "$, $N_{x,{\rm free}}=", self%nxfree, "$, $a=L/W=", self%aratio, "$"
        
    end subroutine
    
    !!**********************************************************************************************
    !! Generate the appropriate file path corresponding to the current "waveguide" type and write
    !! it to the string "outputdir". This string must be a folder, not the full filename.
    !! Example: outputdir = "out/dscat_0.2/dabso_0./wol_25.5/periodic/"
    !!**********************************************************************************************
    subroutine generate_directory_waveguide(self, outputdir)
        class(waveguide_t), intent(in) :: self
        character(len=*), intent(out) :: outputdir
        character(len=30) :: str_scat, str_abso, str_wol, folder_scat, folder_abso, folder_wol
        
        str_scat = "dscat_"
        str_abso = "dabso_"
        str_wol  = "wol_"
        if (self%dscat < 1.) str_scat = trim(str_scat) // "0"   !! Append a zero just before the decimal point.
        if (self%dabso < 1.) str_abso = trim(str_abso) // "0"
        if (self%wol   < 1.) str_wol  = trim(str_wol)  // "0"
        
        write (folder_scat, '(a,f0.6)') trim(str_scat), self%dscat
        write (folder_abso, '(a,f0.6)') trim(str_abso), self%dabso
        write (folder_wol,  '(a,f0.6)') trim(str_wol),  self%wol
        
        !! Prepare the path for the output directory:
        outputdir = "out" // pathsep // trim_zero(folder_scat) // pathsep // &
                trim_zero(folder_abso) // pathsep // &
                trim_zero(folder_wol) // pathsep // trim(self%boundtype) // pathsep 
        
    end subroutine
    
end module
