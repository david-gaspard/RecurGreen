!!**************************************************************************************************
!! Created on 2024-05-15 at 15:00:37 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran module containing the basic utilities of the RecurGreen program including mathematical constants.
!!**************************************************************************************************
module base_utils
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, input_unit, output_unit, error_unit
	implicit none
    
    !! Intrinsic constants:
    integer, parameter :: wp     = real64       !! Current working precision (wp) is double precision.
    integer, parameter :: stdin  = input_unit   !! Standard input.
    integer, parameter :: stdout = output_unit  !! Standard output display.
    integer, parameter :: stderr = error_unit   !! Standard error display.
    
    !! Useful mathematical constants:
    real(wp), parameter :: pi = acos(-1.0_wp)      !! Pi constant at working precision.
    complex(wp), parameter :: iu = (0._wp, 1._wp)  !! Imaginary unit (iu).
    
    !! Global strings regarding the program:
    character(len=*), parameter :: program_shortname = "RecurGreen v1.0"
    character(len=*), parameter :: copyright = "(c) 2024 David GASPARD <david.gaspard@espci.fr>"
    character(len=*), parameter :: program_fullname  = program_shortname // " - Computing the transmission statistics &
        &through a disordered waveguide"
    character(len=*), parameter :: program_copyright = program_shortname // " " // copyright
    
    !! Length of the progress bar (see subroutine below):
    integer, parameter :: progress_bar_length = 50
    
    !! Special characters:
    character, parameter :: char_tab = char(9)       !! Tab character.
    character, parameter :: char_nl = new_line('a')  !! New line character.
    character, parameter :: char_cr = achar(13)      !! Carriage return character (used for progress bar).
    character, parameter :: pathsep = '/'            !! Path separator used on UNIX filesystems.
    
    !! Define colors (for TTY only):
    character(len=*), parameter :: tcolor_nc         = achar(27)//"[0m"
    character(len=*), parameter :: tcolor_bold       = achar(27)//"[1m"
    character(len=*), parameter :: tcolor_red        = achar(27)//"[31m"
    character(len=*), parameter :: tcolor_lightred   = achar(27)//"[91m"
	character(len=*), parameter :: tcolor_lightgreen = achar(27)//"[92m"
	character(len=*), parameter :: tcolor_yellow     = achar(27)//"[93m"
	character(len=*), parameter :: tcolor_blue       = achar(27)//"[94m"
	character(len=*), parameter :: tcolor_purple     = achar(27)//"[95m"
	character(len=*), parameter :: tcolor_cyan       = achar(27)//"[96m"
	
    !! Constant tags:
    character(len=*), parameter :: tag_info  = "[INFO] "
    character(len=*), parameter :: tag_warn  = "["//tcolor_bold//tcolor_yellow//"WARN"//tcolor_nc//"] "
    character(len=*), parameter :: tag_error = "["//tcolor_bold//tcolor_red//"ERROR"//tcolor_nc//"] "
    character(len=*), parameter :: tag_test  = "["//tcolor_bold//tcolor_cyan//"TEST"//tcolor_nc//"] "
    character(len=*), parameter :: tag_fail  = "["//tcolor_bold//tcolor_lightred//"FAIL"//tcolor_nc//"] "
	character(len=*), parameter :: tag_pass  = "["//tcolor_bold//tcolor_lightgreen//" OK "//tcolor_nc//"] "
    
	!! Generic error messages:
	character(len=*), parameter :: errmsg_invalid_arg   = "RecurGreen: Invalid arguments, aborting..."
	character(len=*), parameter :: errmsg_unbound_index = "RecurGreen: Index out of bounds, aborting..."
	
	contains
    
    !!**********************************************************************************************
    !! Returns the n-by-n identity matrix defined on complex numbers.
    !!**********************************************************************************************
    function identity_complex(n) result(identity)
        integer, intent(in) :: n
        complex(wp), allocatable :: identity(:, :)
        integer :: i
        
        allocate(identity(n, n))
        identity = (0.0_wp, 0.0_wp)
        do i = 1, n
            identity(i, i) = (1.0_wp, 0.0_wp)
        end do

    end function
    
    !!**********************************************************************************************
    !! Returns the 2-norm of the given complex matrix.
    !!**********************************************************************************************
    pure function norm2_complex(A) result(norm)
        complex(wp), intent(in) :: A(:, :)
        real(wp) :: norm
        norm = sqrt(sum(real(A*conjg(A))))
    end function
    
    !!**********************************************************************************************
    !! Add the elements "elems" along the main diagonal of "matrix".
    !!**********************************************************************************************
    subroutine add_diag(elems, matrix)
        complex(wp), intent(in) :: elems(:)
        complex(wp), intent(inout) :: matrix(:, :)
        integer :: n, i
        
        !! Check for possible errors:
        n = min(size(matrix, 1), size(matrix, 2))
        if (size(elems) /= n) then
            write (stderr, '(a,i0,a,i0,a)') tag_error // "add_diag: Invalid number of elements &
                &(received ", size(elems), " expected ", n, ")."
            stop errmsg_invalid_arg
        end if
        
        !! Add the band to the main diagonal of "matrix":
        do i = 1, n  !! Loop along the main diagonal.
            matrix(i, i) = matrix(i, i) + elems(i)
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Add the elements "elems" along the diagonals of "matrix" while ensuring the Hermiticity of "matrix".
    !! The pattern of the increment of "matrix" is:
    !! | elem(1)            elem(2)+elem(3)*i  elem(4)+elem(5)*i  elem(6)+elem(7)*i | 
    !! | elem(2)-elem(3)*i  elem(1)            elem(2)+elem(3)*i  elem(4)+elem(5)*i | 
    !! | elem(4)-elem(5)*i  elem(2)-elem(3)*i  elem(1)            elem(2)+elem(3)*i |
    !! | elem(6)-elem(7)*i  elem(4)-elem(5)*i  elem(2)-elem(3)*i  elem(1)           |
    !!**********************************************************************************************
    subroutine add_band_matrix(elems, matrix)
        real(wp), intent(in) :: elems(:)
        complex(wp), intent(inout) :: matrix(:, :)
        integer :: n, k, l
        real(wp) :: re, im
        
        !! Check for possible errors:
        n = size(matrix, 1)
        if (size(matrix, 2) /= n) then
            write (stderr, '(a)') tag_error // "add_band_matrix: Matrix is not square."
            stop errmsg_invalid_arg
        else if (size(elems) /= 2*n - 1) then
            write (stderr, '(a,i0,a,i0,a)') tag_error // "add_band_matrix: Invalid number of elements &
                &(received ", size(elems), " expected ", 2*n - 1, ")."
            stop errmsg_invalid_arg
        end if
        
        !! Add the band to the main diagonal of "matrix":
        do l = 1, n  !! Loop along the main diagonal.
            matrix(l, l) = matrix(l, l) + elems(1)
        end do
        
        !! Add the other bands next to the main diagonal:
        do k = 1, n - 1 
            re = elems(2*k)
            im = elems(2*k + 1)
            do l = 1, n - k  !! Loop along the current diagonal (of shift 'k').
                matrix(l, l+k) = matrix(l, l+k) + re + iu*im
                matrix(l+k, l) = matrix(l+k, l) + re - iu*im
            end do
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Multiply the matrix M by the diagonal matrix D given by its diagonal values to the left.
    !! Returns the product P = D*M, where '*' means the usual matrix product.
    !!**********************************************************************************************
    function diag_mul_left(D, M) result(P)
        complex(wp), intent(in) :: D(:), M(:, :)
        complex(wp) :: P(size(M, 1), size(M, 2))
        integer :: i, j, nrow, ncol
        
        !! Check for possible errors:
        nrow = size(M, 1)  !! Number of rows of M.
        ncol = size(M, 2)  !! Number of columns of M.
        if (size(D) /= nrow) then !! Number of elements in D must be equal to the number of rows of M.
            write (stderr, '(a,i0,a,i0,a)') tag_error // "diag_mul_left: Incorrect number of elements of D &
                &(received ", size(D), ", expected ", nrow, ")."
            stop errmsg_invalid_arg
        end if
        
        !! Multiply the matrices:
        do i = 1, nrow
            do j = 1, ncol
                P(i, j) = D(i) * M(i, j)
            end do
        end do
        
    end function
    
    !!**********************************************************************************************
    !! Multiply the matrix M by the diagonal matrix D given by its diagonal values to the right.
    !! Returns the product P = M*D, where '*' means the usual matrix product.
    !!**********************************************************************************************
    function diag_mul_right(M, D) result(P)
        complex(wp), intent(in) :: D(:), M(:, :)
        complex(wp) :: P(size(M, 1), size(M, 2))
        integer :: i, j, nrow, ncol
        
        !! Check for possible errors:
        nrow = size(M, 1)  !! Number of rows of M.
        ncol = size(M, 2)  !! Number of columns of M.
        if (size(D) /= ncol) then !! Number of elements in D must be equal to the number of rows of M.
            write (stderr, '(a,i0,a,i0,a)') tag_error // "diag_mul_left: Incorrect number of elements of D &
                &(received ", size(D), ", expected ", ncol, ")."
            stop errmsg_invalid_arg
        end if
        
        !! Multiply the matrices:
        do i = 1, nrow
            do j = 1, ncol
                P(i, j) = M(i, j) * D(j)
            end do
        end do
        
    end function
    
    !!**********************************************************************************************
    !! Computes the product of two complex matrices C = A * B using LAPACK's zgemm() subroutine.
    !! It seems that with gfortran, matmul() is faster than zgemm(), so the use of the present
    !! subroutine is not recommended.
    !!**********************************************************************************************
    function product_complex(A, B) result(C)
        complex(wp) :: A(:, :), B(:, :), C(size(A, 1), size(B, 2))
        integer :: m, n, k
        character, parameter :: trans = 'N'   !! Tell LAPACK there is no matrix transpose.
        complex(wp), parameter :: alpha = 1.0, beta = 0.0
        
        external zgemm  !! External procedures provided by LAPACK.
        
        !! Check for possible errors:
        m = size(A, 1)
        n = size(B, 2)
        k = size(A, 2)
        if (k /= size(B, 1)) then
            stop program_shortname // ": zgemm: Matrix size do not match for the product, aborting..."
        end if
        
        call zgemm(trans, trans, m, n, k, alpha, A, m, B, k, beta, C, m)
        
    end function
    
	!!**********************************************************************************************
    !! Computes the inverse of A matrix calculated by finding the LU decomposition, and overwrites A.
    !! This subroutine depends on LAPACK.
    !! Adapted from: https://fortranwiki.org/fortran/show/Matrix+inversion
	!!**********************************************************************************************
    subroutine inverse_complex(A)
        complex(wp), intent(in) :: A(:, :)
        complex(wp), allocatable :: work(:)  !! Work array for LAPACK.
        integer, allocatable :: ipiv(:)   !! Pivot indices.
        integer :: n, info
        character(len=10) :: strinfo  !! String version of "info".
        
        external zgetrf  !! External procedures provided by LAPACK.
        external zgetri
        
        if (size(A, 1) /= size(A, 2)) then
            write (stderr, '(a)') tag_error // "inverse_complex: Matrix A is not square, aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! ZGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
        n = size(A, 1)
        allocate(ipiv(n))
        call zgetrf(n, n, A, n, ipiv, info)
        
        if (info /= 0) then
            write (strinfo, '(i0)') info
            stop program_shortname // ": zgetrf: Matrix is numerically singular &
                &(LAPACK info=" // trim(strinfo) // "), aborting..."
        end if
        
        !! ZGETRI computes the inverse of a matrix using the LU factorization computed by ZGETRF.
        allocate(work(n))
        call zgetri(n, A, n, ipiv, work, n, info)
        
        if (info /= 0) then
            write (strinfo, '(i0)') info
            stop program_shortname // ": zgetri: Matrix inversion failed &
                &(LAPACK info=" // trim(strinfo) // "), aborting..."
        end if
        
        deallocate(ipiv)
        deallocate(work)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Solves the linear system H*X = I, where H is a block-tridiagonal matrix with [A_1, A_2, ..., A_i, ..., A_N]
    !! along the diagonal and [B, B, ..., B] along the subdiagonal and superdiagonal.
    !! Matrix I is the independent term (denoted "imat" in the code).
    !! 
    !! The full system looks like:
    !! | A_1  B                  | | X_1 |   | I_1 | 
    !! | B    A_2  B             | | X_2 |   | I_2 | 
    !! |      B    A_3  ...      |.| X_3 | = | I_3 | 
    !! |           ...  ...  B   | | ... |   | ... | 
    !! |                B    A_N | | X_N |   | I_N | 
    !! 
    !! The computational cost is N block inversions, 2*N full block multiplications, and 4*N multiplications by a diagonal block.
    !! 
    !! Arguments:
    !! amat  = (IN)  List of diagonal blocks amat(i, m, n). i=Major index of the block, (m,n)=Indices of the elements within the blocks.
    !! bdiag = (IN)  List of diagonal elements of the block B.
    !! imat  = (IN)  List of blocks of the independent term imat(i, m, n).
    !! xmat  = (OUT) Solution of the block-tridiagonal system x(i, m, n).
    !!**********************************************************************************************
    subroutine solve_block_tridiag(amat, bdiag, imat, xmat)
        complex(wp), intent(in) :: amat(:, :, :), bdiag(:), imat(:, :, :)
        complex(wp), intent(out) :: xmat(:, :, :)
        complex(wp) :: smat(size(amat, 1), size(amat, 2), size(amat, 3)), zmat(size(imat, 1), size(imat, 2), size(imat, 3))
        !! NB: The allocations of "smat" and "zmat" can be memory-intensive. Typical size is: 200 x 200 x 200 = 8 x 10^6 c-numbers.
        integer :: i, n
        
        !! 1. Check for possible errors:
        n = size(amat, 1)  !! Index of the last block (also number of blocks).
        if (size(amat, 2) /= size(amat, 3)) then
            write (stderr, '(a)') tag_error // "solve_block_tridiag: Matrices A_1, A_2, ..., A_N are not square."
            stop errmsg_invalid_arg
        else if (size(bdiag) /= size(amat, 2)) then
            write (stderr, '(a)') tag_error // "solve_block_tridiag: Incompatible matrix shape of B and A_1."
            stop errmsg_invalid_arg
        else if (size(imat, 1) /= size(amat, 1) .or. size(imat, 2) /= size(amat, 2)) then
            write (stderr, '(a)') tag_error // "solve_block_tridiag: Incompatible matrix shape of I and A."
            stop errmsg_invalid_arg
        else if (any(shape(xmat) /= shape(imat))) then
            write (stderr, '(a)') tag_error // "solve_block_tridiag: Incompatible matrix shape of X and I."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Solves the block-tridiagonal system:
        smat(1, :, :) = amat(1, :, :)  !! Initialize the first recursion.
        call inverse_complex(smat(1, :, :))
        zmat(1, :, :) = matmul(smat(1, :, :), imat(1, :, :))
        
        do i = 2, n   !! Forward loop.
            smat(i, :, :) = amat(i, :, :) - diag_mul_right(diag_mul_left(bdiag, smat(i-1, :, :)), bdiag)
            call inverse_complex(smat(i, :, :))
            zmat(i, :, :) = matmul(smat(i, :, :), imat(i, :, :) - diag_mul_left(bdiag, zmat(i-1, :, :)))
        end do
        
        xmat(n, :, :) = zmat(n, :, :)  !! Initialize the back-recursion.
        do i = n-1, 1, -1   !! Backward loop.
            xmat(i, :, :) = zmat(i, :, :) - matmul(smat(i, :, :), diag_mul_left(bdiag, xmat(i+1, :, :)))
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the singular values of the given matrix "A" and store them in the real array "S".
    !! See LAPACK manual: https://netlib.org/lapack/explore-html-3.6.1/d3/da8/group__complex16_g_esing_gab239ad131cb8f1ab83955894ec7a7b6e.html
    !!**********************************************************************************************
    subroutine singular_values(A, S)
        complex(wp) :: A(:, :)
        real(wp) :: S(:)
        character, parameter :: job = 'N'   !! Tell LAPACK not to calculate the matrices U and V.
        complex(wp) :: U(0, 0), VT(0, 0)    !! Unreferenced matrices. These matrices are not computed by LAPACK.
        complex(wp), allocatable :: work(:) !! LAPACK work space arrays.
        real(wp), allocatable :: rwork(:)   !! LAPACK work space arrays.
        integer :: m, n, lwork, info
        character(len=10) :: strinfo  !! Information string.
        
        external zgesvd  !! External procedures provided by LAPACK.
        
        !! Check for possible errors:
        m = size(A, 1)  !! Number of rows of A.
        n = size(A, 2)  !! Number of columns of A.
        if (size(S) /= min(m, n)) then
            stop program_shortname // ": zgesvd: Invalid number of singular values, aborting..."
        end if
        
        !! Prepare the work space:
        lwork =  max(1, 2*min(m, n) + max(m, n))
        allocate(work(lwork))
        allocate(rwork(5*min(m, n)))
        
        call zgesvd(job, job, m, n, A, m, S, U, m, VT, n, work, lwork, rwork, info)
        
        if (info < 0) then
            write (strinfo, "(a,i0,a,i0,a)") program_shortname // ": zgesvd: The ", -info, "-th argument &
                &had an illegal value (LAPACK info=", info, "), aborting..."
            stop strinfo
        else if (info > 0) then
            write (strinfo, "(a,i0,a)") program_shortname // ": zgesvd: Failed to find the singular values, &
                &zbdsqr() did not converge (LAPACK info=", info, "), aborting..."
            stop strinfo
        end if
        
        deallocate(work)  !! Frees LAPACK's work space.
        deallocate(rwork)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the eigenvalues and eigenvectors of the Hermitian matrix A^H * A, where "^H" denotes the
    !! Hermitian conjugation (transpose conjugation), by resorting to singular value decomposition.
    !! Indeed, building A^H * A explicitly is not recommended since it is numerically badly conditioned.
    !! In practice, this subrutine serves to compute the transmission eigenvalues and transmission eigenvectors.
    !! In other words, this subroutine computes the eigendecomposition: A^H * A = U * T * U^H
    !! 
    !! Arguments:
    !! A = (IN)  Complex rectangular matrix (such as the transmission matrix). This array is destroyed on exit.
    !! T = (OUT) List of eigenvalues of A^H * A in decreasing order. Assumed to be already allocated to "n".
    !! U = (OUT) List of eigenvectors of A^H * A in the same order. Assumed to be already allocated to n-by-n.
    !!**********************************************************************************************
    subroutine diagonalize_ahca(A, T, U)
        complex(wp), intent(inout) :: A(:, :)
        real(wp), intent(out) :: T(:)
        complex(wp), intent(out) :: U(:, :)
        character, parameter :: jobv = 'N', jobwh = 'A'   !! Tell LAPACK to calculate W^H but not V.
        complex(wp) :: V(0, 0)   !! Unreferenced matrix. This matrix will not be computed by LAPACK.
        complex(wp), allocatable :: WH(:, :)  !! The second unitary matrix of the SVD must be computed.
        complex(wp), allocatable :: work(:) !! LAPACK work space arrays.
        real(wp), allocatable :: rwork(:)   !! LAPACK work space arrays.
        integer :: m, n, lwork, info
        character(len=10) :: strinfo  !! Information string.
        
        external zgesvd  !! External procedures provided by LAPACK.
        
        !! Check for possible errors:
        m = size(A, 1)  !! Number of rows of A.
        n = size(A, 2)  !! Number of columns of A.
        if (size(T) /= n) then
            stop program_shortname // ": zgesvd: Invalid number of singular values, aborting..."
        else if (size(U, 1) /= n .or. size(U, 2) /= n) then
            stop program_shortname // ": zgesvd: Invalid size of unitary matrix, aborting..."
        end if
        
        !! Prepare the work space:
        lwork = max(1, 2*min(m, n) + max(m, n))
        allocate(work(lwork))
        allocate(rwork(5*min(m, n)))
        allocate(WH(n, n))
        
        call zgesvd(jobv, jobwh, m, n, A, m, T, V, m, WH, n, work, lwork, rwork, info)
        
        if (info < 0) then
            write (strinfo, "(a,i0,a,i0,a)") program_shortname // ": zgesvd: The ", -info, "-th argument &
                &had an illegal value (LAPACK info=", info, "), aborting..."
            stop strinfo
        else if (info > 0) then
            write (strinfo, "(a,i0,a)") program_shortname // ": zgesvd: Failed to find the singular values, &
                &zbdsqr() did not converge (LAPACK info=", info, "), aborting..."
            stop strinfo
        end if
        
        !! Adapt the results to get the expected eigendecomposition:
        T = T*T  !! The eigenvalues of A^H * A are the squares of the singular values of A.
        U = transpose(conjg(WH))  !! Saves the eigenvectors of A^H * A in U.
        
        !! Free the memory:
        deallocate(work)
        deallocate(rwork)
        deallocate(WH)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Returns true if the given real array is sorted in the increasing order, false otherwise.
    !!**********************************************************************************************
    function is_sorted(array) result(sorted)
        real(wp) :: array(:)
        integer :: i
        logical :: sorted
        
        sorted = .true.
        do i = 2, size(array)
            if (array(i-1) > array(i)) then
                sorted = .false.
                exit
            end if
        end do
        
    end function
    
    !!**********************************************************************************************
    !! Check the validity of the given histogram. This subroutine stops the program when an error
    !! is detected in the given data.
    !!**********************************************************************************************
    subroutine check_histogram(bins, histo)
        real(wp), intent(in) :: bins(:), histo(:)
        integer :: nbins
        
        nbins = size(bins) - 1   !! Number of bins, number of intervals between points given by "bins".
        if (size(histo) /= nbins) then
            write (stderr, '(a)') tag_error // "Invalid histogram length &
                &(received ", size(histo), ", expected ", nbins, "), aborting..."
            stop errmsg_invalid_arg
        else if (.not. is_sorted(bins)) then
            write (stderr, '(a)') tag_error // "Unsorted array of bins, aborting..."
            stop errmsg_invalid_arg
        end if
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the histogram of the values containes in the array "samples" according to the intervals
    !! defined by "bins". The histogram "histo" is incremented accordingly. The number of missing samples
    !! is written in "missing". In principle, "missing" should be zero.
    !! Note that "histo" and "missing" are overwritten whichever their previous value.
    !! Note that the algorithm is written such that: size(samples) = sum(histo) + missing
    !!**********************************************************************************************
    subroutine histogram(samples, bins, histo, missing)
        real(wp), intent(in) :: samples(:), bins(:)
        real(wp), intent(out) :: histo(:), missing
        integer :: i, j, nbins
        
        call check_histogram(bins, histo)  !! Check for possible errors in the data.
        
        nbins = size(bins) - 1   !! Number of bins, or number of intervals between the edges given by "bins".
        missing = 0.0_wp   !! Initialize the histogram and the missing number to zero.
        histo = 0.0_wp
        do i = 1, size(samples)  !! Loop on the samples.
            if (samples(i) < bins(1) .or. bins(nbins+1) < samples(i)) then  !! The sample is out of the histogram, declared as missing.
                missing = missing + 1.0_wp 
            else if (samples(i) == bins(nbins+1)) then  !! The exceptional case where the sample is at the rightmost edge.
                histo(nbins) = histo(nbins) + 1.0_wp
            else
                do j = 1, nbins  !! Loop on the bins - Find the corresponding position in the histogram.
                    if (samples(i) < bins(j+1)) then
                        histo(j) = histo(j) + 1.0_wp
                        exit
                    end if
                end do
            end if 
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Normalize the given histogram "histo" supported by the corresponding "bins" so that the integral
    !! sum( histo * dx ) = 1. The bins are replaced by the centers.
    !! This subroutine overwrites "bins" and "histo" which afterwards can be respectively interpreted
    !! as the abscissas and the ordinates of the density function "rho(x)".
    !!**********************************************************************************************
    subroutine normalize_histogram(bins, histo)
        real(wp), intent(inout) :: bins(:), histo(:)
        real(wp) :: npoints
        integer :: i
        
        call check_histogram(bins, histo)  !! First check for possible errors in the data.
        
        npoints = sum(histo)  !! Compute the total of the given histogram.
        
        do i = 1, size(histo)
            histo(i) = histo(i)/(npoints*(bins(i+1) - bins(i)))
        end do
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute the moving average of the array "signal" over "n" points.
    !! The result is written in the array "res".
    !!**********************************************************************************************
    subroutine moving_average(signal, n, res)
        real(wp), intent(in) :: signal(:)
        integer, intent(in) :: n
        real(wp), intent(out) :: res(size(signal)-n+1)
        real(wp) :: stot
        integer :: i, j, nres
        
        !! Check for possible errors:
        nres = size(signal)-n+1
        if (n <= 1) then
            write (stderr, '(a,i0,a)') tag_error // "Invalid moving average parameter ", n, &
                ", expected >1, aborting..."
            stop errmsg_invalid_arg
        else if (size(res) /= nres) then
            write (stderr, '(a,i0,a)') tag_error // "Invalid moving average array length, &
                &received ", size(res), ", expected ", nres, ", aborting..."
            stop errmsg_invalid_arg
        end if
        
        do i = 1, nres  !! Compute the moving average.
            stot = 0.0_wp  !! First compute the total of local elements of signal().
            do j = 0, n-1
                stot = stot + signal(i+j)
            end do
            res(i) = stot/n  !! Normalize the total to get the moving average.
        end do
        
    end subroutine
    
	!!**********************************************************************************************
	!! Safely compare two real numbers and check if they are equal to within the given tolerance 'tol'.
	!! If not, prints them.
	!! This function returns True if the numbers are equal, False otherwise.
	!!**********************************************************************************************
	function assert_equal_real(x1, x2, tol) result(equal)
		real(wp), intent(in) :: x1, x2, tol
		logical :: equal
		equal = .true.
		if (abs(x1 - x2) > tol*(abs(x1) + abs(x2))/2) then
			print '(a,g0,a,g0,a)', tag_fail // "Numbers are too different: x1=", x1, ", and x2)=", x2, "."
			equal = .false.
		end if
	end function
	
	!!**********************************************************************************************
	!! Safely compare two complex numbers and check if they are equal to within the given tolerance 'tol'.
	!! If not, prints them.
	!! This function returns True if the numbers are equal, False otherwise.
	!!**********************************************************************************************
	function assert_equal_complex(z1, z2, tol) result(equal)
		complex(wp), intent(in) :: z1, z2
		real(wp), intent(in) :: tol
		logical :: equal
		equal = .true.
		if (abs(z1 - z2) > tol*(abs(z1) + abs(z2))/2) then
			print '(a,g0,sp,g0,ss,"i",a,g0,sp,g0,ss,"i",a)', &
                tag_fail // "Numbers are too different: z1=", z1, ", and z2)=", z2, "."
			equal = .false.
		end if
	end function
    
    !!**********************************************************************************************
	!! Subroutine to print a 1D real array to standard output.
	!!**********************************************************************************************
	subroutine print_real_array(strname, fmtreal, a)
		real(wp), intent(in) :: a(:)
		character(len=*), intent(in) :: strname, fmtreal
        character(len=100) :: style
        
        write (style, '(a,i0,a)') "(a,' = [',", size(a), "(1x," // fmtreal // ",1x),']')"
		write (stdout, style) trim(strname), a
        
	end subroutine
    
    !!**********************************************************************************************
    !! Subroutine to print a 2D real matrix to standard output.
    !!**********************************************************************************************
    subroutine print_real_matrix(strname, fmtreal, a)
        real(wp), intent(in) :: a(:, :)
        character(len=*), intent(in) :: strname, fmtreal
        integer :: i, j, nrow, ncol
		nrow = size(a, 1)
		ncol = size(a, 2)
		write (stdout, '(2a)') trim(strname), " = "
		do i = 1, nrow
			write (stdout, '(a)', advance='no') "|"
			do j = 1, ncol
				write (stdout, "(1x,"//fmtreal//",1x)", advance='no') a(i, j)
			end do
			write (stdout, '(a)') "|"
		end do
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to save a 2D real array to the file "fp" using the format "fmtreal" and the separator "sep".
	!!**********************************************************************************************
    subroutine save_real_matrix(fp, fmtreal, sep, a)
        integer, intent(in) :: fp
        character(len=*), intent(in) :: fmtreal, sep
        real(wp), intent(in) :: a(:, :)
        integer :: nrow, ncol
        character(len=100) :: style
        
        nrow = size(a, 1)
        ncol = size(a, 2)
        if (ncol == 1) then
            write (style, '(a,i0,a)') "(", nrow, "(" // trim(fmtreal) // ",/))"
        else
            write (style, '(a,i0,a,i0,a)') "(", nrow, "(", ncol-1, "(" // trim(fmtreal) // &
                ",'" // sep // "')," // trim(fmtreal) // ",/))"
        end if
        !!print '(a)', "Style = " // trim(style)
        write (fp, style) transpose(a)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Read an array of real number from the file "funit" and store it in "array".
    !! Ignore the lines devoid of a valid real number.
    !!**********************************************************************************************
    subroutine read_real_array(funit, array)
        integer, intent(in) :: funit
        real(wp), allocatable, intent(out) :: array(:)
        real(wp) :: buffer
        integer :: cnt, i, iost
        
        !! 1. First count the number of samples:
        cnt = 0  !! Initialize the count.
        do
            read (funit, *, iostat=iost) buffer
            if (iost == 0) then
                cnt = cnt + 1
            else if (iost < 0)    then
                exit   !! End of file detected.
            end if
        end do
        rewind(funit)   !! Rewind the file.
        !!print '(a,i0)', "Count of real numbers: ", cnt
        
        !! 2. Allocate and parse the data:
        allocate(array(cnt))
        
        i = 0
        do
            read (funit, *, iostat=iost) buffer
            if (iost == 0) then
                i = i + 1
                array(i) = buffer
            else if (iost < 0) then
                exit   !! End of file detected.
            end if
        end do
        rewind(funit)  !! Rewind the file to read other data if necessary.
        
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to write a timestamp at the beginning of the file "fp" with the given "prefix"
    !! (typically comment characters).
	!!**********************************************************************************************
    subroutine print_timestamp(fp, prefix)
        integer, intent(in) :: fp
        character(len=*), intent(in) :: prefix
        character(len=40) :: datetime_string, zone
        integer :: dt(8)
        
        call date_and_time(values=dt, zone=zone)
        write (datetime_string, '(i4,a,5(i0.2,a))') &
            dt(1), "-", dt(2), "-", dt(3), " at ", &
            dt(5), ":", dt(6), ":", dt(7), " " // trim(zone)
        
        write (fp, '(a)') trim(prefix) // " Computed on " // trim(datetime_string) // " by " // trim(program_copyright)
        
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to print a 1D complex array to standard output.
	!!**********************************************************************************************
	subroutine print_complex_array(strname, fmtreal, a)
		complex(wp), intent(in) :: a(:)
		character(len=*), intent(in) :: strname, fmtreal
        character(len=100) :: style
        
        write (style, '(a,i0,a)') "(a,' = [',", size(a), "(1x," // fmtreal // ",sp," // fmtreal // ",ss,'i',1x),']')"
		write (stdout, style) trim(strname), a
        
	end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to print a 2D complex array to standard output.
	!!**********************************************************************************************
	subroutine print_complex_matrix(strname, fmtreal, a)
		complex(wp), intent(in) :: a(:, :)
		character(len=*), intent(in) :: strname, fmtreal
		integer :: i, j, nrow, ncol
		nrow = size(a, 1)
		ncol = size(a, 2)
		write (*, '(2a)') trim(strname), " = "
		do i = 1, nrow
			write (*, '(a)', advance='no') "|"
			do j = 1, ncol
				write (*, "(1x,"//fmtreal//",sp,"//fmtreal//",ss,'i',1x)", advance='no') a(i, j)
			end do
			write (*, '(a)') "|"
		end do
	end subroutine
    
	!!**********************************************************************************************
	!! If the IO status flag 'istat' is zero, then tries to rewind the file.
	!! Otherwise, displays the invalid entry in the namelist.
	!! In principle, this routine should be called after each namelist reading.
	!!**********************************************************************************************
    subroutine check_and_rewind_nml(funit, istat)
		integer, intent(in) :: funit
		integer, intent(in) :: istat
		integer :: rewind_stat
		character(len=150) :: line
		if (istat == 0) then
			!!print '(a)', tag_info // "Rewinding the namelist file..."
			rewind (unit=funit, iostat=rewind_stat)
			if (rewind_stat /= 0) then
				print '(a)', tag_error // "Error in rewinding."
				stop errmsg_invalid_arg
			end if
		else
			backspace(funit)  !! Go back one line.
			read (funit, '(a)') line
			write (stderr, '(a)') tag_error // "Invalid entry in namelist: " // trim(line) // "..."
			if (line == "/") then
				write (stderr, '(a)') tag_warn // "Try checking the spelling of the namelist groups..."
			end if
			stop errmsg_invalid_arg
		end if
	end subroutine
    
    !!**********************************************************************************************
    !! Trim the spaces and zeros at the end of the string, and the spaces at the beginning of the string.
    !! This function does not need to be passed to trim().
    !!**********************************************************************************************
    function trim_zero(string) result(res)
        character(len=*), intent(in) :: string
        character(len=:), allocatable :: res
        integer :: imin, imax
        imin = 1
        do while (string(imin:imin) == " ")
            imin = imin + 1
        end do
        imax = len(string)
        do while (string(imax:imax) == "0" .or. string(imax:imax) == " ")
            imax = imax - 1
        end do
        res = string(imin:imax)
    end function
    
    !!**********************************************************************************************
    !! Separate the path from the basename of the given path name.
    !! The directory is written to the string "directory", and the basename to "basename".
    !! This subroutine assumes the path separator given by "pathsep" (see global variables above).
    !!**********************************************************************************************
    subroutine path_split(pathname, directory, basename)
        character(len=*), intent(in) :: pathname
        character(len=*), intent(out) :: directory, basename
        integer :: i
        
        i = len(pathname)
        do while (pathname(i:i) /= pathsep)
            i = i - 1
        end do
        directory = pathname(:i)
        basename = pathname(i+1:)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to generate a new unique filename by checking the availability of the filename.
	!! The template is "base_x.ext" where "x" is a number from 1 to MAX_INT.
	!! If the file already exists, then increments the index "x" until the filename does not exist.
	!! The purpose is of course to avoid overwriting existing data.
	!! base     = Basename of the file, including the folder.
	!! ext      = File extension, including the dot prefix.
	!! filename = Resulting output filename, guaranteed to exist if no error is sent.
	!!**********************************************************************************************
	subroutine unique_filename(base, ext, filename)
		character(len=*), intent(in) :: base, ext
		character(len=*), intent(out) :: filename
		integer, parameter :: nmax = 1000000000  !! Large value still representable on int32.
		integer :: i
		logical :: fexist
		do i = 1, nmax
			write (filename, '(a,i0,a)') trim(base), i, trim(ext)
			!!print '(a)', tag_info // "Testing filename '" // filename // "'..."
			inquire(file=filename, exist=fexist)
			if (.not. fexist) return
		end do
		!! The loop should never reach this point:
		write (stderr, '(a,i0,a)') tag_error // "Too many files (nmax=", nmax, "), no valid filename !"
		stop errmsg_invalid_arg
	end subroutine
	
	!!**********************************************************************************************
	!! Check if the given output directory "dirname" exists. If not, it will be created.
    !! WARNING : This subroutine requires the "mkdir" command !
	!!**********************************************************************************************
	subroutine create_directory(dirname)
		character(len=*), intent(in) :: dirname
		call execute_command_line("mkdir -p " // trim(dirname))
	end subroutine
    
    !!**********************************************************************************************
    !! Copy the preamble (namely the first commented lines) of the input file "inputfp" to the given
    !! output file "outputfp". These files are supposed to the opened in reading and writing modes,
    !! respectively.
    !!**********************************************************************************************
    subroutine copy_preamble(inputfp, outputfp)
        integer, intent(in) :: inputfp, outputfp
        character(len=300) :: line  !! Only 300 characters per line are copied.
        integer :: iost
        
        do
            read (inputfp, '(a)', iostat=iost) line
            !!print '(3a)', "Line = '", trim(line), "'..."
            if (iost == 0 .and. (line(1:1) == '%' .or. line(1:1) == '#')) then
                write (outputfp, '(a)') trim(line)
            else
                exit   !! End of file detected.
            end if
        end do
        
    end subroutine
    
    !!**********************************************************************************************
	!! Returns the given time duration "time" (in seconds) in a string of the format "DDDd HH:MM:SS".
	!!**********************************************************************************************
    function format_time(time) result(str)
        real(wp), intent(in) :: time   !! Time duration in seconds.
        character(len=30) :: str
        real(wp) :: t
        integer :: day, hour, minute, second
        
        day = floor(time/86400)
        t = time - 86400*day
        hour = floor(t/3600)
        t = t - 3600*hour
        minute = floor(t/60)
        second = floor(t - 60*minute)
        
        if (day == 0) then  !! If day=0, then print a short version:
            write (str, '(i0.2,":",i0.2,":",i0.2)') hour, minute, second
        else  !! If day is not zero, then print a longer version:
            write (str, '(i0,"d ",i0,"h ",i0,"m ",i0,"s")') day, hour, minute, second
        end if
        
    end function
	
    !!**********************************************************************************************
    !! Prints the progress bar to standard output with the expected time of arrival.
    !! cjob  = Index of the current completed job [1, ..., njob].
    !! njob  = Total number of jobs.
    !! start = Date and time array given by the subroutine: date_and_time(values=start).
    !! msg   = Short message printed before the progress bar.
    !!**********************************************************************************************
    subroutine print_progress_bar(cjob, njob, start, msg)
        integer, intent(in) :: cjob, njob, start(8)
        character(len=*), intent(in) :: msg
        integer :: nchar, now(8)
        real(wp) :: delta(8), time, speed, eta
        
        !! 1. Compute the expected time of arrival:
        call date_and_time(values=now)
        delta = real(now - start, kind=wp)
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Elapsed time in seconds.
        speed = cjob/time           !! Compute the mean speed in job per second.
        eta = (njob - cjob)/speed   !! Compute the expected time of arrival (in seconds).
        
        !! 2. Print the progress bar:
        nchar = (progress_bar_length*cjob)/njob   !! Integer division.
        write(stdout, '(a,f0.2,a,f0.1,a)', advance='no') &
            tag_info // trim(msg) // " [" // repeat("#", nchar) // &
            repeat(" ", progress_bar_length - nchar) // "] ", 100.0*cjob/njob, "% ETA " // &
            trim(format_time(eta)) // " (", speed, " job/s)" // char_cr
        
    end subroutine
    
    !!**********************************************************************************************
    !! Ends the progress bar and computes the total computation time.
    !!**********************************************************************************************
    subroutine end_progress_bar(start, time)
        integer, intent(in) :: start(8)
        real(wp), intent(out) :: time
        integer :: time_end(8)
        real(wp) :: delta(8)
        
        call date_and_time(values=time_end)    !! End the time measurement.
        delta = real(time_end - start, kind=wp)
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Computation time in seconds.
        print '(/,a,f0.2,a)', tag_info // "Done, computation time is ", time, " s."
        
    end subroutine
    
end module
