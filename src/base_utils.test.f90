!!**************************************************************************************************
!! Created on 2024-05-16 at 10:13:45 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran program to test the procedures of the base_utils module.
!!**************************************************************************************************
program base_utils_test
    use base_utils
	implicit none
	
	!call test_add_band_matrix()
	!call test_diag_mul_left()
	!call test_diag_mul_right()
    !call test_inverse_complex()
    !call test_product_complex()
    !call perf_matmul_vs_zgemm()
    !call test_matmul_vector()
    !call test_singular_values()
    !call test_diagonalize_ahca()
    !call test_norm2_complex()
    !call test_histogram()
    !call test_moving_average()
    !call test_save_real_matrix()
    !call test_system_clock()
    !call test_format_time()
    !call test_progress_bar()
    !call test_solve_block_tridiag()
    !call test_read_real_array()
    !call test_path_split()
    
	contains
    
    !!**********************************************************************************************
    !! Test the norm2_complex() function.
    !!**********************************************************************************************
    subroutine test_norm2_complex()
        complex(wp) :: A(4, 3)
        
        A = reshape([   ( 0.56822371494419420_wp,  0.25427911592041810_wp), &
                        ( 0.92842629910986170_wp,  0.41627873454887654_wp), &
                        ( 0.28171163164123714_wp, -0.45353926993390250_wp), &
                        (-0.24135095684205732_wp, -0.32814750945176563_wp), &
                        ( 0.34256543681552910_wp, -0.45025306988077496_wp), &
                        (-0.30110268435099340_wp,  0.03998748001092878_wp), &
                        ( 0.63042919841714260_wp, -0.46963146547886137_wp), &
                        ( 0.90685286456350810_wp,  0.12375768017873767_wp), &
                        ( 0.10100461913328873_wp, -0.93995673371119360_wp), &
                        (-0.13628841781638990_wp,  0.94900245310786870_wp), &
                        (-0.52576497526710990_wp,  0.97769098014051450_wp), &
                        (-0.05815140897705184_wp,  0.45950148578633376_wp) ], shape(A), order=[2, 1])
        
        call print_complex_matrix("A", "f11.6", A)
        print '(a,g0.16)', "Norm_calc = ", norm2_complex(A)
        print '(a,g0.16)', "Norm_expc = 2.646045307275178"
        
        !!print '(a,g0.16,"   ",g0.16,"i")', "Calc = ", sum(A)
        !!print '(a)', "Expc = 2.4965553213711593 + 0.57896988123718 I"
        
    end subroutine
	
	!!**********************************************************************************************
	!! Subroutine to test the add_band_matrix() subroutine.
	!!**********************************************************************************************
	subroutine test_add_band_matrix()
		complex(wp) :: matrix(4, 4), newmatrix(4, 4), diff(4, 4)
        real(wp) :: elems(7)
        
        print '(a)', "====== TEST ADD_BAND_MATRIX ======"
        
        matrix = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], shape(matrix), order=[2, 1])
        newmatrix = matrix 
        
        elems = [1, 2, 3, 4, 5, 6, 7]
        
        call add_band_matrix(elems, newmatrix)
        !!call print_complex_matrix(strname, "f9.6", matrix)
        !!call print_complex_matrix(strname, "f9.6", newmatrix)
        
        diff = newmatrix - matrix
        
        call print_complex_matrix("Diff", "f9.6", diff)
        
	end subroutine
	
	!!**********************************************************************************************
	!! Subroutine to test the diag_mul_left() subroutine.
	!!**********************************************************************************************
    subroutine test_diag_mul_left()
        complex(wp) :: D(3), M(3, 4), P(3, 4)
        
        print '(a)', "====== TEST DIAG_MUL_LEFT ======"
        
        D = [3, 2, 6]
        M = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], shape(M), order=[2, 1])
        
        P = diag_mul_left(D, M)
        
        call print_complex_array("D",  "f11.6", D)
        call print_complex_matrix("M", "f11.6", M)
        call print_complex_matrix("P", "f11.6", P)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the diag_mul_right() subroutine.
	!!**********************************************************************************************
    subroutine test_diag_mul_right()
        complex(wp) :: D(4), M(3, 4), P(3, 4)
        
        print '(a)', "====== TEST DIAG_MUL_RIGHT ======"
        
        D = [3, 2, 6, -1]
        M = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], shape(M), order=[2, 1])
        
        P = diag_mul_right(M, D)
        
        call print_complex_array("D",  "f11.6", D)
        call print_complex_matrix("M", "f11.6", M)
        call print_complex_matrix("P", "f11.6", P)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the matrix product product_complex() based on LAPACK's zgemm().
	!!**********************************************************************************************
    subroutine test_product_complex()
        complex(wp) :: A(3, 5), B(5, 4), C(3, 4), Cexpc(3, 4)
        
        print '(a)', "====== TEST PRODUCT_COMPLEX ======"
        
        A = reshape([   (-0.6915716510885130_wp,  0.0724952835238364_wp), &
                        (-0.4611786855999176_wp,  0.9638526711679933_wp), &
                        ( 0.9903086181688163_wp, -0.1368014057397304_wp), &
                        (-0.6602774660990369_wp,  0.3259633376548603_wp), &
                        (-0.4566913352355613_wp,  0.7892413430657932_wp), &
                        (-0.3243387663117661_wp, -0.2342522243843050_wp), &
                        ( 0.9635405020462309_wp,  0.5250309478303370_wp), &
                        ( 0.6323863158906908_wp, -0.7102752820607163_wp), &
                        ( 0.8236099856672578_wp,  0.0542580642984696_wp), &
                        ( 0.6309210359052750_wp,  0.4360210588154851_wp), &
                        ( 0.5117899727591046_wp, -0.6694183019815863_wp), &
                        ( 0.7973832988008160_wp,  0.8866004390044604_wp), &
                        ( 0.6478776288753254_wp, -0.7247529625863822_wp), &
                        ( 0.6910756943080449_wp, -0.1917028464836923_wp), &
                        ( 0.9785817857483745_wp,  0.6165605000093395_wp) ], shape(A), order=[2, 1])
        
        B = reshape([   (-0.7475807692479659_wp,  0.1051578325333438_wp), &
                        (-0.9931569706498324_wp,  0.2366981710525344_wp), &
                        ( 0.8510930676672674_wp,  0.5659441221064094_wp), &
                        (-0.7780243784557905_wp, -0.2023095236700652_wp), &
                        ( 0.3415283086102203_wp,  0.5516632147067315_wp), &
                        (-0.1912587576667275_wp,  0.8783103842511961_wp), &
                        ( 0.7484072572523148_wp,  0.8166301415504860_wp), &
                        (-0.8471326576677565_wp, -0.6580094874907068_wp), &
                        (-0.0554930340606696_wp, -0.1940850004385273_wp), &
                        ( 0.3731799812737075_wp, -0.8561060886847955_wp), &
                        (-0.3940467289994212_wp,  0.7986349768734584_wp), &
                        ( 0.7822440709948602_wp,  0.9165646189069192_wp), &
                        (-0.9950226371623558_wp, -0.7444521472006715_wp), &
                        ( 0.0599025932743462_wp, -0.9386063632869499_wp), &
                        (-0.7189472021908725_wp,  0.3224448474332870_wp), &
                        ( 0.9816387897189496_wp,  0.0752499638876500_wp), &
                        (-0.8813070978758768_wp, -0.3858252959621060_wp), &
                        (-0.7301025200135407_wp, -0.7859279309089686_wp), &
                        ( 0.3043604777091882_wp, -0.5490757112232876_wp), &
                        (-0.3644357864265317_wp,  0.7560407680115855_wp) ], shape(B), order=[2, 1])
        
        C = product_complex(A, B)
        Cexpc = matmul(A, B)
        
        !!call print_complex_matrix("A", "f11.6", A)
        !!call print_complex_matrix("B", "f11.6", B)
        
        call print_complex_matrix("C", "f20.16", C)
        call print_complex_matrix("C (expected)", "f20.16", Cexpc)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to compare the performance of the matrix products given by matmul() and LAPACK's zgemm().
    !! It seems that matmul() is faster.
	!!**********************************************************************************************
    subroutine perf_matmul_vs_zgemm()
        integer, allocatable :: seedarr(:)
        complex(wp), allocatable :: A(:, :), B(:, :), C(:, :)
        real(wp), allocatable :: randseq(:)
        integer, parameter :: n = 200      !! Size of the matrix that will be multiplied.
        integer, parameter :: nloop = 500  !! Number of times the matrix product is done to get an accurate time estimate.
        integer :: i, nseedarr
        real(wp) :: matmul_start, matmul_end, zgemm_start, zgemm_end
        
        print '(a)', "====== PERF MATMUL VS ZGEMM ======"
        
        !! Initialize the random generator to get reproducible results:
        call random_seed(size=nseedarr)  !! Get the expected length of the seed (8 * 32-bit integers = 256 bits for xoshiro256).
        allocate(seedarr(nseedarr))
        seedarr = 1
        call random_seed(put=seedarr)  !! Initialize the seed of the random generator.
        
        !! Prepare the matrices:
        allocate(randseq(2*n*n))
        call random_number(randseq)  !! Generate all the random numbers at once (xoshiro256 for GNU Fortran).
        allocate(A(n, n))
        allocate(B(n, n))
        allocate(C(n, n))
        A = reshape(randseq(1:(n*n)), shape(A))
        B = reshape(randseq((n*n+1):(2*n*n)), shape(B))
        
        !call print_complex_matrix("A", "f11.6", A)
        !call print_complex_matrix("B", "f11.6", B)
        
        call cpu_time(matmul_start)
        do i = 1, nloop
            C = matmul(A, B)
        end do
        call cpu_time(matmul_end)
        
        call cpu_time(zgemm_start)
        do i = 1, nloop
            C = product_complex(A, B)
        end do
        call cpu_time(zgemm_end)
        
        print '(a,f12.6,a)', tag_info // "matmul time  = ", matmul_end-matmul_start, " s"
        print '(a,f12.6,a)', tag_info // "zgemm  time  = ", zgemm_end-zgemm_start, " s"
        print '(a,f12.6,a)', tag_info // "matmul/zgemm = ", (matmul_end-matmul_start)/(zgemm_end-zgemm_start), " s"
        
        !! Deallocate memory:
        deallocate(randseq)
        deallocate(A)
        deallocate(B)
        deallocate(C)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the matrix multiplication over a vector.
    !!**********************************************************************************************
    subroutine test_matmul_vector()
        real(wp) :: A(3, 3), B(3), E(3)
        
        print '(a)', "====== TEST MATMUL_VECTOR ======"
        
        A = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp], shape(A), order=[2, 1])
        
        B = [2.0_wp, 0.0_wp, -0.5_wp]
        
        E = [0.5_wp, 5.0_wp, 9.5_wp]
        
        call print_real_matrix("A", "f0.3", A)
        call print_real_array("B", "f0.3", B)
        call print_real_array("Product A*B", "f0.3", matmul(A, B))
        call print_real_array("Expected___", "f0.3", E)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the inversion of complex matrix.
	!!**********************************************************************************************
    subroutine test_inverse_complex()
        complex(wp) :: A(4, 4), Ainv(4, 4)
        
        print '(a)', "====== TEST INVERSE_COMPLEX ======"
        
        A = reshape([   (-0.28059215604538010_wp, -0.3461434824696745_wp), &
                        ( 0.18473331224680223_wp,  0.6318812387085075_wp), &
                        (-0.40746040147950290_wp,  0.1382727026518729_wp), &
                        ( 0.87389608120820170_wp, -0.2661422854657700_wp), &
                        ( 0.95491046318197360_wp, -0.2105937134153009_wp), &
                        ( 0.16813499164823710_wp,  0.0418637468832905_wp), &
                        ( 0.28781284790753325_wp, -0.4760898733399750_wp), &
                        ( 0.85653516925394910_wp,  0.4349539392523778_wp), &
                        ( 0.01216264260162925_wp,  0.3239628778542527_wp), &
                        (-0.19292937387052378_wp, -0.9451826160100225_wp), &
                        (-0.28724113783282945_wp, -0.7265218262414725_wp), &
                        ( 0.03149998011019228_wp,  0.4644297650159581_wp), &
                        (-0.82541754199947800_wp,  0.9077152529671042_wp), &
                        (-0.83414745798494080_wp, -0.3596418431929087_wp), &
                        ( 0.51455340226931460_wp, -0.6296860149451851_wp), &
                        ( 0.14143895745785207_wp,  0.9369393180454373_wp) ], shape(A), order=[2, 1])
        
        Ainv = A
        call inverse_complex(Ainv)  !! Call the inversion subroutine (based on LAPACK).
        call print_complex_matrix("A^-1 * A", "f20.16", matmul(Ainv, A))
        
    end subroutine
    
    !!**********************************************************************************************
	!! Subroutine to test the solution of a block-tridiagonal system returned by the solve_block_tridiag() subroutine.
	!!**********************************************************************************************
    subroutine test_solve_block_tridiag()
        complex(wp) :: amat(3, 2, 2), bdiag(2), imat(3, 2, 1)
        
        print '(a)', "====== TEST SOLVE_BLOCK_TRIDIAG ======"
        
        amat = reshape([(-0.3450721935254197_wp,  0.5392304997496717_wp), &
                        ( 0.6897836443224770_wp,  0.5412262828699856_wp), &
                        ( 0.3231770298367444_wp,  0.0785420140112380_wp), &
                        ( 0.3128024256862249_wp,  0.0550528538949643_wp), &
                        ( 0.3334441501416232_wp,  0.8645183075981557_wp), &
                        (-0.7888526093217032_wp, -0.8561361200407882_wp), &
                        ( 0.9215704435800176_wp, -0.1313107808488722_wp), &
                        (-0.9883756365893341_wp,  0.7297537414237123_wp), &
                        (-0.8048043239128400_wp,  0.9277252313541604_wp), &
                        ( 0.9113402180562988_wp, -0.2620342537740332_wp), &
                        ( 0.5449654877189336_wp,  0.1398942224192785_wp), &
                        (-0.2826158268607601_wp, -0.4251008049348050_wp) ], shape(amat))
        
        imat = reshape([(-0.1832660166300050_wp, -0.9263008028015078_wp), &
                        (-0.9377593996701266_wp,  0.9927737116226631_wp), &
                        ( 0.8876721334230475_wp, -0.0677127726588384_wp), &
                        ( 0.9551250020205959_wp,  0.9161932570770022_wp), &
                        ( 0.4131234542792588_wp, -0.3505411051103899_wp), &
                        ( 0.0751291445949542_wp,  0.8756862438708004_wp) ], shape(imat))
        
        !imat = reshape([(-0.1832660166300050_wp, -0.9263008028015078_wp), &
        !                (-0.9377593996701266_wp,  0.9927737116226631_wp), &
        !                ( 0.8876721334230475_wp, -0.0677127726588384_wp), &
        !                ( 0.9551250020205959_wp,  0.9161932570770022_wp), &
        !                ( 0.4131234542792588_wp, -0.3505411051103899_wp), &
        !                ( 0.0751291445949542_wp,  0.8756862438708004_wp), &
        !                ( 0.0842531410118012_wp, -0.8494624480371211_wp), &
        !                ( 0.9596497570040681_wp, -0.1558778671431713_wp), &
        !                (-0.1424153566210044_wp,  0.6073310899414994_wp), &
        !                ( 0.5869594774563454_wp,  0.0299935357770002_wp), &
        !                ( 0.1310804335743585_wp,  0.4016149175953827_wp), &
        !                ( 0.7323816310625264_wp,  0.0283754403548432_wp) ], shape(imat))
        
        bdiag = [2.0_wp, 1.0_wp]
        
        call test_solve_block_tridiag_with(amat, bdiag, imat)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the solve_block_tridiag() subroutine with specific arguments.
    !! This subroutine uses LAPACK to compare with the solution of the system.
    !!**********************************************************************************************
    subroutine test_solve_block_tridiag_with(amat, bdiag, imat)
        complex(wp), intent(in) :: amat(:, :, :), bdiag(size(amat, 2)), imat(:, :, :)
        complex(wp) :: xmat(size(imat, 1), size(imat, 2), size(imat, 3))
        complex(wp) :: x(size(imat, 1)*size(imat, 2), size(imat, 3)), x_expc(size(imat, 1)*size(imat, 2), size(imat, 3))
        complex(wp) :: h(size(amat, 1)*size(amat, 2), size(amat, 1)*size(amat, 2)), ind(size(imat, 1)*size(imat, 2), size(imat, 3))
        integer :: i, j, imin, imax, nb, bs
        
        nb = size(amat, 1)  !! Number of blocks.
        bs = size(amat, 2)  !! Size of blocks.
        
        !! 1. First solve with the custom subroutine:
        call solve_block_tridiag(amat, bdiag, imat, xmat)
        do i = 1, nb  !! Reshape "xmat" into a matrix to compare with LAPACK's result.
            imin = (i-1)*bs + 1  !! Minimum index.
            imax = i*bs          !! Maximum index.
            x(imin:imax, :) = xmat(i, :, :)  !! Prepare the solution for comparison with LAPACK.
            ind(imin:imax, :) = imat(i, :, :)  !! Prepare the independent term for LAPACK.
            h(imin:imax, imin:imax) = amat(i, :, :)  !! Prepare the tridiagonal matrix for LAPACK.
            if (i < nb) then   !! Exclude the last step (out-of-bound memory access).
                do j = 1, bs   !! Fill the superdiagonal and subdiagonal.
                    h((i-1)*bs + j, i*bs + j) = bdiag(j)
                    h(i*bs + j, (i-1)*bs + j) = bdiag(j)
                end do
            end if
        end do
        
        call print_complex_matrix("I_1", "f0.3", imat(1, :, :))
        call print_complex_matrix("I_2", "f0.3", imat(2, :, :))
        call print_complex_matrix("I_3", "f0.3", imat(3, :, :))
        call print_complex_matrix("I", "f0.3", ind)
        
        call print_complex_matrix("A_1", "f0.3", amat(1, :, :))
        call print_complex_matrix("A_2", "f0.3", amat(2, :, :))
        call print_complex_matrix("A_3", "f0.3", amat(3, :, :))
        call print_complex_matrix("H", "f0.3", h)
        
        !! 2. Solve the system with LAPACK:
        call inverse_complex(h)
        x_expc = matmul(h, ind)
        
        call print_complex_matrix("X computed", "f0.3", x)
        call print_complex_matrix("X expected", "f0.3", x_expc)
        print '(a,g0.3)', "Relative error = ", norm2_complex(x - x_expc)/norm2_complex(x_expc)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the singular values computed by LAPACK's zgesvd() subroutine.
	!!**********************************************************************************************
    subroutine test_singular_values()
        complex(wp) :: A(4, 3)
        real(wp) :: S(3), S_expc(3)
        
        print '(a)', "====== TEST SINGULAR_VALUES ======"
        
        A = reshape([   (-0.01283288021173989_wp,  0.9377239102900901_wp), &
                        (-0.95242616135770230_wp, -0.8994963544167294_wp), &
                        (-0.06886775632719067_wp,  0.6636594333362238_wp), &
                        ( 0.70466932341279700_wp,  0.7651115855726442_wp), &
                        (-0.89922794465691420_wp, -0.5377498780142429_wp), &
                        ( 0.13723225079543510_wp, -0.5763274633716735_wp), &
                        (-0.22856244103566636_wp,  0.1867127015776933_wp), &
                        (-0.70086940334971310_wp, -0.8956086270489085_wp), &
                        (-0.98866683961216850_wp,  0.9166435318614940_wp), &
                        ( 0.89492167606273610_wp, -0.4782826178810580_wp), &
                        ( 0.82361640475279610_wp,  0.3462092784006585_wp), &
                        (-0.84349089333020190_wp, -0.8929737044010584_wp) ], shape(A), order=[2, 1])
        
        S_expc = [2.9024423080155306_wp, 1.860396504848104_wp, 0.47034257779202826_wp]
        
        call print_complex_matrix("A", "f20.16", A)
        call singular_values(A, S)
        call print_real_array("S_calc", "f20.16", S)
        call print_real_array("S_expc", "f20.16", S_expc)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Subroutine to test the diagonalization (eigendecomposition) of the matrix A^H * A resorting to
    !! singular value decomposition in LAPACK.
    !!**********************************************************************************************
    subroutine test_diagonalize_ahca()
        complex(wp) :: A(4, 3), A_copy(4, 3), U(3, 3)
        real(wp) :: T(3), T_expc(3)
        
        print '(a)', "====== TEST DIAGONALIZE_AHA ======"
        
        A = reshape([   (-0.01283288021173989_wp,  0.9377239102900901_wp), &
                        (-0.95242616135770230_wp, -0.8994963544167294_wp), &
                        (-0.06886775632719067_wp,  0.6636594333362238_wp), &
                        ( 0.70466932341279700_wp,  0.7651115855726442_wp), &
                        (-0.89922794465691420_wp, -0.5377498780142429_wp), &
                        ( 0.13723225079543510_wp, -0.5763274633716735_wp), &
                        (-0.22856244103566636_wp,  0.1867127015776933_wp), &
                        (-0.70086940334971310_wp, -0.8956086270489085_wp), &
                        (-0.98866683961216850_wp,  0.9166435318614940_wp), &
                        ( 0.89492167606273610_wp, -0.4782826178810580_wp), &
                        ( 0.82361640475279610_wp,  0.3462092784006585_wp), &
                        (-0.84349089333020190_wp, -0.8929737044010584_wp) ], shape(A), order=[2, 1])
        
        T_expc = [8.42417135135852_wp, 3.461075155251041_wp, 0.22122214048405014_wp]
        
        call print_complex_matrix("A", "f20.16", A)
        
        A_copy = A
        call diagonalize_ahca(A_copy, T, U)
        
        call print_real_array("T_calc", "f20.16", T)
        call print_real_array("T_expc", "f20.16", T_expc)
        
        call print_complex_matrix("Error, A^H * A - U * T * U^H", "g0.3", &
            matmul(transpose(conjg(A)), A) - matmul(U, diag_mul_left(cmplx(T, kind=wp), transpose(conjg(U)))))
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the histogram() subroutine.
	!!**********************************************************************************************
    subroutine test_histogram()
        real(wp) :: samples(12), bins(5)
        real(wp) :: histo(4), histo_expc(4), missing
        
        print '(a)', "====== TEST ADD_HISTOGRAM ======"
        
        samples = [  0.21_wp,  0.00_wp,  0.26_wp, -0.22_wp, &
                    -1.00_wp,  0.77_wp,  1.00_wp, -0.37_wp, &
                    -0.43_wp, -0.80_wp, -0.17_wp,  0.50_wp  ]
        
        bins = [ -1.0_wp, -0.5_wp, 0.0_wp, 0.5_wp, 1.0_wp ]
        
        histo_expc = [ 2.0_wp, 4.0_wp, 3.0_wp, 3.0_wp ]
        
        call histogram(samples, bins, histo, missing)
        call print_real_array("Histo_calc", "f7.2", histo)
        call print_real_array("Histo_expc", "f7.2", histo_expc)
        print '(a,f7.2)', "Missing  = ", missing 
        print '(a,f7.2,a,i0)', "Checksum = ", sum(histo) + missing, " should be ", size(samples)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the moving_average() subroutine.
	!!**********************************************************************************************
    subroutine test_moving_average()
        real(wp) :: signal(6), res(5), res_expc(5)
        integer, parameter :: n = 2
        
        signal = [ 1.0_wp, 5.0_wp, -3.0_wp, 0.0_wp, 10.0_wp, -4.0_wp ]
        res_expc = [ 3.0_wp, 1.0_wp, -1.5_wp, 5.0_wp, 3.0_wp ]
        
        call moving_average(signal, n, res)
        
        call print_real_array("Res_calc", "f7.2", res)
        call print_real_array("Res_expc", "f7.2", res_expc)
        
    end subroutine
    
	!!**********************************************************************************************
	!! Subroutine to test the save_real_matrix() subroutine.
	!!**********************************************************************************************
    subroutine test_save_real_matrix()
        real(wp) :: A(3, 4)
        character(len=*), parameter :: fmtreal = "f0.1", sep = ", "
        
        A = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], shape(A), order=[2, 1])
        
        call print_complex_matrix("A", fmtreal, cmplx(A, kind=wp))
        call save_real_matrix(stdout, fmtreal, sep, A)
        
    end subroutine
    
	!!**********************************************************************************************
    !! Test the measurement of the wall time in Fortran. This is useful to compare the true
    !! performance of various programs and algorithms.
	!!**********************************************************************************************
    subroutine test_system_clock()
        integer :: values_start(8), values_end(8), i
        real(wp) :: delta(8), time, x(100000000)
        
        call date_and_time(values=values_start)
        write (stdout, '(a,"[",8(1x,i0,1x),"]")') "Start = ", values_start
        
        do i = 1, 100000000
            x(i) = sin(real(i))
        end do
        
        call date_and_time(values=values_end)
        delta = real(values_end - values_start)
        write (stdout, '(a,"[",8(1x,f0.1,1x),"]")') "Delta = ", delta
        
        !! Compute the time in seconds:
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000
        write (stdout, '(a,f0.3,a)') "Time = ", time, "s"
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the format_time() subroutine.
    !!**********************************************************************************************
    subroutine test_format_time()
        real(wp) :: time
        time = 578965.5225_wp
        print '(a,f12.2,a)', "Time = ", time, " -> " // trim(format_time(time))
        time = 25418.0_wp
        print '(a,f12.2,a)', "Time = ", time, " -> " // trim(format_time(time))
        time = 569.0_wp
        print '(a,f12.2,a)', "Time = ", time, " -> " // trim(format_time(time))
    end subroutine
    
    !!**********************************************************************************************
    !! Test the progress bar.
    !!**********************************************************************************************
    subroutine test_progress_bar()
        integer :: cjob, start_time(8), end_time(8), i
        integer, parameter :: njob = 85120, joblen = 75000
        character(len=*), parameter :: msg = "Test"
        real(wp) :: job(joblen), delta(8), time
        
        job = [( real(i, kind=wp), i = 1, joblen )]
        
        call date_and_time(values=start_time)
        
        do cjob = 1, njob
            job = sin(job)
            if (mod(cjob, 1000) == 0 .or. cjob == njob) then
                call print_progress_bar(cjob, njob, start_time, msg)
            end if
        end do
        print *
        
        call date_and_time(values=end_time)
        delta = real(end_time - start_time, kind=wp)
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Elapsed time in seconds.
        print '(a,f0.2,a)', tag_info // "Final elapsed time: ", time, " s"
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the reading of a real array.
    !!**********************************************************************************************
    subroutine test_read_real_array()
        integer :: funit
        character(len=*), parameter :: fname = "test/test_data.csv"
        real(wp), allocatable :: array(:)
        
        print '(a)', "====== TEST READ_REAL_ARRAY ======"
        
        open(newunit=funit, file=fname, action='read')
        print '(a)', tag_info // "Reading test file: '" // trim(fname) // "'..."
        call read_real_array(funit, array)
        call print_real_array("array", "g0", array)
        
        close(funit)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Test the path_split() subroutine.
    !!**********************************************************************************************
    subroutine test_path_split()
        character(len=*), parameter :: pathname = "/home/user/Documents/path/to/my/file.dot.ext"
        character(len=50) :: directory, basename
        
        print '(a)', "====== TEST PATH_SPLIT ======"
        
        call path_split(pathname, directory, basename)
        print '(a)', tag_info // "pathname  = '" // trim(pathname)  // "'"
        print '(a)', tag_info // "directory = '" // trim(directory) // "'"
        print '(a)', tag_info // "basename  = '" // trim(basename)  // "'"
        
    end subroutine
    
end program
