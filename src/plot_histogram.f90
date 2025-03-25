!!**************************************************************************************************
!! Created on 2024-07-10 at 11:19:33 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the MIT License.
!! Fortran program to plot or replot the histogram corresponding to the list of samples given by the
!! data file in argument. One sample is given by line.
!!**************************************************************************************************
program plot_histogram
    use base_utils
	implicit none
	character(len=150) :: configfname, datafname  !! This assumes that filenames are not too long.
	integer :: configfp, datafp
	logical :: fexist
    
	write (stdout, '(a)') "====== This is PlotHisto (from RecurGreen) " // trim(copyright) // " ======"
	if (command_argument_count() /= 2) then
		write (stdout, '(a)') tag_info // "Usage: $ ./plothisto plothisto_settings.nml path/to/samples.csv"
        write (stderr, '(a)') tag_error // "Two files expected, doing nothing..."
	end if
    
    !! Check the configuration file:
	call get_command_argument(1, configfname)
	inquire(file=configfname, exist=fexist)
    if (.not. fexist) then
        write (stderr, '(a)') tag_error // "File not found: " // trim(configfname)
        stop errmsg_invalid_arg
    end if
    
    !! Check the data file:
    call get_command_argument(2, datafname)
    inquire(file=datafname, exist=fexist)
    if (.not. fexist) then
        write (stderr, '(a)') tag_error // "File not found: " // trim(datafname)
        stop errmsg_invalid_arg
    end if
    
    open(newunit=configfp, file=configfname, action='read')
    open(newunit=datafp, file=datafname, action='read')
    
    call plot_histo_from_file(configfp, datafp)
    
    close(configfp)
    close(datafp)
    
	contains
	
	!!**********************************************************************************************
	!! Subroutine to plot (or replot) the histogram corresponding to the data in the given file "datafp".
    !! The settings are given in the file "configfp".
	!!**********************************************************************************************
	subroutine plot_histo_from_file(configfp, datafp)
		integer, intent(in) :: configfp, datafp
        real(wp), allocatable :: samples(:), bins(:), histo(:, :)
        character(len=50) :: bintype
        real(wp) :: xmin, xmax, missing
        integer :: nsamp, nbins, i, istat
        namelist /settings/ bintype, nbins, xmin, xmax
        
        !! 1. Read the data from the files:
        read(unit=configfp, nml=settings, iostat=istat)
		call check_and_rewind_nml(configfp, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        call read_real_array(datafp, samples)  !! Read thes samples from the file.
        nsamp = size(samples)   !! Determine the number of samples.
        print '(a,i0,a)', tag_info // "Found ", nsamp, " samples in the data file..."
        
        !! 2. Check for possible errors:
        if (nbins < 1) then
            write (stderr, '(a,i0,a)') tag_error // &
                "Invalid number of bins (nbins=", nbins, " must be positive), aborting..."
            stop errmsg_invalid_arg
        else if (xmax <= xmin .or. isnan(xmin) .or. isnan(xmax)) then
            write (stderr, '(a)') tag_error // "Invalid binning range, aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 3. Prepare the bins:
        allocate(bins(nbins + 1))  !! This array contains the counting intervals, the "bins": [x(1), x(2), x(3), ..., x(n), x(n+1)].
        allocate(histo(nbins, 2))  !! This matrix contains the histogram data as a list of couples [x, y] where "x" is the center of each bins.
        select case (bintype)
            case ("linear")  !! Linear lattice of bins.
                bins = [( xmin + (xmax - xmin)*(i-1)/nbins, i = 1, nbins+1 )]
            case ("chebyshev")  !! Bins at the Chebyshev nodes.
                bins = [( xmin + (xmax - xmin)*(1.0 - cos(pi*(i-1)/nbins))/2, i = 1, nbins+1 )]  
            case default
                write (stderr, '(a)') tag_error // "Unknown binning type '" // trim(bintype) // "', aborting..."
                stop errmsg_invalid_arg
        end select
        
        !! 4. Compute the histogram:
        call histogram(samples, bins, histo(:, 2), missing)   !! Compute the histogram of the transmission eigenvalues.
        
        if (missing /= 0.0_wp) then  !! Check for possible missing eigenvalues on the histogram:
            write (stdout, '(a,f0.6,a,f0.6,a,/,a,f0.0,a,f0.2,a)') tag_warn // "Some samples are &
                &out of the bounds [", bins(1), ", ", bins(nbins+1), "].", &
                tag_warn // "There are ", missing, " missing samples in the histogram &
                &(", 100.0*missing/nsamp, "%), continuing anyway..."
        end if
        
        histo(:, 2) = histo(:, 2) + epsilon(histo(1,2))/10.  !! Perturb the data to avoid pure zero density in log scale (points are removed by TikZ/PGFPlots).
        call normalize_histogram(bins, histo(:, 2))   !! Normalize the histogram to get rho(t) such that integral(rho(T)*dT, [0, 1]) = 1.
        call moving_average(bins, 2, histo(:, 1))    !! Calculate the centers of the bins using 2-points sliding average.
        
        !! 5. Save the histogram to a file:
        call save_histo(histo, datafp)
        
        deallocate(bins)  !! Free the memory.
        deallocate(histo)
        deallocate(samples)
        
	end subroutine
    
    !!**********************************************************************************************
    !! Save the histogram in a file. The destination file will be located at the same place as the data file "datafp".
    !! The header of "datafp" is copied to the header of the output file.
    !!**********************************************************************************************
    subroutine save_histo(histo, datafp)
        real(wp), intent(in) :: histo(:, :)
        integer, intent(in) :: datafp
        integer :: outfp
        logical :: openedq
        character(len=150) :: datafname, outputdir, basename, outfname
        
        !! 1. Get the data filename, and check for errors:
        inquire(unit=datafp, opened=openedq, name=datafname)
        if (.not. openedq) then
            write (stderr, '(a,i0,a)') tag_error // "No file connected to unit ", datafp, ", aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 2. Generate a unique filename for the output:
        call path_split(datafname, outputdir, basename)
        call unique_filename(trim(outputdir) // "replot_", ".distrib.csv", outfname)
        
        !! 3. Open the file and write the data:
        open(newunit=outfp, file=outfname, action='write')
        call copy_preamble(datafp, outfp)
        write (outfp, '(a)') "tm, rho"
        call save_real_matrix(outfp, "g0.16", ", ", histo)
        
        close(outfp)
        
        print '(a)', tag_info // "Replot written in file: " // trim(outfname)
        
    end subroutine
	
end program
