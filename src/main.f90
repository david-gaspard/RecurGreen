!!**************************************************************************************************
!! Created on 2024-05-14 at 16:58:43 CEST by David Gaspard <david.gaspard@espci.fr>
!! This program is distributed under the Creative Commons (CC) BY-NC-SA license.
!! Fortran program to compute the transmission eigenvalues in a disordered waveguide
!! using the recursive Green method in the modal representation.
!!**************************************************************************************************
program main
    use omp_lib       !! OpenMP library for parallelization (compile with -fopenmp).
    use waveguide     !! Module to compute the transmission eigenvalues.
    implicit none
    character(len=150) :: arg  !! This assumes that filenames are not too long.
	integer :: i, funit, rstat
	logical :: fexist
    
	write (stdout, '(a)') "====== This is " // trim(program_copyright) // " ======"
	if (command_argument_count() == 0) then
		write (stderr, '(a)') tag_error // "No settings file provided, doing nothing..."
	end if
	do i = 1, command_argument_count()
		call get_command_argument(i, arg)
		inquire(file=arg, exist=fexist)
		if (fexist) then
			open(newunit=funit, file=arg, action='read', iostat=rstat)
			if (rstat == 0) then
				write (stdout, '(a)') tag_info // "Reading file: " // trim(arg)
				call execute_file(funit)
				close(funit)
			else
				write (stderr, '(a)') tag_error // "Cannot open file: " // trim(arg)
			end if
		else
			write (stderr, '(a)') tag_error // "File not found: " // trim(arg)
		end if
	end do
    
    contains
    
	!!**********************************************************************************************
	!! Executes the instructions given in the given namelist file, assuming this file exists and can be opened.
	!!**********************************************************************************************
    subroutine execute_file(funit)
        integer, intent(in) :: funit
        character(len=50) :: task
        integer :: istat
        namelist /settings/ task
        
        read(unit=funit, nml=settings, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        
        select case (task)
            case ("tdistrib")
                call task_tdistrib(funit)
            case ("wavefunction")
                call task_wavefunction(funit)
            case ("tprofile")
                call task_tprofile(funit)
            case ("tprofile_group")
                call task_tprofile_group(funit)
            case ("tpattern")
                call task_tpattern(funit)
            case default
                write (stderr, '(a)') tag_error // "Unknown task '" // trim(task) // "', aborting..."
                stop errmsg_invalid_arg
        end select
        
    end subroutine
    
    !!**********************************************************************************************
    !! Execute the task "tdistrib" with the settings given in the file "funit".
    !! Compute the distribution of transmission eigenvalues using the recursive Green method
    !! either in serial or in parallel (with OpenMP).
    !!**********************************************************************************************
    subroutine task_tdistrib(funit)
        integer, intent(in) :: funit
        type(waveguide_t) :: wguide
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb
        integer :: nxdiso, nxfree, nseed, nbins, nthreads, istat
        real(wp) :: time
        logical :: save_eigvals
        real(wp), allocatable :: tallvals(:), rho_data(:, :)
        
        namelist /waveguide/ boundtype, wol, exfac, dscat, dabso, aratio, nxdiso, nxfree, napera, naperb
        namelist /tdistrib/ nseed, nbins, nthreads, save_eigvals
        
        !! 1. Parse the settings from a namelist and check for errors:
		read(unit=funit, nml=waveguide, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        read(unit=funit, nml=tdistrib, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        
        !! 2. Initialize the waveguide:
        call wguide%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide%print_parameters(stdout, tag_info)
        print '(3(a,i0))', tag_info // "Nseed=", nseed, ", Nbins=", nbins, ", Nthreads=", nthreads
        
        !! 3. Compute the transmission eigenvalues for all the realizations of the disorder:
        allocate(tallvals(nseed * wguide%get_neigv()))  !! Allocate space for all the transmission eigenvalues.
        tallvals = 0.0_wp   !! Initialize to zero (not needed but avoid random memory).
        
        if (nthreads >= 2) then
            call compute_eigvals_omp(wguide, nseed, nthreads, tallvals, time)
        else if (nthreads == 1) then
            call compute_eigvals_serial(wguide, nseed, tallvals, time)
        else
            write (stderr, '(a,i0,a)') tag_error // "Invalid number of threads, &
                &received ", nthreads, " expected >0, aborting..."
            stop errmsg_invalid_arg
        end if
        
        !! 4. Prepare and save the histogram:
        allocate(rho_data(nbins, 2))   !! Allocate space for the transmisison density: rho(T).
        call prepare_distrib(tallvals, nbins, rho_data)
        call save_distrib(rho_data, wguide, nseed, nbins, nthreads, time)
        if (save_eigvals) then   !! Optionally returns all the transmission eigenvalues (to regenerate the histogram).
            call save_raw_eigvals(tallvals, wguide, nseed, nbins, nthreads, time)
        end if
        
        deallocate(rho_data)
        deallocate(tallvals)
        call wguide%del()
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute all the required transmission eigenvalues using OpenMP (parallel computation) and
    !! store them in the array "tallvals".
    !! wguide_shared = Waveguide type containg the simulation parameters (OpenMP-shared variable).
    !! nseed         = Number of random realizations of the disorder.
    !! nthreads      = Number of threads used for the parallelization.
    !! tallvals      = Real array where the transmission eigenvalues are written (output).
    !! time          = Total computation time in seconds (output).
    !!**********************************************************************************************
    subroutine compute_eigvals_omp(wguide_shared, nseed, nthreads, tallvals, time)
        class(waveguide_t), intent(inout) :: wguide_shared
        integer, intent(in) :: nseed, nthreads
        real(wp), intent(out) :: tallvals(:), time
        type(waveguide_t) :: wguide  !! Declare local (OMP-Private) objects of type Waveguide.
        integer :: iseed, cjob, time_start(8), time_end(8)
        real(wp) :: delta(8)
        character(len=30) :: msg
        
        write (msg, '(a,i0,a)') "RecurGreen tdistrib ", nthreads, " thr"
        cjob = 0    !! Initialize the number of completed jobs (OMP-Shared).
        call date_and_time(values=time_start)   !! Start the time measurement. date_and_time() is the true wall time even with parallelization.
        
        !! Main computational loop to compute the transmission eigenvalues:
        !$omp parallel num_threads(nthreads) default(shared) private(iseed, wguide)
            
            call wguide%init_copy(wguide_shared)  !! Create a local copy of the waveguide, one per thread.
                                                  !! Needed because this object stores the disorder realizations, and this is OMP-Private.
            
            !$omp do schedule(dynamic, 1)
            do iseed = 1, nseed   !! Loop on the random realization of the disorder.
                
                call wguide%compute_tdistrib(iseed, tallvals)    !! Compute the transmission eigenvalues for the realization "iseed" of the disorder.
                
                !! Critical section to deal with the progress bar:
                !$omp critical
                cjob = cjob + 1      !! Use a critical section to print a progress bar.
                !!if (mod(cjob, 10) == 0) then  !! Reduce the number of execution.
                call print_progress_bar(cjob, nseed, time_start, msg)  !! Print the progress bar with the expected time of arrival.
                !!end if
                !$omp end critical
                
            end do
            !$omp end do
            
            call wguide%del()
            
        !$omp end parallel
        
        !! Determine the computation time:
        call date_and_time(values=time_end)    !! End the time measurement.
        delta = real(time_end - time_start, kind=wp)
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Computation time in seconds.
        print '(/,a,f0.2,a)', tag_info // "Done, computation time is ", time, " s."
        
    end subroutine
    
    !!**********************************************************************************************
    !! Compute all the required transmission eigenvalues without OpenMP (serial computation) and
    !! store them in the array "tallvals".
    !! wguide   = Waveguide type containg the simulation parameters and where the computations are done.
    !! nseed    = Number of random realizations of the disorder.
    !! tallvals = Real array where the transmission eigenvalues are written (output).
    !! time     = Total computation time in seconds (output).
    !!**********************************************************************************************
    subroutine compute_eigvals_serial(wguide, nseed, tallvals, time)
        class(waveguide_t), intent(inout) :: wguide
        integer, intent(in) :: nseed
        real(wp), intent(out) :: tallvals(:), time
        integer :: iseed, time_start(8), time_end(8)
        real(wp) :: delta(8)
        character(len=30) :: msg
        
        write (msg, '(a)') "RecurGreen tdistrib serial"
        
        call date_and_time(values=time_start)   !! Start the time measurement. date_and_time() is the true wall time even with parallelization.
        
        do iseed = 1, nseed   !! Loop on the random realization of the disorder.
            
            call wguide%compute_tdistrib(iseed, tallvals)    !! Compute a save the transmission eigenvalues for the realization "iseed" of the disorder.
            
            call print_progress_bar(iseed, nseed, time_start, msg)  !! Print the progress bar with the expected time of arrival.
            
        end do
        
        !! Determine the computation time:
        call date_and_time(values=time_end)    !! End the time measurement.
        delta = real(time_end - time_start, kind=wp)
        time = 86400*delta(3) + 3600*delta(5) + 60*delta(6) + delta(7) + delta(8)/1000  !! Computation time in seconds.
        print '(/,a,f0.2,a)', tag_info // "Done, computation time is ", time, " s."
        
    end subroutine
    
    !!**********************************************************************************************
    !! Prepare a histogram from the given list of transmission eigenvalues "tallvals" using bins located
    !! at the Chebyshev nodes. The number of bins is given by "nbins".
    !! The resulting histogram is stored in the 2D array "rho_data". The columns of this array are (tm, rho),
    !! i.e., the transmission eigenvalue, and the density, respectively.
    !!**********************************************************************************************
    subroutine prepare_distrib(tallvals, nbins, rho_data)
        real(wp), intent(in) :: tallvals(:)
        integer, intent(in) :: nbins
        real(wp), intent(out) :: rho_data(:, :)
        real(wp), allocatable :: bins(:), histo(:)
        real(wp) :: missing
        
        allocate(bins(nbins + 1))
        allocate(histo(nbins))
        
        bins = [( (1.0 - cos(pi*(i-1)/nbins))/2, i = 1, nbins+1 )]   !! Bins at the Chebyshev nodes.
        
        call histogram(tallvals, bins, histo, missing)   !! Compute the histogram of the transmission eigenvalues.
        
        if (missing /= 0.0_wp) then  !! Check for possible missing eigenvalues on the histogram:
            write (stdout, '(a,f0.6,a,f0.6,a,/,a,f0.1,a)') tag_warn // "Some transmission eigenvalues are &
                &out of the bounds [", bins(1), ", ", bins(nbins+1), "].", &
                tag_warn // "There are ", missing, " missing eigenvalues in the histogram, continuing anyway..."
        end if
        
        histo = histo + epsilon(histo(1))/10.   !! Perturb the data to avoid pure zero density in log scale (points are removed by TikZ/PGFPlots).
        call normalize_histogram(bins, histo)   !! Normalize the histogram to get rho(t) such that integral(rho(T)*dT, [0, 1]) = 1.
        call moving_average(bins, 2, rho_data(:, 1))    !! Calculate the centers of the bins using 2-points sliding average.
        rho_data(:, 2) = histo   !! Saves the density "histo" in the "rho_data" matrix.
        
        deallocate(bins)
        deallocate(histo)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the 2D distribution stored in "rho_data" to a file. This subroutine also generates a
    !! PGFPlots/TikZ file (a LaTeX file) which is finally compiled.
    !!**********************************************************************************************
    subroutine save_distrib(rho_data, wguide, nseed, nbins, nthreads, time)
        real(wp), intent(in) :: rho_data(:, :), time
        class(waveguide_t), intent(in) :: wguide
        integer, intent(in) :: nseed, nbins, nthreads
        integer :: outfp 
        character(len=120) :: outputdir, filename
        character(len=*), parameter :: script_name = "./tikz/distrib.py"  !! Location of script to be executed (relative to the root directory of the program).
        character(len=800) :: title_tikz
        character(len=len(script_name)+len(filename)+len(title_tikz)+20) :: cmd
        real(wp) :: speed
        
        speed = nseed/time   !! Deduce the mean computation speed, i.e., the number of seeds computed per seconds.
        
        !! 1. Generate a filename and create the appropriate directory:
        call wguide%generate_directory(outputdir)
        outputdir = trim(outputdir) // "distrib" // pathsep
        call create_directory(outputdir)    !! Ensure the output directory exists.
        call unique_filename(trim(outputdir) // "result_", ".distrib.csv", filename)
        
        !! 2. Save the histogram in a CSV file:
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (outfp, '(3(a,i0),2(a,g0.6),a)') "%% Nseed=", nseed, ", Nbins=", nbins, ", Nthreads=", nthreads, &
            ", Time=", time, " s, Speed=", speed, " seed/s"
        write (outfp, '(a)') "tm, rho"
        call save_real_matrix(outfp, "g0.16", ", ", rho_data)
        close(outfp)
        print '(a)', tag_info // "Result written in file: " // trim(filename)
        
        !! 3. Generate the corresponding TikZ file using templates:
        call wguide%latex_parameters(title_tikz)
        write (title_tikz, '(3(a,i0),2(a,f0.2),a)') trim(title_tikz) // ", $N_{\rm seed}=", nseed, &
            "$, $N_{\rm bins}=", nbins, "$, $N_{\rm thr}=", nthreads, "$, Time=", time, " s, Speed=", speed, " seed/s"
        cmd = script_name // " " // trim(filename) // " '" // trim(title_tikz) // "' "
        call execute_command_line(cmd)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the given transmission eigenvalues "tallvals" in a file.
    !!**********************************************************************************************
    subroutine save_raw_eigvals(tallvals, wguide, nseed, nbins, nthreads, time)
        real(wp), intent(in) :: tallvals(:), time
        class(waveguide_t), intent(in) :: wguide
        integer, intent(in) :: nseed, nbins, nthreads
        integer :: outfp 
        character(len=120) :: outputdir, filename
        character(len=30) :: style
        real(wp) :: speed
        
        speed = nseed/time   !! Deduce the mean computation speed, i.e., the number of seeds computed per seconds.
        
        !! 1. Generate a filename and create the appropriate directory:
        call wguide%generate_directory(outputdir)
        outputdir = trim(outputdir) // "distrib" // pathsep
        call create_directory(outputdir)    !! Ensure the output directory exists.
        call unique_filename(trim(outputdir) // "result_", ".eigvals.csv", filename)
        
        !! 2. Save the transmission eigenvalues in a file:
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (outfp, '(3(a,i0),2(a,g0.6),a)') "%% Nseed=", nseed, ", Nbins=", nbins, ", Nthreads=", nthreads, &
            ", Time=", time, " s, Speed=", speed, " seed/s"
        write (outfp, '(a,i0,a)') "%% Transmission eigenvalues (", size(tallvals), " values)"
        write (style, '(a,i0,a)')  "(", size(tallvals), "(g0.16,/))"
        write (outfp, style) tallvals
        close(outfp)
        print '(a)', tag_info // "Eigenvalues saved in file: " // trim(filename)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Execute the task "wavefunction" with the settings given in the file "funit".
    !! Compute the square modulus of the wavefunction using the recursive Green method
    !! either in serial or in parallel (with OpenMP).
    !!**********************************************************************************************
    subroutine task_wavefunction(funit)
        integer, intent(in) :: funit
        type(waveguide_t) :: wguide
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb
        integer :: nxdiso, nxfree, iseed, ipos, imode, istat
        real(wp), allocatable :: wfundata(:, :)
        
        namelist /waveguide/ boundtype, wol, exfac, dscat, dabso, aratio, nxdiso, nxfree, napera, naperb
        namelist /wavefunction/ iseed, ipos, imode
        
        !! 1. Parse the settings from a namelist and check for errors:
		read(unit=funit, nml=waveguide, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        read(unit=funit, nml=wavefunction, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        
        !! 2. Initialize the waveguide:
        call wguide%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide%print_parameters(stdout, tag_info)
        print '(3(a,i0))', tag_info // "Wavefunction parameters: iseed=", iseed, ", ipos=", ipos, ", imode=", imode
        
        !! 3. Compute the wavefunction:
        call wguide%compute_wavefunction(iseed, ipos, imode, wfundata)
        
        !! 4. Saves the wavefunction to a CSV file:
        print '(a)', tag_info // "Saving the wavefunction to a CSV file..."
        call save_wavefunction(wfundata, wguide, iseed, ipos, imode)
        
        call wguide%del()      !! Frees the memory allcoated by the waveguide.
        deallocate(wfundata)   !! Deallocate "wfundata". Allocation took place in the subroutine: compute_wavefunction().
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the 2D wavefunction stored in "wfundata" to a file.
    !!**********************************************************************************************
    subroutine save_wavefunction(wfundata, wguide, iseed, ipos, imode)
        real(wp), intent(in) :: wfundata(:, :)
        class(waveguide_t), intent(in) :: wguide
        integer, intent(in) :: iseed, ipos, imode
        character(len=120) :: outputdir, filename
        integer :: outfp 
        
        !! 1. Generate a filename and create the appropriate directory:
        call wguide%generate_directory(outputdir)
        outputdir = trim(outputdir) // "wavefun" // pathsep
        call create_directory(outputdir)    !! Ensure the output directory exists.
        call unique_filename(trim(outputdir) // "result_", ".wavefun.csv", filename)
        
        !! 2. Save the histogram in a CSV file:
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (outfp, '(3(a,i0))') "%% Wavefunction parameters: iseed=", iseed, ", ipos=", ipos, ", imode=", imode
        write (outfp, '(a)') "x, y, repsi, impsi, upot"
        call save_real_matrix(outfp, "g0.16", ", ", wfundata)
        close(outfp)
        print '(a)', tag_info // "Result written in file: " // trim(filename)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Execute the task "tprofile" with the settings given in the file "funit".
    !! Compute the transmission-eigenstate intensity profile <|psi_T(x)|^2> using the recursive Green method
    !! and OpenMP to parallelize the disorder average.
    !!**********************************************************************************************
    subroutine task_tprofile(funit)
        integer, intent(in) :: funit
        type(waveguide_t) :: wguide, wguide_shared
        real(wp), allocatable :: tprofiledata(:, :), tprofileavg(:)
        integer, allocatable :: ntdata(:)
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb, tmin, tmax, time
        integer :: nxdiso, nxfree, iseed, nseed, nx, nsample, nthreads, istat, start(8), cjob
        character(len=30) :: msg
        
        namelist /waveguide/ boundtype, wol, exfac, dscat, dabso, aratio, nxdiso, nxfree, napera, naperb
        namelist /tprofile/ tmin, tmax, nseed, nthreads
        
        !! 1. Parse the settings from a namelist and check for errors:
		read(unit=funit, nml=waveguide, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        read(unit=funit, nml=tprofile, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        
        !! 2. Initialize the waveguide:
        call wguide_shared%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide_shared%print_parameters(stdout, tag_info)
        
        !! 3. Check for possible errors before allocating memory:
        if (nseed <= 0) then
            write (stderr, '(a,i0,a)') tag_error // "Invalid number of disorder realizations (received ", nseed, "), aborting..."
            stop errmsg_invalid_arg
        end if
        
        nx = wguide_shared%get_nx()
        allocate(tprofiledata(nseed, nx))
        allocate(ntdata(nseed))
        
        write (msg, '(a,i0,a)') "RecurGreen tprofile ", nthreads, " thr"
        cjob = 0    !! Initialize the number of completed jobs (OMP-Shared).
        call date_and_time(values=start)   !! Start the time measurement. date_and_time() is the true wall time even with parallelization.
        
        !! 4. Loop on the realization of disorder:
        !$omp parallel num_threads(nthreads) default(shared) private(iseed, wguide)
            
            call wguide%init_copy(wguide_shared)  !! Create a local copy of the waveguide, one per thread.
                                                  !! Needed because this object stores the disorder realizations, and this is OMP-Private.
            
            !$omp do schedule(dynamic, 1)
            do iseed = 1, nseed   !! Loop on the random realization of the disorder.
                
                call wguide%compute_tprofile(iseed, tmin, tmax, tprofiledata(iseed, :), ntdata(iseed))  !! Compute the transmission-eigenstate profile.
                
                !! Critical section to deal with the progress bar:
                !$omp critical
                cjob = cjob + 1      !! Use a critical section to print a progress bar.
                call print_progress_bar(cjob, nseed, start, msg)  !! Print the progress bar with the expected time of arrival.
                !$omp end critical
                
            end do
            !$omp end do
            
            call wguide%del()
            
        !$omp end parallel
        
        call end_progress_bar(start, time)  !! Stops the progress bar.
        
        !! 5. Compute the disorder average here in serial:
        nsample = count(ntdata /= 0)   !! Counts the number of nonzero elements. Actual number of samples for the disorder average. 
        
        if (nsample == 0) then
            write (stdout, '(a,i0,a)') tag_warn // &
                "Found no transmission eigenfunction (nsample=", nsample, "), results will be undetermined..."
        else
            write (stdout, '(2(a,i0),a)') tag_info // &
                "Found ", nsample, "/", nseed, " transmission eigenfunctions, now averaging..."
        end if
        
        allocate(tprofileavg(nx))   !! Allocate space for the average profile function, <|psi_T(x)|^2>.
        tprofileavg = 0.0_wp
        do iseed = 1, nseed
            if (ntdata(iseed) /= 0) then
                tprofileavg = tprofileavg + tprofiledata(iseed, :)
            end if
        end do
        tprofileavg = tprofileavg/nsample
        
        !! 6. Save the transmission-eigenfunction profile to a file:
        call save_tprofile(tprofiledata, tprofileavg, nsample, wguide_shared, tmin, tmax, nseed, time)
        
        call wguide_shared%del()   !! Frees the memory allcoated by the waveguide.
        deallocate(tprofiledata)
        deallocate(tprofileavg)
        deallocate(ntdata)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the simulation data for the transmission-eigenstate intensity profile into files.
    !! - One file for the raw data, i.e., the function |psi_T(x)|^2 for individual realizations of the disorder.
    !! - Another file for the average function <|psi_T(x)|^2> along with the corresponding abscissas ("x" values).
    !!**********************************************************************************************
    subroutine save_tprofile(tprofiledata, tprofileavg, nsample, wguide, tmin, tmax, nseed, time)
        real(wp), intent(in) :: tprofiledata(:, :), tprofileavg(:)
        class(waveguide_t), intent(in) :: wguide
        real(wp), intent(in) :: tmin, tmax, time
        integer, intent(in) :: nsample, nseed
        real(wp), allocatable :: xmesh(:)
        character(len=*), parameter :: script_name = "./tikz/tprofile.py"  !! Location of script to be executed (relative to the root directory of the program).
        character(len=120) :: outputdir, filename
        character(len=800) :: title_tikz
        character(len=len(script_name)+len(filename)+len(title_tikz)+20) :: cmd
        integer :: outfp, ix, nx
        
        !! 1. Extract the xmesh from the waveguide:
        nx = wguide%get_nx()
        allocate(xmesh(nx))
        call wguide%get_xmesh(xmesh)  !! Extracts the xmesh from the waveguide.
        
        !! 2. Generate a filename and create the appropriate directory:
        call wguide%generate_directory(outputdir)
        outputdir = trim(outputdir) // "tprofile" // pathsep
        call create_directory(outputdir)    !! Ensure the output directory exists.
        call unique_filename(trim(outputdir) // "result_", ".raw.csv", filename)
        
        !! 3. Save the raw intensity profile data in a CSV file:
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (outfp, '(3(a,i0),2(a,g0.6),a,f0.3,a)') "%% Transmission profile parameters: Found nsample=", &
            nsample, " / nseed=", nseed, ", nx=", nx, ", tmin=", tmin, ", tmax=", tmax, ", time=", time, " s."
        write (outfp, '(a)') "%% Here follows, on each line, the function |psi_T(x)|^2 for a given realization of the disorder."
        call save_real_matrix(outfp, "g0.16", ", ", tprofiledata)
        close(outfp)
        print '(a)', tag_info // "Raw profile data written in file: " // trim(filename)
        
        !! 4. Saves the average profile <|psi_T(x)|^2>:
        call unique_filename(trim(outputdir) // "result_", ".avg.csv", filename)
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (outfp, '(3(a,i0),2(a,g0.6),a,f0.3,a)') "%% Transmission profile parameters: Found nsample=", &
            nsample, " / nseed=", nseed, ", nx=", nx, ", tmin=", tmin, ", tmax=", tmax, ", time=", time, " s."
        write (outfp, '(a)') "x, tprofile"
        do ix = 1, nx
            write (outfp, '(g0.16,", ",g0.16)') xmesh(ix), tprofileavg(ix)
        end do
        close(outfp)
        print '(a)', tag_info // "Average profile written in file: " // trim(filename)
        
        !! 5. Generate the corresponding TikZ file using templates:
        call wguide%latex_parameters(title_tikz)
        write (title_tikz, '(2(a,i0),2(a,f0.6),a,f0.2,a)') trim(title_tikz) // ", $N_{\rm sample}=", nsample, &
            "$, $N_{\rm seed}=", nseed, "$, $T_{\rm min}=", tmin, "$, $T_{\rm max}=", tmax, "$, Time=", time, " s"
        cmd = script_name // " " // trim(filename) // " '" // trim(title_tikz) // "' "
        call execute_command_line(cmd)
        
        deallocate(xmesh)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Execute the task "tprofile_group" with the settings given in the file "funit".
    !! Compute a group of 4 transmission-eigenstate profiles, <|psi_T(x)|^2>, using the recursive Green method
    !! and OpenMP to parallelize the disorder average.
    !!**********************************************************************************************
    subroutine task_tprofile_group(funit)
        integer, intent(in) :: funit
        type(waveguide_t) :: wguide, wguide_shared
        integer, parameter :: nprofile = 4  !! Number of transmission profiles computed at the same time.
        real(wp), allocatable :: tprofiledata(:, :, :)  !! Transmission profiles. Dimensions=(nx, nprofile, nseed).
        real(wp), allocatable :: tprofileavg(:, :)  !! Transmission profile averages. Dimensions=(nx, nprofile).
        integer, allocatable :: ntdata(:, :)  !! Number of found transmission values. Dimension=(nprofile, nseed).
        real(wp), allocatable :: tavgdata(:)  !! List of transmission probability averages. Dimensions=(nseed).
        character(len=50) :: boundtype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb, tm(nprofile), dt(nprofile), tavg, uctavg, time
        integer :: nxdiso, nxfree, iseed, nseed, nx, nsample(nprofile), nthreads, istat, start(8), cjob, ip
        character(len=50) :: msg
        
        namelist /waveguide/ boundtype, wol, exfac, dscat, dabso, aratio, nxdiso, nxfree, napera, naperb
        namelist /tprofile_group/ tm, dt, nseed, nthreads
        
        !! 1. Parse the settings from a namelist and check for errors:
		read(unit=funit, nml=waveguide, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        read(unit=funit, nml=tprofile_group, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        
        !! 2. Initialize the waveguide:
        call wguide_shared%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        call wguide_shared%print_parameters(stdout, tag_info)
        
        !! 3. Check for possible errors before allocating memory:
        if (nseed <= 0) then
            write (stderr, '(a,i0,a)') tag_error // "Invalid number of disorder realizations (received ", nseed, "), aborting..."
            stop errmsg_invalid_arg
        end if
        
        nx = wguide_shared%get_nx()
        allocate(tprofiledata(nx, nprofile, nseed))
        allocate(tprofileavg(nx, nprofile))
        allocate(ntdata(nprofile, nseed))
        allocate(tavgdata(nseed))
        
        write (msg, '(a,i0,a)') "RecurGreen tprofile_group ", nthreads, " thr"
        cjob = 0    !! Initialize the number of completed jobs (OMP-Shared).
        call date_and_time(values=start)   !! Start the time measurement. date_and_time() is the true wall time even with parallelization.
        
        !! 4. Loop on the realization of disorder:
        !$omp parallel num_threads(nthreads) default(shared) private(iseed, wguide)
            
            call wguide%init_copy(wguide_shared)  !! Create a local copy of the waveguide, one per thread.
                                                  !! Needed because this object stores the disorder realizations, and this is OMP-Private.
            
            !$omp do schedule(dynamic, 1)
            do iseed = 1, nseed   !! Loop on the random realization of the disorder.
                
                call wguide%compute_tprofile_group(iseed, tm, dt, tprofiledata(:, :, iseed), ntdata(:, iseed), tavgdata(iseed))  !! Compute the transmission-eigenstate profiles.
                
                !! Critical section to deal with the progress bar:
                !$omp critical
                cjob = cjob + 1      !! Use a critical section to print a progress bar.
                call print_progress_bar(cjob, nseed, start, msg)  !! Print the progress bar with the expected time of arrival.
                !$omp end critical
                
            end do
            !$omp end do
            
            call wguide%del()
            
        !$omp end parallel
        
        call end_progress_bar(start, time)  !! Stops the progress bar.
        
        !! 5. Compute the disorder average here in serial:
        do ip = 1, nprofile  !! Loop on the transmission profiles to average.
            tprofileavg(:, ip) = 0.0_wp  !! Initialize the average profile to zero because it will be incremented.
            nsample(ip) = 0  !! Initialize the number of nonzero transmission profiles to zero because it will be incremented.
            do iseed = 1, nseed   !! Loop on the realizations of disorder to compute the average.
                if (ntdata(ip, iseed) /= 0) then
                    tprofileavg(:, ip) = tprofileavg(:, ip) + tprofiledata(:, ip, iseed)
                    nsample(ip) = nsample(ip) + 1
                end if
            end do
            if (nsample(ip) /= 0) then
                write (stdout, '(4(a,i0),a)') tag_info // &
                    "Profile ", ip, "/", nprofile, " : Found ", nsample(ip), "/", nseed, &
                    " transmission eigenfunctions, now averaging..."
                tprofileavg(:, ip) = tprofileavg(:, ip)/nsample(ip)
            else
                write (stdout, '(3(a,i0),a)') tag_warn // &
                    "Profile ", ip, "/", nprofile, " : Found no transmission eigenstate (nsample=", &
                    nsample(ip), "), results will be undetermined..."
            end if
        end do
        
        !! 6. Compute the average transmission and estimate the uncertainty:
        tavg = sum(tavgdata)/nseed   !! Average transmission probability.
        uctavg = sqrt(sum((tavgdata - tavg)**2))/nseed   !! Uncertainty over the average value "tavg".
        
        !! 7. Save the transmission-eigenfunction profile to a file:
        call save_tprofile_group(tprofileavg, nsample, tavg, uctavg, wguide_shared, tm, dt, nseed, time)
        
        call wguide_shared%del()   !! Frees the memory allcoated by the waveguide.
        deallocate(tprofiledata)
        deallocate(tprofileavg)
        deallocate(ntdata)
        deallocate(tavgdata)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the simulation data for the transmission-eigenstate intensity profile into files.
    !! - One file for the raw data, i.e., the function |psi_T(x)|^2 for individual realizations of the disorder.
    !! - Another file for the average function <|psi_T(x)|^2> along with the corresponding abscissas ("x" values).
    !!**********************************************************************************************
    subroutine save_tprofile_group(tprofileavg, nsample, tavg, uctavg, wguide, tm, dt, nseed, time)
        real(wp), intent(in) :: tprofileavg(:, :)  !! Average transmission profiles. Dimensions=(nx, nprofile).
        integer, intent(in) :: nsample(:), nseed
        class(waveguide_t), intent(in) :: wguide
        real(wp), intent(in) :: tavg, uctavg, tm(:), dt(:), time
        real(wp), allocatable :: xmesh(:)
        character(len=*), parameter :: script_name = "./tikz/tprofile_group.py"  !! Location of script to be executed (relative to the root directory of the program).
        character(len=120) :: outputdir, filename, strfmt
        character(len=800) :: title_tikz
        character(len=len(script_name)+len(filename)+len(title_tikz)+20) :: cmd
        integer :: outfp, nprofile, ix, nx
        
        !! 1. Extract the xmesh from the waveguide:
        nprofile = size(tprofileavg, 2)
        nx = wguide%get_nx()
        allocate(xmesh(nx))
        call wguide%get_xmesh(xmesh)  !! Extracts the xmesh from the waveguide.
        
        !! 2. Create the appropriate directory:
        call wguide%generate_directory(outputdir)
        outputdir = trim(outputdir) // "tprofile_group" // pathsep
        call create_directory(outputdir)   !! Ensure the output directory exists.
        
        !! 3. Saves the average profiles <|psi_T(x)|^2>:
        call unique_filename(trim(outputdir) // "result_", ".csv", filename)
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (strfmt, '(a,i0,a)') "(a,i0,2(a,g0.6,", nprofile-1, "(', ',g0.6)),a)"
        write (outfp, strfmt) "%% Transmission profile parameters: nx=", nx, ", tm=[", tm, "], dt=[", dt, "]."
        write (strfmt, '(a,i0,a)') "(a,i0,", nprofile-1, "(', ',i0),a,i0,2(a,g0.6),a,f0.3,a)"
        write (outfp, strfmt) "%% Results: Found nsample=[", nsample, "] over nseed=", nseed, &
            ", tavg=", tavg, "+-", uctavg, ", time=", time, " s."
        write (strfmt, '(a,i0,a)') "('x',", nprofile, "(', tprofile',f0.4))"
        write (outfp, strfmt) tm
        write (strfmt, '(a,i0,a)') "(g0.16,", nprofile, "(', ',g0.16))"
        do ix = 1, nx
            write (outfp, strfmt) xmesh(ix), tprofileavg(ix, :)
        end do
        close(outfp)
        print '(a)', tag_info // "Average profiles written in file: " // trim(filename)
        
        !! 4. Generate the corresponding TikZ file using templates:
        call wguide%latex_parameters(title_tikz)
        write (strfmt, '(2(a,i0),a)') "(a,i0,", nprofile-1, "(', ',i0),a,i0,2(a,f0.6,", &
            nprofile-1, "(', ',f0.6)),2(a,f0.6),a,f0.3,a)"
        write (title_tikz, strfmt) trim(title_tikz) // ", $N_{\rm sample}=\{", nsample, &
            "\}$, $N_{\rm seed}=", nseed, "$, $T_{\rm m}=\{", tm, "\}$, $\delta T=\{", dt, &
            "\}$, $\bar{T}_{\rm obs}=", tavg, "\pm", uctavg, "$, Time=", time, " s"
        cmd = script_name // " " // trim(filename) // " '" // trim(title_tikz) // "' "
        call execute_command_line(cmd)
        
        deallocate(xmesh)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Execute the task "tpattern" with the settings given in the file "funit".
    !! Compute the radiation pattern of the transmitted probabilities per channel using the recursive Green method
    !! and OpenMP to parallelize the disorder average.
    !!**********************************************************************************************
    subroutine task_tpattern(funit)
        integer, intent(in) :: funit
        type(waveguide_t) :: wguide, wguide_shared
        complex(wp), allocatable :: tmat(:,:), rmat(:,:), depmat(:,:)
        real(wp), allocatable :: tprob(:), rprob(:), tprob_array(:,:), rprob_array(:,:)
        character(len=50) :: boundtype, probatype
        real(wp) :: wol, aratio, exfac, dscat, dabso, napera, naperb, time
        integer :: nxdiso, nxfree, iseed, imode, nmu, imu, nx, ipos, fpos, nseed, istat, start(8), cjob, nthreads
        character(len=50) :: msg
        
        namelist /waveguide/ boundtype, wol, exfac, dscat, dabso, aratio, nxdiso, nxfree, napera, naperb
        namelist /tpattern/ nseed, imode, nthreads, probatype
        
        !! 1. Parse the settings from a namelist and check for errors:
		read(unit=funit, nml=waveguide, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        read(unit=funit, nml=tpattern, iostat=istat)
		call check_and_rewind_nml(funit, istat)  !! Check for reading error and rewind the namelist file (for further reading).
        
        if (nseed <= 0) then
            write (stderr, '(a,i0,a)') tag_error // "Invalid number of disorder realizations (received ", nseed, "), aborting..."
            stop errmsg_invalid_arg
        else if (probatype /= "incident" .and. probatype /= "total" .and. probatype /= "deposition") then
            write (stderr, '(a,i0,a)') tag_error // &
                "Invalid probatype, received '" // trim(probatype) // "', expected incident/total/deposition, aborting..."
            stop errmsg_invalid_arg
        end if
        
        call wguide_shared%init(boundtype, nxdiso, nxfree, wol, exfac, aratio, dscat, dabso, napera, naperb)
        
        !! 2. Shift the mode index in the periodic case to get imode=1 correspond to the ground state.
        nx = wguide_shared%get_nx()  !! Extract the number of x points.
        nmu = wguide_shared%get_nmu()  !! Extract the dimensions of the transmission/reflection matrices.
        
        ipos = 1  !! Initial position is the first x point at the leftmost side (only used for probatype="deposition").
        fpos = nxfree + nxdiso/2 + 1  !! The final position is the middle point in the center of the medium (only used for probatype="deposition").
        
        if (boundtype == "periodic") then
            imu = mod(imode + nmu/2 - 1, nmu) + 1
        else
            imu = imode
        end if
        
        !! 3. Allocate space for the transmission/reflection matrices:
        allocate(tprob(nmu))
        allocate(rprob(nmu))
        allocate(tprob_array(nmu, nseed))
        allocate(rprob_array(nmu, nseed))
        
        write (msg, '(a,i0,a)') "RecurGreen tpattern ", nthreads, " thr"
        cjob = 0    !! Initialize the number of completed jobs (OMP-Shared).
        call date_and_time(values=start)   !! Start the time measurement. date_and_time() is the true wall time even with parallelization.
        
        !! 4. Computes the transmission/reflection matrices for different realizations of the disorder:
        !$omp parallel num_threads(nthreads) default(shared) private(iseed, wguide, tmat, rmat, depmat)
            
            call wguide%init_copy(wguide_shared)  !! Create a local copy of the waveguide, one per thread.
                                                  !! Needed because this object stores the disorder realizations, and this is OMP-Private.
            
            if (probatype == "deposition") then
                allocate(depmat(nmu, nmu))
            else
                allocate(tmat(nmu, nmu))
                allocate(rmat(nmu, nmu))
            end if
            
            !$omp do schedule(dynamic, 1)
            do iseed = 1, nseed   !! Loop on the random realization of the disorder.
                
                if (probatype == "deposition") then
                    
                    call wguide%compute_deposition(iseed, fpos, ipos, depmat)
                    
                    tprob_array(:, iseed) = abs(depmat(:, imu))**2
                    rprob_array(:, iseed) = 0.0_wp  !! Transmission and reflection are the same in this case (cannot be separated easily).
                    
                else
                    
                    call wguide%compute_refl_trans_v1(iseed, rmat, tmat)
                    
                    if (probatype == "incident") then !! Compute the reflection/transmission probability per channel starting from channel "imu".
                        tprob_array(:, iseed) = abs(tmat(:, imu))**2
                        rprob_array(:, iseed) = abs(rmat(:, imu))**2
                    else !! Compute the total reflection/transmission probability per channel.
                        tprob_array(:, iseed) = sum(abs(tmat)**2, dim=1)
                        rprob_array(:, iseed) = sum(abs(rmat)**2, dim=1)
                    end if
                    
                end if
                
                !! Critical section to deal with the progress bar:
                !$omp critical
                cjob = cjob + 1      !! Use a critical section to print a progress bar.
                call print_progress_bar(cjob, nseed, start, msg)  !! Print the progress bar with the expected time of arrival.
                !$omp end critical
                
            end do
            !$omp end do
            
            call wguide%del()
            if (allocated(tmat)) deallocate(tmat)
            if (allocated(rmat)) deallocate(rmat)
            if (allocated(depmat)) deallocate(depmat)
            
        !$omp end parallel
        
        call end_progress_bar(start, time)  !! Stops the progress bar.
        
        !! 5. Computes the average over the realizations of the disorder:
        tprob = sum(tprob_array, dim=2)/nseed
        rprob = sum(rprob_array, dim=2)/nseed
        
        !! 6. Save the transmission-eigenfunction profile to a file:
        call save_tpattern(tprob, rprob, wguide_shared, probatype, imode, nseed, time)
        
        call wguide_shared%del()   !! Frees the memory allcoated by the waveguide.
        deallocate(tprob)
        deallocate(rprob)
        deallocate(tprob_array)
        deallocate(rprob_array)
        
    end subroutine
    
    !!**********************************************************************************************
    !! Save the simulation data for the radiation pattern of transmitted/reflected probabilities into files.
    !!**********************************************************************************************
    subroutine save_tpattern(tprob, rprob, wguide, probatype, imode, nseed, time)
        class(waveguide_t), intent(in) :: wguide
        real(wp), intent(in) :: tprob(:), rprob(:), time
        character(len=*), intent(in) :: probatype
        integer, intent(in) :: imode, nseed
        integer :: imu, nmu, outfp
        character(len=*), parameter :: script_name = "./tikz/tpattern.py"  !! Location of script to be executed (relative to the root directory of the program).
        character(len=120) :: outputdir, filename
        character(len=800) :: title_tikz
        character(len=len(script_name)+len(filename)+len(title_tikz)+20) :: cmd
        
        nmu = size(tprob)
        
        !! 1. Create the appropriate directory:
        call wguide%generate_directory(outputdir)
        outputdir = trim(outputdir) // "tpattern" // pathsep
        call create_directory(outputdir)   !! Ensure the output directory exists.
        
        !! 2. Saves the average transmission/reflection pattern:
        call unique_filename(trim(outputdir) // "result_", ".csv", filename)
        open(newunit=outfp, file=filename, action='write')
        call print_timestamp(outfp, "%%")
        call wguide%print_parameters(outfp, "%%")
        write (outfp, '(2(a,i0),a,f0.3,a)') &
            "%% Radiation pattern parameters: probatype='" // trim(probatype) // &
            "', nseed=", nseed, ", imode=", imode, ", time=", time, " s."
        write (outfp, '(a)') "imu, tprob, rprob"
        do imu = 1, nmu
            write (outfp, '(i0,", ",g0.16,", ",g0.16)') imu, tprob(imu), rprob(imu)
        end do
        close(outfp)
        print '(a)', tag_info // "Average profiles written in file: " // trim(filename)
        
        !! 3. Generate the corresponding TikZ file using templates:
        call wguide%latex_parameters(title_tikz)
        write (title_tikz, '(2(a,i0),a,f0.3,a)') trim(title_tikz) // ", Probatype=" // trim(probatype) // &
            ", $i_{\rm mode}=", imode, "$, $N_{\rm seed}=", nseed, "$, Time=", time, " s"
        cmd = script_name // " " // trim(filename) // " '" // trim(title_tikz) // "' "
        call execute_command_line(cmd)
        
    end subroutine
    
end program
