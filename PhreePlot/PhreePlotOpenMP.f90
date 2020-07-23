   program PhreePlotOpenMP
   ! demonstrates the use of PhreeqcRM with OpenMP to calculate the predominant species over a grid

   USE IPhreeqc
   USE PhreeqcRM
   USE omp_lib

   implicit none

   integer, parameter :: dp = kind(1.d0)
   integer :: nres = 50
   integer :: nthreads = 0
   integer :: nxyz
   integer :: status_rm
   integer :: status_ip
   real(dp) :: t0, t1, t_start

   integer :: id
   integer :: iphreeqc_id1
   integer :: w
   character(1024) :: line
   character(:), allocatable :: output(:)
   character(:), allocatable :: string
   integer :: i, j
   integer, allocatable :: print_chemistry_mask(:)
   integer, allocatable :: cell_numbers(:)

   integer :: cell = 1
   integer, allocatable :: startcell(:), endcell(:)
   real(dp) :: ph
   real(dp) :: o2
   real(dp) :: min_ph = -4.0_dp
   real(dp) :: max_ph = 10.0_dp
   real(dp) :: d_ph
   real(dp) :: min_o2 = -80._dp
   real(dp) :: max_o2 = 0.0_dp
   real(dp) :: d_o2

   integer, allocatable :: ic1(:,:)
   integer, allocatable :: ic2(:,:)
   real(dp), allocatable :: f1(:,:)

   integer :: type
   real(dp) :: dvalue
   character(:), allocatable :: svalue
   integer :: svalue_size = 4 !dlp
   character(1024) :: out
   integer :: slength
   integer :: tc, nrow, ncol, ir
   integer :: n
   integer :: unit
   integer :: ios
   character(1), parameter :: tab = CHAR(09)
   character(80) :: fn
   character(16) :: arg

   character(1), parameter :: eol = ';'
   character(32) :: cval(3)
   character(32), parameter :: str(3) = ["EQUILIBRIUM_PHASES_MODIFY",&
      &"  -component Fix_H+;  -si",&
      &"  -component O2(g);   -si"] !dlp, gfortran required equal lengths?

   ! get resolution and nthreads from command line
   if (COMMAND_ARGUMENT_COUNT() .ge. 1) then
      call GET_COMMAND_ARGUMENT(1,arg)
      read (arg,*, iostat = ios) nres
      if (ios .ne. 0  .or. nres .le. 1) stop 'Bad nres.'
      if (COMMAND_ARGUMENT_COUNT() .ge. 2) then
         call GET_COMMAND_ARGUMENT(2,arg)
         read (arg,*, iostat = ios) nthreads
         if (ios .ne. 0 .or. nthreads .lt. 0) stop 'Bad nthreads.'
      else
         nthreads = 1
      endif
   else
      nres = 50
   endif
   
   if (nthreads .eq. 0) nthreads = omp_get_max_threads()

   write (*,'("nres = ", I0, 1x, ", nthreads = ", I0)')  nres, nthreads

   allocate (startcell(0:nthreads-1), endcell(0:nthreads-1), source = 0)

   nxyz = nres * nres

   d_ph = (max_ph - min_ph) / real(nres - 1, dp)
   d_o2 = (max_o2 - min_o2) / real(nres - 1, dp)

   t0 = omp_get_wtime()
   t_start = t0

   ! initialize
   id = RM_Create(nxyz, nthreads)
   status_rm = RM_loaddatabase(id, "wateq4f.dat")  !dlp
   status_rm = RM_RunFile(id, 1, 1, 0, "preloop2.pqi") !dlp
   status_rm = RM_SetFilePrefix(id, "PhreePlotTest")
   status_rm = RM_SetErrorHandlerMode(id, 3)   !OpenMP only
   status_rm = RM_SetErrorOn(id, 0) ! 0 = off; 1 = on ! dlp
   status_rm = RM_SetScreenOn(id, 0)  ! 0 = off; 1 = on ! dlp
   !$OMP PARALLEL DO PRIVATE(w, status_rm, status_ip)
   do w = 0, nthreads - 1
      iphreeqc_id1 = RM_GetIPhreeqcId(id, w)
      !status_ip = SetErrorFileOn(iphreeqc_id1, .false.) !dlp
      status_ip = SetOutputStringOn(iphreeqc_id1, .true.)
   enddo
   !$OMP END PARALLEL DO 
   
   status_rm = RM_OpenFiles(id)

   ! only one RunCells so no need to rebalance the load for the next RunCells
   status_rm = RM_SetRebalanceFraction(id, 0._dp)
   ! delete all reactants in worker, selected_output still defined
   !status_rm = RM_RunString(id, 1, 0, 0, "DELETE; -all")
   ! timing
   t1 = omp_get_wtime()
   !write (*, '("Before equilibrium phases definition", t55, F8.3, " sec")') t1 - t0
   t0 = t1

   ! Make equilibrium_phases grid using initial iphreeqc instance

   call omp_set_num_threads(nthreads)

   status_rm = RM_GetStartCell(id,startcell)
   status_rm = RM_GetEndCell(id,endcell)

   allocate(character(1024) :: output(0:nthreads-1))
   
   !$OMP PARALLEL DO PRIVATE(w, iphreeqc_id1, line, ph, o2, cval, string, cell, output)
   do w = 0, nthreads - 1
      output(w) = ''
      iphreeqc_id1 = RM_GetIPhreeqcId(id, w)
      write (line, '("COPY solution 1 ",I0,"-",I0,";")') startcell(w), endcell(w)
      output(w) = trim(line)
      write (line, '("COPY surface 1 ",I0,"-",I0,";")') startcell(w), endcell(w)
      output(w) = trim(output(w))//trim(line)
      write (line, '("COPY equilibrium_phases 1 ",I0,"-",I0,";")') startcell(w), endcell(w)
      output(w) = trim(output(w))//trim(line)
      !write (*,'("RunString(",I0,"): ", A)') w, trim(output(w))
      status_ip = RunString (iphreeqc_id1, output(w))

      t1 = omp_get_wtime()
      !write (*, '("After COPY block 1", t55, F8.3, " sec")') t1 - t0
      t0 = t1

      !write (*,'("Process ", I0,"-",I0)') startcell(w), endcell(w)

      ph = min_ph - d_ph

      cell = 0

      ph_loop: do
         ph = ph + d_ph
         if (ph > max_ph + 1.e-8_dp) exit ph_loop
         write (cval(2), '(1x, F0.8)') -ph
         o2 = min_o2 - d_o2
         o2_loop: do
            o2 = o2 + d_o2
            if (o2 > max_o2 + 1.e-8_dp) exit o2_loop
            if (cell >= startcell(w) .and. cell <= endcell(w)) then
               write (cval(3), '(1x, F0.8)') o2
               write (cval(1) , '(1x, I0)') cell
               string = trim(str(1))//trim(cval(1))//eol//trim(str(2))//trim(cval(2))//eol//trim(str(3))//trim(cval(3))//eol
               status_ip = AccumulateLine(iphreeqc_id1, string)
               !write (*,'("Cell = ", I0,": ",A)') cell, string
               !!!!!status_ip = RunString(iphreeqc_id1, string)
            endif
            cell = cell + 1
         enddo o2_loop
      enddo ph_loop
      status_ip = RunAccumulated(iphreeqc_id1)
      !write (*,'(/"RunAccumated ", I0)') status_ip
   enddo
   !$OMP END PARALLEL DO

   t1 = omp_get_wtime()
   write (*, '("After accumulating EQUILIBRIUM_PHASES_MODIFY string", t55, F8.3, " sec")') t1 - t0
   t0 = t1

   !! modify all equilibrium_phases in initial iphreeqc instance
   !!status_ip = RunAccumulated(iphreeqc_id1)
   !!status_ip = RunString(iphreeqc_id1, string)
   !
   !! timing
   !!t1 = omp_get_wtime()
   !!write (*,'("After RunAccumulated", t55, F8.3," sec")') t1 - t0
   !!t0 = t1
   !
   !! Set array of initial conditions
   !!status_rm = RM_SetUnitsSurface(id, 0)
   !!status_rm = RM_SetUnitsPPassemblage(id, 0)
   !
   !! transfer all cell definitions to workers
   !!status_rm = RM_InitialPhreeqc2Module(id, ic1)
   !
   !!allocate (cell_numbers(nxyz))
   !!cell_numbers = [1:nxyz]
   !!status_rm = RM_InitialPhreeqcCell2Module(id, 1, cell_numbers, nxyz)

   ! timing
   t1 = omp_get_wtime()
   write (*, '("After initial conditions", t55, F8.3," sec" )') t1 - t0
   t0 = t1

   status_rm = RM_SetSelectedOutputOn(id, 1)
   ! for output file
   status_rm = RM_SetPrintChemistryOn(id, 0, 0, 0)   ! workers, initial_phreeqc, utility
   allocate (print_chemistry_mask(nxyz), source = 1)
   status_rm = RM_SetPrintChemistryMask(id, print_chemistry_mask)

   ! run nxyz cells in parallel
   status_rm = RM_RunCells(id)

   if (status_rm .ne. 0) stop 'Failed in RunCells'

   ! timing
   t1 = omp_get_wtime()
   write (*, '("Time for RunCells", t55, F8.3," sec" )') t1 - t0
   write (*, '("Total without processing results", t55, F8.3," sec" )') t1 - t_start
   t0 = t1;

   ! process results
   write (fn, '("output-",I0,"-",I0,".dat")') nres, nthreads
   open (newunit = unit, file = fn, status = 'UNKNOWN')

   allocate(character(svalue_size) :: svalue) !dlp

   n = -1

   nrow = 0
   ncol  = 0

   tc = RM_GetThreadCount(id)

   do w = 0, tc - 1
      iphreeqc_id1 = RM_GetIPhreeqcId(id, w)
      status_ip = SetCurrentSelectedOutputUserNumber(iphreeqc_id1, 1)
      nrow = GetSelectedOutputRowCount(iphreeqc_id1)
      ncol = GetSelectedOutputColumnCount(iphreeqc_id1)
      
      if (nrow .le. 0 .or. ncol .le. 0) write (*,'("Phreeqc failed for iteration: ",I0)') n

      write (*, '("Worker number: ", I0)') w

      do i = 0, nrow
         out = ''
         if (w > 0 .and. i == 0) cycle
         do j = 1, ncol
            status_ip = GetSelectedOutputValue(iphreeqc_id1, i, j, type, dvalue, svalue, slength)
			
			if(type .eq. 0) then
				svalue(:) = "NA"
				slength = 2
			else
				do while (slength .ge. svalue_size)
					deallocate(svalue)
					svalue_size = svalue_size*2
					write(*,*) "Doubling to ", svalue_size
					allocate(character(svalue_size) :: svalue)
					status_ip = GetSelectedOutputValue(iphreeqc_id1, i, j, type, dvalue, svalue, slength)
				enddo
			endif
			!write (*,'(2I5, " slength = ", I0, " len = ",I0,1x,A)') w, n, slength, len(svalue), trim(svalue)
            out = trim(out) // tab //trim(adjustl(svalue(1:slength)))
         enddo
         n = n + 1
         !write (unit,'(I0, A, I0, A, A)')  n, tab, w, trim(out)
		 write(*,*) trim(out) !dlp
      enddo
   enddo
   ! timing
   t1 = omp_get_wtime()
   write (*, '("Time for processing output", t55, F8.3," sec" )') t1 - t0
   t0 = t1;
   
!$OMP PARALLEL DO
   do w = 0, tc - 1
      iphreeqc_id1 = RM_GetIPhreeqcId(id, w)
      !!status_ip = RunString(iphreeqc_id1, "DELETE; -all") ! throws an exception
   enddo
!$OMP END PARALLEL DO

   status_RM = RM_CloseFiles(id)

   close (unit)

   ! timing
   t1 = omp_get_wtime()
   write (*, '("Time for deleting worker definitions", t55, F8.3," sec" )') t1 - t0
   write (*, '("Total time before C++ cleanup ", t55, F8.3," sec" )') t1 - t_start
   t0 = t1;
   Write (*,'(/"Done."/)')
   

   end program PhreePlotOpenMP
