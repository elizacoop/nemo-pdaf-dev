!> Utility Routines
!!
!! This module contains several routines useful for common
!! model tasks. The initial routines included output configuration
!! information about the PDAF library, and configuration information
!! about the assimilation parameters.
!!
module utils_pdaf
   use mod_kind_pdaf
   implicit none
   save
contains
  !> This routine performs a model-sided screen output about
  !! the configuration of the data assimilation system.
  !!
  !! **Calling Sequence**
  !!
  !! - Called from: `init_pdaf`
  !!
   subroutine init_info_pdaf()

      use config_pdaf, only: filtertype, subtype, dim_ens, forget
      ! *****************************
      ! *** Initial Screen output ***
      ! *****************************
      if (filtertype == 1) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: SEIK'
         if (subtype == 2) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed error-space basis'
         else if (subtype == 3) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed state covariance matrix'
         else if (subtype == 4) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- use ensemble transformation'
         else if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      else if (filtertype == 2) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: EnKF'
         if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      else if (filtertype == 3) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LSEIK'
         if (subtype == 2) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed error-space basis'
         else if (subtype == 3) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- fixed state covariance matrix'
         else if (subtype == 4) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- use ensemble transformation'
         else if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      else if (filtertype == 4) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: ETKF'
         if (subtype == 0) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant using T-matrix'
         else if (subtype == 1) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant following Hunt et al. (2007)'
         else if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      else if (filtertype == 5) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LETKF'
         if (subtype == 0) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant using T-matrix'
         else if (subtype == 1) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Variant following Hunt et al. (2007)'
         else if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      else if (filtertype == 6) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: ESTKF'
         if (subtype == 0) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Standard mode'
         else if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      else if (filtertype == 7) then
         write (*, '(a,21x, a)') 'NEMO-PDAF', 'Filter: LESTKF'
         if (subtype == 0) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Standard mode'
         else if (subtype == 5) then
            write (*, '(a,6x, a)') 'NEMO-PDAF', '-- Offline mode'
         end if
         write (*, '(a,14x, a, i5)') 'NEMO-PDAF', 'ensemble size:', dim_ens
         write (*, '(a,10x, a, f5.2)') 'NEMO-PDAF', 'forgetting factor:', forget
      end if
   end subroutine init_info_pdaf

   !----------------------------------------------------------------------------
   !> This routine reads the namelist file with parameters
   !! controlling data assimilation with PDAF and outputs to
   !! screen.
   !!
   !! **Calling Sequence**
   !!
   !! - Called from: `init_pdaf`
   !!
   subroutine read_config_pdaf()
      use config_pdaf, only: pdaf_nml
      use io_pdaf, only: io_nml, add_slash
      use obs_ssh_mgrid_pdafomi, only: mgrid_ssh_nml
      use obs_sit_mgrid_pdafomi, only: mgrid_sit_nml
      use obs_fbd_mgrid_pdafomi, only: mgrid_fbd_nml
      use statevector_pdaf, only: sv_nml
      ! ****************************************************
      ! ***   Initialize PDAF parameters from namelist   ***
      ! ****************************************************
      open (20, file='namelist_cfg.pdaf')
      rewind(20)
      read (20, NML=pdaf_nml)
      rewind(20)
      read (20, NML=io_nml)
      rewind(20)
      read (20, NML=sv_nml)
      rewind(20)
      read(20, NML=mgrid_ssh_nml)
      rewind (20)
      read(20, NML=mgrid_sit_nml)
      rewind (20)
      read(20, NML=mgrid_fbd_nml)
      close (20)


   end subroutine read_config_pdaf

   subroutine print_config()
      use parallel_pdaf, only: mype_ens
      use config_pdaf, only: print_pdaf_configuration
      use io_pdaf, only: print_io_configuration
      use obs_ssh_mgrid_pdafomi, only: print_ssh_mgrid_configuration
      use obs_sit_mgrid_pdafomi, only: print_sit_mgrid_configuration
      use obs_fbd_mgrid_pdafomi, only: print_fbd_mgrid_configuration
      ! Print PDAF parameters to screen
      showconf: if (mype_ens == 0) then
         write (*, '(/a,1x,a)') 'NEMO-PDAF','-- Overview of PDAF configuration --'
         call print_pdaf_configuration()
         call print_io_configuration()
         call print_ssh_mgrid_configuration()
         call print_sit_mgrid_configuration()
         call print_fbd_mgrid_configuration()
         write (*, '(a,1x,a/)') 'NEMO-PDAF','-- End of PDAF configuration overview --'
      end if showconf
   end subroutine print_config


   !-------------------------------------------------------------------------------

   !> Timing and clean-up of PDAF
   !!
   !! The routine prints timing and memory information.
   !! It further deallocates PDAF internal arrays
   !! and the ASMINC increment array for BGC variables
   !!
   !! - Called from: `nemogcm`
   !!
   subroutine finalize_pdaf()
      use PDAF, only: PDAF_print_info, PDAF_deallocate
      use parallel_pdaf, only: mype_ens, npes_ens, comm_ensemble, mpierr
      use timer, only: timeit, time_tot
      implicit none
      ! Show allocated memory for PDAF
      if (mype_ens==0) call PDAF_print_info(10)
      if (npes_ens>1) call PDAF_print_info(11)
      ! Print PDAF timings onto screen
      if (mype_ens==0) call PDAF_print_info(3)
      ! Deallocate PDAF arrays
      call PDAF_deallocate()

      call timeit(4,'old')
      call timeit(5,'old')

      if (mype_ens==0) then
         write (*, '(/a,10x,a)') 'NEMO-PDAF', 'Model-sided timing overview'
         write (*, '(a,2x,a)') 'NEMO-PDAF', '-----------------------------------'
         write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize MPI:  ', time_tot(1), 's'
         write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize model:', time_tot(2), 's'
         write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'initialize PDAF :', time_tot(3), 's'
         write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'main part:       ', time_tot(4), 's'
         write (*, '(a,8x,a,F11.3,1x,a)') 'NEMO-PDAF', 'total:         ', time_tot(5), 's'
      end if

      call mpi_barrier(comm_ensemble, mpierr)

   end subroutine finalize_pdaf

end module utils_pdaf
