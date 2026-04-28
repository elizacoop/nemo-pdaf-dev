!> Module holding IO operations for NEMO-PDAF
!!
!! This code bases in wide parts on the implementation
!! by Wibke Duesterhoeft-Wriggers, BSH, Germany for the
!! CMEMS Baltic Monitoring and forecasting center
!!
module io_pdaf
   use nemo_pdaf, only: nn_time0
   use mod_kind_pdaf
   use mpi

   implicit none
   save

   integer :: verbose_io=0   ! Set verbosity of IO routines (0,1,2,3)

   ! Control of IO
   logical :: save_ens_sngl=.false.           ! write set of files holding ensemble of selected field
   logical :: do_deflate=.false.              ! Deflate variables in NC files (this seems to fail for parallel nc)
   character(len=3) :: sgldbl_io='sgl'        ! Write PDAF output in single (sgl) or double (dbl) precision

   character(len=lc) :: fname_dom               ! Name of domain file
   character(len=lc) :: path_dom                ! Path for NEMO file holding dimensions
   character(len=lc) :: f_basename_rst          ! Name of restart file
   ! the parth to restart files is constructed by the following variables
   ! path_rst_root//ens_prefix[1-N]//path_rst_suffix
   character(len=lc) :: path_rst_root           ! Path for NEMO file holding dimensions
   character(len=lc) :: ens_prefix = 'ens_'
   character(len=lc) :: path_rst_suffix
   ! the parth to asm files is constructed by the following variables
   ! path_asm_root//ens_prefix[1-N]/assim_background_increments_xxxx.nc
   ! path_asm_root//ens_prefix[1-N]/assim_background_state_DI_xxxx.nc
   character(len=lc) :: path_asm_root              ! Path for increment files
   ! *** temporary variables for IO
   real(pwp), allocatable :: tmp_4d(:,:,:,:)     ! 4D array used to represent full NEMO grid box
                                                 ! the dimension should be (x,y,z,time)
   real(4),   allocatable :: stmp_4d(:,:,:,:)    ! 4D array used to represent full NEMO grid box
                                                 ! the dimension should be (x,y,z,time)

   namelist /io_nml/ verbose_io, sgldbl_io, path_dom, fname_dom, path_asm_root, &
                     path_rst_root, ens_prefix, path_rst_suffix,f_basename_rst, &
                     nn_time0

contains
   !> Print configuration of IO module
   !!
   SUBROUTINE print_io_configuration()
      implicit none
      character(len=256) :: path_rst ! path to restart files
      ! Construct restart file path
      call add_slash(path_rst_root)
      write(path_rst, '(a,i0,"/",a)') TRIM(path_rst_root)//TRIM(ens_prefix),1,TRIM(path_rst_suffix)
      call add_slash(path_rst)
      ! Print PDAF IO configuration to screen
      write (*, '(a,3x,a)')     'NEMO-PDAF','[io_nml]:'
      write (*, '(a,5x,a,i10)') 'NEMO-PDAF','verbose_io                   ', verbose_io
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','sgldbl_io                    ', trim(sgldbl_io)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path to model domain file    ', trim(path_dom)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','model domain filename        ', trim(fname_dom)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path to restart files        ', trim(path_rst)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','basename of restart files    ', trim(f_basename_rst)
      write (*, '(a,5x,a,6x,a)')'NEMO-PDAF','path to write increment files', trim(path_asm_root)
   END SUBROUTINE print_io_configuration

   !> Read domain local information from restart files
   !!
   SUBROUTINE read_local_domain()
      use mpi
      USE netcdf
      use config_pdaf, only: screen
      use mod_memcount_pdaf, only: memcount
      use nemo_pdaf, only: i0, j0, ni_p, nj_p, nk_p, nav_lat, nav_lon, nav_lev,&
                           halo0, halo1, time_counter, ndastp
      use parallel_pdaf, only: mype_model, npes_model, comm_model,MPIerr
      IMPLICIT NONE
      ! Local variables
      character(len=lc) :: fname ! file name
      character(len=lc) :: path_rst ! path to restart files
      INTEGER :: ncid             ! netCDF file identifier
      INTEGER :: varid            ! variable identifier

      integer :: dom_size_local(2)
      integer :: dom_pos_first(2)

      ! Construct restart file path
      call add_slash(path_rst_root)
      write(path_rst, '(a,i0,"/",a)') TRIM(path_rst_root)//TRIM(ens_prefix),1,TRIM(path_rst_suffix)
      call add_slash(path_rst)
      write(fname, '(a,i4.4)') TRIM(f_basename_rst)//'_', mype_model
      !!----------------------------------------------------------------------
      if (mype_model == 0 .and. verbose_io>0) then
         WRITE(*,*)
         WRITE(*,'(/a,1x,a)') 'read_local_domain : Reading local domain information from file'
         WRITE(*,'(a,1x,a/)') '   Input file: ', TRIM(fname)//'.nc'
      end if
      ! Open the NetCDF file
      call check(nf90_open( trim(path_rst)//trim(fname)//'.nc', NF90_NOWRITE, ncid ))
      ! Read 2D grid variables for local domain
      ! dom_size_local
      call check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_size_local', dom_size_local ))
      ! dom_pos_first
      call check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_position_first', dom_pos_first ))
      ! dom_pos_first
      call check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', halo0 ))
      ! dom_pos_first
      call check(nf90_get_att( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', halo1 ))
      !ndastp
      call check(nf90_inq_varid( ncid, 'ndastp', varid ))
      call check(nf90_get_var( ncid, varid, ndastp))

      i0 = dom_pos_first(1)
      j0 = dom_pos_first(2)
      ni_p = dom_size_local(1)
      nj_p = dom_size_local(2)
      ! jpk (z dimension)
      call check(nf90_inq_dimid( ncid, 'nav_lev', varid ))
      call check(nf90_inquire_dimension( ncid, varid, len=nk_p ))
      ! nav_lon
      allocate( nav_lon(ni_p, nj_p) )
      allocate( nav_lat(ni_p, nj_p) )
      ALLOCATE( nav_lev(nk_p) )
      call memcount(1, 'r', 2*ni_p*nj_p + nk_p )
      call check(nf90_inq_varid( ncid, 'nav_lon', varid ))
      call check(nf90_get_var( ncid, varid, nav_lon, [1, 1], [ni_p, nj_p] ))
      ! nav_lat
      call check(nf90_inq_varid( ncid, 'nav_lat', varid ))
      call check(nf90_get_var( ncid, varid, nav_lat, [1, 1], [ni_p, nj_p] ))
      ! nav_lev
      call check(nf90_inq_varid( ncid, 'nav_lev', varid ))
      call check(nf90_get_var( ncid, varid, nav_lev) )
      ! time_counter
      call check(nf90_inq_varid( ncid, 'time_counter', varid ))
      call check(nf90_get_var( ncid, varid, time_counter, [1], [1] ))
      ! Close the NetCDF files
      call check (nf90_close( ncid ))
      ! Screen output
      if (npes_model>1 .and. screen>0) then
         if (mype_model == 0) then
            write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:'
            write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
               'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
         end if
         call MPI_Barrier(comm_model, MPIerr)
         write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7/)') 'NEMO-PDAF', 'RANK', mype_model, i0, i0+ni_p-1, j0, j0+nj_p-1, ni_p, nj_p
      end if
      call MPI_Barrier(comm_model, MPIerr)
   END SUBROUTINE read_local_domain

   !> Read global domain information from restart files
   !!
   SUBROUTINE read_global_domain()
      use mpi
      USE netcdf
      use config_pdaf, only: screen
      use mod_memcount_pdaf, only: memcount
      use nemo_pdaf, only: jpiglo, jpjglo, jpk, i0, j0, ni_p, nj_p, nk_p, &
                           glamt, glamu, glamv, gphit, gphiu, gphiv, tmask
      use parallel_pdaf, only: mype_model, npes_model, comm_model, MPIerr
      IMPLICIT NONE
      ! Local variables
      INTEGER :: ncid       ! netCDF file identifier
      INTEGER :: varid      ! variable identifier
      INTEGER :: i, j, iktop, ikbot ! counter
      integer, allocatable :: k_top(:, :), k_bot(:, :) ! top and bottom wet levels
      !!----------------------------------------------------------------------
      if (mype_model == 0 .and. verbose_io>0) then
         WRITE(*,*)
         WRITE(*,*) 'NEMO-PDAF', 'read_grid_variables : Reading grid variables from file'
         WRITE(*,*) 'NEMO-PDAF', '   Input file: ', TRIM(fname_dom)
      end if
      call add_slash(path_dom)
      ! Open the NetCDF file
      call check(nf90_open( trim(path_dom)//trim(fname_dom), NF90_NOWRITE, ncid ))
      ! Read dimension sizes
      ! jpiglo (x dimension)
      call check(nf90_inq_dimid( ncid, 'x', varid ))
      call check(nf90_inquire_dimension( ncid, varid, len=jpiglo ))
      ! jpjglo (y dimension)
      call check(nf90_inq_dimid( ncid, 'y', varid ))
      call check(nf90_inquire_dimension( ncid, varid, len=jpjglo ))
      ! jpk (z dimension)
      call check(nf90_inq_dimid( ncid, 'z', varid ))
      call check(nf90_inquire_dimension( ncid, varid, len=jpk ))
      nk_p = jpk
      ! Allocate arrays with dimensions (time, y, x)
      ALLOCATE( glamt(ni_p, nj_p) )
      ALLOCATE( glamu(ni_p, nj_p) )
      ALLOCATE( glamv(ni_p, nj_p) )
      ALLOCATE( gphit(ni_p, nj_p) )
      ALLOCATE( gphiu(ni_p, nj_p) )
      ALLOCATE( gphiv(ni_p, nj_p) )
      call memcount(1, 'r', 6*ni_p*nj_p)
      ! Read 3D grid variables
      ! glamt
      call check(nf90_inq_varid( ncid, 'glamt', varid ))
      call check(nf90_get_var( ncid, varid, glamt, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! glamu
      call check(nf90_inq_varid( ncid, 'glamu', varid ))
      call check(nf90_get_var( ncid, varid, glamu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! glamv
      call check(nf90_inq_varid( ncid, 'glamv', varid ))
      call check(nf90_get_var( ncid, varid, glamv, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphit
      call check(nf90_inq_varid( ncid, 'gphit', varid ))
      call check(nf90_get_var( ncid, varid, gphit, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphiu
      call check(nf90_inq_varid( ncid, 'gphiu', varid ))
      call check(nf90_get_var( ncid, varid, gphiu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      call check(nf90_get_var( ncid, varid, gphiu, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! gphiv
      call check(nf90_inq_varid( ncid, 'gphiv', varid ))
      call check(nf90_get_var( ncid, varid, gphiv, [i0, j0, 1], [ni_p, nj_p, 1] ))
      ! calculate t_mask
      allocate( k_top(ni_p, nj_p) )
      allocate( k_bot(ni_p, nj_p) )
      allocate( tmask(ni_p, nj_p, nk_p) )
      call memcount(1, 'r', ni_p*nj_p*nk_p )
      ! k_top
      call check(nf90_inq_varid( ncid, 'top_level', varid ))
      call check(nf90_get_var( ncid, varid, k_top, [i0, j0], [ni_p, nj_p] ))
      ! k_bot
      call check(nf90_inq_varid( ncid, 'bottom_level', varid ) )
      call check(nf90_get_var( ncid, varid, k_bot, [i0, j0], [ni_p, nj_p] ))
      ! k_top and k_bot
      tmask(:,:,:) = 0._pwp
      DO j = 1, nj_p
         DO i = 1, ni_p
            iktop = k_top(i,j)
            ikbot = k_bot(i,j)
            IF( iktop /= 0 ) THEN       ! water in the column
               tmask(i, j, iktop:ikbot  ) = 1._pwp
            ENDIF
         END DO
      END DO
      ! deallocate k_top and k_bot
      deallocate( k_top, k_bot )
      ! Close the NetCDF file
      call check (nf90_close( ncid ))
      ! *** Screen output ***
      if (mype_model==0 .and. screen>0) then
         write (*,'(/a,5x,a)') 'NEMO-PDAF', '*** NEMO: grid dimensions ***'
         write(*,'(a,3x,2(6x,a),9x,a)') 'NEMO-PDAF', 'jpiglo','jpjglo','jpk'
         write(*,'(a,3x,3i12)') 'NEMO-PDAF', jpiglo, jpjglo, jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Dimension of global 3D grid box', jpiglo*jpjglo*jpk
         write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Number of global surface points', jpiglo*jpjglo
      end if
      !
      call MPI_Barrier(comm_model, MPIerr)
      if (npes_model>1 .and. screen>1) then
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
               'NEMO-PDAF', 'PE', mype_model, 'Dimension of local 3D grid box', ni_p*nj_p*nk_p
         write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
               'NEMO-PDAF', 'PE', mype_model, 'Number of local surface points', ni_p * nj_p
      end if
   END SUBROUTINE read_global_domain

   !> Read restart file to form state vector
   !!
   SUBROUTINE read_restart(ens_member, state_p)
      USE netcdf
      use mod_memcount_pdaf, only: memcount
      use nemo_pdaf, only: ni_p, nj_p, nk_p, numcat
      use parallel_pdaf, only: mype_model
      use statevector_pdaf, only: sfields, n_fields
      use transforms_pdaf, only: field2state
      IMPLICIT NONE
      !*** Arguments ***
      integer,   intent(in)    :: ens_member !< Ensemble member index
      real(pwp), intent(inout) :: state_p(:) !< State vector
      ! Local variables
      character(len=lc) :: fname ! file name
      character(len=lc) :: path_rst ! path to restart files
      INTEGER :: ncid              ! netCDF file identifier
      INTEGER :: varid             ! variable identifier
      integer :: i, k              ! counter
      integer :: nk                ! number of levels for given variable
      integer :: start_k           ! starting level for reading variable (for handling sea ice categories)
      integer :: ndims             ! number of dimensions of variable in restart file
      ! Construct restart file path
      call add_slash(path_rst_root)
      write(path_rst, '(a,i0,"/",a)') TRIM(path_rst_root)//TRIM(ens_prefix), &
                                      ens_member,TRIM(path_rst_suffix)
      call add_slash(path_rst)

      if (verbose_io>0 .and. mype_model==0) &
            write(*,'(a,4x,a)') 'NEMO-PDAF','*** Ensemble: Reading model restart file'

      if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

      ! Initialize state
      state_p = 0.0_pwp

      do i = 1, n_fields
         write(fname, '(a,i4.4)') TRIM(sfields(i)%rst_file)//'_', mype_model
         if (verbose_io>1 .and. mype_model==0) then
            write(*,'(a,2x,a)') 'NEMO-PDAF', 'Reading: '//trim(path_rst)//trim(fname)//'.nc'
            write (*,'(a,i5,1x,a,a,a,i10)') &
                  'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
         end if
         ! Open the file
         call check( nf90_open(trim(path_rst)//trim(fname)//'.nc', nf90_nowrite, ncid) )
         !  Read field
         call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_n), varid) )
         call check( nf90_inquire_variable(ncid, varid, ndims=ndims) )
         if (ndims == 4) then
            start_k = 1
            if (sfields(i)%k_name == 'cat') then
               nk = numcat
               ! specific operation for selecting sea ice category to reduce repeating reading
               if (trim(adjustl(sfields(i)%operation)) == 'select_cat_ice') then
                  if (sfields(i)%ice_cat > 0 .and. sfields(i)%ice_cat <= numcat) then
                     nk = 1
                     start_k = sfields(i)%ice_cat
                  else
                     write(*,'(a,2x,a,2x,i0,2x,a)') 'NEMO-PDAF', &
                        'Invalid sea ice category (', sfields(i)%ice_cat, ') selected for variable: '// &
                        trim(sfields(i)%variable)
                  end if
               end if
            else
               nk = nk_p
            end if
            call check( nf90_get_var(ncid, varid, tmp_4d(:, :, :nk, 1), &
                  start=[1, 1, start_k, 1], count=[ni_p, nj_p, nk, 1]) )
         else
            call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
                  start=[1, 1, 1], count=[ni_p, nj_p, 1]) )
         end if
         ! operations to form state vector
         if (trim(adjustl(sfields(i)%operation)) == 'sum_over_cat_ice') then
            do k = 2, nk
               tmp_4d(:,:,1,1) = tmp_4d(:,:,1,1) + tmp_4d(:,:,k,1)
            end do
         else if (trim(adjustl(sfields(i)%operation)) == 'select_cat_ice') then
         else if (trim(adjustl(sfields(i)%operation)) == '') then
         else
            write(*,'(a,2x,a)') 'NEMO-PDAF', 'Unknown operation for combining variables in restart file: '// &
                  trim(sfields(i)%operation)//'. Using first names_rest_n only.'
         end if
         ! Close the file
         call check( nf90_close(ncid) )
         ! Convert field to state vector
         call field2state(tmp_4d, state_p, sfields(i)%off, sfields(i)%ndims)
      end do

      if (verbose_io>2) then
         do i = 1, n_fields
            write(*,*) 'Min and max for ',trim(sfields(i)%variable),' :     ', &
                  minval(state_p(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
                  maxval(state_p(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
         enddo
      end if
   END SUBROUTINE read_restart
   !===========================================================================
   !> Write NEMO ASM increment compatible for IAU
   !!
   subroutine write_asminc_mv(ens_member, state, state_f)
      use netcdf
      use config_pdaf, only: screen
      use nemo_pdaf, only: ni_p, nj_p, nk_p, jpiglo, jpjglo, halo0, halo1, &
                           i0, j0, nav_lat, nav_lon, time_counter, nav_lev, ndastp
      use parallel_pdaf, only: mype_model, npes_model
      use statevector_pdaf, only: n_fields, sfields
      use transforms_pdaf, only: state2field, transform_field_mv
      implicit none
      ! *** Arguments ***
      integer,          intent(in) :: ens_member    ! Ensemble member index
      real(pwp),        intent(inout) :: state(:)   ! Analysis state vector
      real(pwp),        intent(in) :: state_f(:) ! Forecast state vector
      ! *** Local variables ***
      integer :: ncid
      integer :: dimids_field(4)
      integer :: i
      integer :: dimid_time, dimid_lvls, dimid_lat, dimid_lon
      integer :: id_dateb, id_datef
      integer :: id_lat, id_lon, id_lev, id_time, id_time_counter, id_incr
      integer :: startC(2), countC(2)
      integer :: startt(4), countt(4)
      integer :: nf_prec      ! Precision for netcdf output of model fields
      integer :: verbose
      real(pwp) :: time
      character(len=lc) :: filename
      character(len=lc) :: path_asm ! path to increment file
      ! **********************
      ! *** Initialization ***
      ! **********************
      if (verbose_io>0 .and. mype_model==0) write (*,'(8x,a)') '--- Write increment file'
      ! *** Set increment times ***
      time = ndastp + real(nn_time0, pwp)*0.0001_pwp
      ! Prepare file writing
      if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))
      nf_prec = NF90_DOUBLE
      ! *****************************
      ! *** Create and write file ***
      ! *****************************
      ! *** Create file ***
      call add_slash(path_asm_root)
      write(path_asm, '(a,i0,"/",a)') TRIM(path_asm_root)//TRIM(ens_prefix), &
                                      ens_member, 'assim_background_increments'
      write(filename, '(a,i4.4,a)') TRIM(path_asm)//'_', mype_model, '.nc'

      if (verbose_io>0 .and. mype_model==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Create file: ', trim(filename)

      call check( NF90_CREATE(trim(filename), NF90_NETCDF4, ncid))
      ! define global attributes
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'title', &
                               'Increment for NEMO-PDAF data assimilation'))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_number_total', &
                               npes_model))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_number', &
                               mype_model))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_dimensions_ids', &
                               [1, 2]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_size_global', &
                               [jpiglo, jpjglo]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_size_local', &
                               [ni_p, nj_p]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_position_first', &
                               [i0, j0]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_position_last', &
                               [i0 + ni_p - 1, j0 + nj_p - 1]))
      CALL check(NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start',&
                               halo0))
      CALL check(NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end'  ,&
                               halo1))
      CALL check(NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_type'           ,&
                               'BOX'))

      ! define dimensions for NEMO-input file
      ! define spatial dimensions
      call check( NF90_DEF_DIM(ncid, 'x', ni_p, dimid_lon) )
      call check( NF90_DEF_DIM(ncid, 'y', nj_p, dimid_lat) )
      call check( NF90_DEF_DIM(ncid, 'nav_lev', nk_p, dimid_lvls))
      call check( NF90_DEF_DIM(ncid,'time_counter', NF90_UNLIMITED, dimid_time))

      dimids_field=[dimid_lon, dimid_lat, dimid_lvls, dimid_time]
      ! define variables
      call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
      call check( NF90_DEF_VAR(ncid, 'z_inc_dateb', NF90_DOUBLE, id_dateb))
      call check( NF90_DEF_VAR(ncid, 'z_inc_datef', NF90_DOUBLE, id_datef))
      call check( NF90_DEF_VAR(ncid, 'time_counter', NF90_DOUBLE, dimids_field(4), id_time_counter))
      call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
      call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
      call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
      if (do_deflate) then
         call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
         call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
         call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
      end if

      do i = 1, n_fields
         if (sfields(i)%ndims==3) then
            dimids_field(3)=dimid_lvls
            call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), nf_prec, &
                                     dimids_field, id_incr) )
         else
            dimids_field(3)=dimid_time
            call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), &
                                     nf_prec, dimids_field(1:3), id_incr) )
         end if
         if (do_deflate) &
               call check( NF90_def_var_deflate(ncid, id_incr, 0, 1, 1) )
      end do

      ! End define mode
      call check( NF90_ENDDEF(ncid) )

      ! write coordinates
      startC = [1, 1]
      countC = [ni_p, nj_p]

      call check( nf90_put_var(ncid, id_lon, nav_lon, startC, countC))
      call check( nf90_put_var(ncid, id_lat, nav_lat, startC, countC))
      call check( nf90_put_var(ncid, id_lev, nav_lev, [1], [nk_p]))
      call check( nf90_put_var(ncid, id_time_counter, time_counter, start=[1], count=[1]))
      ! keep all time attributes identical as NEMO assigns dateb=datef=time
      call check( nf90_put_var(ncid, id_time, time))
      call check( nf90_put_var(ncid, id_dateb, time))
      call check( nf90_put_var(ncid, id_datef, time))
      ! *** Write fields
      ! Backwards transformation of state fields
      if (mype_model==0) then
         verbose = screen
      else
         verbose = 0
      end if
      call transform_field_mv(2, state, 0, verbose)

      ! Compute increment
      state = state - state_f
      ! Write each updated field
      startt = [1, 1, 1, 1]
      countt = [ni_p, nj_p, nk_p, 1]
      do i = 1, n_fields
         tmp_4d = 0._pwp
         call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)
         if (verbose_io>1 .and. mype_model==0) &
         write (*,'(a,1x,a,a)') 'NEMO-PDAF', '--- write variable: ', trim(sfields(i)%variable)
         call check( nf90_inq_varid(ncid, trim(sfields(i)%name_incr), id_incr) )
         if (sfields(i)%ndims==3) then
            countt(3) = nk_p
            call check( nf90_put_var(ncid, id_incr, tmp_4d, startt, countt))
         else
            countt(3) = 1
            call check( nf90_put_var(ncid, id_incr, tmp_4d, startt(1:3), countt(1:3)))
         end if
      end do
      ! *** close file with state sequence ***
      call check( NF90_CLOSE(ncid) )
   end subroutine write_asminc_mv
   !===========================================================================
   !> Write NEMO ASM compatible for IAU
   !!
   subroutine write_asmdin_mv(ens_member, state_f)
      use netcdf
      use config_pdaf, only: screen
      use nemo_pdaf, only: ni_p, nj_p, nk_p, jpiglo, jpjglo, halo0, halo1, &
                           i0, j0, nav_lat, nav_lon, time_counter, nav_lev, ndastp
      use parallel_pdaf, only: mype_model, npes_model
      use statevector_pdaf, only: n_fields, sfields
      use transforms_pdaf, only: state2field, transform_field_mv
      implicit none
      ! *** Arguments ***
      integer,          intent(in) :: ens_member    ! Ensemble member index
      real(pwp),        intent(in) :: state_f(:) ! Forecast state vector
      ! *** Local variables ***
      integer :: ncid
      integer :: dimids_field(4)
      integer :: i
      integer :: dimid_time, dimid_lvls, dimid_lat, dimid_lon
      integer :: id_rdastp
      integer :: id_lat, id_lon, id_lev, id_time_counter, id_bkg
      integer :: startC(2), countC(2)
      integer :: startt(4), countt(4)
      integer :: nf_prec      ! Precision for netcdf output of model fields
      character(len=lc) :: filename
      character(len=lc) :: path_asm ! path to increment file
      ! **********************
      ! *** Initialization ***
      ! **********************
      if (verbose_io>0 .and. mype_model==0) write (*,'(8x,a)') '--- Write increment file'
      ! Prepare file writing
      if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))
      nf_prec = NF90_DOUBLE
      ! *****************************
      ! *** Create and write file ***
      ! *****************************
      ! *** Create file ***
      call add_slash(path_asm_root)
      write(path_asm, '(a,i0,"/",a)') TRIM(path_asm_root)//TRIM(ens_prefix), &
                                      ens_member, 'assim_background_state_DI'
      write(filename, '(a,i4.4,a)') TRIM(path_asm)//'_', mype_model, '.nc'

      if (verbose_io>0 .and. mype_model==0) &
            write (*,'(a,1x,a,a)') 'NEMO-PDAF', 'Create file: ', trim(filename)

      call check( NF90_CREATE(trim(filename), NF90_NETCDF4, ncid))
      ! define global attributes
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'title', &
                               'Increment for NEMO-PDAF data assimilation'))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_number_total', &
                               npes_model))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_number', &
                               mype_model))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_dimensions_ids', &
                               [1, 2]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_size_global', &
                               [jpiglo, jpjglo]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_size_local', &
                               [ni_p, nj_p]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_position_first', &
                               [i0, j0]))
      call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'DOMAIN_position_last', &
                               [i0 + ni_p - 1, j0 + nj_p - 1]))
      CALL check(NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start',&
                               halo0))
      CALL check(NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end'  ,&
                               halo1))
      CALL check(NF90_PUT_ATT( ncid, NF90_GLOBAL, 'DOMAIN_type'           ,&
                               'BOX'))

      ! define dimensions for NEMO-input file
      ! define spatial dimensions
      call check( NF90_DEF_DIM(ncid, 'x', ni_p, dimid_lon) )
      call check( NF90_DEF_DIM(ncid, 'y', nj_p, dimid_lat) )
      call check( NF90_DEF_DIM(ncid, 'nav_lev', nk_p, dimid_lvls))
      call check( NF90_DEF_DIM(ncid,'time_counter', NF90_UNLIMITED, dimid_time))

      dimids_field=[dimid_lon, dimid_lat, dimid_lvls, dimid_time]
      ! define variables
      call check( NF90_DEF_VAR(ncid, 'rdastp', NF90_DOUBLE, id_rdastp))
      call check( NF90_DEF_VAR(ncid, 'time_counter', NF90_DOUBLE, dimid_time, id_time_counter))
      call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
      call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
      call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
      if (do_deflate) then
         call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
         call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
         call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
      end if

      do i = 1, n_fields
         if (sfields(i)%ndims==3) then
            dimids_field(3)=dimid_lvls
            call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_bkg_din), &
                                     nf_prec, dimids_field(1:4), id_bkg) )
         else
            dimids_field(3)=dimid_time
            call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_bkg_din), &
                                     nf_prec, dimids_field(1:3), id_bkg) )
         end if
         if (do_deflate) &
               call check( NF90_def_var_deflate(ncid, id_bkg, 0, 1, 1) )
      end do

      ! End define mode
      call check( NF90_ENDDEF(ncid) )

      ! write coordinates
      startC = [1, 1]
      countC = [ni_p, nj_p]

      call check( nf90_put_var(ncid, id_lon, nav_lon, startC, countC))
      call check( nf90_put_var(ncid, id_lat, nav_lat, startC, countC))
      call check( nf90_put_var(ncid, id_lev, nav_lev, [1], [nk_p]))
      call check( nf90_put_var(ncid, id_time_counter, time_counter, start=[1], count=[1]))
      ! keep all time attributes identical as NEMO assigns dateb=datef=time
      call check( nf90_put_var(ncid, id_rdastp, ndastp))
      ! *** Write fields
      ! Write each updated field
      startt = [1, 1, 1, 1]
      countt = [ni_p, nj_p, nk_p, 1]
      do i = 1, n_fields
         tmp_4d = 0._pwp
         call state2field(state_f, tmp_4d, sfields(i)%off, sfields(i)%ndims)
         if (verbose_io>1 .and. mype_model==0) &
         write (*,'(a,1x,a,a)') 'NEMO-PDAF', '--- write bkg variable: ', trim(sfields(i)%variable)
         call check( nf90_inq_varid(ncid, trim(sfields(i)%name_bkg_din), id_bkg) )
         if (sfields(i)%ndims==3) then
            countt(3) = nk_p
            call check( nf90_put_var(ncid, id_bkg, tmp_4d, startt, countt))
         else
            countt(3) = 1
            call check( nf90_put_var(ncid, id_bkg, tmp_4d, startt(1:3), countt(1:3)))
         end if
      end do
      ! *** close file with state sequence ***
      call check( NF90_CLOSE(ncid) )
   end subroutine write_asmdin_mv
   !============================================================================
   !> Check status of NC operation
   !!
   subroutine check(status)
      use netcdf
      use parallel_pdaf, only: abort_parallel
      ! *** Arguments ***
      ! Reading status
      integer, intent ( in) :: status
      ! end program with error message if status is not nf90_noerr
      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         call abort_parallel()
      end if
   end subroutine check
   ! ===========================================================================
   !> Add a trailing slash to a path string
   !!
   !! This routine ensures that a string defining a path
   !! has a trailing slash.
   !!
   subroutine add_slash(path)
      implicit none
      ! *** Arguments ***
      !< String holding the path
      character(len=100) :: path
      ! *** Local variables ***
      integer :: strlength
      ! *** Add trailing slash ***
      strlength = len_trim(path)
      if (path(strlength:strlength) /= '/') then
         path = trim(path) // '/'
      end if
   end subroutine add_slash
   !============================================================================
   !> Convert an integer to a strong of length 4
   !!
   character(len=4) function str(k)
      implicit none
      !< number
      integer, intent(in) :: k
      ! string representation
      write (str, '(i4.4)') k
   end function str
   !============================================================================
   !> Check whether a file exists
   !!
   function file_exists(filename) result(res)
      implicit none
      !< File name
      character(len=*),intent(in) :: filename
      !< Status of file
      logical                     :: res
      ! Check if the file exists
      inquire( file=trim(filename), exist=res )
   end function file_exists

end module io_pdaf
