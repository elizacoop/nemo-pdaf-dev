!> PDAF-OMI observation module for ssh observations (on model grid)
!!
!! Observation type: SSH on model grid
!!
!! The subroutines in this module are for the particular handling of
!! ssh observations available on the model grid. The observation module
!! also allows to perform a twin experiment in which model output
!! is read and used as observations after adding noise to the values.
!!
!! The routines are called by the different call-back routines of PDAF.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initialized by the generic routines from `PDAFomi`.
!!
!! Author: Nicholas Byrne, NCEO & University of Reading, UK
!!
module obs_ssh_mgrid_pdafomi
   use mod_kind_pdaf
   use parallel_pdaf, only: mype_filter, abort_parallel
   use PDAF, only: obs_f, obs_l
   use nemo_pdaf, only: i0, j0
   use netcdf
   implicit none
   save

   !> Whether to assimilate this data type
   logical :: assim_ssh_mgrid = .false.
   !> (1) use global obs.; (0) use domain-reduced full obs.
   integer :: use_global_obs = 1
   !> index in sfields array
   integer :: id_sfields = 1
   !> Observation error standard deviation (for constant errors)
   real(pwp) :: rms = 0.1 !1
   !> Localization cut-off radius
   real(pwp) :: lradius = 1.0
   !> Support radius for weight function
   real(pwp) :: sradius = 1.0
   !> Standard deviation for Gaussian noise in twin experiment
   real(pwp) :: noise_amp = 1
   !> path to obs file
   character(lc) :: file_ssh_mgrid = 'ssh_obs.nc'
   !> Name of SIC variable in the observation file
   character(lc) :: varname_ssh_mgrid = 'sic'
   !> Time step index to read from the observation file
   integer :: nc_step = 1 !WILL NEED TO BE ABLE TO LINE THESE UP

   namelist /mgrid_ssh_nml/ assim_ssh_mgrid, use_global_obs, id_sfields, rms, &
                            lradius, sradius, noise_amp, file_ssh_mgrid, &
                            nc_step, varname_ssh_mgrid

   !> Instance of full observation data type - see `PDAFomi` for details.
   type(obs_f), target, public :: thisobs
   !> Instance of local observation data type - see `PDAFomi` for details.
   type(obs_l), target, public :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)

contains
   !> print configuration of this observation type to screen
   subroutine print_ssh_mgrid_configuration()
      write (*, '(a,3x,a)') 'NEMO-PDAF', '[obs_ssh_mgrid_nml]:'
      write (*, '(a,5x,a,5x,l)') 'NEMO-PDAF',      'assim_ssh_mgrid', assim_ssh_mgrid
      if (assim_ssh_mgrid) then
          write (*, '(a,5x,a,f12.4)')  'NEMO-PDAF','rms            ', rms
          write (*, '(a,5x,a,i0)')     'NEMO-PDAF','use_global_obs ', use_global_obs
          write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF', 'noise_amp      ', noise_amp
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','lradius        ', lradius
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','sradius        ', sradius
          write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'obs_file   ', trim(file_ssh_mgrid)
          write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'varname', trim(varname_ssh_mgrid)
          write (*, '(a,5x,a,i0)')   'NEMO-PDAF', 'nc_step          ', nc_step
      end if
   end subroutine print_ssh_mgrid_configuration

   !>Initialize information on the observation
   !!
   !! The routine is called by each filter process.
   !! at the beginning of the analysis step before
   !! the loop through all local analysis domains.
   !!
   !! It has to count the number of process-local and full
   !! observations, initialize the vector of observations
   !! and their inverse variances, initialize the coordinate
   !! array and index array for indices of observed elements
   !! of the state vector.
   !!
   !! The following four variables have to be initialized in this routine:
   !!
   !! - **thisobs%doassim** - Whether to assimilate ssh
   !! - **thisobs%disttype** - type of distance computation for localization
   !! with ssh
   !! - **thisobs%ncoord** - number of coordinates used for distance
   !! computation
   !! - **thisobs%id_obs_p** - index of module-type observation in PE-local state
   !! vector
   !!
   !!
   !! Optional is the use of:
   !! - **thisobs%icoeff_p**       - Interpolation coefficients for obs. operator (only if interpolation is used)
   !! - **thisobs%domainsize**     - Size of domain for periodicity for *disttype=1* (<0 for no periodicity)
   !! - **thisobs%obs_err_type**   - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
   !! - **thisobs%use_global_obs** - Whether to use global observations or restrict the observations
   !!                               to the relevant ones (default: *.true.* i.e use global full observations)
   !!
   !! Further variables are set when the routine PDAFomi_gather_obs is called.
   !!
   subroutine init_dim_obs_ssh_mgrid(step, dim_obs)
      use PDAF, only: PDAFomi_gather_obs
      use statevector_pdaf, only: sfields
      use nemo_pdaf, only: nwet, wet_pts, glamt, gphit,jpiglo,jpjglo
      use io_pdaf, only: check
      use, intrinsic :: ieee_arithmetic ! to detect Nans

      integer, intent(in)    :: step    !< Current time step
      integer, intent(inout) :: dim_obs !< Dimension of full observation vector

      integer :: i                             !> Counters
      integer :: dim_obs_p                     !> Number of process-local observations
      integer :: ncid, id_var
      integer :: pos(3), cnt(3)         
      real(pwp), allocatable :: obs_p(:)       !> PE-local observation vector

      real(pwp), allocatable :: ivar_obs_p(:)  !> PE-local inverse observation error variance
      real(pwp), allocatable :: ocoord_p(:, :) !> PE-local observation coordinates
      real(pwp) :: rad_conv = 3.141592653589793/180.0 !> Degree to radian conversion
      real(pwp), allocatable :: obs_global(:,:) !< Global obs field read from file
      integer :: ndims, natts, xtype
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: ik, dimlen
      character(len=NF90_MAX_NAME) :: varname, dimname

      ! *****************************
      ! *** Global setting config ***
      ! *****************************
      if (mype_filter == 0) &
         write (*, '(a,4x,a)') 'NEMO-PDAF', 'Assimilate observations - obs_ssh_mgrid'
      ! Store whether to assimilate this observation type (used in routines
      ! below)
      if (assim_ssh_mgrid) thisobs%doassim = 1
      ! Specify type of distance computation
      thisobs%disttype = 3   ! 3=Haversine
      ! Number of coordinates used for distance computation.
      ! The distance compution starts from the first row
      thisobs%ncoord = 2
      ! Set to use limited full observations
      thisobs%use_global_obs = use_global_obs
      ! Omit large innovation to omit nwet = 0 points
      thisobs%inno_omit = 1e6
      
      ! **********************************
      ! *** Read observations from file ***
      ! **********************************
      if (mype_filter == 0) then
         write (*, '(a,4x,a,a)')  'NEMO-PDAF', '--- obs file: ', trim(file_ssh_mgrid)
         write (*, '(a,4x,a,a)')  'NEMO-PDAF', '--- variable: ', trim(varname_ssh_mgrid)
         write (*, '(a,4x,a,i0)') 'NEMO-PDAF', '--- reading nc_step: ', nc_step
      end if

      allocate(obs_global(jpjglo, jpiglo)) !148 x 180
      obs_global = -1
      call check( nf90_open(trim(file_ssh_mgrid), NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, trim(varname_ssh_mgrid), id_var) )
      call check( nf90_inq_varid(ncid, trim(varname_ssh_mgrid), id_var) )
      
      call check( nf90_inquire_variable(ncid, id_var, name=varname, xtype=xtype, &
                                 ndims=ndims, dimids=dimids, natts=natts) )

      print *, 'Variable:', trim(varname)
      print *, 'Number of dimensions:', ndims

      ! Loop over dimensions in the order Fortran sees them to check same as ens
      do ik = 1, ndims
         call check( nf90_inquire_dimension(ncid, dimids(ik), name=dimname, len=dimlen) )
         print *, 'Dim', ik, ':', trim(dimname), ' size =', dimlen
      end do

      pos = (/1, 1, nc_step/)
      cnt = (/jpiglo, jpjglo, 1/)
      call check( nf90_get_var(ncid, id_var, obs_global, start=pos, count=cnt) )
      call check( nf90_close(ncid) )
      where (ieee_is_nan(obs_global))
          obs_global = 0.0             !set nan values to zero
      end where

      obs_global = obs_global / 100.0          !divide by 10 to convert from % to fraction
      write(*,*) minval(obs_global), maxval(obs_global)
      


      ! ***********************************************************
      ! *** Count available observations for the process domain ***
      ! *** and initialize index and coordinate arrays.         ***
      ! ***********************************************************
      ! Set number of local observations
      dim_obs_p = nwet
      if (dim_obs_p > 0) then
         ! Vector of observations on the process sub-domain
         allocate (obs_p(dim_obs_p))
         ! Coordinate array of observations on the process sub-domain
         allocate (ocoord_p(2, dim_obs_p))
         ! Coordinate array for observation operator
         allocate (thisobs%id_obs_p(1, dim_obs_p))
         allocate (ivar_obs_p(dim_obs_p))
         obs_p = -1
         do i = 1, nwet
            ! State vector index counter for observation operator.
            !obs_p(i) = 0.0_pwp
            !*** READING FROM GLOBAL_OBS - BUT IS THIS CORRECT??
            obs_p(i) = obs_global(wet_pts(2, i), wet_pts(1, i))
            !if (i == 1) then

            !    write(*,*) obs_global(wet_pts(6, i),wet_pts(7, i))
            !    write(*,*) obs_global(wet_pts(2, i),wet_pts(1, i))
            !endif
            ! if (i == nwet) then
            !      write(*,*) obs_global(wet_pts(6, i),wet_pts(7, i))
            !      write(*,*) obs_global(wet_pts(2, i),wet_pts(1, i))

             !endif
           !write(*,*)
            ! Observation coordinates - must be in radians for PDAFOMI
            ocoord_p(1, i) = glamt(wet_pts(6, i), wet_pts(7, i))*rad_conv
            ocoord_p(2, i) = gphit(wet_pts(6, i), wet_pts(7, i))*rad_conv
            ! Coordinates for observation operator (gridpoint)
            thisobs%id_obs_p(1, i) = i + sfields(id_sfields)%off
         end do
         write (*,*) minval(obs_p),maxval(obs_p)


         ! *********************************************************
         ! *** For twin experiment: Read synthetic observations  ***
         ! *********************************************************
         !call add_noise(dim_obs_p, obs_p) DON'T WANT TO DO THIS - need to read in WHOLE OBS
      else
         ! Set empty arrays if no local observations
         allocate (obs_p(1))
         obs_p(1) = thisobs%inno_omit
         allocate (ocoord_p(2, 1))
         allocate (thisobs%id_obs_p(1, 1))
         allocate (ivar_obs_p(1))
      end if
      deallocate(obs_global)

      ! ****************************************************************
      ! *** Define observation errors for process-local observations ***
      ! ****************************************************************
      ! Set inverse observation error variances
      ivar_obs_p(:) = 1.0/(rms*rms)

      ! ****************************************
      ! *** Gather global observation arrays ***
      ! ****************************************
      call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs%ncoord, lradius, dim_obs)

      ! ********************
      ! *** Finishing up ***
      ! ********************
      ! Deallocate all local arrays
      deallocate (obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs_ssh_mgrid

   !>###Implementation of observation operator
   !>
   !>This routine applies the full observation operator
   !>for the ssh observations.
   !>
   !>The routine is called by all filter processes.
   !>
   subroutine obs_op_ssh_mgrid(dim_p, dim_obs, state_p, ostate)

      use PDAF, only: PDAFomi_obs_op_gridpoint
      !> PE-local state dimension
      integer, intent(in) :: dim_p
      !> Dimension of full observed state (all observed fields)
      integer, intent(in) :: dim_obs
      !> PE-local model state
      real(pwp), intent(in) :: state_p(dim_p)
      !> Full observed state
      real(pwp), intent(inout) :: ostate(dim_obs)

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      if (thisobs%doassim == 1) then
         call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
      end if

   end subroutine obs_op_ssh_mgrid

   !>###Initialize local information on the module-type observation
   !>
   !>The routine is called during the loop over all local
   !>analysis domains. It has to initialize the information
   !>about local ssh observations.
   !>
   !>This routine calls the routine `PDAFomi_init_dim_obs_l`
   !>for each observation type. The call allows to specify a
   !>different localization radius and localization functions
   !>for each observation type and local analysis domain.
   !>
   subroutine init_dim_obs_l_ssh_mgrid(domain_p, step, dim_obs, dim_obs_l)

      use PDAF, only: PDAFomi_init_dim_obs_l
      use assimilation_pdaf, only: domain_coords
      use config_pdaf, only: locweight

      !> Index of current local analysis domain
      integer, intent(in)  :: domain_p
      !> Current time step
      integer, intent(in)  :: step
      !> Full dimension of observation vector
      integer, intent(in)  :: dim_obs
      !> Local dimension of observation vector
      integer, intent(out) :: dim_obs_l

      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, domain_coords, &
           locweight, lradius, sradius, dim_obs_l)

   end subroutine init_dim_obs_l_ssh_mgrid

   !>###Routine to add model error.
   !>
   subroutine add_noise(dim_obs_p, obs)

      !> Number of process-local observations
      integer, intent(in) :: dim_obs_p
      !> Process-local observations
      real, intent(inout) :: obs(dim_obs_p)

      !> Random noise
      real, allocatable :: noise(:)
      !> Seed for random number generator
      integer, save :: iseed(4)
      !> Flag for first call
      logical, save :: firststep = .true.

      ! Seeds taken from PDAF Lorenz96 routine
      if (firststep) then
         write (*, '(9x, a)') '--- Initialize seed for ssh_mgrid noise'
         iseed(1) = 2*220 + 1
         iseed(2) = 2*100 + 5
         iseed(3) = 2*10 + 7
         iseed(4) = 2*30 + 9
         firststep = .false.
      end if

      ! Generate random Gaussian noise
      allocate (noise(dim_obs_p))
      call dlarnv(3, iseed, dim_obs_p, noise)

      obs = obs + (noise_amp*noise)

      deallocate (noise)

   end subroutine add_noise

 end module obs_ssh_mgrid_pdafomi
