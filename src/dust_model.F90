!===============================================================================
! Dust for Bulk Aerosol Model
!===============================================================================
module dust_model
  use shr_kind_mod,    only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun
  use shr_dust_emis_mod,only: is_dust_emis_zender, is_zender_soil_erod_from_atm

  use aerosol_properties_mod, only: aerosol_properties
  use sectional_aerosol_properties_mod, only: sectional_aerosol_properties

  implicit none
  private

  ! Public data (TODO: move to object?) also, we have index-mania now!
  ! TODO: use internal indices for the aerosol scheme, connect only with get_transported/set_transported
  public :: dust_names   ! names of dust tracers (dust_nrange)
  public :: dust_nbin    ! nr bins containing dust
  public :: dust_nrange  ! nr of ranges containing dust
  public :: dust_bin_tracer_ndx ! indices of num_ tracers containing dust (from const_get_ind)
  public :: dust_range_tracer_ndx ! indices of DU_ tracers (from const_get_ind)
  public :: dust_species_ndx ! index in the species object array
  public :: dust_active

  ! Public procedures
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init

  ! private routines (previously in soil_erod_mod in CAM)
  private :: soil_erod_init
  private :: dust_emis_fraction_bin
  private :: dust_dist

  ! TODO: move to object?
  integer :: dust_nbin = 0
  integer :: dust_nrange = 0
  integer, allocatable :: dust_bin_ndx(:) ! object internal bin index (array index for bins containing dust)
  integer, allocatable :: dust_range_ndx(:) ! object internal bin index (array index for bins containing dust)
  integer :: dust_species_ndx = 0 ! object internal index
  character(len=6), protected, allocatable :: dust_names(:)
  character(len=10), allocatable :: dust_bin_names(:)

  ! TODO: move to obj somehow, currently the format in the base obj. is too rigid?
  integer, protected, allocatable :: dust_bin_tracer_ndx(:)
  integer, protected, allocatable :: dust_range_tracer_ndx(:)

  ! TODO: get proper distribution & map onto bins, 11 dust bins currently (e.g. Kok et al 2011) see dust_emis
  real(r8), allocatable :: emis_fraction_in_bin(:)

  logical :: dust_active = .false.
  class(aerosol_properties), pointer :: aero_props=>null()

  ! soil parameters from oslo_aero
  real(r8)          :: dust_emis_fact = 0._r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'none'       ! full pathname for soil erodibility dataset

  real(r8), allocatable ::  soil_erodibility(:,:)        ! soil erodibility factor

contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,    only: find_group_name
    use spmd_utils,        only: mpicom, mstrid=>masterprocid, mpi_character, mpi_real8, MPI_SUCCESS
    use shr_dust_emis_mod, only: shr_dust_emis_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer                     :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call MPI_Bcast(dust_emis_fact, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr/=MPI_SUCCESS) then
        call endrun(subname//' ERROR: Broadcasting dust_emis_fact')
    end if
    call MPI_Bcast(soil_erod_file, len(soil_erod_file), mpi_character, mstrid, mpicom, ierr)
    if (ierr/=MPI_SUCCESS) then
        call endrun(subname//' ERROR: Broadcasting soil_erod_file')
    end if

    call shr_dust_emis_readnl(mpicom, 'drv_flds_in')

    if ((soil_erod_file /= 'none') .and. (.not.is_zender_soil_erod_from_atm())) then
       call endrun(subname//': should not specify soil_erod_file if Zender soil erosion is not in CAM')
    end if

    ! Report
    if (masterproc) then
        if (is_dust_emis_zender()) then
          write(iulog,*) subname,': Zender_2003 dust emission method is being used.'
        end if
        if (is_zender_soil_erod_from_atm()) then
          write(iulog,*) subname,': Zender soil erod file is handled in atm'
          write(iulog,*) subname,': soil_erod_file = ',trim(soil_erod_file)
          write(iulog,*) subname,': dust_emis_fact = ',dust_emis_fact
        end if
        write(iulog, *) 'Dust namelist'
        write(iulog, *) 'dust_emis_fact: ', dust_emis_fact
        write(iulog, *) 'soil_erod_file: ', soil_erod_file
    end if

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================

  subroutine dust_init(aero_props)

    use constituents,  only: cnst_get_ind
    use chem_mods,     only: gas_pcnst !TODO,currently not used..
    use mo_tracname,   only: solsym !TODO, currently not used..
    use string_utils,      only: int2str

   ! use soil_erod_mod, only: soil_erod_init ! TODO, QUESTION: Oslo_aero has its own nearly identical version, why?

    type(sectional_aerosol_properties), intent(in) :: aero_props

    ! local variables
    integer               :: dust_nspecies
    integer               :: ispec, ind, istat, ibin, irange
    character(len=6)      :: name
    real(r8), allocatable :: dust_bin_bounds(:,:)
    real(r8), allocatable :: bin_bounds(:,:)

    character(len=*), parameter :: subname = 'dust_init'

    dust_nspecies = 0
    dust_nrange = 0
    dust_nbin = 0

    do ispec = 1, aero_props%nspecies_tot() !TODO: currently this just works with one dust species
        name = aero_props%specname(ispec)
        if (trim(name) == 'DU') then
            dust_species_ndx = ispec
            dust_nspecies = dust_nspecies + 1
            dust_nrange = aero_props%spec_nrange(ispec)
            dust_nbin = aero_props%spec_nbin(ispec)

            allocate(dust_bin_ndx(dust_nbin), dust_range_ndx(dust_nrange), dust_names(dust_nrange))

            dust_bin_ndx = aero_props%spec_bin_ndx(ispec, dust_nbin)
            dust_range_ndx = aero_props%spec_range_ndx(ispec, dust_nrange)
            dust_names = aero_props%spec_tracernames(ispec, dust_nrange)
            exit ! TODO: Change this to allow for more dust species/compositions
        end if
    end do

    ! find ndst from aero_props species props or nr of DU_ tracers -> move these to sectional_aerosol_properties?
    !allocate( dust_names(dust_nrange), stat=istat )
    !if ( istat /= 0 ) then
    !    call endrun(subname//":: ERROR could not allocate 'dust_names'")
    !end if
    allocate( dust_bin_names(dust_nbin), stat=istat )
    if (istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'dust_bin_names'")
    end if
    allocate(dust_bin_tracer_ndx(dust_nbin), stat=istat )
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'dust_bin_tracer_ndx'")
    end if
    allocate(dust_range_tracer_ndx(dust_nrange), stat=istat )
    if ( istat /= 0) then
        call endrun(subname//":: ERROR could not allocate 'dust_range_tracer_ndx'")
    end if
    allocate(bin_bounds(aero_props%nbins(), 2), stat=istat)
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'bin_bounds'")
    end if
    allocate(dust_bin_bounds(dust_nbin, 2), stat=istat)
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'dust_bin_bounds'")
    end if
    allocate(emis_fraction_in_bin(dust_nbin), stat=istat)
    if (istat/=0) then
        call endrun(subname//":: Error could not allocate 'emis_fraction_in_bin'")
    end if

! TODO: use model-internal indices instead! no "cnst_get_ind" anymore
    do ibin = 1, dust_nbin
        dust_bin_names(ibin) = 'num_'//int2str(dust_bin_ndx(ibin))
        call cnst_get_ind(dust_bin_names(ibin), dust_bin_tracer_ndx(ibin))
    end do
    do irange = 1, dust_nrange
        call cnst_get_ind(dust_names(irange), dust_range_tracer_ndx(irange))
    end do

    dust_active = dust_nrange > 0
    if (.not.dust_active) return

    if (is_zender_soil_erod_from_atm()) then
        call soil_erod_init( dust_emis_fact, soil_erod_file )
    end if

    ! calculate emission fraction per bin

    bin_bounds = aero_props%bin_bounds(aero_props%nbins())
    dust_bin_bounds(:,1) = bin_bounds(dust_bin_ndx(:dust_nbin),1)
    dust_bin_bounds(:,2) = bin_bounds(dust_bin_ndx(:dust_nbin),2)

    call dust_emis_fraction_bin(dust_nbin, dust_bin_bounds, emis_fraction_in_bin)

    deallocate(bin_bounds)
    deallocate(dust_bin_bounds)
  end subroutine dust_init

  !==============================================================================
  !==============================================================================

  subroutine dust_emis( lchnk, ncol, dust_flux_in, cflx, aero_props )
! TODO: move dust_flux_in and cflx out of here somehow?

    !-----------------------------------------------------------------------
    ! Purpose: Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model and
    ! the soil erodibility factor is applied here.
    ! Copied and adapted from oslo_aero
    !-----------------------------------------------------------------------
    use ppgrid,          only: pcols
    use constituents,    only: pcnst

    ! Arguments:
    type(sectional_aerosol_properties), intent(in) :: aero_props
    integer  , intent(in)    :: lchnk
    integer  , intent(in)    :: ncol
    real(r8) , intent(in)    :: dust_flux_in(:,:)   ! Leung emissions (?)
    real(r8) , intent(inout) :: cflx(pcols,pcnst) ! Surface fluxes

    ! Local variables
    integer  :: icol, ibin, irange
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: totalEmissionFlux(pcols)    ! sum emission flux over all sizes
    real(r8) :: cflx_tmp(pcols,dust_nbin)
    character(len=*), parameter :: subname = 'dust_emis'


    ! Note that following CESM use of "dust_emis_fact", the emissions are
    ! scaled by the INVERSE of the factor!!
    ! There is another random scale factor of 1.15 there. Adapting the exact
    ! same formulation as MAM now and tune later
    ! As of NE-380: Oslo dust emissions are 2/3 of CAM emissions
    ! gives better AOD close to dust sources

    totalEmissionFlux(:) = 0.0_r8
    totalEmissionFlux(:ncol) = sum(dust_flux_in(:ncol,:), dim=2)

if (masterproc) then
   write(6,*)" DEBUG: dust_flux_in: ", maxval(abs(sum(dust_flux_in(:ncol,:), dim=2)))
end if

    if (is_zender_soil_erod_from_atm()) then
        ! Filter away unreasonable values for soil erodibility
        ! (using low values e.g. gives emissions in greenland..)
        where(soil_erodibility(:,lchnk) < 0.1_r8)
            soil_erod_tmp(:)=0.0_r8
        elsewhere
            soil_erod_tmp(:)=soil_erodibility(:,lchnk)
        end where

    ! Sectional model: dust is emitted to the bins, then transferred to ranges
    ! TODO: check compatability with bins! this needs to be number concentration, mass to ranges
! TODO: use aerosol model internal indices instead

        do ibin = 1, dust_nbin
            cflx_tmp(:ncol, ibin) = -1.0_r8*emis_fraction_in_bin(ibin) & ! calculate dust flux kg/m2/s
                *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8
            cflx(:ncol, dust_bin_tracer_ndx(ibin)) = cflx_tmp(:ncol, ibin) / aero_props%density(dust_species_ndx) / aero_props%particle_volume(dust_bin_ndx(ibin)) ! emission in nr/m2/s
            do irange = 1, dust_nrange
                ! emissions in kg/m2/s to ranges
                if (aero_props%bins2ranges(dust_bin_ndx(ibin)) == dust_range_ndx(irange)) then
                    cflx(:ncol, dust_range_tracer_ndx(irange)) = cflx(:ncol, dust_range_tracer_ndx(irange)) + cflx_tmp(:ncol, ibin)
                end if
            end do
        end do

    else ! Leung emissions

        do ibin = 1, dust_nbin
            cflx_tmp(:ncol, ibin) = -1.0_r8*emis_fraction_in_bin(ibin) & ! calculate dust flux kg/m2/s
                *totalEmissionFlux(:ncol) / dust_emis_fact
            cflx(:ncol, dust_bin_tracer_ndx(ibin)) = cflx_tmp(:ncol, ibin) / aero_props%density(dust_species_ndx) / aero_props%particle_volume(dust_bin_ndx(ibin)) ! emission in nr/m2/s

            do irange = 1, dust_nrange
                ! emissions in kg/m2/s to ranges
                if (aero_props%bins2ranges(dust_bin_ndx(ibin)) == dust_range_ndx(irange)) then
                    cflx(:ncol, dust_range_tracer_ndx(irange)) = cflx(:ncol, dust_range_tracer_ndx(irange)) + cflx_tmp(:ncol, ibin)
                end if
            end do
        end do
    end if

  end subroutine dust_emis
  !=============================================================================
  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file )

! TODO, QUESTION: This is the oslo_aero version but it hardly differs from soil_erod_mod, which one to use?

  use ppgrid,           only: pcols, begchunk, endchunk
  use cam_pio_utils,    only: cam_pio_openfile
  use ioFileMod,        only: getfil
  use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE
  use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
  use mo_constants,     only: pi, d2r

    ! arguments
    real(r8),         intent(in) :: dust_emis_fact
    character(len=*), intent(in) :: soil_erod_file

    ! localvaraibles
    real(r8), allocatable :: soil_erodibility_in(:,:)
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
    character(len=cl)     :: infile
    integer               :: did, vid, nlat, nlon
    type(file_desc_t)     :: ncid
    type(interp_type)     :: lon_wgts, lat_wgts
    real(r8)              :: to_lats(pcols), to_lons(pcols)
    integer               :: c, ncols, ierr, istat
    real(r8), parameter   :: zero=0._r8
    real(r8), parameter   :: twopi=2._r8*pi
    character(len=*), parameter :: subname = 'soil_erod_init'

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: dust_emis_fact = ', dust_emis_fact
    end if

    ! read in soil erodibility factors, similar to Zender's boundary conditions

    ! Get file name.
    call getfil(soil_erod_file, infile, 0)
    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlon )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlat )

    allocate(dst_lons(nlon), stat=istat )
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could allocate 'dst_lons'")
    end if
    allocate(dst_lats(nlat), stat=istat )
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could allocate 'dst_lats'")
    end if
    allocate(soil_erodibility_in(nlon,nlat), stat=istat )
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could allocate 'soil_erodibility_in'")
    end if

    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, dst_lons  )

    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, dst_lats  )

    ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
    ierr = pio_get_var( ncid, vid, soil_erodibility_in )

    ! convert to radians and setup regridding
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

    allocate( soil_erodibility(pcols,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to allocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to allocate soil_erodibility_in')
    end if

    soil_erodibility(:,:)=0._r8

    ! regrid
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon, nlat, soil_erodibility(:,c), ncols, lon_wgts, lat_wgts)

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( soil_erodibility_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to deallocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to deallocate soil_erodibility_in')
    end if

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine soil_erod_init

!=============================================================================
! Helper functions for dust emis chrisvbr@uio.no
!=============================================================================

  subroutine dust_emis_fraction_bin(nbin, bin_bounds, vol_frac)
    use string_utils,      only: int2str

    ! input
    integer, intent(in)  :: nbin
    real(r8), intent(in) :: bin_bounds(nbin,2)

    ! local variables
    real(r8)             :: vol(nbin)                 ! volume in each subinterval
    real(r8)             :: D1, D2                    ! diameters at bin bounds
    real(r8)             :: subint_width              ! width of each subinterval
    real(r8)             :: vol_old                   ! volume of previous iteration
    real(r8)             :: err=0.0001_r8, rel_err    ! TODO: better error measure?
    real(r8)             :: subint_number             ! number of subintervals in iteration
    integer              :: ibin, iteration, isubint

    ! output
    real(r8), intent(out):: vol_frac(nbin)            ! volume fraction of emitted dust in each bin

    do ibin = 1, nbin
        ! integrate function in log space for each bin diameter
        subint_number = 5._r8      ! start with small subinterval number
        vol_old = 0.0_r8
        vol(ibin) = 0._r8

        D1 = log( bin_bounds(ibin,1) *1.e6_r8 * 2._r8 )                         ! transform to diameter and um
        D2 = log( bin_bounds(ibin,2) *1.e6_r8 * 2._r8 )

        ! initialize "vol_old" to a coarse distribution
        subint_width = (D2 - D1) / subint_number                            ! with of each subinterval
        vol_old = 0.5d0 * (dust_dist(D1) + dust_dist(D2))                 ! at bounds -> only half of the trapezoid counts

        do isubint = 1, int(subint_number)-1                                     ! calc volume for each sub-increment
            vol_old = vol_old + dust_dist(D1 + isubint*subint_width)
        end do
        vol_old = subint_width*vol_old

        ! start with 10 subintervals for the actual calculation
        subint_number = subint_number*2._r8

        do iteration = 1, 1000                                                  ! loop to some large nr to make smaller and smaller increments

            vol(ibin) = 0._r8

            subint_width = (D2 - D1) / subint_number                            ! with of each subinterval
            vol(ibin) = 0.5d0 * (dust_dist(D1) + dust_dist(D2))                 ! at bounds -> only half of the trapezoid counts

            do isubint = 1, int(subint_number)-1                                ! calc volume for each sub-increment
                vol(ibin) = vol(ibin) + dust_dist(D1 + isubint*subint_width)
            end do

            vol(ibin) = subint_width*vol(ibin)                                  ! times the width of the subinterval

            ! check for convergence
            rel_err = abs(vol(ibin) - vol_old) / max( vol(ibin), 1.e-30_r8)     ! avoid divide by 0 and tiny values
            if ( rel_err < err ) then
                exit
            end if

            subint_number = subint_number*2._r8
            vol_old = vol(ibin)

        end do

        if (iteration == 1000) then
            call endrun('ERROR: bin'//int2str(ibin)//' did not converge')
        end if

    end do

    vol_frac = vol / sum(vol)

    end subroutine dust_emis_fraction_bin

!=============================================================================
!=============================================================================

  real(r8) function dust_dist(logD_d)
    ! Size distribution of emitted dust
    ! local parameters and volume size distribution from Kok et al. (2011) eq. 6 (https://doi.org/10.1073/pnas.1014798108)
        real(r8), intent(in) :: logD_d           ! (um) size (within one bin) (natural log)
        real(r8)             :: D_d
        real(r8), parameter  :: c_V = 12.62_r8   ! (um) normalization constant volume
        real(r8), parameter  :: D_s = 3.5_r8     ! (um) median diameter by volume
        real(r8), parameter  :: lambda = 12.0_r8 ! (um) crack propagation length
        real(r8), parameter  :: sigma = 3.0_r8   ! geometric stndard deviation

        D_d = exp(logD_d)
        dust_dist = (D_d / c_V) * (1._r8 + erf(log(D_d/D_s) / (sqrt(2._r8)*log(sigma)))) &
                    * exp(-1._r8 * (D_d / lambda)**3)

  end function dust_dist

end module dust_model
