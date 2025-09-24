!===============================================================================
! Dust for Bulk Aerosol Model
!===============================================================================
module dust_model
  use shr_kind_mod,    only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun

  use aerosol_properties_mod, only: aerosol_properties
  use sectional_aerosol_properties_mod, only: sectional_aerosol_properties

  implicit none
  private

  ! Public data (TODO: move to object?) also, we have index-mania now!
  public :: dust_names   ! names of dust tracers (dust_nrange)
  public :: dust_nbin    ! nr bins containing dust
  public :: dust_nrange  ! nr of ranges containing dust
  public :: dust_bin_tracer_idx ! indices of num_ tracers containing dust (from const_get_ind)
  public :: dust_range_tracer_idx ! indices of DU_ tracers (from const_get_ind)
  public :: dust_species_idx ! index in the species object array
  public :: dust_active

  ! Public procedures
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init

  ! private routines (previously in soil_erod_mod in CAM)
  private :: soil_erod_init

  ! TODO: move to object?
  integer :: dust_nbin = 0
  integer :: dust_nrange = 0
  integer :: dust_bin_idx(100) ! object internal bin index (array index for bins containing dust)
  integer :: dust_range_idx(100) ! object internal bin index (array index for bins containing dust)
  integer :: dust_species_idx = 0 ! object internal index
  character(len=6), protected, allocatable :: dust_names(:)
  character(len=10), allocatable :: dust_bin_names(:)

  ! TODO: move to obj somehow, currently the format in the base obj. is too rigid?
  integer, protected, allocatable :: dust_bin_tracer_idx(:)
  integer, protected, allocatable :: dust_range_tracer_idx(:)

  ! TODO: get proper distribution & map onto bins, 11 dust bins currently (e.g. Kok et al 2011) see dust_emis
  real(r8), allocatable :: emis_fraction_in_bin(:)

  logical :: dust_active = .false.
  class(aerosol_properties), pointer :: aero_props=>null()

  ! soil parameters from oslo_aero
  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  real(r8), allocatable ::  soil_erodibility(:,:)        ! soil erodibility factor
  real(r8) :: soil_erod_fact                             ! tuning parameter for dust emissions

contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,    only: find_group_name
    use spmd_utils,        only: mpicom, mstrid=>masterprocid, mpi_character, mpi_real8, MPI_SUCCESS

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

    ! Report
    if (masterproc) then
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
    character(len=6)      :: dust_name_list(100)
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
        call aero_props%get(1, ispec, specname=name)
        if (trim(name) == 'DU') then
            dust_species_idx = ispec
            dust_nspecies = dust_nspecies + 1
            dust_nrange = aero_props%spec_nrange(ispec)
            dust_nbin = aero_props%spec_nbin(ispec)
            dust_bin_idx = aero_props%spec_bin_idx(ispec)
            dust_range_idx = aero_props%spec_range_idx(ispec)
            dust_name_list = aero_props%spec_tracernames(ispec) !TODO, Question: make part of "get" routine?
        end if
    end do

    ! find ndst from aero_props species props or nr of DU_ tracers -> move these to sectional_aerosol_properties?
    allocate( dust_names(dust_nrange), stat=istat )
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'dust_names'")
    end if
    allocate( dust_bin_names(dust_nbin), stat=istat )
    if (istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'dust_bin_names'")
    end if
    allocate(dust_bin_tracer_idx(dust_nbin), stat=istat )
    if ( istat /= 0 ) then
        call endrun(subname//":: ERROR could not allocate 'dust_bin_tracer_idx'")
    end if
    allocate(dust_range_tracer_idx(dust_nrange), stat=istat )
    if ( istat /= 0) then
        call endrun(subname//":: ERROR could not allocate 'dust_range_tracer_idx'")
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

    dust_names = dust_name_list(:dust_nrange)
    do ibin = 1, dust_nbin
        dust_bin_names(ibin) = 'num_'//int2str(dust_bin_idx(ibin))
        call cnst_get_ind(dust_bin_names(ibin), dust_bin_tracer_idx(ibin))
    end do
    do irange = 1, dust_nrange
        call cnst_get_ind(dust_names(irange), dust_range_tracer_idx(irange))
    end do

    dust_active = dust_nrange > 0
    if (.not.dust_active) return

    call soil_erod_init( dust_emis_fact, soil_erod_file )

    ! calculate emission fraction per bin

    bin_bounds = aero_props%bin_bounds(aero_props%nbins())
    dust_bin_bounds(:,1) = bin_bounds(dust_bin_idx(:dust_nbin),1)
    dust_bin_bounds(:,2) = bin_bounds(dust_bin_idx(:dust_nbin),2)

    call dust_emis_fraction_bin(dust_nbin, dust_bin_bounds, emis_fraction_in_bin)

    deallocate(bin_bounds)
    deallocate(dust_bin_bounds)
  end subroutine dust_init

  !==============================================================================
  !==============================================================================

  subroutine dust_emis( lchnk, ncol, dstflx, cflx, aero_props )
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
    real(r8) , intent(in)    :: dstflx(pcols,4)
    real(r8) , intent(inout) :: cflx(pcols,pcnst) ! Surface fluxes

    ! Local variables
    integer  :: icol, ibin, irange
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: totalEmissionFlux(pcols)    ! sum emission flux over all sizes
    real(r8) :: cflx_tmp(pcols,dust_nbin)
    integer  :: bins2ranges(aero_props%nbins())
    character(len=*), parameter :: subname = 'dust_emis'

    bins2ranges = aero_props%bins2ranges(aero_props%nbins())
    ! Filter away unreasonable values for soil erodibility
    ! (using low values e.g. gives emissions in greenland..)
    where(soil_erodibility(:,lchnk) < 0.1_r8)
       soil_erod_tmp(:)=0.0_r8
    elsewhere
       soil_erod_tmp(:)=soil_erodibility(:,lchnk)
    end where

    totalEmissionFlux(:) = 0.0_r8
    do icol=1,ncol
       totalEmissionFlux(icol) = totalEmissionFlux(icol) + sum(dstflx(icol,:))
    end do

    ! Note that following CESM use of "dust_emis_fact", the emissions are
    ! scaled by the INVERSE of the factor!!
    ! There is another random scale factor of 1.15 there. Adapting the exact
    ! same formulation as MAM now and tune later
    ! As of NE-380: Oslo dust emissions are 2/3 of CAM emissions
    ! gives better AOD close to dust sources

    ! Sectional model: dust is emitted to the bins, then transferred to ranges
    ! TODO: check compatability with bins! this needs to be number concentration, mass to ranges

    do ibin = 1, dust_nbin
        cflx_tmp(:ncol, ibin) = -1.0_r8*emis_fraction_in_bin(ibin) & ! calculate dust flux
            *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8
        cflx(:ncol, dust_bin_tracer_idx(ibin)) = cflx_tmp(:ncol, ibin) / aero_props%density(dust_species_idx) / aero_props%particle_volume(ibin) ! emission in nr/m2/s
        do irange = 1, dust_nrange
            ! emissions in kg/m2/s to ranges
            if (bins2ranges(dust_bin_idx(ibin)) == dust_range_idx(irange)) then
                cflx(:ncol, dust_range_tracer_idx(irange)) = cflx(:ncol, dust_range_tracer_idx(irange)) + cflx_tmp(:ncol, ibin)
            end if
        end do
    end do

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

    soil_erod_fact = dust_emis_fact

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: soil_erod_fact = ', soil_erod_fact
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
! Helper functions for dust emis (C. Brodowsky)
!=============================================================================

  subroutine dust_emis_fraction_bin(nbin, bin_bounds, vol_frac)
    use string_utils,      only: int2str

    implicit none

    ! input
    integer, intent(in)  :: nbin
    real(r8), intent(in) :: bin_bounds(nbin,2)

    ! local variables
    real(r8)             :: vol(nbin)
    real(r8)             :: D1, D2
    real(r8)             :: h, vol_old
    real(r8)             :: err=0.0001_r8, rel_err ! TODO: better error measure?
    integer              :: ibin, i, j, n

    ! output
    real(r8), intent(out):: vol_frac(nbin)

    do ibin = 1, nbin
        ! integrate function in log space for each bin diameter
        n = 10
        vol_old = 0.0_r8

        D1 = log( bin_bounds(ibin,1) /500._r8 ) ! transform to diameter and um
        D2 = log( bin_bounds(ibin,2) /500._r8 )

        do i = 1, 1000
            h = (D2 - D1) / n ! with of each subinterval
            vol(ibin) = 0.5d0 * (dust_dist(D1) + dust_dist(D2)) ! at bounds -> only half of the trapezoid counts

            do j = 1, n-1 ! calc volume for each sub-increment
                vol(ibin) = vol(ibin) + dust_dist(D1 + i*h)
            end do
            vol(ibin) = h*vol(ibin)

            ! check for convergence
            rel_err = abs(vol(ibin) - vol_old) / (abs(vol(ibin)) )
            if ( rel_err < err ) then
                exit
            end if

            n = n*2
            vol_old = vol(ibin)
        end do

        if (j == 1000) write(iulog,*) 'bin', int2str(ibin), ' did not converge'

    end do

    vol_frac = vol / sum(vol)

  end subroutine dust_emis_fraction_bin

!=============================================================================
!=============================================================================

  real(r8) function dust_dist(logD_d)
    ! Size distribution of emitted dust
    ! local parameters and volume size distribution from Kok et al. (2011) eq. 6
        implicit none
        real(r8), intent(in) :: logD_d           ! (um) size (within one bin)
        real(r8), parameter  :: c_V = 12.62_r8   ! (um) normalization constant volume
        real(r8), parameter  :: D_s = 3.5_r8     ! (um) median diameter by volume
        real(r8), parameter  :: lambda = 12.0_r8 ! (um) crack propagation length
        real(r8), parameter  :: sigma = 3.0_r8   ! geometric stndard deviation

        dust_dist = (exp(logD_d) / c_V) * (1._r8 + erf(log(exp(logD_d)/D_s) / (sqrt(2._r8)*log(sigma)))) &
                    * exp(-1._r8 * (exp(logD_d) / lambda)**3)

        ! TODO: cut-off to set very small values to 0?
  end function dust_dist

end module dust_model
