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

  ! Public data (TODO: move to object?)
  public :: dust_names   ! names of dust tracers (dust_nrange)
  public :: dust_nbin    ! nr bins containing dust
  public :: dust_nrange  ! nr of ranges containing dust
  public :: dust_indices ! indices of dust tracers

  ! Public procedures
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

  ! private routines (previously in soil_erod_mod in CAM)
  private :: soil_erod_init

  ! TODO: move to object?
  integer :: dust_nbin = 0
  integer :: dust_nrange = 0
  integer :: dust_bin_idx(100)
  character(len=6), protected, allocatable :: dust_names(:)

  real(r8), allocatable :: dust_dmt_grd(:) ! TODO: ?? diameter?

  integer, protected, allocatable :: dust_indices(:)
!  real(r8) :: dust_dmt_vwr(dust_nbin) !TODO: wet diameter??
!  real(r8) :: dust_stk_crc(dust_nbin)

  ! TODO: get proper distribution & map onto bins, 11 dust bins currently (e.g. Kok et al 2011) see dust_emis
  real(r8), parameter :: emis_fraction_in_bin(11) = (/0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.09_r8,0.1_r8/)

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
    integer :: unitn, ierr
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
    call mpibcast(dust_emis_fact, 1,                   mpi_real8,     mstrid, mpicom, ierr)
    if (ierr/=MPI_SUCCESS) then
        call endrun(subname//' ERROR: Broadcasting dust_emis_fact')
    end if
    call mpibcast(soil_erod_file, len(soil_erod_file), mpi_character, mstrid, mpicom, ierr)
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
    use chem_mods,     only: gas_pcnst
    use mo_tracname,   only: solsym
   ! use soil_erod_mod, only: soil_erod_init ! TODO, QUESTION: Oslo_aero has its own nearly identical version, why?

    type(sectional_aerosol_properties), intent(in) :: aero_props

    ! local variables
    character(len=6) :: dust_name_list(100)
    integer :: n, dust_nspecies
    integer :: ispec, ibin, ndst, nbin, ind
    character(len=6) :: name

    character(len=*), parameter :: subname = 'dust_init'

    !aero_props => sectional_aerosol_properties()

    dust_nspecies = 0
    dust_nrange = 0
    dust_nbin = 0

    do ispec = 1, aero_props%nspecies_tot()
        call aero_props%get(1, ispec, specname=name)
        if (trim(name) == 'DU') then
            dust_nspecies = dust_nspecies + 1
            dust_nrange = aero_props%spec_nrange(ispec)
            dust_nbin = aero_props%spec_nbin(ispec)
            dust_bin_idx = aero_props%spec_bin_idx(ispec)
            dust_name_list = aero_props%spec_tracernames(ispec) !TODO, Question: make part of "get" routine?
        end if
    end do

    ! find ndst from aero_props species props or nr of DU_ tracers
    allocate( dust_names(dust_nrange) )

    dust_names = dust_name_list(:dust_nrange)

    allocate(dust_indices(dust_nrange))

!    do ind = 1, dust_nrange
!        call cnst_get_ind(trim(dust_names(ind)), dust_indices(ind), abort=.false.)
!    end do
!    dust_active = any(dust_indices(:) > 0)
 !   deallocate(dust_indices)
    dust_active = dust_nrange > 0

    if (.not.dust_active) return

call soil_erod_init( dust_emis_fact, soil_erod_file )

end subroutine dust_init

  !==============================================================================
  !==============================================================================

  subroutine dust_emis( lchnk, ncol, dstflx, cflx )
    !-----------------------------------------------------------------------
    ! Purpose: Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model and
    ! the soil erodibility factor is applied here.
    ! Copied and adapted from oslo_aero
    !-----------------------------------------------------------------------
    use ppgrid,          only: pcols
    use constituents,    only: pcnst

    ! Arguments:
    integer  , intent(in)    :: lchnk
    integer  , intent(in)    :: ncol
    real(r8) , intent(in)    :: dstflx(pcols,4)
    real(r8) , intent(inout) :: cflx(pcols,pcnst) ! Surface fluxes

    ! Local variables
    integer  :: icol, ibin, irange
    integer  :: dust_ind
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: totalEmissionFlux(pcols)    ! sum emission flux over all sizes
    character(len=*), parameter :: subname = 'dust_emis'

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
    ! TODO: units??
    do ibin = 1, dust_nbin ! TODO: change to dust_nbin and num_ -> add to dust_nrange after
        dust_ind = dust_bin_idx(ibin)
        cflx(:ncol, dust_ind) = -1.0_r8*emis_fraction_in_bin(ibin) & !TODO: fix emis_fraction_in_bin
            *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8
    end do

    contains

   ! function dust_emis_fraction_bin()
    ! TODO: fix initial dust distribution
    ! QUESTION: how to integrate nicely?
        !type(sectional_aerosol_properties), intent(in) :: aero_props
   !     dust_total_volume
   !     dust_bin_volume
   !     dust_bin_fraction
   ! end function dust_emis_fraction_bin

  !  function dust_emis_distribution(D_d) result(dV_dlogD)
    ! see Kok et al 2011, equation 6

   !     implicit none
   !     real(r8), intent(in) :: D_d =
   !     real(r8), parameter  :: c_N = 12620.0d0      ! normalization constant volume (nm)
   !     real(r8), parameter  :: D_s = 3400.0d0       ! (nm)
   !     real(r8), parameter  :: sigma = 3.0d0
   !     real(r8), parameter  :: lambda = 12000.0d0   ! (nm)
   !     real(8)              :: dV_dlogD

    !    dV_dlogD = (D_d / c_V) * (1.0d0 + erf(log(D_d/D_s) / (sqrt(2.0d0)*log(sigma)))) * exp(-1.0d0 * (D_d / lambda)**3)

    !end function dust_distribution
  end subroutine dust_emis

  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file )

! TODO, QUESTION: This is the oslo_aero version but it hardly differs from soil_erod_mod, which one to use?

  use ppgrid,           only: pcols, begchunk, endchunk
  use cam_pio_utils,    only: cam_pio_openfile
  use ioFileMod,        only: getfil
  use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE
  use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type

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
    integer               :: c, ncols, ierr
    real(r8), parameter   :: zero=0._r8
    real(r8), parameter   :: twopi=2._r8*pi

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

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    allocate(soil_erodibility_in(nlon,nlat))

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

end module dust_model
