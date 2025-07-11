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

  public :: dust_names
  public :: dust_nbin
  public :: dust_nrange
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

  public :: dust_depvel

  integer :: dust_nbin = 0
  integer :: dust_nrange = 0
  character(len=6), protected, allocatable :: dust_names(:)

  real(r8), allocatable :: dust_dmt_grd(:) ! TODO: ?? diameter?

  integer, protected, allocatable :: dust_indices(:)
!  real(r8) :: dust_dmt_vwr(dust_nbin) !TODO: wet diameter??
!  real(r8) :: dust_stk_crc(dust_nbin)

  real(r8)          :: dust_emis_fact = 0._r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'none' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.
  class(aerosol_properties), pointer :: aero_props=>null()

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
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use dust_common,   only: dust_set_params
    use chem_mods,     only: gas_pcnst
    use constituents,  only: cnst_name
    use mo_tracname, only: solsym

    type(sectional_aerosol_properties), intent(in) :: aero_props

    ! local variables
    character(len=6), dimension(100) :: dust_name_list
    integer :: n, dust_nspecies
    integer :: ispec, ibin, ndst
    integer,allocatable :: bins2ranges(:)
    character(len=6) :: name

    character(len=*), parameter :: subname = 'dust_init'

    !aero_props => sectional_aerosol_properties()
    ! ndst = dust_nbin what bins the dust is in should be dust_nrange here
    ! dust_nnum also ndst
    ! dust_indices(2*ndst) -> (dust_nbin + dust_nnum) -> one per tracer here num_1, dust_r3, ... = nbin(where du) + nrange(where du)

    allocate(bins2ranges(aero_props%nbins()))

    bins2ranges = aero_props%bins2ranges(aero_props%nbins())

    ! TODO: generalize and move to obj?
    dust_nspecies = 0
    dust_nbin = 0 ! TODO: make allocatable with dust_nspecies
    dust_nrange = 0 ! TODO: make allocatable with dust_nspecies
    do ispec = 1, aero_props%nspecies_tot()
        call aero_props%get(1, ispec, specname=name)
        if (trim(name) == 'DU') then
            dust_nspecies = dust_nspecies + 1
            dust_nrange = aero_props%spec_range_idx(ispec,'upper') - aero_props%spec_range_idx(ispec,'lower') + 1 !TODO: move to second loop to allow for several dust species
            do ibin = 1, aero_props%nbins()
                if (bins2ranges(ibin) >= aero_props%spec_range_idx(ispec,'lower') .and. &
                    bins2ranges(ibin) <= aero_props%spec_range_idx(ispec,'upper')) then
                    dust_nbin = dust_nbin + 1 ! TODO: move to second loop to allow for several dust species
                end if
            end do
        end if
    end do

    allocate(dust_indices(dust_nbin + dust_nrange)) ! TODO fix with dust_nspecies, sth sum(dust_nbin) + sum(dust_nrange)

    ! TODO: get dust_names from aero_props?
    ! inspired by chemistry.F90 chem_implements_cnst
    ndst = 0

    do n = 1,gas_pcnst
        name = solsym(n)
        if( name(:3) == 'DU_' ) then
            ndst=ndst+1
            dust_name_list(ndst) = name
        endif
    enddo


    ! find ndst from aero_props species props or nr of DU_ tracers
    allocate( dust_names(dust_nrange) )

    dust_names = dust_name_list(:dust_nrange)

    if (masterproc) then
        write(iulog,*) ' dust names: ', dust_names
        write(iulog,*) ' dust_nspecies: ', dust_nspecies
        write(iulog,*) ' dust_nrange: ', dust_nrange
        write(iulog,*) ' dust_nbin: ', dust_nbin
    end if

    call  soil_erod_init( dust_emis_fact, soil_erod_file )

  end subroutine dust_init

  !==============================================================================
  !==============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use mo_constants,  only : dust_density
    use physconst,     only : pi

   ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

   ! local vars
    integer :: i, m, idst
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8

    character(len=*), parameter :: subname = 'dust_emis'

  !  real(r8), parameter :: dust_emis_sclfctr(dust_nbin) &
  !       = (/ 0.011_r8/0.032456_r8, 0.087_r8/0.174216_r8, 0.277_r8/0.4085517_r8, 0.625_r8/0.384811_r8 /)

    ! set dust emissions

    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       ! adjust emissions based on soil erosion
  !     do m = 1,dust_nbin

  !        idst = dust_indices(m)
  !        cflx(i,idst) = -dust_flux_in(i,m) &
  !             * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

   !    enddo

    end do col_loop

  end subroutine dust_emis

  !===============================================================================
  !===============================================================================
  subroutine dust_depvel( temp, pmid, ram1, fv, ncol,  vlc_dry,vlc_trb,vlc_grv )
    ! TODO: Is this routine deprecated?
    use aerosol_depvel, only: aerosol_depvel_compute
    use mo_constants,   only: dust_density
    use ppgrid,         only: pver

    real(r8), intent(in) :: temp(:,:)  ! temperature
    real(r8), intent(in) :: pmid(:,:)  ! mid point pressure
    real(r8), intent(in) :: ram1(:)    ! aerodynamical resistance (s/m)
    real(r8), intent(in) :: fv(:)      ! friction velocity (m/s)
    integer,  intent(in) :: ncol

    real(r8), intent(out) :: vlc_trb(:,:)    !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(:,:,:)  !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(:,:,:)  !dry deposn velocity (m/s)

  !  real(r8) :: diam(ncol,pver,dust_nbin)
    integer :: m
    character(len=*), parameter :: subname = 'dust_depvel'

   ! do m=1,dust_nbin
   !    diam(:,:,m) = dust_dmt_vwr(m)
   ! enddo
   ! call aerosol_depvel_compute( ncol, pver, dust_nbin, temp, pmid, ram1, fv, diam, dust_stk_crc, dust_density, &
   !                              vlc_dry,vlc_trb,vlc_grv)
    return
  endsubroutine dust_depvel

end module dust_model
