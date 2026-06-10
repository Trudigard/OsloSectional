!===============================================================================
! Seasalt for Bulk Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,       only: pcols, pver

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

  public :: seasalt_depvel

  logical :: seasalt_active = .false.

  integer, parameter :: seasalt_nbin = 0
  integer, parameter :: seasalt_nnum = 0

  integer :: seasalt_indices(seasalt_nbin)

 contains

   !=============================================================================
   !=============================================================================
   subroutine seasalt_init
     use cam_history,   only: addfld, fieldname_len
     use constituents,  only: cnst_get_ind
     use string_utils,      only: int2str

     character(len=fieldname_len) :: dummy
    character(len=*), parameter :: subname = 'seasalt_init'
    call endrun(subname//" is not yet implemented")

   end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx )


    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    ! local vars
    character(len=*), parameter :: subname = 'seasalt_emis'
    call endrun(subname//" is not yet implemented")

  end subroutine seasalt_emis

  !=============================================================================
  !=============================================================================
 ! subroutine seasalt_emis_frac(emis_frac)

    ! CoefA
    ! CoefB
    ! CoefC
    ! CoefD

    ! From OsloAero:
    ! New whitecap area fraction / air entrainment flux from eqn. 6 in Salter et al. (2015)
    ! JCA & MS Using Hanson & Phillips 99 air entrainment vs. wind speed
    ! (Note the uncertainty in the factor 2, written as 2 pluss/minus 1 in Eq. 6 -> possible tuning factor)
 !   whitecapAreaFraction(:ncol) = (2.0_r8*10.0_r8**(-8.0_r8))*(u10m(:ncol)**3.74_r8)
 !   whitecapAreaFraction(:ncol) = ocnfrc(:ncol) * (1._r8-icefrc(:ncol)) * whitecapAreaFraction(:ncol)
 ! end subroutine seasalt_emis_frac


end module seasalt_model
