module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o, rair
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use wv_saturation,    only: qsat_water
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props
use cam_history,      only: addfld, add_default, outfld, horiz_only
use cam_logfile,      only: iulog
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,   only: endrun
  use spmd_utils,           only: masterproc
  use cam_logfile,          only: iulog
implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr,   &
   modal_aero_wateruptake_sub,  &
   modal_aero_kohler

public :: modal_aero_wateruptake_reg

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8


! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: sulfeq_idx     = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0
integer :: hygro_idx      = 0
integer :: dryvol_idx     = 0
integer :: dryrad_idx     = 0
integer :: drymass_idx    = 0
integer :: so4dryvol_idx  = 0
integer :: naer_idx       = 0


logical, public :: modal_strat_sulfate = .false.   ! If .true. then MAM sulfate surface area density used in stratospheric heterogeneous chemistry

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_reg()

  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

   integer :: nmodes


end subroutine modal_aero_wateruptake_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_wateruptake_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field
   use infnan,       only : nan, assignment(=)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   real(r8) :: real_nan

   integer :: m, nmodes
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------


end subroutine modal_aero_wateruptake_init

!===============================================================================


subroutine modal_aero_wateruptake_dr(state, pbuf, list_idx_in, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m, hygro_m, dryvol_m, dryrad_m, drymass_m,&
                                     so4dryvol_m, naer_m)
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
! *** N.B. *** The calculation has been enabled for diagnostic mode lists
!              via optional arguments.  If the list_idx arg is present then
!              all the optional args must be present.
!
!-----------------------------------------------------------------------

   use time_manager,  only: is_first_step
   use cam_history,   only: outfld, fieldname_len
   use tropopause,    only: tropopause_find_cam, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer

   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional,          pointer       :: dgnumdry_m(:,:,:)
   real(r8), optional,          pointer       :: dgnumwet_m(:,:,:)
   real(r8), optional,          pointer       :: qaerwat_m(:,:,:)
   real(r8), optional,          pointer       :: wetdens_m(:,:,:)
   real(r8), optional,          pointer       :: hygro_m(:,:,:)
   real(r8), optional,          pointer       :: dryvol_m(:,:,:)
   real(r8), optional,          pointer       :: dryrad_m(:,:,:)
   real(r8), optional,          pointer       :: drymass_m(:,:,:)
   real(r8), optional,          pointer       :: so4dryvol_m(:,:,:)
   real(r8), optional,          pointer       :: naer_m(:,:,:)

   ! local variables

   integer  :: lchnk              ! chunk index
   integer  :: ncol               ! number of columns
   integer  :: list_idx           ! radiative constituents list index
   integer  :: stat

   integer :: i, k, l, m
   integer :: itim_old
   integer :: nmodes
   integer :: nspec
   integer :: tropLev(pcols)

   character(len=fieldname_len+3) :: fieldname

   real(r8), pointer :: h2ommr(:,:) ! specific humidity
   real(r8), pointer :: t(:,:)      ! temperatures (K)
   real(r8), pointer :: pmid(:,:)   ! layer pressure (Pa)

   real(r8), pointer :: cldn(:,:)      ! layer cloud fraction (0-1)
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: dgncur_awet(:,:,:)
   real(r8), pointer :: wetdens(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)

   real(r8), pointer :: raer(:,:)        ! aerosol species MRs (kg/kg)
   real(r8), pointer :: maer(:,:,:)      ! accumulated aerosol mode MRs
   real(r8), pointer :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), pointer :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), pointer :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), pointer :: so4dryvol(:,:,:) ! single-particle-mean so4 dry volume (m3)
   real(r8), pointer :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), pointer :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)
   real(r8), allocatable :: specdens_1(:)

   real(r8), pointer :: sulfeq(:,:,:) ! H2SO4 equilibrium mixing ratios over particles (mol/mol)
   real(r8), allocatable :: wtpct(:,:,:)  ! sulfate aerosol composition, weight % H2SO4
   real(r8), allocatable :: sulden(:,:,:) ! sulfate aerosol mass density (g/cm3)

   real(r8) :: specdens, so4specdens
   real(r8) :: sigmag
   real(r8), allocatable :: alnsg(:)
   real(r8) :: rh(pcols,pver)        ! relative humidity (0-1)
   real(r8) :: dmean, qh2so4_equilib, wtpct_mode, sulden_mode

   real(r8) :: es(pcols)             ! saturation vapor pressure
   real(r8) :: qs(pcols)             ! saturation specific humidity

   real(r8) :: pm25(pcols,pver)      ! PM2.5 diagnostics
   real(r8) :: rhoair(pcols,pver)

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   character(len=32) :: spectype
   !-----------------------------------------------------------------------

end subroutine modal_aero_wateruptake_dr

!===============================================================================

subroutine modal_aero_wateruptake_sub( &
   ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
   hygro, rh, dryvol, so4dryvol, so4specdens, troplev, &
   wetrad, wetvol, wtrvol, sulden, wtpct)

!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------


   ! Arguments
   integer, intent(in)  :: ncol                    ! number of columns
   integer, intent(in)  :: nmodes
   integer, intent(in)  :: troplev(:)

   real(r8), intent(in) :: rhcrystal(:)
   real(r8), intent(in) :: rhdeliques(:)
   real(r8), intent(in) :: dryrad(:,:,:)         ! dry volume mean radius of aerosol (m)
   real(r8), intent(in) :: hygro(:,:,:)          ! volume-weighted mean hygroscopicity (--)
   real(r8), intent(in) :: rh(:,:)               ! relative humidity (0-1)
   real(r8), intent(in) :: dryvol(:,:,:)         ! dry volume of single aerosol (m3)
   real(r8), intent(in) :: so4dryvol(:,:,:)      ! dry volume of sulfate in single aerosol (m3)
   real(r8), intent(in) :: so4specdens           ! mass density sulfate in single aerosol (kg/m3)
   real(r8), intent(in) :: wtpct(:,:,:)          ! sulfate aerosol composition, weight % H2SO4
   real(r8), intent(in) :: sulden(:,:,:)         ! sulfate aerosol mass density (g/cm3)

   real(r8), intent(out) :: wetrad(:,:,:)        ! wet radius of aerosol (m)
   real(r8), intent(out) :: wetvol(:,:,:)        ! single-particle-mean wet volume (m3)
   real(r8), intent(out) :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol (m3)

   ! local variables

   integer :: i, k, m

   real(r8) :: hystfac                ! working variable for hysteresis
   !-----------------------------------------------------------------------



end subroutine modal_aero_wateruptake_sub

!-----------------------------------------------------------------------
      subroutine modal_aero_kohler(   &
          rdry_in, hygro, s, rwet_out, im )

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols

      implicit none

! arguments
      integer :: im         ! number of grid points to be processed
      real(r8) :: rdry_in(:)    ! aerosol dry radius (m)
      real(r8) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
      real(r8) :: s(:)          ! relative humidity (1 = saturated)
      real(r8) :: rwet_out(:)   ! aerosol wet radius (m)

! local variables
      integer, parameter :: imax=200
      integer :: i, n, nsol

      real(r8) :: a, b
      real(r8) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
      real(r8) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
      real(r8) :: p
      real(r8) :: r3, r4
      real(r8) :: r(im)         ! wet radius (microns)
      real(r8) :: rdry(imax)    ! radius of dry particle (microns)
      real(r8) :: ss            ! relative humidity (1 = saturated)
      real(r8) :: slog(imax)    ! log relative humidity
      real(r8) :: vol(imax)     ! total volume of particle (microns**3)
      real(r8) :: xi, xr

      complex(r8) :: cx4(4,imax),cx3(3,imax)

      real(r8), parameter :: eps = 1.e-4_r8
      real(r8), parameter :: mw = 18._r8
      real(r8), parameter :: pi = 3.14159_r8
      real(r8), parameter :: rhow = 1._r8
      real(r8), parameter :: surften = 76._r8
      real(r8), parameter :: tair = 273._r8
      real(r8), parameter :: third = 1._r8/3._r8
      real(r8), parameter :: ugascon = 8.3e7_r8



      end subroutine modal_aero_kohler


!-----------------------------------------------------------------------
      subroutine makoh_cubic( cx, p2, p1, p0, im )
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx)
      complex(r8) :: cx(3,imx)

      integer :: i
      real(r8) :: eps, q(imx), r(imx), sqrt3, third
      complex(r8) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      save eps


      end subroutine makoh_cubic


!-----------------------------------------------------------------------
      subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex(r8) :: cx(4,imx)

      integer :: i
      real(r8) :: third, q(imx), r(imx)
      complex(r8) :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero



      end subroutine makoh_quartic

!----------------------------------------------------------------------
      subroutine calc_h2so4_equilib_mixrat( temp, pres, qh2o, dmean, &
                                            qh2so4_equilib, wtpct, sulden )

      implicit none

      real(r8), intent(in)  :: temp            ! temperature (K)
      real(r8), intent(in)  :: pres            ! pressure (Pa)
      real(r8), intent(in)  :: qh2o            ! water vapor specific humidity (kg/kg)
      real(r8), intent(in)  :: dmean           ! mean diameter of particles in each mode
      real(r8), intent(out) :: qh2so4_equilib  ! h2so4 saturation mixing ratios over the particles (mol/mol)
      real(r8), intent(out) :: wtpct           ! sulfate composition, weight % H2SO4
      real(r8), intent(out) :: sulden          ! sulfate aerosol mass density (g/cm3)

      ! Local declarations
      real(r8)            :: qh2o_kelvin ! water vapor specific humidity adjusted for Kelvin effect (kg/kg)
      real(r8)            :: wtpct_flat  ! sulfate composition over a flat surface, weight % H2SO4
      real(r8)            :: fk1, fk4, fk4_1, fk4_2
      real(r8)            :: factor_kulm                     ! Kulmala correction terms
      real(r8)            :: en, t, sig1, sig2, frac, surf_tens, surf_tens_mode, dsigma_dwt
      real(r8)            :: den1, den2, sulfate_density, drho_dwt
      real(r8)            :: akelvin, expon, akas, rkelvinH2O, rkelvinH2O_a, rkelvinH2O_b
      real(r8)            :: sulfequil, r
      real(r8), parameter :: t0_kulm     = 340._r8           !  T0 set in the low end of the Ayers measurement range (338-445K)
      real(r8), parameter :: t_crit_kulm = 905._r8           !  Critical temperature = 1.5 * Boiling point
      real(r8), parameter :: fk0         = -10156._r8 / t0_kulm + 16.259_r8   !  Log(Kulmala correction factor)
      real(r8), parameter :: fk2         = 1._r8 / t0_kulm
      real(r8), parameter :: fk3         = 0.38_r8 / (t_crit_kulm - t0_kulm)
      real(r8), parameter :: RGAS = 8.31430e+07_r8           ! ideal gas constant (erg/mol/K)
      real(r8), parameter :: wtmol_h2so4 =  98.078479_r8     ! molecular weight of sulphuric acid
      real(r8), parameter :: wtmol_h2o   =  18.015280_r8     ! molecular weight of water vapor
      real(r8), parameter :: wtmol_air   =  28.9644_r8       ! molecular weight of dry air
      real(r8)            :: gv_wt_abcd_en(6,95), gvh2ovp(95)
      real(r8)            :: dnwtp(46), dnc0(46), dnc1(46)
      real(r8)            :: stwtp(15), stc0(15), stc1(15)
      integer             :: i, k

      end subroutine calc_h2so4_equilib_mixrat


!----------------------------------------------------------------------
      subroutine calc_h2so4_wtpct( temp, pres, qh2o, wtpct )

  !!  This function calculates the weight % H2SO4 composition of
  !!  sulfate aerosol, using Tabazadeh et. al. (GRL, 1931, 1997).
  !!  Rated for T=185-260K, activity=0.01-1.0
  !!
  !!  Argument list input:
  !!    temp = temperature (K)
  !!    pres = atmospheric pressure (Pa)
  !!    qh2o = water specific humidity (kg/kg)
  !!
  !!  Output:
  !!    wtpct = weight % H2SO4 in H2O/H2SO4 particle (0-100)
  !!
  !! @author Mike Mills
  !! @ version October 2013

      use wv_saturation, only: qsat_water

      implicit none

      real(r8), intent(in)  :: temp  ! temperature (K)
      real(r8), intent(in)  :: pres  ! pressure (Pa)
      real(r8), intent(in)  :: qh2o  ! water vapor specific humidity (kg/kg)
      real(r8), intent(out) :: wtpct ! sulfate weight % H2SO4 composition

      ! Local declarations
      real(r8)            :: atab1,btab1,ctab1,dtab1,atab2,btab2,ctab2,dtab2
      real(r8)            :: contl, conth, contt, conwtp
      real(r8)            :: activ
      real(r8)            :: es ! saturation vapor pressure over water (Pa) (dummy)
      real(r8)            :: qs ! saturation specific humidity over water (kg/kg)



      end subroutine calc_h2so4_wtpct


!----------------------------------------------------------------------

   end module modal_aero_wateruptake


