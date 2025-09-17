module modal_aero_calcsize
! TODO: call to "modal_aero" in main cam code
!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use physconst,        only: pi, rhoh2o, gravit

use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state, physics_ptend
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use phys_control,     only: phys_getopts
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num
  use spmd_utils,       only : masterproc
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use cam_history,      only: addfld, add_default, fieldname_len, horiz_only, outfld
use constituents,     only: pcnst, cnst_name
use modal_aero_wateruptake, only: modal_strat_sulfate

use ref_pres,         only: top_lev => clim_modal_aero_top_lev


implicit none
private
save

public modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
public :: modal_aero_calcsize_reg

logical :: do_adjust_default
logical :: do_aitacc_transfer_default

integer :: dgnum_idx     = -1
integer :: hygro_idx     = -1
integer :: dryvol_idx    = -1
integer :: dryrad_idx    = -1
integer :: drymass_idx   = -1
integer :: so4dryvol_idx = -1
integer :: naer_idx      = -1
integer :: sulfeq_idx    = -1


!===============================================================================
contains
!===============================================================================

subroutine modal_aero_calcsize_reg()
  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

  integer :: nmodes


end subroutine modal_aero_calcsize_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   !    set do_adjust_default and do_aitacc_transfer_default flags
   !    create history fields for column tendencies associated with
   !       modal_aero_calcsize
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local
   integer  :: ipair, iq
   integer  :: jac
   integer  :: lsfrm, lstoo
   integer  :: n, nacc, nait
   logical  :: history_aerosol

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
   !-----------------------------------------------------------------------


end subroutine modal_aero_calcsize_init

!===============================================================================

subroutine modal_aero_calcsize_sub(state, ptend, deltat, pbuf, do_adjust_in, &
   do_aitacc_transfer_in)

   !-----------------------------------------------------------------------
   !
   ! Calculates aerosol size distribution parameters
   !    mprognum_amode >  0
   !       calculate Dgnum from mass, number, and fixed sigmag
   !    mprognum_amode <= 0
   !       calculate number from mass, fixed Dgnum, and fixed sigmag
   !
   ! Also (optionally) adjusts prognostic number to
   !    be within bounds determined by mass, Dgnum bounds, and sigma bounds
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), target, intent(in)    :: state       ! Physics state variables

   type(physics_ptend), target, intent(inout)   :: ptend       ! indivdual parameterization tendencies

   real(r8),                    intent(in)    :: deltat      ! model time-step size (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer

   logical, optional :: do_adjust_in
   logical, optional :: do_aitacc_transfer_in

   ! process-specific column tracer tendencies
   ! 3rd index --
   !    1="standard" number adjust gain;
   !    2="standard" number adjust loss;
   !    3=aitken-->accum renaming; 4=accum-->aitken)
   ! 4th index --
   !    1="a" species; 2="c" species
   !-----------------------------------------------------------------------


end subroutine modal_aero_calcsize_sub


!----------------------------------------------------------------------


subroutine modal_aero_calcsize_diag(state, pbuf, list_idx_in, dgnum_m, &
                                    hygro_m, dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)

   !-----------------------------------------------------------------------
   !
   ! Calculate aerosol size distribution parameters
   !
   ! ***N.B.*** DGNUM for the modes in the climate list are put directly into
   !            the physics buffer.  For diagnostic list calculations use the
   !            optional list_idx and dgnum args.
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), intent(in), target :: state   ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)      ! physics buffer

   integer,  optional, intent(in)   :: list_idx_in    ! diagnostic list index
   real(r8), optional, pointer      :: dgnum_m(:,:,:) ! interstital aerosol dry number mode radius (m)
   real(r8), optional, pointer      :: hygro_m(:,:,:)
   real(r8), optional, pointer      :: dryvol_m(:,:,:)
   real(r8), optional, pointer      :: dryrad_m(:,:,:)
   real(r8), optional, pointer      :: drymass_m(:,:,:)
   real(r8), optional, pointer      :: so4dryvol_m(:,:,:)
   real(r8), optional, pointer      :: naer_m(:,:,:)


   ! local
   integer  :: i, k, l1, n
   integer  :: lchnk, ncol
   integer  :: list_idx, stat
   integer  :: nmodes
   integer  :: nspec

   real(r8), pointer :: dgncur_a(:,:) ! (pcols,pver)


   real(r8), parameter :: third = 1.0_r8/3.0_r8

   real(r8), pointer :: mode_num(:,:) ! mode number mixing ratio
   real(r8), pointer :: specmmr(:,:)  ! specie mmr
   real(r8)          :: specdens      ! specie density

   real(r8) :: dryvol_a(pcols,pver)   ! interstital aerosol dry volume (cm^3/mol_air)

   real(r8) :: dgnum, dgnumhi, dgnumlo
   real(r8) :: dgnyy, dgnxx           ! dgnumlo/hi of current mode
   real(r8) :: drv_a                  ! dry volume (cm3/mol_air)
   real(r8) :: dumfac, dummwdens      ! work variables
   real(r8) :: num_a0                 ! initial number (#/mol_air)
   real(r8) :: num_a                  ! final number (#/mol_air)
   real(r8) :: voltonumbhi, voltonumblo
   real(r8) :: v2nyy, v2nxx           ! voltonumblo/hi of current mode
   real(r8) :: sigmag, alnsg
   !-----------------------------------------------------------------------


end subroutine modal_aero_calcsize_diag

subroutine modal_aero_calcdry(state, pbuf, list_idx_in, dgnumdry_m, hygro_m, dryvol_m, dryrad_m, drymass_m, so4dryvol_m, naer_m)

   type(physics_state), target, intent(in)    :: state       ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer
   integer,  optional, intent(in)             :: list_idx_in ! diagnostic list index
   real(r8), optional,          pointer       :: dgnumdry_m(:,:,:)
   real(r8), optional,          pointer       :: hygro_m(:,:,:)
   real(r8), optional,          pointer       :: dryvol_m(:,:,:)
   real(r8), optional,          pointer       :: dryrad_m(:,:,:)
   real(r8), optional,          pointer       :: drymass_m(:,:,:)
   real(r8), optional,          pointer       :: so4dryvol_m(:,:,:)
   real(r8), optional,          pointer       :: naer_m(:,:,:)

   real(r8), parameter :: third = 1._r8/3._r8
   real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8

   real(r8), pointer :: maer(:,:)        ! aerosol wet mass MR (including water) (kg/kg-air)
   real(r8), pointer :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), pointer :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), pointer :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)
   real(r8), pointer :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), pointer :: so4dryvol(:,:,:) ! single-particle-mean so4 dry volume (m3)
   real(r8), pointer :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)

   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: raer(:,:)   ! aerosol species MRs (kg/kg and #/kg)

   real(r8), pointer :: sulfeq(:,:,:) ! H2SO4 equilibrium mixing ratios over particles (mol/mol)

   real(r8) :: dryvolmr(pcols,pver)          ! volume MR for aerosol mode (m3/kg)
   real(r8) :: so4dryvolmr(pcols,pver)       ! volume MR for sulfate aerosol in mode (m3/kg)

   real(r8) :: specdens
   real(r8) :: spechygro, spechygro_1
   real(r8) :: sigmag
   real(r8) :: duma, dumb
   real(r8) :: alnsg

   real(r8) :: v2ncur_a
   real(r8) :: drydens               ! dry particle density  (kg/m^3)

   character(len=fieldname_len+3) :: fieldname
   character(len=32) :: spectype

   integer :: nmodes, lchnk, ncol, list_idx, i, k, l, m
   integer :: nspec



end subroutine modal_aero_calcdry
!----------------------------------------------------------------------

end module modal_aero_calcsize
