!===============================================================================
! Sectional Aerosol Model
!===============================================================================
module aero_model
  use shr_kind_mod,      only: r8 => shr_kind_r8
  use constituents,      only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,            only: pcols, pver, pverp
  use cam_abortutils,    only: endrun
  use cam_logfile,       only: iulog
  use perf_mod,          only: t_startf, t_stopf
  use camsrfexch,        only: cam_in_t, cam_out_t
  use aerodep_flx,       only: aerodep_flx_prescribed
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,    only: physics_buffer_desc
  use physconst,         only: gravit, rair
  use dust_model,        only: dust_active, dust_names, dust_nbin, dust_nrange
  use seasalt_model,     only: sslt_active=>seasalt_active, seasalt_names, seasalt_nbin
  use spmd_utils,        only: masterproc
  use physics_buffer,    only: pbuf_get_field, pbuf_get_index
  use cam_history,       only: outfld
  use infnan,            only: nan, assignment(=)
  use sectional_aerosol_properties_mod, only: sectional_aerosol_properties

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep     ! aerosol dry deposition and sediment
  public :: aero_model_wetdep     ! aerosol wet removal
  public :: aero_model_emissions  ! aerosol emissions
  public :: aero_model_surfarea    ! tropospheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea   ! stub

 ! Misc private data

  integer :: so4_ndx, cb2_ndx, oc2_ndx, nit_ndx
  integer :: soa_ndx, soai_ndx, soam_ndx, soab_ndx, soat_ndx, soax_ndx

  ! aerosol_nl Namelist variables
  character(len=16), allocatable :: wetdep_list(:)
  character(len=16), allocatable :: drydep_list(:)

  integer :: ndrydep = 0
  integer :: nwetdep = 0
  logical :: drydep_lq(pcnst)
  logical :: wetdep_lq(pcnst)

  real(r8) :: aer_sol_facti(pcnst) ! in-cloud solubility factor
  real(r8) :: aer_sol_factb(pcnst) ! below-cloud solubility factor
  real(r8) :: aer_scav_coef(pcnst)

  integer :: fracis_idx = 0

  type(sectional_aerosol_properties), pointer :: aero_props=>null()


contains

  !=============================================================================
  ! reads aerosol namelist options
  !=============================================================================
  subroutine aero_model_readnl(nlfile)
    use dust_model,      only: dust_readnl

    use oslo_aero_control, only: oslo_aero_ctl_readnl ! TODO: use the file in oslo_aero directly? currently this is a copy in the local chemistry.F90

    ! filepath for file containing namelist input
    character(len=*), intent(in) :: nlfile

    ! Local variables
    integer                     :: unitn, ierr, ind, pos
    character(len=50)           :: tmp

    character(len=*), parameter :: subname = 'aero_model_readnl'

    ! modal_aero: read aerosol_nl: aer_drydep_list, modal_strat_sulfate, modal_accum_coarse_exch, seasalt_emis_scale
    ! oslo aero: read aerosol_nl: sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial

    call oslo_aero_ctl_readnl(nlfile)

    call dust_readnl(nlfile)

  end subroutine aero_model_readnl

  !=============================================================================
  !=============================================================================
  subroutine aero_model_register()
    character(len=*), parameter :: subname = 'aero_model_register'

    ! modal_aero: modal_aero_data_reg: allocate all kind of stuff, etc with number modes
    ! oslo aero: aero_register: lots of cnst_get_ind calls for all tracers, set aerosol types

    ! TODO: find out how to get tracer indices -> part of constructor?

    if (masterproc) then
        write(iulog,*) subname//' nothing to do here yet..'
    end if

  end subroutine aero_model_register

  !=============================================================================
  !=============================================================================
  subroutine aero_model_init( pbuf2d, nlfile )

    use mo_chem_utls,   only: get_inv_ndx, get_spc_ndx
    use cam_history,    only: addfld, add_default, horiz_only
    use phys_control,   only: phys_getopts
    use dust_model,     only: dust_init
    use mo_setsox,   only : setsox, has_sox
    use shr_mem_mod,      only: shr_mem_getusage
    use mpi,              only: MPI_REAL8, MPI_MAX
      use spmd_utils,      only: mpicom, masterprocid, masterproc
   use shr_mpi_mod,       only: shr_mpi_barrier
    !use aer_drydep_mod, only: inidrydep
    !use wetdep,         only: wetdep_init

    !use oslo_aero_ocean, only: oslo_aero_ocean_init ! TODO: DMS, add to build-namelist and chemistry.F90 and as well
   logical :: calc_memory_increase = .true.
    real(r8)                            :: mem_hw_beg, mem_hw_end
    real(r8)                            :: mem_beg, mem_end
      real(r8)                            :: temp ! For MPI

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    character(len=12), parameter :: subname = 'aero_model_init'
    integer :: m, id, ierr, ibin, ispec
    character(len=20) :: dummy
    logical  :: history_aerosol ! Output MAM or SECT aerosol tendencies
    logical  :: history_dust    ! Output dust
    logical  :: history_chemistry ! Output Chemistry

    character(len=*), intent(in) :: nlfile
    mem_hw_beg = 0.0
    mem_hw_end = 0.0
    mem_beg = 0.0
    mem_end = 0.0
    !call oslo_aero_ocean_init() ! DMS
!    if (calc_memory_increase) then
       call shr_mem_getusage(mem_hw_beg, mem_beg)
!    end if
    ! TODO: fix this :)
    call phys_getopts( history_aerosol_out   = history_aerosol, &
                       history_dust_out      = history_dust,    &
                       history_chemistry_out = history_chemistry   )

    aero_props => sectional_aerosol_properties(nlfile) ! calls constructor function in sectional_aerosol_properties

    call dust_init(aero_props)
    fracis_idx = pbuf_get_index('FRACIS')

    if (dust_active) then
        do m = 1, dust_nrange !TODO CHECK!!
          dummy = trim(dust_names(m)) // 'SF'
          call addfld (dummy,horiz_only, 'A','kg/m2/s',trim(dust_names(m))//' dust surface emission')
          if (history_aerosol.or.history_chemistry) then
             call add_default (dummy, 1, ' ')
          endif
       enddo

       dummy = 'DSTSFMBL'
       call addfld (dummy,horiz_only, 'A','kg/m2/s','Mobilization flux at surface')
       if (history_aerosol .or. history_dust) then
          call add_default (dummy, 1, ' ')
       endif

       dummy = 'LND_MBL'
       call addfld (dummy,horiz_only, 'A','frac','Soil erodibility factor')
       if (history_aerosol) then
          call add_default (dummy, 1, ' ')
       endif
    endif

 !     if (calc_memory_increase) then
    temp = 0.0
    call shr_mem_getusage(mem_hw_end, mem_end)
         temp = mem_end - mem_beg
         call MPI_barrier(mpicom, ierr)
         call MPI_reduce(temp, mem_end, 1, MPI_REAL8, MPI_MAX, masterprocid,  &
              mpicom, ierr)
         if (masterproc) then
            write(iulog, *) subname, ': Increase in memory usage = ',    &
                 mem_end, ' (MB)'
         end if
         temp = mem_hw_end - mem_hw_beg
         call MPI_barrier(mpicom, ierr)
         call MPI_reduce(temp, mem_hw_end, 1, MPI_REAL8, MPI_MAX,             &
              masterprocid, mpicom, ierr)
         if (masterproc) then
            write(iulog, *) subname, 'Increase in memory highwater = ',       &
                 mem_hw_end, ' (MB)'
         end if
 !     end if
    ! TODO add aq chem (if has_sox ...)

    ! call aero_wetdep_init()
  end subroutine aero_model_init

  !=============================================================================
  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

  !  use dust_sediment_mod, only: dust_sediment_tend
    use aer_drydep_mod,    only: d3ddflux, calcram
    use dust_model,        only: dust_names, dust_nbin
    use seasalt_model,     only: sslt_depvel=>seasalt_depvel, sslt_nbin=>seasalt_nbin, sslt_names=>seasalt_names

    ! args
    type(physics_state),    intent(in)    :: state     ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:)  ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in    ! import state
    real(r8),               intent(in)    :: dt        ! time step
    type(cam_out_t),        intent(inout) :: cam_out   ! export state
    type(physics_ptend),    intent(out)   :: ptend     ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

  ! local vars
    real(r8), pointer :: landfrac(:) ! land fraction
    real(r8), pointer :: icefrac(:)  ! ice fraction
    real(r8), pointer :: ocnfrac(:)  ! ocean fraction
    real(r8), pointer :: fvin(:)     !
    real(r8), pointer :: ram1in(:)   ! for dry dep velocities from land model for progseasalts

    real(r8) :: fv(pcols)            ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)          ! for dry dep velocities, from land modified over ocean & ice

     ! local decarations

    integer, parameter :: begdst = 1 ! TODO now: index in aeronames where dust names start (from bulk aero)

    integer :: ncol, lchnk


    real(r8) :: tsflx_dst(pcols)
    real(r8) :: tsflx_slt(pcols)
    real(r8) :: pvaeros(pcols,pverp)    ! sedimentation velocity in Pa
    real(r8) :: sflx(pcols)

    real(r8) :: tvs(pcols,pver)
    real(r8) :: rho(pcols,pver)      ! air density in kg/m3

    integer :: m,mm, i, im

    character(len=*), parameter :: subname = 'aero_model_drydep'

    if (ndrydep<1) return

    call endrun(subname//":: is not yet implemented")

    landfrac => cam_in%landfrac(:)
    icefrac  => cam_in%icefrac(:)
    ocnfrac  => cam_in%ocnfrac(:)
    fvin     => cam_in%fv(:)
    ram1in   => cam_in%ram1(:)

    lchnk = state%lchnk
    ncol  = state%ncol

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac,icefrac,ocnfrac,obklen,&
                  ustar,ram1in,ram1,state%t(:,pver),state%pmid(:,pver),&
                  state%pdel(:,pver),fvin,fv)

    !call outfld( 'airFV', fv(:), pcols, lchnk )
    !call outfld( 'RAM1', ram1(:), pcols, lchnk )

    ! note that tendencies are not only in sfc layer (because of sedimentation)
    ! and that ptend is updated within each subroutine for different species

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_drydep', lq=drydep_lq)

    lchnk = state%lchnk
    ncol  = state%ncol

    tvs(:ncol,:) = state%t(:ncol,:)
    rho(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

    tsflx_dst(:)=0._r8
    tsflx_slt(:)=0._r8

    ! do drydep for each of the bins of dust and seasalt
    do m=1,ndrydep

       pvaeros(:ncol,1)=0._r8

       call outfld( trim(cnst_name(mm))//'DV', pvaeros(:,2:pverp), pcols, lchnk )

       if(.true.) then ! use phil's method
          !      convert from meters/sec to pascals/sec
          pvaeros(:ncol,2:pverp) = pvaeros(:ncol,2:pverp) * rho(:ncol,:)*gravit

       endif

       if ( any( dust_names(:)==trim(cnst_name(mm)) ) ) &
            tsflx_dst(:ncol)=tsflx_dst(:ncol)+sflx(:ncol)

       ! if the user has specified prescribed aerosol dep fluxes then
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not. aerodep_flx_prescribed()) then
          ! set deposition in export state
          if (im==begdst) then
             cam_out%dstdry1(:ncol) = max(sflx(:ncol), 0._r8)
          elseif(im==begdst+1) then
             cam_out%dstdry2(:ncol) = max(sflx(:ncol), 0._r8)
          elseif(im==begdst+2) then
             cam_out%dstdry3(:ncol) = max(sflx(:ncol), 0._r8)
          elseif(im==begdst+3) then
             cam_out%dstdry4(:ncol) = max(sflx(:ncol), 0._r8)
          endif
       endif
    end do

  endsubroutine aero_model_drydep

  !=============================================================================
  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    use wetdep,        only : wetdepa_v1, wetdep_inputs_set, wetdep_inputs_t
    use dust_model,    only : dust_names
  !  use seasalt_model, only : sslt_names=>seasalt_names

    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars

    integer  :: ncol                     ! number of atmospheric columns
    integer  :: lchnk                    ! chunk identifier
    integer  :: m,mm, i,k

    real(r8) :: sflx_tot_dst(pcols)
    real(r8) :: sflx_tot_slt(pcols)

    real(r8) :: iscavt(pcols, pver)
    real(r8) :: scavt(pcols, pver)
    real(r8) :: scavcoef(pcols,pver)     ! Dana and Hales coefficient (/mm) (0.1)
    real(r8) :: sflx(pcols)              ! deposition flux

    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)

    real(r8) :: sol_factb, sol_facti

    real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
    real(r8) :: cldv(pcols,pver)         ! cloudy volume undergoing scavenging
    real(r8) :: cldvcu(pcols,pver)       ! Convective precipitation area at the top interface of current layer
    real(r8) :: cldvst(pcols,pver)       ! Stratiform precipitation area at the top interface of current layer

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    type(wetdep_inputs_t) :: dep_inputs  ! obj that contains inputs to wetdepa routine

    character(len=*), parameter :: subname = 'aero_model_wetdep'


    if (nwetdep<1) return

    call endrun(subname//":: is not yet implemented")


    call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq)

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    lchnk = state%lchnk
    ncol  = state%ncol

    sflx_tot_dst(:) = 0._r8
    sflx_tot_slt(:) = 0._r8

    do m = 1, nwetdep

       sol_factb = aer_sol_factb(m)
       sol_facti = aer_sol_facti(m)

       scavcoef(:ncol,:) = aer_scav_coef(m)

       call wetdepa_v1( state%t, state%pmid, state%q(:,:,1), state%pdel, &
            dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
            dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
            dep_inputs%evapr, dep_inputs%totcond, state%q(:,:,mm), dt, &
            scavt, iscavt, dep_inputs%cldv, &
            fracis(:,:,mm), sol_factb, ncol, &
            scavcoef, &
            sol_facti_in=sol_facti, &
            icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt )

       ptend%q(:ncol,:,mm)=scavt(:ncol,:)

       call outfld( trim(cnst_name(mm))//'WET', ptend%q(:,:,mm), pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SIC', icscavt , pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
       call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

       sflx(:)=0._r8

       do k=1,pver
          do i=1,ncol
             sflx(i)=sflx(i)+ptend%q(i,k,mm)*state%pdel(i,k)/gravit
          enddo
       enddo
       call outfld( trim(cnst_name(mm))//'SFWET', sflx, pcols, lchnk)

     !  if ( any( sslt_names(:)==trim(cnst_name(mm)) ) ) &
     !       sflx_tot_slt(:ncol) = sflx_tot_slt(:ncol) + sflx(:ncol)
       if ( any( dust_names(:)==trim(cnst_name(mm)) ) ) &
            sflx_tot_dst(:ncol) = sflx_tot_dst(:ncol) + sflx(:ncol)

       ! if the user has specified prescribed aerosol dep fluxes then
       ! do not set cam_out dep fluxes according to the prognostic aerosols
       if (.not.aerodep_flx_prescribed()) then
          ! export deposition fluxes to coupler ??? why "-" sign ???
          if (trim(cnst_name(mm))=='CB2') then
             cam_out%bcphiwet(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))=='OC2') then
             cam_out%ocphiwet(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(1))) then
             cam_out%dstwet1(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(2))) then
             cam_out%dstwet2(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(3))) then
             cam_out%dstwet3(:) = max(-sflx(:), 0._r8)
          elseif (trim(cnst_name(mm))==trim(dust_names(4))) then
             cam_out%dstwet4(:) = max(-sflx(:), 0._r8)
          endif
       endif

    enddo

  !  if (sslt_active) then
  !     call outfld( 'SSTSFWET', sflx_tot_slt, pcols, lchnk)
  !  endif
    if (dust_active) then
       call outfld( 'DSTSFWET', sflx_tot_dst, pcols, lchnk)
    endif

  endsubroutine aero_model_wetdep

  !-------------------------------------------------------------------------
  ! provides aerosol surface area info for sectional aerosols
  ! called from mo_usrrxt
  !-------------------------------------------------------------------------
  subroutine aero_model_surfarea( &
                  state, mmr, radmean, relhum, pmid, temp, strato_sad, sulfate,  m, ltrop, &
                  dlat, het1_ndx, pbuf, ncol, sfc, dm_aer, sad_total, reff_trop )

    use mo_constants, only : pi, avo => avogadro


    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: m(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_total(:,:)
    real(r8), intent(out)   :: reff_trop(:,:)

    ! local vars

    integer  :: i,k
    real(r8) :: rho_air
    real(r8) :: v, n, n_exp, r_rd, r_sd
    real(r8) :: dm_sulf, dm_sulf_wet, log_sd_sulf, sfc_sulf, sfc_nit
    real(r8) :: dm_orgc, dm_orgc_wet, log_sd_orgc, sfc_oc, sfc_soa
    real(r8) :: sfc_soai, sfc_soam, sfc_soab, sfc_soat, sfc_soax
    real(r8) :: dm_bc, dm_bc_wet, log_sd_bc, sfc_bc
    real(r8) :: rxt_sulf, rxt_nit, rxt_oc, rxt_soa
    real(r8) :: c_n2o5, c_ho2, c_no2, c_no3
    real(r8) :: s_exp

    !-----------------------------------------------------------------
    ! 	... parameters for log-normal distribution by number
    ! references:
    !   Chin et al., JAS, 59, 461, 2003
    !   Liao et al., JGR, 108(D1), 4001, 2003
    !   Martin et al., JGR, 108(D3), 4097, 2003
    !-----------------------------------------------------------------
    real(r8), parameter :: rm_sulf  = 6.95e-6_r8        ! mean radius of sulfate particles (cm) (Chin)
    real(r8), parameter :: sd_sulf  = 2.03_r8           ! standard deviation of radius for sulfate (Chin)
    real(r8), parameter :: rho_sulf = 1.7e3_r8          ! density of sulfate aerosols (kg/m3) (Chin)

    real(r8), parameter :: rm_orgc  = 2.12e-6_r8        ! mean radius of organic carbon particles (cm) (Chin)
    real(r8), parameter :: sd_orgc  = 2.20_r8           ! standard deviation of radius for OC (Chin)
    real(r8), parameter :: rho_orgc = 1.8e3_r8          ! density of OC aerosols (kg/m3) (Chin)

    real(r8), parameter :: rm_bc    = 1.18e-6_r8        ! mean radius of soot/BC particles (cm) (Chin)
    real(r8), parameter :: sd_bc    = 2.00_r8           ! standard deviation of radius for BC (Chin)
    real(r8), parameter :: rho_bc   = 1.0e3_r8          ! density of BC aerosols (kg/m3) (Chin)

    real(r8), parameter :: mw_so4 = 98.e-3_r8     ! so4 molecular wt (kg/mole)

    integer  ::  irh, rh_l, rh_u
    real(r8) ::  factor, rfac_sulf, rfac_oc, rfac_bc, rfac_ss
    logical :: zero_aerosols

    !-----------------------------------------------------------------
    ! 	... table for hygroscopic growth effect on radius (Chin et al)
    !           (no growth effect for mineral dust)
    !-----------------------------------------------------------------
    real(r8), dimension(7) :: table_rh, table_rfac_sulf, table_rfac_bc, table_rfac_oc, table_rfac_ss

    character(len=*), parameter :: subname = 'aero_model_surfarea'

    call endrun(subname//":: is not yet implemented")

  end subroutine aero_model_surfarea

  !-------------------------------------------------------------------------
  ! stub
  !-------------------------------------------------------------------------
  subroutine aero_model_strat_surfarea( state, ncol, mmr, pmid, temp, ltrop, pbuf, strato_sad, reff_strat )

    !-------------------------------------------------------------------------
    ! provides WET stratospheric aerosol surface area info for modal aerosols
    ! if modal_strat_sulfate = TRUE -- called from mo_gas_phase_chemdr
    ! Copied from bulk_aero
    !-------------------------------------------------------------------------

    ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:)
    real(r8), intent(out)   :: reff_strat(:,:)

    character(len=*), parameter :: subname = 'aero_model_strat_surfarea'

    strato_sad(:,:) = 0._r8
    reff_strat(:,:) = 0._r8

  end subroutine aero_model_strat_surfarea

  !=============================================================================
  !=============================================================================
  subroutine aero_model_gasaerexch( state, loffset, ncol, lchnk, troplev, delt, reaction_rates, &
                                    tfld, pmid, pdel, mbar, relhum, &
                                    zm,  qh2o, cwat, cldfr, cldnum, &
                                    airdens, invariants, del_h2so4_gasprod,  &
                                    vmr0, vmr, pbuf )


    use chem_mods,   only : gas_pcnst
    use mo_aerosols, only : aerosols_formation, has_aerosols
    use mo_setsox,   only : setsox, has_sox
    use mo_setsoa,   only : setsoa, has_soa

    !-----------------------------------------------------------------------
    !      ... dummy arguments
    !-----------------------------------------------------------------------
        ! dummy args
    type(physics_state), intent(in) :: state           ! Physics state variables
    integer,  intent(in) :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in) :: ncol                   ! number columns in chunk
    integer,  intent(in) :: lchnk                  ! chunk index
    integer,  intent(in) :: troplev(:)
    real(r8), intent(in) :: delt                   ! time step size (sec)
    real(r8), intent(in) :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in) :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in) :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in) :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in) :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: relhum(:,:)            ! relative humidity
    real(r8), intent(in) :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:)
    real(r8), intent(in) :: zm(:,:)
    real(r8), intent(in) :: qh2o(:,:)
    real(r8), intent(in) :: cwat(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), intent(in) :: cldfr(:,:)
    real(r8), intent(in) :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8), intent(in) :: vmr0(:,:,:)       ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)         ! mixing ratios ( vmr )

    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars

    real(r8) :: vmrcw(ncol,pver,gas_pcnst)            ! cloud-borne aerosol (vmr)

    real(r8) ::  aqso4(ncol,1)               ! aqueous phase chemistry
    real(r8) ::  aqh2so4(ncol,1)             ! aqueous phase chemistry
    real(r8) ::  aqso4_h2o2(ncol)            ! SO4 aqueous phase chemistry due to H2O2
    real(r8) ::  aqso4_o3(ncol)              ! SO4 aqueous phase chemistry due to O3
    real(r8) ::  xphlwc(ncol,pver)           ! pH value multiplied by lwc

    character(len=*), parameter :: subname = 'aero_model_gasaerexch'

   ! nstep = get_nstep()

    ! Get height of boundary layer (needed for boundary layer nucleation)
   ! call pbuf_get_field(pbuf, pblh_idx, pblh)

    ! calculate tendency due to gas phase chemistry and processes
   ! dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
   ! do icnst = 1, gas_pcnst
   !    wrk(:) = 0._r8
   !    do ilev = 1,pver
   !       wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
   !    end do
   !    name = 'GS_'//trim(solsym(icnst))
   !    call outfld( name, wrk(:ncol), ncol, lchnk )
   ! enddo

! vmr2mmr (oslo_aero and carma)
! call to qqcw2tvmr (oslo_aero and mam)
! dvmrdt and dvmrcwdt (all)

!
    ! save h2so4 change by gas phase chem (for later new particle nucleation)
 !   if (ndx_h2so4 > 0) then
 !      del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - vmr0(1:ncol,:,ndx_h2so4)
 !   endif

! TODO: aq chem setsox (all)

!

    !call endrun(subname//":: is not yet implemented")

    if (masterproc) then
        write(iulog,*) subname, ":: is not yet implemented"
    end if

  end subroutine aero_model_gasaerexch

  !=============================================================================
  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
     !use oslo_aero_control, only: dms_from_ocn ! DMS
     use constituents,      only: cnst_get_ind, sflxnam
     !use oslo_aero_ocean,   only: oslo_aero_dms_emis ! DMS
     use dust_model,        only: dust_active, dust_emis, dust_names, dust_nbin, dust_nrange

     ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

    ! local vars

    integer  :: lchnk, ncol
    integer  :: m, mm
    real(r8) :: soil_erod_tmp(pcols)
    real(r8) :: sflx(pcols)   ! accumulate over all bins for output
    !integer  :: pndx_fdms  ! DMS surface flux physics index

    character(len=*), parameter :: subname = 'aero_model_emissions'

    lchnk = state%lchnk
    ncol  = state%ncol

    if (dust_active) then

        call dust_emis( lchnk, ncol, cam_in%dstflx, cam_in%cflx, aero_props )

       ! some dust emis diagnostics ...
       sflx(:)=0._r8
        do m=1,dust_nrange
          if (m<=dust_nrange) sflx(:ncol)=sflx(:ncol)+cam_in%cflx(:ncol,m) ! TODO: check indices
          call outfld(trim(dust_names(m))//'SF',cam_in%cflx(:,m),pcols, lchnk)
       enddo
       call outfld('DSTSFMBL',sflx(:),pcols,lchnk)
       call outfld('LND_MBL',soil_erod_tmp(:),pcols, lchnk )
    endif
    call endrun(subname//":: is not yet implemented")
  end subroutine aero_model_emissions

end module aero_model
