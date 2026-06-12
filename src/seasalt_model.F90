!===============================================================================
! Seasalt for Bulk Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,       only: pcols, pver
  use spmd_utils,      only: masterproc
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun
  use aerosol_properties_mod, only: aerosol_properties
  use sectional_aerosol_properties_mod, only: sectional_aerosol_properties

  implicit none
  private

  ! public variables
  public :: seasalt_active
  public :: seasalt_nspecies

  ! public procedures
  public :: seasalt_init
  public :: seasalt_emis

  ! initialization of public variables
  logical :: seasalt_active = .false.
  integer :: seasalt_nspecies = 0

  ! module variables
  integer :: seasalt_specprop_ndx(10)

  ! positions in q arrays
  integer, allocatable :: seasalt_bin_tracer_ndx(:,:)
  integer, allocatable :: seasalt_range_tracer_ndx(:,:)

  contains

   !=============================================================================
   !=============================================================================
   subroutine seasalt_init(aero_props)
     use cam_history,   only: addfld, fieldname_len
     use constituents,  only: cnst_get_ind
     use string_utils,  only: int2str

     type(sectional_aerosol_properties), intent(in) :: aero_props

     ! local variables
     integer              :: ispecprop, ibin, irange, nname, isaltspec
     integer              :: specprop_ndx
     integer, allocatable :: bin_ndx(:,:)
     character(len=10)    :: type
     character(len=10)    :: names_tmp(200)

     character(len=fieldname_len) :: dummy
     character(len=*), parameter :: subname = 'seasalt_init'

    ! Initalize seasalt vars
     nname = 0

     do ispecprop = 1, aero_props%nspecies_tot()
        type = aero_props%spectype(ispecprop)

        if (trim(type) == 'seasalt') then
            seasalt_nspecies = seasalt_nspecies + 1
            seasalt_specprop_ndx(seasalt_nspecies) = ispecprop
        end if
    end do

    seasalt_active = seasalt_nspecies > 0

    if ( .not. seasalt_active ) return

    allocate(seasalt_bin_tracer_ndx(seasalt_nspecies, aero_props%nbins())) ! use all bins
    allocate(seasalt_range_tracer_ndx(seasalt_nspecies, aero_props%nranges()))
    allocate(bin_ndx(seasalt_nspecies, aero_props%nbins()))

    ! loop through different seasalt species, if there is more than one
    do isaltspec = 1, seasalt_nspecies
        specprop_ndx = seasalt_specprop_ndx(isaltspec)

        ! get indices of bins with the correct seasalt species in them
        bin_ndx(isaltspec,:aero_props%spec_nbin(specprop_ndx)) = aero_props%spec_bin_ndx(specprop_ndx, aero_props%spec_nbin(specprop_ndx))

        ! Get the indices for the bins containing sea salt
        do ibin = 1, aero_props%spec_nbin(specprop_ndx)
            dummy = 'num_'//int2str(bin_ndx(isaltspec,ibin))
            call cnst_get_ind(dummy, seasalt_bin_tracer_ndx(isaltspec,ibin))
        end do

        ! Get the indices for the range/mmr tracers
        do irange = 1, aero_props%spec_nrange(specprop_ndx)
            dummy = aero_props%spec_tracernames(specprop_ndx, irange)
            call cnst_get_ind(dummy, seasalt_range_tracer_ndx(isaltspec,irange))
        end do
    end do

   end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx, aero_props )

    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)
    type(sectional_aerosol_properties), intent(in) :: aero_props

    ! local vars
    character(len=*), parameter :: subname = 'seasalt_emis'

 ! from CARMA: carma_model_mod.F90 /src/physics/carma/models/sea_salt
  !--------CMS (Clarke, Monahan, and Smith source function)-------

    ! ------------------------------------------------------------
    ! ----  Clarke Source Function. Coefficients for Ai    -------
    ! ------------------------------------------------------------
    real(r8), parameter :: beta01 =-5.001e3_r8
    real(r8), parameter :: beta11 = 0.808e6_r8
    real(r8), parameter :: beta21 =-1.980e7_r8
    real(r8), parameter :: beta31 = 2.188e8_r8
    real(r8), parameter :: beta41 =-1.144e9_r8
    real(r8), parameter :: beta51 = 2.290e9_r8
    real(r8), parameter :: beta02 = 3.854e3_r8
    real(r8), parameter :: beta12 = 1.168e4_r8
    real(r8), parameter :: beta22 =-6.572e4_r8
    real(r8), parameter :: beta32 = 1.003e5_r8
    real(r8), parameter :: beta42 =-6.407e4_r8
    real(r8), parameter :: beta52 = 1.493e4_r8
    real(r8), parameter :: beta03 = 4.498e2_r8
    real(r8), parameter :: beta13 = 0.839e3_r8
    real(r8), parameter :: beta23 =-5.394e2_r8
    real(r8), parameter :: beta33 = 1.218e2_r8
    real(r8), parameter :: beta43 =-1.213e1_r8
    real(r8), parameter :: beta53 = 4.514e-1_r8
    real(r8)            :: A1                                ! Coefficient Ak in Clarkes's source function
    real(r8)            :: A2
    real(r8)            :: A3
    real(r8)            :: wcap                         ! whitecap coverage
    real(r8)            :: rpdry,rpdry_cm                              ! dry radius
    real(r8)            :: Monahan, Clarke, Smith             ! dF/dr [#/m2/s/um]
    real(r8), allocatable :: bin_ndx(:)
    real(r8), allocatable :: bin_centers(:)
    integer             :: ibin, icol
    real(r8)            :: B_mona                            ! the parameter used in Monahan
    real(r8)            :: r80, r80cm
    real(r8)            :: ncflx                              ! dF/dr [#/m2/s/um]
    real(r8), parameter :: xkar = 0.4_r8                      ! Von Karman constant

    ! --------------------------------------------------------------------
    ! ---- constants in calculating the particle wet radius [Gerber, 1985]
    ! --------------------------------------------------------------------
    real(r8), parameter :: c1   = 0.7674_r8        ! .
    real(r8), parameter :: c2   = 3.079_r8         ! .
    real(r8), parameter :: c3   = 2.573e-11_r8     ! .
    real(r8), parameter :: c4   = -1.424_r8        ! constants in calculating the particle wet radius

    real(r8)            :: u14, ustar_smith, cd_smith         ! 14m wind velocity, friction velocity and drag
    ! ---------------------------------------------
    ! coefficient A1, A2 in Andreas's Source funcion
    ! ---------------------------------------------
    real(r8)            ::A1A92
    real(r8)            ::A2A92
                                                              ! coefficient as desired by Andreas source function

    ! ---------------------------------------------
    ! coefficient in Smith's Source funcion
    ! ---------------------------------------------
    real(r8), parameter ::  f1 = 3.1_r8
    real(r8), parameter ::  f2 = 3.3_r8
    real(r8), parameter ::  r1 = 2.1_r8
    real(r8), parameter ::  r2 = 9.2_r8

    if (.not. seasalt_active) return

    allocate(bin_ndx(aero_props%spec_nbin(seasalt_specprop_ndx(1))))
    allocate(bin_centers(aero_props%nbins()))
          ! Add any surface flux here.
          ncflx       = 0.0_r8
          Monahan     = 0.0_r8
          Clarke      = 0.0_r8
          Smith       = 0.0_r8

    bin_ndx = aero_props%spec_bin_ndx(seasalt_specprop_ndx(1), aero_props%spec_nbin(seasalt_specprop_ndx(1)))
    bin_centers = aero_props%bin_centers(aero_props%nbins())

   ! do ibin = seasalt_bin_start, seasalt_bin_end
    do ibin = bin_ndx(1), bin_ndx(aero_props%spec_nbin(seasalt_specprop_ndx(1)))
! TODO: check unit of rpdry

        do icol = 1, ncol
            wcap = 3.84e-6_r8 * u10cubed(icol)      ! in percent, ie., 75%, wcap = 0.75
!    W(:ncol)=3.84e-6_r8*u10cubed(:ncol)*0.1_r8 ! whitecap area

            rpdry = bin_centers(ibin) ! in m
            rpdry_cm = rpdry * 100._r8

            r80cm   = (c1 *  (rpdry_cm) ** c2 / (c3 * rpdry_cm ** c4 - log10(0.8_r8)) + (rpdry_cm)**3._r8) ** (1._r8/3._r8) ! [cm]
            r80     = r80cm *1.e4_r8    ! [um]

            A1 = beta01 + beta11*(2._r8*rpdry) + beta21*(2._r8*rpdry)**2 + &
                beta31*(2._r8*rpdry)**3 + beta41*(2._r8*rpdry)**4 + beta51*(2._r8*rpdry)**5
            A2 = beta02 + beta12*(2._r8*rpdry) + beta22*(2._r8*rpdry)**2 + &
                 beta32*(2._r8*rpdry)**3 + beta42*(2._r8*rpdry)**4 + beta52*(2._r8*rpdry)**5
            A3 = beta03 + beta13*(2._r8*rpdry) + beta23*(2._r8*rpdry)**2 + &
                 beta33*(2._r8*rpdry)**3 + beta43*(2._r8*rpdry)**4 + beta53*(2._r8*rpdry)**5

! TODO: replaced rdry with rpdry here
            !Clarke
            if (rpdry .lt. 0.066_r8) then
              Clarke = A1 * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2]
              Clarke = Clarke / (2.30258509_r8 * rpdry)               ! dF/dr    [#/s/m2/um]
            elseif ((rpdry .ge. 0.066_r8) .and. (rpdry .lt. 0.6_r8)) then
              Clarke = A2 * 1.e4_r8 * wcap                     ! dF/dlogr [#/s/m2]
              Clarke = Clarke / (2.30258509_r8 * rpdry)               ! dF/dr    [#/s/m2/um]
            elseif ((rpdry .ge. 0.6_r8) .and. (rpdry .lt. 4.0_r8)) then
              Clarke = A3 * 1.e4_r8 * wcap                      ! dF/dlogr [#/s/m2]
              Clarke = Clarke / (2.30258509_r8 * rpdry)                 ! dF/dr    [#/s/m2/um]
            end if

            !Monahan
            B_Mona = (0.38_r8 - log10(r80)) / 0.65_r8
            Monahan = 1.373_r8 * u10cubed(icol) * r80 ** (-3._r8) * &
                 (1._r8 + 0.057_r8 * r80**1.05_r8) * 10._r8 ** (1.19_r8 * exp(- B_Mona **2))

            !Smith
            u14 = u10cubed(icol) ** (1._r8/3.41_r8) * (1._r8 + cd_smith**0.5_r8 / xkar*log(14._r8 / 10._r8))  ! 14 meter wind
            A1A92 = 10._r8 ** (0.0676_r8 * u14 + 2.430_r8)
            A2A92 = 10._r8 ** (0.9590_r8 * u14**0.5_r8 - 1.476_r8)
            Smith = A1A92*exp(-f1 *(log(r80 / r1))**2) + A2A92*exp(-f2 * (log(r80 / r2))**2)     ! dF/dr   [#/m2/s/um]

            !%%%%%%%%%%%%%%%%%%%%%%%%%
            !     CMS1 or CMS2
            !%%%%%%%%%%%%%%%%%%%%%%%%%
  !          if (rdry .lt. 0.1_r8) then   ! originally cut at 0.1 um
            ! ***CMS1*****
            if (rpdry .lt. 1e-6_r8) then    ! cut at 1.0 um, rpdry in [m]
            ! ***CMS2*****
  !          if (rdry .lt. 2._r8) then    ! cut at 2.0 um
                ncflx = Clarke
            else
                if (u10cubed(icol)**(1._r8/3.41_r8) .lt. 9._r8) then
                    ncflx = Monahan
                else
                    if (Monahan .gt. Smith) then
                        ncflx = Monahan
                    else
                        ncflx = Smith
                    end if
                end if
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Apply Hoppel correction
            !%%%%%%%%%%%%%%%%%%%%%%%%%
!            ncflx = ncflx * fref
!          end if

          ! convert ncflx [#/m^2/s/um] to surfaceFlx [kg/m^2/s]
!          surfaceFlux(icol) = ncflx * dr(ibin) * rmass(ibin) * 10._r8      ! *1e4[um/cm] * 1.e-3[kg/g]

          ! weighted by the ocean fraction
!          surfaceFlux(icol) = surfaceFlux(icol) * cam_in%ocnfrac(icol)

        end do

    end do


  end subroutine seasalt_emis

end module seasalt_model
