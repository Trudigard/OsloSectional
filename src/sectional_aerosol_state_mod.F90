module sectional_aerosol_state_mod

! TODO: make object to internally hold bin number concentrations and range bulk masses
! TODO: Range density function/array + set_density and get_density?
! TODO: Update range function
! TODO: range_state object: density, mass
! TODO: bin_state: mmr/number, surface area, hygroscopicity, ...

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_spfn_mod, only: erf => shr_spfn_erf
  use aerosol_state_mod, only: aerosol_state, ptr2d_t
  use physics_types, only: physics_state
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use sectional_aerosol_properties_mod, only: sectional_aerosol_properties
  use physconst,  only: rhoh2o, mwh2o

  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use ppgrid,         only: pcols, pver

  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index

  implicit none

  private

  public :: sectional_aerosol_state

  type aerosol_range_state ! one instance per range
!TODO: chunks?
     real(r8), allocatable :: dry_density(:,:)        ! density of the species mixture in a range without water
     real(r8), allocatable :: hygroscopicity(:,:)     ! hygroscopicity of the species mixture
! TODO: mass as separate array?
     real(r8), allocatable :: mass(:, :, :)            ! mass of each species in this range len = maxval(range_nspecies) so that index is same in each range
     ! ...
  end type aerosol_range_state

  type, extends(aerosol_state) :: sectional_aerosol_state
     private
     type(aerosol_range_state), allocatable :: aer_range_state(:)
     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     real(r8), allocatable :: bin_numconc(:,:,:)
   contains

     procedure :: get_transported
     procedure :: set_transported
     procedure :: ambient_total_bin_mmr
     procedure :: get_ambient_mmr_0list
     procedure :: get_ambient_mmr_rlist
     procedure :: get_cldbrne_mmr
     procedure :: get_ambient_num
     procedure :: get_cldbrne_num
     procedure :: get_states
     procedure :: icenuc_size_wght_arr
     procedure :: icenuc_size_wght_val
     procedure :: icenuc_type_wght
     procedure :: update_bin
     procedure :: hetfrz_size_wght
     procedure :: hygroscopicity
     procedure :: water_uptake
     procedure :: dry_volume
     procedure :: wet_volume
     procedure :: water_volume
     procedure :: wet_diameter
     procedure :: convcld_actfrac
     procedure :: wgtpct
     procedure :: dry_density
     procedure :: update_range

     final :: destructor

  end type sectional_aerosol_state

  interface
    function sas_state_obj_initialize(state, pbuf, aero_props) result(newobj)
      import :: physics_state, physics_buffer_desc, sectional_aerosol_state, sectional_aerosol_properties
      type(physics_state), target :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      type(sectional_aerosol_state), pointer :: newobj
      type(sectional_aerosol_properties), intent(in) :: aero_props
    end function sas_state_obj_initialize
  end interface
  
  interface sectional_aerosol_state
     procedure :: constructor(sas_state_obj_initialize)
  end interface sectional_aerosol_state

  real(r8), parameter :: rh2odens = 1._r8/rhoh2o

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf, aero_props) result(newobj)
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(sectional_aerosol_state), pointer :: newobj
    type(sectional_aerosol_properties), intent(in) :: aero_props

    integer :: ierr, maxspec, irange

    character(len=*), parameter :: subname = 'constructor'

    call endrun(subname//' is not yet implemented')

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    newobj%state => state
    newobj%pbuf => pbuf
    maxspec = maxval(aero_props%range_nspecies(aero_props%nranges()))

    allocate(newobj%aer_range_state(aero_props%nranges()), stat=ierr)
    if( ierr /= 0 ) then
        nullify(newobj)
        return
    end if
    do irange = 1, aero_props%nranges()
        allocate(newobj%aer_range_state(irange)%mass(maxspec, pcols, pver), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if
        allocate(newobj%aer_range_state(irange)%dry_density(pcols, pver), stat=ierr)
        if( ierr /= 0 )then
            nullify(newobj)
            return
        end if
        allocate(newobj%aer_range_state(irange)%hygroscopicity(pcols,pver), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

       ! newobj%aer_range_state%dry_density = 0._r8
       ! newobj%aer_range_state%hygroscopicity = 0._r8
       ! newobj%aer_range_state%mass = 0._r8
    end do

    allocate(newobj%bin_numconc(aero_props%nbins(),pcols, pver), stat=ierr)
    if( ierr /= 0 ) then
        nullify(newobj)
        return
    end if

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(sectional_aerosol_state), intent(inout) :: self

    character(len=*), parameter :: subname = 'destructor'

    call endrun(subname//' is not yet implemented')

    nullify(self%state)
    nullify(self%pbuf)
    if (allocated(self%aer_range_state)) then
        deallocate(self%aer_range_state)
    end if
    if (allocated(self%bin_numconc)) then
        deallocate(self%bin_numconc)
    end if

end subroutine destructor

  !------------------------------------------------------------------------------
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(sectional_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)

    character(len=*), parameter :: subname = 'set_transported'

    call endrun(subname//' is not yet implemented')
    ! BEFORE advection time step
    ! Connect internal arrays for bin_number concentrations to num_1, num_2, ...
    ! and internal masses for ranges to DU_R3, DU_R4, etc
    ! internal bin_numconc -> cnst_get_ind for num_1, num_2, ...
    ! internal aer_range_state%mass -> cnst_get_ind for DU_R3, DU_R4, ...
    ! internal bin_numconc and aer_range_state%mass = 0._r8
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(sectional_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)

    character(len=*), parameter :: subname = 'get_transported'

    call endrun(subname//' is not yet implemented')
    ! AFTER advection time step
    ! Retrieve new values for number and masses and put them back into internal array
    ! cnst_get_ind num_1, num_2, ... -> bin_numconc
    !
  end subroutine get_transported

  !------------------------------------------------------------------------
  ! Total aerosol mass mixing ratio for a bin in a given grid box location (column and layer)
  !------------------------------------------------------------------------
  function ambient_total_bin_mmr(self, aero_props, bin_ndx, col_ndx, lyr_ndx) result(mmr_tot)
    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: bin_ndx      ! bin index
    integer, intent(in) :: col_ndx      ! column index
    integer, intent(in) :: lyr_ndx      ! vertical layer index

    real(r8) :: mmr_tot                 ! mass mixing ratios totaled for all species

    character(len=*), parameter :: subname = 'ambient_total_bin_mmr'

    call endrun(subname//' is not yet implemented')

    ! bin_num -> convert to mmr using range%density and airdensity

  end function ambient_total_bin_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_0list(self, species_ndx, bin_ndx, mmr)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    character(len=*), parameter :: subname = 'get_ambient_mmr_0list'

    call endrun(subname//' is not yet implemented')

    ! use range index instead of bin index here?
    ! use range species mass and airdens (+num for a single bin)

  end subroutine get_ambient_mmr_0list

  !------------------------------------------------------------------------------
  ! returns ambient aerosol mass mixing ratio for a given radiation diagnostics
  ! list index, species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_mmr_rlist(self, list_ndx, species_ndx, bin_ndx, mmr)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: list_ndx     ! rad climate list index
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    character(len=*), parameter :: subname = 'get_ambient_mmr_rlist'

    call endrun(subname//' is not yet implemented')

  end subroutine get_ambient_mmr_rlist

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_mmr(self, species_ndx, bin_ndx, mmr)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: species_ndx  ! species index
    integer, intent(in) :: bin_ndx      ! bin index
    real(r8), pointer :: mmr(:,:)       ! mass mixing ratios (ncol,nlev)

    character(len=*), parameter :: subname = 'get_cldbrne_mmr'

    call endrun(subname//' is not yet implemented')

  end subroutine get_cldbrne_mmr

  !------------------------------------------------------------------------------
  ! returns ambient aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_ambient_num(self, bin_ndx, num)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx     ! bin index
    real(r8), pointer   :: num(:,:)    ! number densities

    character(len=*), parameter :: subname = 'get_ambient_num'

    call endrun(subname//' is not yet implemented')

    ! use bin_num, range_density and range_mass
  end subroutine get_ambient_num

  !------------------------------------------------------------------------------
  ! returns cloud-borne aerosol number mixing ratio for a given species index and bin index
  !------------------------------------------------------------------------------
  subroutine get_cldbrne_num(self, bin_ndx, num)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    real(r8), pointer :: num(:,:)

    character(len=*), parameter :: subname = 'get_cldbrne_num'

    call endrun(subname//' is not yet implemented')

  end subroutine get_cldbrne_num

  !------------------------------------------------------------------------------
  ! returns interstitial and cloud-borne aerosol states
  !------------------------------------------------------------------------------
  subroutine get_states( self, aero_props, raer, qqcw )
    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    type(ptr2d_t), intent(out) :: raer(:)
    type(ptr2d_t), intent(out) :: qqcw(:)

    integer :: ibin,ispc, iidx

    character(len=*), parameter :: subname = 'get_states'

    call endrun(subname//' is not yet implemented')


  end subroutine get_states

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_arr(self, bin_ndx, ncol, nlev, species_type, use_preexisting_ice, wght)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice ! pre-existing ice flag
    real(r8), intent(out) :: wght(:,:)

    character(len=*), parameter :: subname = 'icenuc_size_wght_arr'

    call endrun(subname//' is not yet implemented')
! ??
  end subroutine icenuc_size_wght_arr

  !------------------------------------------------------------------------------
  ! return aerosol bin size weights for a given bin, column and vertical layer
  !------------------------------------------------------------------------------
  subroutine icenuc_size_wght_val(self, bin_ndx, col_ndx, lyr_ndx, species_type, use_preexisting_ice, wght)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    character(len=*), intent(in) :: species_type  ! species type
    logical, intent(in) :: use_preexisting_ice    ! pre-existing ice flag
    real(r8), intent(out) :: wght

    character(len=*), parameter :: subname = 'icenuc_size_wght_val'

    call endrun(subname//' is not yet implemented')

  end subroutine icenuc_size_wght_val

  !------------------------------------------------------------------------------
  ! returns aerosol type weights for a given aerosol type and bin
  !------------------------------------------------------------------------------
  subroutine icenuc_type_wght(self, bin_ndx, ncol, nlev, species_type, aero_props, rho, wght, cloud_borne)

    use aerosol_properties_mod, only: aerosol_properties

    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: ncol                   ! number of columns
    integer, intent(in) :: nlev                   ! number of vertical levels
    character(len=*), intent(in) :: species_type  ! species type
    class(aerosol_properties), intent(in) :: aero_props ! aerosol properties object
    real(r8), intent(in) :: rho(:,:)              ! air density (kg m-3)
    real(r8), intent(out) :: wght(:,:)            ! type weights
    logical, optional, intent(in) :: cloud_borne  ! if TRUE cloud-borne aerosols are used
                                                  ! otherwise ambient aerosols are used

    character(len=*), parameter :: subname = 'icenuc_type_wght'

    call endrun(subname//' is not yet implemented')

  end subroutine icenuc_type_wght

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine update_bin( self, bin_ndx, col_ndx, lyr_ndx, delmmr_sum, delnum_sum, tnd_ndx, dtime, tend )
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx                ! bin number
    integer, intent(in) :: col_ndx                ! column index
    integer, intent(in) :: lyr_ndx                ! vertical layer index
    real(r8),intent(in) :: delmmr_sum             ! mass mixing ratio change summed over all species in bin
    real(r8),intent(in) :: delnum_sum             ! number mixing ratio change summed over all species in bin
    integer, intent(in) :: tnd_ndx                ! tendency index
    real(r8),intent(in) :: dtime                  ! time step size (sec)
    real(r8),intent(inout) :: tend(:,:,:)         ! tendency

    character(len=*), parameter :: subname = 'update_bin'

    call endrun(subname//' is not yet implemented')

    ! bin_num(bin_ndx) = bin_num + tendency
  end subroutine update_bin

  !------------------------------------------------------------------------------
  ! returns the volume-weighted fractions of aerosol subset `bin_ndx` that can act
  ! as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_size_wght(self, bin_ndx, ncol, nlev) result(wght)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin number
    integer, intent(in) :: ncol                ! number of columns
    integer, intent(in) :: nlev                ! number of vertical levels

    real(r8) :: wght(ncol,nlev)                 !

    character(len=*), parameter :: subname = 'hetfrz_size_wght'

    call endrun(subname//' is not yet implemented')

  end function hetfrz_size_wght

  !------------------------------------------------------------------------------
  ! returns hygroscopicity for a given radiation diagnostic list number and
  ! bin number
  !------------------------------------------------------------------------------
  subroutine hygroscopicity(self, list_ndx, bin_ndx, kappa)
    class(sectional_aerosol_state), intent(in) :: self
 !   class(aerosol_properties), intent(in) :: aero_props
! TODO: what is list_ndx?
    integer, intent(in) :: list_ndx        ! rad climate list number
    integer, intent(in) :: bin_ndx         ! bin number
    real(r8), intent(out) :: kappa(:,:)                 !
    real(r8), allocatable :: bins2ranges(:)
    integer :: irange

    character(len=*), parameter :: subname = 'hygroscopicity'

 !   allocate(bins2ranges(aero_props%nbins()))
    !bins2ranges = aero_props%bins2ranges(aero_props%nbins())
    !irange = bins2ranges(bin_ndx)
    !kappa = self%aer_range_state%hygroscopicity(irange)

  end subroutine hygroscopicity

  !------------------------------------------------------------------------------
  ! returns aerosol wet diameter and aerosol water concentration for a given
  ! radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  subroutine water_uptake(self, aero_props, list_ndx, bin_ndx, ncol, nlev, dgnumwet, qaerwat)

    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in) :: list_ndx             ! rad climate/diags list number
    integer, intent(in) :: bin_ndx              ! bin number
    integer, intent(in) :: ncol                 ! number of columns
    integer, intent(in) :: nlev                 ! number of levels
    real(r8),intent(out) :: dgnumwet(ncol,nlev) ! aerosol wet diameter (m)
    real(r8),intent(out) :: qaerwat(ncol,nlev)  ! aerosol water concentration (g/g)

    character(len=*), parameter :: subname = 'water_uptake'

    call endrun(subname//' is not yet implemented')

  end subroutine water_uptake

  !------------------------------------------------------------------------------
  ! aerosol dry volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function dry_volume(self, aero_props, list_ndx, bin_ndx, ncol, nlev) result(vol)

    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_ndx  ! rad climate/diags list number
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8), pointer :: mmr(:,:)
    real(r8) :: specdens              ! species density (kg/m3)

    integer :: ispec

    character(len=*), parameter :: subname = 'dry_volume'

    call endrun(subname//' is not yet implemented')
! bin_num * aero_props%volume

  end function dry_volume

  !------------------------------------------------------------------------------
  ! aerosol wet volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function wet_volume(self, aero_props, list_ndx, bin_ndx, ncol, nlev) result(vol)

    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_ndx  ! rad climate/diags list number
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: dryvol(ncol,nlev)
    real(r8) :: watervol(ncol,nlev)

    character(len=*), parameter :: subname = 'wet_volume'

    call endrun(subname//' is not yet implemented')
! dry_volume + water volume

  end function wet_volume

  !------------------------------------------------------------------------------
  ! aerosol water volume (m3/kg) for given radiation diagnostic list number and bin number
  !------------------------------------------------------------------------------
  function water_volume(self, aero_props, list_ndx, bin_ndx, ncol, nlev) result(vol)

    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props

    integer, intent(in) :: list_ndx  ! rad climate/diags list number
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: vol(ncol,nlev)       ! m3/kg

    real(r8) :: dgnumwet(ncol,nlev)
    real(r8) :: qaerwat(ncol,nlev)

    character(len=*), parameter :: subname = 'water_volume'

    call endrun(subname//' is not yet implemented')

  end function water_volume

  !------------------------------------------------------------------------------
  ! aerosol wet diameter
  !------------------------------------------------------------------------------
  function wet_diameter(self, bin_ndx, ncol, nlev) result(diam)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin number
    integer, intent(in) :: ncol      ! number of columns
    integer, intent(in) :: nlev      ! number of levels

    real(r8) :: diam(ncol,nlev)

    real(r8), pointer :: dgnumwet(:,:,:)

    character(len=*), parameter :: subname = 'wet_diameter'

    call endrun(subname//' is not yet implemented')
! wet_volume/bin_num

  end function wet_diameter

  !------------------------------------------------------------------------------
  ! prescribed aerosol activation fraction for convective cloud
  !------------------------------------------------------------------------------
  function convcld_actfrac(self, ibin, ispc, ncol, nlev) result(frac)

    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) :: ibin   ! bin index
    integer, intent(in) :: ispc   ! species index
    integer, intent(in) :: ncol   ! number of columns
    integer, intent(in) :: nlev   ! number of vertical levels

    real(r8) :: frac(ncol,nlev)

    character(len=*), parameter :: subname = 'convcld_actfrac'

    call endrun(subname//' is not yet implemented')

  end function convcld_actfrac

  !------------------------------------------------------------------------------
  ! aerosol weight precent of H2SO4/H2O solution
  !------------------------------------------------------------------------------
  function wgtpct(self, ncol, nlev) result(wtp)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in) ::  ncol, nlev
    real(r8) :: wtp(ncol,nlev)  ! weight precent of H2SO4/H2O solution for given icol, ilev
    character(len=*), parameter :: subname = 'wgtpct'

    call endrun(subname//' is not yet implemented')

  end function wgtpct

  real(r8) function dry_density(self, bin_ndx, aero_props)
    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    integer, intent(in)   :: bin_ndx
    real(r8), allocatable :: bins2ranges(:)
    integer               :: irange

    character(len=*), parameter :: subname = 'dry_density'

 !   allocate(bins2ranges(aero_props%nbins()))
 !   bins2ranges = aero_props%bins2ranges(aero_props%nbins())

!    irange = bins2ranges(bin_ndx)
!    aero_state%aer_range_state%density(irange)

  end function dry_density

  subroutine update_range(self, mass_tend, nranges, aero_props)
    class(sectional_aerosol_state), intent(in) :: self
    class(aerosol_properties), intent(in) :: aero_props
    real(r8), intent(in) :: mass_tend(:,:) ! shape aero_props%range_nspecies (range, (max(nspecies in range))
    integer              :: irange, nranges
    real(r8)             :: range_dry_volume
    real(r8)             :: range_bounds(nranges,2)

! update state of range
    ! update mass
 !   self%aer_range_state%mass = self%aer_range_state%mass + mass_tend
    ! reset density and hygroscopicity
 !   self%aer_range_state%dry_density = 0._r8
 !   self%aer_range_statehygroscopicity = 0._r8
 !   range_bounds = aero_props%range_bounds()
 !   do irange = 1, aero_props%nranges()
 !       range_dry_volume = sum(dry_volume(range_bounds(irange,1), range_bounds(irange, 2)))
 !       self%aer_range_statedry_density(irange) = self%aer_range_state%mass(:,irange)/range_dry_volume
 !       do ispec = 1, max(aero_props%range_nspecies)
! TODO: source
 !           self%aer_range_state%hygroscopicity(irange) = self%aer_range_state%hygroscopicity(irange) + &
 !                    mass_species/density_species/range_dry_volume * aero_props%kappa(ispec)
  !      end do
  !  end do

  end subroutine update_range

end module sectional_aerosol_state_mod
