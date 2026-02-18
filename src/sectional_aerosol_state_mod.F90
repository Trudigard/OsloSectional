module sectional_aerosol_state_mod

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
  use ppgrid,         only: pver

  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index

  implicit none

  private

  public :: sectional_aerosol_state

  type aerosol_range_state ! one instance per range
     character(len=16), allocatable :: range_name(:)
     integer, allocatable  :: transport_ndx(:)
     integer, allocatable  :: spec_ndx(:)             ! same length as transport_ndx or third dimension of mass(:,:,range_nspecies), indices corresponding to species properties object
     real(r8), allocatable :: dry_density(:,:)        ! density of the species mixture in a range without water, ncol, pver
     real(r8), allocatable :: hygroscopicity(:,:)     ! hygroscopicity of the species mixture
     real(r8), allocatable :: mass(:, :, :)           ! (ncol, pver, range_nspecies)
     ! ...

  end type aerosol_range_state

  type, extends(aerosol_state) :: sectional_aerosol_state
     !private

     type(physics_state), pointer :: state => null()
     type(physics_buffer_desc), pointer :: pbuf(:) => null()
     type(sectional_aerosol_properties), pointer :: sec_aero_props => null()
     real(r8), allocatable :: bin_numconc(:,:,:)
     integer, allocatable :: num_transport_ndx(:)
     type(aerosol_range_state), allocatable :: aero_range_state(:)
     integer :: ncol

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
     procedure :: bin_dry_density
     procedure :: update_range

     final :: destructor

  end type sectional_aerosol_state

  interface sectional_aerosol_state
     procedure :: constructor
  end interface sectional_aerosol_state

  real(r8), parameter :: rh2odens = 1._r8/rhoh2o

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(state,pbuf) result(newobj)
    use mo_tracname, only: solsym
    use constituents, only: cnst_get_ind
    use string_utils, only: int2str
    type(physics_state), target :: state
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(sectional_aerosol_state), pointer :: newobj
    type(sectional_aerosol_properties), target :: aero_props
    integer :: ierr, irange, solsym_ndx, ispec, ubar_ndx, ibin, ispecprops
    character(len=:), allocatable :: num_name, specname, specname_props
    logical :: solsym_found
    character(len=*), parameter :: subname = 'constructor'

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    newobj%state => state
    newobj%pbuf => pbuf
    newobj%sec_aero_props => sectional_aerosol_properties()

    newobj%ncol = state%ncol

    allocate(newobj%bin_numconc(newobj%ncol, pver, aero_props%nbins()), stat=ierr)
    if( ierr /= 0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%aero_range_state(aero_props%nranges()), stat=ierr)
    if( ierr /= 0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%num_transport_ndx(aero_props%nbins()), stat=ierr)
    if( ierr /= 0 ) then
        nullify(newobj)
        return
    end if

    do irange = 1, aero_props%nranges()
        allocate(newobj%aero_range_state(irange)%mass(newobj%ncol, pver, newobj%sec_aero_props%range_nspecies(irange) ), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if
        allocate(newobj%aero_range_state(irange)%dry_density(newobj%ncol, pver), stat=ierr)
        if( ierr /= 0 )then
            nullify(newobj)
            return
        end if
        allocate(newobj%aero_range_state(irange)%hygroscopicity(newobj%ncol, pver), stat=ierr)
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if
        allocate(newobj%aero_range_state(irange)%range_name(newobj%sec_aero_props%range_nspecies(irange)))
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if
        allocate(newobj%aero_range_state(irange)%transport_ndx(newobj%sec_aero_props%range_nspecies(irange)))
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        allocate(newobj%aero_range_state(irange)%spec_ndx(newobj%sec_aero_props%range_nspecies(irange)))
        if( ierr /= 0 ) then
            nullify(newobj)
            return
        end if

        newobj%aero_range_state(irange)%dry_density = 0._r8
        newobj%aero_range_state(irange)%hygroscopicity = 0._r8
        newobj%aero_range_state(irange)%mass = 0._r8
        newobj%aero_range_state(irange)%range_name = ''
        newobj%aero_range_state(irange)%transport_ndx = 0
        newobj%aero_range_state(irange)%spec_ndx = 0

        ispec=0
        do solsym_ndx = 1, size(solsym)
            ubar_ndx = index(solsym(solsym_ndx), '_R'//int2str(irange))
            if ( ubar_ndx > 0 ) then
                ispec = ispec + 1
                if ( ispec > newobj%sec_aero_props%range_nspecies(irange) ) then
                    call endrun(subname//':: ERROR : number of species is larger than number of species in range')
                end if
                newobj%aero_range_state(irange)%range_name(ispec) = solsym(solsym_ndx)
                call cnst_get_ind(solsym(solsym_ndx), newobj%aero_range_state(irange)%transport_ndx(ispec), abort=.false.)
                if ( newobj%aero_range_state(irange)%transport_ndx(ispec) < 0 ) then
                    call endrun(subname//":: ERROR: transport array index for"//trim(solsym(solsym_ndx))//" not found")
                end if

                ! add index to connect to the species objects in aero_props
                specname = solsym(solsym_ndx)
                do ispecprops = 1, newobj%sec_aero_props%nspecies_tot()
                    call newobj%sec_aero_props%get(bin_ndx=1,species_ndx=ispecprops, specname=specname_props)
                    if ( specname(1:ubar_ndx-1) == trim( specname_props ) ) then ! before ubar_ndx -> species name
                        newobj%aero_range_state(irange)%spec_ndx(ispec) = ispecprops
                    end if
                end do

            end if
        end do

    end do

    do ibin = 1, newobj%sec_aero_props%nbins()
        num_name = 'num_'//int2str(ibin)
        solsym_found = .false.
        do solsym_ndx = 1, size(solsym)
            if ( trim(solsym(solsym_ndx)) == trim(num_name) ) then
                solsym_found = .true.
                exit
            end if
        end do
        if ( .not. solsym_found ) then
            call endrun(subname//':: ERROR: bin '//trim(num_name)//' not found')
        end if
        call cnst_get_ind(num_name, newobj%num_transport_ndx(ibin), abort = .false.)
        if ( newobj%num_transport_ndx(ibin) < 0 ) then
            call endrun(subname//" :: ERROR: transport array index for "//trim(num_name)//' not found')
        end if
    end do

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(sectional_aerosol_state), intent(inout) :: self

    character(len=*), parameter :: subname = 'destructor'

    call endrun(subname//' is not yet implemented')

    nullify(self%state)
    nullify(self%pbuf)
    if (allocated(self%aero_range_state)) then
        deallocate(self%aero_range_state)
    end if
    if (allocated(self%bin_numconc)) then
        deallocate(self%bin_numconc)
    end if

end subroutine destructor

  !------------------------------------------------------------------------------
  ! sets transported components
  ! This aerosol model with the state of the transported aerosol     ! This updates the transported aerosol constituent array to match the aerosol model state.
! constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine set_transported( self, transported_array )
    class(sectional_aerosol_state), intent(inout) :: self
    real(r8), intent(in) :: transported_array(:,:,:)
    integer              :: irange, ispec, ibin

    character(len=*), parameter :: subname = 'set_transported'

    do irange = 1, self%sec_aero_props%nranges()
        do ispec = 1, self%sec_aero_props%range_nspecies(irange)
            self%aero_range_state(irange)%mass(:,:,ispec) = self%state%q(:self%ncol,:,self%aero_range_state(irange)%transport_ndx(ispec))
        end do
    end do

    do ibin = 1, self%sec_aero_props%nbins()
        self%bin_numconc(:,:,ibin) = self%state%q(:self%ncol,:,self%num_transport_ndx(ibin))
    end do
  end subroutine set_transported

  !------------------------------------------------------------------------------
  ! returns transported components
  ! This returns to current state of the transported aerosol constituents
  ! (mass mixing ratios or number mixing ratios)
  !------------------------------------------------------------------------------
  subroutine get_transported( self, transported_array )
    class(sectional_aerosol_state), intent(in) :: self
    real(r8), intent(out) :: transported_array(:,:,:)
    integer               :: irange, ispec, ibin

    character(len=*), parameter :: subname = 'get_transported'

    do irange = 1, self%sec_aero_props%nranges()
        do ispec = 1, self%sec_aero_props%range_nspecies(irange)
            self%state%q(:self%ncol,:,self%aero_range_state(irange)%transport_ndx(ispec)) = self%aero_range_state(irange)%mass(:,:,ispec)
        end do
    end do

    do ibin = 1, self%sec_aero_props%nbins()
        self%state%q(:self%ncol,:,self%num_transport_ndx(ibin)) = self%bin_numconc(:,:,ibin)
    end do

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
! TODO: what is list_ndx?
    integer, intent(in) :: list_ndx        ! rad climate list number
    integer, intent(in) :: bin_ndx         ! bin number
    real(r8), intent(out) :: kappa(:,:)                 !
    integer, allocatable :: bins2ranges(:)
    integer :: irange

    character(len=*), parameter :: subname = 'hygroscopicity'
    allocate(bins2ranges(self%sec_aero_props%nbins()))
    bins2ranges = self%sec_aero_props%bins2ranges(self%sec_aero_props%nbins())
    irange = bins2ranges(bin_ndx)
    kappa = self%aero_range_state(irange)%hygroscopicity

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

    vol = self%bin_numconc(:, :, bin_ndx) * self%sec_aero_props%particle_volume(bin_ndx)

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

  function bin_dry_density(self, bin_ndx, ncol) result(ddens)
    class(sectional_aerosol_state), intent(in) :: self
    integer, intent(in)   :: bin_ndx
    integer, intent(in)   :: ncol                 ! number of columns
    integer, allocatable  :: bins2ranges(:)
    integer               :: irange
    real(r8)              :: ddens(ncol,pver)

    character(len=*), parameter :: subname = 'dry_density'
    allocate(bins2ranges(self%sec_aero_props%nbins()))
    bins2ranges = self%sec_aero_props%bins2ranges(self%sec_aero_props%nbins())

    irange = bins2ranges(bin_ndx)
    ddens = self%aero_range_state(irange)%dry_density

  end function bin_dry_density

  subroutine update_range(self, mass_tend, irange, range_bounds, ncol)
    class(sectional_aerosol_state), intent(inout) :: self
    real(r8), optional, intent(in) :: mass_tend(:,:,:) ! shape aero_props%range_nspecies (ncol, pver, range_nspecies) -> one array for one range
    integer, intent(in)  :: irange
    integer, intent(in)  :: range_bounds(:,:)
    integer, intent(in) :: ncol                 ! number of columns
    integer              :: ispec, ibin, ispecprop
    real(r8)             :: range_dry_volume(ncol, pver)
    real(r8)             :: range_total_mass(ncol, pver)
    real(r8)             :: test(ncol, pver)

    ! get range bounds
  !  range_bounds = self%sec_aero_props%range_bounds(nranges)
!    bins2ranges = self%sec_aero_props%bins2ranges(nbins)

    range_dry_volume = 0._r8
    range_total_mass = 0._r8

    ! update mass of each component if mass tendency has been passed as an argument
    ! else: update other state variables with mass from before
    if ( present(mass_tend) ) then
        self%aero_range_state(irange)%mass = self%aero_range_state(irange)%mass + mass_tend(:,:,:)
    end if

    range_total_mass = sum(self%aero_range_state(irange)%mass, dim=3) ! sum over species

    ! reset density and hygroscopicity
    self%aero_range_state(irange)%dry_density = 0._r8
    self%aero_range_state(irange)%hygroscopicity = 0._r8

    ! volume of entire mass in a range
    do ibin = range_bounds(irange, 1), range_bounds(irange, 2)
        range_dry_volume = range_dry_volume + self%dry_volume(self%sec_aero_props, 1, ibin, 1, 1) !TODO: change input parameters
    end do

    self%aero_range_state(irange)%dry_density = range_total_mass/range_dry_volume
    do ispec = 1,self%sec_aero_props%range_nspecies(irange)
        ispecprop = self%aero_range_state(irange)%spec_ndx(ispec)
! TODO: source, total hygroscopicity parameter kappa_tot = SUM_OVER_ALL_SPECIES(volume_i/volume_tot * kappa_i)
        self%aero_range_state(irange)%hygroscopicity = self%aero_range_state(irange)%hygroscopicity &
            + self%aero_range_state(irange)%mass(ispec,:,:) / range_dry_volume / &
            self%sec_aero_props%density(ispecprop) * self%sec_aero_props%kappa(ispecprop) ! TODO: probably not the least ugly way to do this
    end do

  end subroutine update_range

end module sectional_aerosol_state_mod
