module sectional_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len

  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog

  implicit none

  private

  public :: sectional_aerosol_properties

  type aerosol_species_properties
     character(len=10)       :: specname
     integer, dimension(2)   :: range_idx
     logical                 :: mixed
  end type aerosol_species_properties

  type, extends(aerosol_properties) :: sectional_aerosol_properties
     private
     integer                                                     :: nranges_ = 0
     integer, dimension(:), allocatable                          :: range_nspecies_
     real(r8), dimension(:), allocatable                         :: bin_centers_ ! radii at bin center (nm)
     real(r8), dimension(:,:), allocatable                       :: bin_bounds_ ! radii at bin bounds (nm)
     integer, dimension(:,:), allocatable                        :: range_bounds_ ! index of bins at range bounds
     integer, dimension(:), allocatable                          :: bins2ranges_ ! range index for each bin
     type(aerosol_species_properties), dimension(:), allocatable :: aer_spec_prop

     integer,  allocatable :: sulfate_mode_ndxs_(:)
     integer,  allocatable :: dust_mode_ndxs_(:)
     integer,  allocatable :: ssalt_mode_ndxs_(:)
     integer,  allocatable :: ammon_mode_ndxs_(:)
     integer,  allocatable :: nitrate_mode_ndxs_(:)
     integer,  allocatable :: msa_mode_ndxs_(:)
     integer,  allocatable :: bcarbon_mode_ndxs_(:,:)
     integer,  allocatable :: porganic_mode_ndxs_(:,:)
     integer,  allocatable :: sorganic_mode_ndxs_(:,:)
     integer :: num_soa_ = 0
     integer :: num_poa_ = 0
     integer :: num_bc_ = 0
   contains
     procedure :: number_transported
     procedure :: get
     procedure :: amcube
     procedure :: actfracs
     procedure :: num_names
     procedure :: mmr_names
     procedure :: amb_num_name
     procedure :: amb_mmr_name
     procedure :: species_type
     procedure :: icenuc_updates_num
     procedure :: icenuc_updates_mmr
     procedure :: apply_number_limits
     procedure :: hetfrz_species
     procedure :: optics_params
     procedure :: nbins_rlist
     procedure :: nspecies_per_bin_rlist
     procedure :: alogsig_rlist
     procedure :: soluble
     procedure :: min_mass_mean_rad
     procedure :: bin_name
     procedure :: scav_diam
     procedure :: resuspension_resize
     procedure :: rebin_bulk_fluxes
     procedure :: hydrophilic

     final :: destructor
  end type sectional_aerosol_properties

  interface sectional_aerosol_properties
     procedure :: constructor
  end interface sectional_aerosol_properties

  logical, parameter :: debug = .false.

contains
  !------------------------------------------------------------------------------
  function constructor(nlfile) result(newobj)

    use mpi,               only: mpi_integer, mpi_real8, mpi_character, mpi_logical, MPI_SUCCESS
    use spmd_utils,        only: mstrid=>masterprocid, mpicom
    use string_utils,      only: int2str
    use namelist_utils,    only: find_group_name
    use infnan, only: nan, assignment(=)


    type(sectional_aerosol_properties), pointer :: newobj

    character(len=*), intent(in) :: nlfile
    integer                      :: nbins=0, nranges=0, ncnst_tot=0, nspecies_tot=0
    ! TODO: ncnst_tot = tracers ??
    integer,allocatable          :: nspecies(:) ! nspecies per bin (given by base_object)
    integer,allocatable          :: range_nspecies(:) ! nspecies per range
    real(r8),allocatable         :: alogsig(:) ! given by base obj
    real(r8),allocatable         :: f1(:) ! given by base obj
    real(r8),allocatable         :: f2(:) ! given by base obj
    real(r8),dimension(:),allocatable  :: bin_centers
    real(r8),dimension(:,:),allocatable:: bin_bounds
    integer,dimension(:,:),allocatable :: range_bounds
    integer,dimension(:),allocatable   :: bins2ranges
    integer                      :: ierr, unitn, ind, pos, ibin, irange
    character(len=50)            :: tmp

    ! namelist variables
    integer                               :: oslo_sectional_nbins
    integer                               :: oslo_sectional_nranges
    integer                               :: oslo_sectional_nspecies_tot
    character(len=10), dimension(500)     :: oslo_sectional_nspecies ! range_nspecies

    character(len=50), dimension(500)     :: oslo_sectional_bin_centers
    character(len=50), dimension(500)     :: oslo_sectional_bin_bounds
    character(len=50), dimension(500)     :: oslo_sectional_range_bounds

    ! namelist aerosol species variables
    type(aerosol_species_properties), dimension(:), allocatable :: oslo_sectional_species_properties
    character(len=10)                     :: oslo_sectional_aerosol_name
    character(len=10)                     :: oslo_sectional_aerosol_range
    logical                               :: oslo_sectional_aerosol_mixed

    character(len=aero_name_len) :: spectype

    character(len=*), parameter :: subname = 'constructor'


    ! Namelists (constructed in bin_config.py)
    namelist /oslo_sectional_properties_nl/ oslo_sectional_nspecies_tot, &
                                            oslo_sectional_nspecies, &
                                            oslo_sectional_nbins, &
                                            oslo_sectional_nranges, &
                                            oslo_sectional_bin_bounds, &
                                            oslo_sectional_bin_centers, &
                                            oslo_sectional_range_bounds

    namelist /oslo_sectional_properties_aerosol_nl/ oslo_sectional_aerosol_name, &
                                            oslo_sectional_aerosol_range, &
                                            oslo_sectional_aerosol_mixed

    ! initialize variables
    oslo_sectional_nspecies_tot = 0
    oslo_sectional_nspecies = ''
    oslo_sectional_nbins = 0
    oslo_sectional_nranges = 0
    oslo_sectional_bin_centers = ''
    oslo_sectional_bin_bounds = ''
    oslo_sectional_range_bounds = ''

    oslo_sectional_aerosol_name = ''
    oslo_sectional_aerosol_range = ''
    oslo_sectional_aerosol_mixed = .false.

    ! read aerosol properties namelist
    if (masterproc) then
        open(newunit=unitn, file=trim(nlfile), status='old')
        call find_group_name(unitn, 'oslo_sectional_properties_nl', ierr)
        if ( ierr /= 0 ) then
            close(unitn)
            call endrun(subname//":: ERROR could not find group 'oslo_sectional_properties_nl'")
        end if

        read(unitn, oslo_sectional_properties_nl, iostat=ierr)
        if ( ierr /= 0 ) then
            close(unitn)
            call endrun(subname // ':: ERROR reading oslo_sectional_properties_nl namelist')
        end if
    end if

    ! allocate array for species objects
    allocate(oslo_sectional_species_properties(oslo_sectional_nspecies_tot), stat=ierr)
    if(ierr/=0) then
        if (masterproc) close(unitn)
        call endrun(subname// ": ERROR "//int2str(ierr)//" allocating oslo_sectional_species_properties")
    end if

    ! read aerosol species properties namelists in sequence
    if (masterproc) then
        do ind=1,oslo_sectional_nspecies_tot
            oslo_sectional_species_properties(ind)%range_idx = 0
            oslo_sectional_species_properties(ind)%specname = ''
            oslo_sectional_species_properties(ind)%mixed = .false.

            call find_group_name(unitn, 'oslo_sectional_properties_aerosol_nl', ierr)
            if (ierr /= 0) then
                close(unitn)
                call endrun(subname//":: ERROR could not find group 'oslo_sectional_properties_aerosol_nl'")
            end if

            read(unitn, oslo_sectional_properties_aerosol_nl, iostat=ierr)
            if (ierr /= 0) then
                close(unitn)
                call endrun(subname // ':: ERROR reading oslo_sectional_properties_aerosol_nl')
            end if

            pos = index(oslo_sectional_aerosol_range, ':')
            if ( pos == 0 ) then
                close(unitn)
                call endrun(subname//":: ERROR invalid format for 'oslo_sectional_aerosol_range'")
            end if

            read(oslo_sectional_aerosol_range(1:pos-1), *) oslo_sectional_species_properties(ind)%range_idx(1)
            read(oslo_sectional_aerosol_range(pos+1:), *) oslo_sectional_species_properties(ind)%range_idx(2)

            oslo_sectional_species_properties(ind)%specname = oslo_sectional_aerosol_name
            oslo_sectional_species_properties(ind)%mixed = oslo_sectional_aerosol_mixed

            oslo_sectional_aerosol_name = ''
            oslo_sectional_aerosol_range = ''
            oslo_sectional_aerosol_mixed = .false.
        end do
        close(unitn)
    end if

    nbins = oslo_sectional_nbins
    nranges = oslo_sectional_nranges
    nspecies_tot = oslo_sectional_nspecies_tot

    ! allocate with size nbins
    allocate( nspecies(nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating nspecies")
            return
        end if

    allocate( bins2ranges(nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating bins2ranges'")
            return
        end if

    allocate( bin_centers(nbins), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating bin_centers'")
           return
       end if
    allocate( bin_bounds(nbins, 2), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating bin_bounds'")
           return
       end if

    ! allocate with size nranges
    allocate( range_nspecies(nranges), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating nspecies'")
           return
       end if

    allocate( range_bounds(nranges, 2), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating range_bounds'")
           return
       end if

    ! TODO: find out what to do with these three
    allocate( alogsig(nbins), stat=ierr)
        if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating alogsig'")
           return
       end if
    allocate( f1(nbins), stat=ierr)
           if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating f1'")
           return
       end if
    allocate( f2(nbins), stat=ierr)
           if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating f2'")
           return
       end if
    alogsig = nan
    f1 = nan
    f2 = nan

    ! TODO: check indexer_ var in aerosol_properties mod. same as the indices from chemical pp?
    do ibin=1,nbins
        read(oslo_sectional_bin_centers(ibin),'(D)') bin_centers(ibin)
        tmp = oslo_sectional_bin_bounds(ibin)
        pos = index(tmp, ':')
        read(tmp(1:pos-1), '(D)') bin_bounds(ibin,1)
        read(tmp(pos+1:), '(D)') bin_bounds(ibin,2)
    end do

    ! parse range bounds
    do irange=1,nranges
        tmp = oslo_sectional_range_bounds(irange)
        pos = index(tmp, ':')
        read(tmp(1:pos-1), *) range_bounds(irange,1)
        read(tmp(pos+1:), *) range_bounds(irange,2)
    end do
    read(oslo_sectional_nspecies, *) range_nspecies

    ncnst_tot = nbins + sum(range_nspecies) ! nbins + sum(nspecies)

    ! initialize bins2ranges array

    do ibin=1,nbins
        do irange=1,nranges
            if (ibin >= range_bounds(irange,1) .and. ibin <= range_bounds(irange,2)) then
                bins2ranges(ibin) = irange
                nspecies(ibin) = range_nspecies(irange)
            end if
        end do
    end do


       ! Report
    if (masterproc) then
        write(iulog ,*) 'sectional aerosol properties namelist: '
        write(iulog ,*) 'nspecies_tot = ', oslo_sectional_nspecies_tot
        write(iulog ,*) 'nbins = ', oslo_sectional_nbins
        write(iulog ,*) 'nranges = ', oslo_sectional_nranges
        write(iulog ,*) 'range_nspecies = '
        do ind = 1, nranges, 5
            write(iulog,*) oslo_sectional_nspecies(ind:min(ind+4, oslo_sectional_nranges))
        end do
        write(iulog ,*) 'bin_centers: '
        do ind = 1, oslo_sectional_nbins, 5
            write(iulog, *) oslo_sectional_bin_centers(ind:min(ind+4, oslo_sectional_nbins))
        end do
        write(iulog ,*) 'bin_bounds: '
        do ind = 1, oslo_sectional_nbins, 5
            write(iulog, *) oslo_sectional_bin_bounds(ind:min(ind+4, oslo_sectional_nbins))
        end do
        write(iulog ,*) 'range_bounds: '
        do ind = 1, oslo_sectional_nranges, 5
            write(iulog, *) oslo_sectional_range_bounds(ind:min(ind+4, oslo_sectional_nranges))
        end do

        do ind = 1,oslo_sectional_nspecies_tot
            write(iulog ,*) 'sectional aerosol species properties namelist: '
            write(iulog ,*) 'species name = ', oslo_sectional_species_properties(ind)%specname
            write(iulog ,*) 'range bound idces = ', oslo_sectional_species_properties(ind)%range_idx
            write(iulog ,*) 'mixed = ', oslo_sectional_species_properties(ind)%mixed
            write(iulog ,*) 'bins2ranges = ', bins2ranges
       end do
    end if

    ! allocate the aerosol_properties object
    allocate(newobj,stat=ierr)
    if ( ierr/=0 ) then
        nullify(newobj)
        return
    end if


    allocate(newobj%bins2ranges_(nbins), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%range_nspecies_(nranges), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%bin_centers_(nbins), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%bin_bounds_(nbins, 2), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%range_bounds_(nranges, 2), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%aer_spec_prop(nspecies_tot), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    if (masterproc) then

        newobj%bins2ranges_ = bins2ranges
        newobj%nranges_ = nranges
        newobj%range_nspecies_ = range_nspecies
        newobj%bin_centers_ = bin_centers(:nbins)
        newobj%bin_bounds_ = bin_bounds(:nbins, :)
        newobj%range_bounds_ = range_bounds(:nranges, :)
        newobj%aer_spec_prop = oslo_sectional_species_properties(:nspecies_tot)

    end if

    ! Broadcast the namelist properties
    call MPI_Bcast(newobj%bins2ranges_, nbins, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) then
        nullify(newobj)
        if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
        call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'bins2ranges'")
    end if

    call MPI_Bcast(newobj%nranges_, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) then
        nullify(newobj)
        if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
        call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'nranges'")
    end if

    call MPI_Bcast(newobj%range_nspecies_, size(newobj%range_nspecies_), mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) then
        nullify(newobj)
        if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
        call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'range_nspecies_'")
    end if

    call MPI_Bcast(newobj%bin_centers_, nbins, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) then
        nullify(newobj)
        if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
        call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'bin_centers'")
    end if

    call MPI_Bcast(newobj%bin_bounds_, size(newobj%bin_bounds_), mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) then
        nullify(newobj)
        if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
        call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'bin_bounds'")
    end if


    call MPI_Bcast(newobj%range_bounds_, size(newobj%range_bounds_), mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= MPI_SUCCESS) then
        nullify(newobj)
        if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
        call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'range_bounds'")
    end if

    ! Broadcast the species properties
    do ind=1,nspecies_tot

        call MPI_Bcast(newobj%aer_spec_prop(ind)%specname, 1, mpi_character, mstrid, mpicom, ierr)
        if (ierr /= MPI_SUCCESS) then
            nullify(newobj)
            if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
            call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'aer_spec_prop specname'")
        end if

        call MPI_Bcast(newobj%aer_spec_prop(ind)%range_idx, 2, mpi_integer, mstrid, mpicom, ierr)
        if (ierr /= MPI_SUCCESS) then
            nullify(newobj)
            if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
            call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'aer_spec_prop range_idx'")
        end if

        ! TODO: fix this broadcasting somehow
        call MPI_Bcast(newobj%aer_spec_prop(ind)%mixed, 1, mpi_logical, mstrid, mpicom, ierr)
        if (ierr/= MPI_SUCCESS) then
            nullify(newobj)
            if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)
            call endrun(subname// ": Error "//int2str(ierr)//" broadcasting 'aer_spec_prop mixed'")
        end if

    end do

    call newobj%initialize(nbins, ncnst_tot, nspecies, nspecies, alogsig, f1, f2, ierr)

    ! deallocate local variables
    if (allocated(nspecies)) deallocate(nspecies)
    if (allocated(range_nspecies)) deallocate(range_nspecies)
    if (allocated(alogsig)) deallocate(alogsig)
    if (allocated(f1)) deallocate(f1)
    if (allocated(f2)) deallocate(f2)
    if (allocated(bin_centers)) deallocate(bin_centers)
    if (allocated(bin_bounds)) deallocate(bin_bounds)
    if (allocated(range_bounds)) deallocate(range_bounds)
    if (allocated(bins2ranges)) deallocate(bins2ranges)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(sectional_aerosol_properties), intent(inout) :: self

    character(len=*), parameter :: subname = 'destructor'

    if (allocated(self%range_nspecies_)) then
        deallocate(self%range_nspecies_)
    end if

    if (allocated(self%bin_centers_)) then
        deallocate(self%bin_centers_)
    end if

    if (allocated(self%bin_bounds_)) then
        deallocate(self%bin_bounds_)
    end if

    if (allocated(self%range_bounds_)) then
        deallocate(self%range_bounds_)
    end if

    if (allocated(self%bins2ranges_)) then
        deallocate(self%bins2ranges_)
    end if

    if (allocated(self%aer_spec_prop)) then
        deallocate(self%aer_spec_prop)
    end if
    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! returns number of transported aerosol constituents
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(sectional_aerosol_properties), intent(in) :: self
    character(len=*), parameter :: subname = 'number_transported'

    call endrun(subname//' is not yet implemented')

  end function number_transported

  !------------------------------------------------------------------------
  ! returns aerosol properties:
  !  density
  !  hygroscopicity
  !  species type
  !  species name
  !  short wave species refractive indices
  !  long wave species refractive indices
  !  species morphology
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, list_ndx, density, hygro, &
                 spectype, specname, specmorph, refindex_sw, refindex_lw)

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    integer, optional, intent(in) :: list_ndx  ! climate or a diagnostic list number
    real(r8), optional, intent(out) :: density ! density (kg/m3)
    real(r8), optional, intent(out) :: hygro   ! hygroscopicity
    character(len=*), optional, intent(out) :: spectype  ! species type
    character(len=*), optional, intent(out) :: specname  ! species name
    character(len=*), optional, intent(out) :: specmorph ! species morphology
    complex(r8), pointer, optional, intent(out) :: refindex_sw(:) ! short wave species refractive indices
    complex(r8), pointer, optional, intent(out) :: refindex_lw(:) ! long wave species refractive indices

    integer :: ilist
    character(len=*), parameter :: subname = 'get'

    call endrun(subname//' is not yet implemented')

  end subroutine get

  !------------------------------------------------------------------------
  ! returns optics type and table parameters
  !------------------------------------------------------------------------
  subroutine optics_params(self, list_ndx, bin_ndx, opticstype, extpsw, abspsw, asmpsw, absplw, &
       refrtabsw, refitabsw, refrtablw, refitablw, ncoef, prefr, prefi, sw_hygro_ext_wtp, &
       sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, wgtpct, nwtp, &
       sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, &
       corefrac, bcdust, kap, relh, nfrac, nbcdust, nkap, nrelh )

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: list_ndx            ! rad climate/diags list

    character(len=*), optional, intent(out) :: opticstype

    ! refactive index table parameters
    real(r8),  optional, pointer     :: extpsw(:,:,:,:) ! short wave specific extinction
    real(r8),  optional, pointer     :: abspsw(:,:,:,:) ! short wave specific absorption
    real(r8),  optional, pointer     :: asmpsw(:,:,:,:) ! short wave asymmetry factor
    real(r8),  optional, pointer     :: absplw(:,:,:,:) ! long wave specific absorption
    real(r8),  optional, pointer     :: refrtabsw(:,:)  ! table of short wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitabsw(:,:)  ! table of short wave imaginary refractive indices for aerosols
    real(r8),  optional, pointer     :: refrtablw(:,:)  ! table of long wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitablw(:,:)  ! table of long wave imaginary refractive indices for aerosols
    integer,   optional, intent(out) :: ncoef  ! number of chebychev polynomials
    integer,   optional, intent(out) :: prefr  ! number of real refractive indices in table
    integer,   optional, intent(out) :: prefi  ! number of imaginary refractive indices in table

    ! hygrowghtpct table parameters
    real(r8),  optional, pointer     :: sw_hygro_ext_wtp(:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_ssa_wtp(:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_asm_wtp(:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_ext_wtp(:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: wgtpct(:)   ! weight precent of H2SO4/H2O solution
    integer,   optional, intent(out) :: nwtp        ! number of weight precent values

    ! hygrocoreshell table parameters
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ext(:,:,:,:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ssa(:,:,:,:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_asm(:,:,:,:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_coreshell_ext(:,:,:,:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: corefrac(:) ! core fraction dimension values
    real(r8),  optional, pointer     :: bcdust(:)   ! bc/(bc + dust) fraction dimension values
    real(r8),  optional, pointer     :: kap(:)      ! hygroscopicity dimension values
    real(r8),  optional, pointer     :: relh(:)     ! relative humidity dimension values
    integer,   optional, intent(out) :: nfrac       ! core fraction dimension size
    integer,   optional, intent(out) :: nbcdust     ! bc/(bc + dust) fraction dimension size
    integer,   optional, intent(out) :: nkap        ! hygroscopicity dimension size
    integer,   optional, intent(out) :: nrelh       ! relative humidity dimension size

    character(len=*), parameter :: subname = 'optics_params'

    call endrun(subname//' is not yet implemented')

  end subroutine optics_params

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    character(len=*), parameter :: subname = 'amcube'

    amcube = -1.0_r8
    ! TODO: do we need this? cannot call endrun, due to "pure elemental"


  end function amcube

  !------------------------------------------------------------------------------
  ! returns mass and number activation fractions
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    use shr_spfn_mod, only: erf => shr_spfn_erf
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin index
    real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
    real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
    real(r8),intent(out) :: fn       ! activation fraction for aerosol number
    real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

    character(len=*), parameter :: subname = 'actfracs'

    call endrun(subname//' is not yet implemented')

  end subroutine actfracs

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol number dens
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens

    character(len=*), parameter :: subname = 'num_names'

    call endrun(subname//' is not yet implemented')

  end subroutine num_names

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol MMR
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR

    character(len=*), parameter :: subname = 'mmr_names'

    call endrun(subname//' is not yet implemented')

  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol number dens
    character(len=*), parameter :: subname = 'amb_num_name'

    call endrun(subname//' is not yet implemented')

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR
    character(len=*), parameter :: subname = 'amb_mmr_name'

    call endrun(subname//' is not yet implemented')

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: spectype ! species type
    character(len=*), parameter :: subname = 'species_type'

    call endrun(subname//' is not yet implemented')

  end subroutine species_type

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
  !------------------------------------------------------------------------------
  function icenuc_updates_num(self, bin_ndx) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    logical :: res

    character(len=aero_name_len) :: spectype
    character(len=aero_name_len) :: modetype
    integer :: spc_ndx
    character(len=*), parameter :: subname = 'icenuc_updates_num'

    res = .false.

    call endrun(subname//' is not yet implemented')

  end function icenuc_updates_num

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
  !------------------------------------------------------------------------------
  function icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    logical :: res

    character(len=32) :: spectype
    character(len=32) :: modetype
    character(len=*), parameter :: subname = 'icenuc_updates_mmr'

    res = .false.

    call endrun(subname//' is not yet implemented')

  end function icenuc_updates_mmr

  !------------------------------------------------------------------------------
  ! apply max / min to number concentration
  !------------------------------------------------------------------------------
  subroutine apply_number_limits( self, naerosol, vaerosol, istart, istop, m )
    class(sectional_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
    integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop           ! stop column index
    integer,  intent(in) :: m               ! mode or bin index
    character(len=*), parameter :: subname = 'apply_number_limits'

    call endrun(subname//' is not yet implemented')

  end subroutine apply_number_limits

  !------------------------------------------------------------------------------
  ! returns TRUE if species `spc_ndx` in aerosol subset `bin_ndx` contributes to
  ! the particles' ability to act as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_species(self, bin_ndx, spc_ndx) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    integer, intent(in) :: spc_ndx  ! species number

    logical :: res
    character(len=*), parameter :: subname = 'hetfrz_species'

    call endrun(subname//' is not yet implemented')

  end function hetfrz_species

  !------------------------------------------------------------------------------
  ! returns TRUE if soluble
  !------------------------------------------------------------------------------
  logical function soluble(self,bin_ndx)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), parameter :: subname = 'soluble'

    call endrun(subname//' is not yet implemented')

  end function soluble

  !------------------------------------------------------------------------------
  ! returns minimum mass mean radius (meters)
  !------------------------------------------------------------------------------
  function min_mass_mean_rad(self,bin_ndx,species_ndx) result(minrad)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    real(r8) :: minrad  ! meters

    integer :: nbins
    character(len=*), parameter :: subname = 'min_mass_mean_rad'

    call endrun(subname//' is not yet implemented')

  end function min_mass_mean_rad

  !------------------------------------------------------------------------------
  ! returns the total number of bins for a given radiation list index
  !------------------------------------------------------------------------------
  function nbins_rlist(self, list_ndx)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx  ! radiation list number

    integer :: res
    character(len=*), parameter :: subname = 'nbins_rlist'

    call endrun(subname//' is not yet implemented')

  end function nbins_rlist

  !------------------------------------------------------------------------------
  ! returns number of species in a bin for a given radiation list index
  !------------------------------------------------------------------------------
  function nspecies_per_bin_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    integer :: res
    character(len=*), parameter :: subname = 'nspecies_per_bin_rlist'

    call endrun(subname//' is not yet implemented')

  end function nspecies_per_bin_rlist

  !------------------------------------------------------------------------------
  ! returns the natural log of geometric standard deviation of the number
  ! distribution for radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function alogsig_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: res
    character(len=*), parameter :: subname = 'alogsig_rlist'

    call endrun(subname//' is not yet implemented')

  end function alogsig_rlist

  !------------------------------------------------------------------------------
  ! returns name for a given radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function bin_name(self, list_ndx,  bin_ndx) result(name)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    character(len=32) name
    character(len=*), parameter :: subname = 'bin_name'

    call endrun(subname//' is not yet implemented')

  end function bin_name

  !------------------------------------------------------------------------------
  ! returns scavenging diameter (cm) for a given aerosol bin number
  !------------------------------------------------------------------------------
  function scav_diam(self, bin_ndx) result(diam)
    use modal_aero_data, only: dgnum_amode

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: diam
    character(len=*), parameter :: subname = 'scav_diam'

    call endrun(subname//' is not yet implemented')

  end function scav_diam

  !------------------------------------------------------------------------------
  ! adjust aerosol concentration tendencies to create larger sizes of aerosols
  ! during resuspension
  !------------------------------------------------------------------------------
  subroutine resuspension_resize(self, dcondt)

    use modal_aero_data, only:  mode_size_order

    class(sectional_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: dcondt(:)
    character(len=*), parameter :: subname = 'resuspension_resize'

    call endrun(subname//' is not yet implemented')

  end subroutine resuspension_resize

  !------------------------------------------------------------------------------
  ! returns bulk deposition fluxes of the specified species type
  ! rebinned to specified diameter limits
  !------------------------------------------------------------------------------
  subroutine rebin_bulk_fluxes(self, bulk_type, dep_fluxes, diam_edges, bulk_fluxes, &
                               error_code, error_string)
    use infnan, only: nan, assignment(=)

    class(sectional_aerosol_properties), intent(in) :: self
    character(len=*),intent(in) :: bulk_type       ! aerosol type to rebin
    real(r8), intent(in) :: dep_fluxes(:)          ! kg/m2
    real(r8), intent(in) :: diam_edges(:)          ! meters
    real(r8), intent(out) :: bulk_fluxes(:)        ! kg/m2
    integer,  intent(out) :: error_code            ! error code (0 if no error)
    character(len=*), intent(out) :: error_string  ! error string
    character(len=*), parameter :: subname = 'rebin_bulk_fluxes'

    call endrun(subname//' is not yet implemented')

  end subroutine rebin_bulk_fluxes

  !------------------------------------------------------------------------------
  ! Returns TRUE if bin is hydrophilic, otherwise FALSE
  !------------------------------------------------------------------------------
  logical function hydrophilic(self, bin_ndx)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx ! bin number

    character(len=aero_name_len) :: modetype
    character(len=*), parameter :: subname = 'hydrophilic'

    call endrun(subname//' is not yet implemented')

  end function hydrophilic

end module sectional_aerosol_properties_mod
