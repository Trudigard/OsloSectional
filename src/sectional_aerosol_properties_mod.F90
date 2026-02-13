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
     character(len=:), allocatable :: specname        ! e.g. DU
     integer                 :: nbin            ! number of bins containing species
     integer                 :: nrange          ! nr of ranges containing species
     integer, allocatable    :: range_ndx(:)    ! indices of ranges containing species
     integer, allocatable    :: bin_ndx(:)      ! bin indices containing species
     real(r8)                :: density
     real(r8)                :: molecular_weight ! kg/kmol
     real(r8)                :: kappa
     logical                 :: mixed           ! true if internally mixed
     !TODO: add hydrophilic/phobic
     ! integer               :: idx             ! indices of species tracers
     character(len=10), allocatable :: tracernames(:) ! e.g. DU_R3
     end type aerosol_species_properties

  type, extends(aerosol_properties) :: sectional_aerosol_properties
     private
     integer               :: nranges_ = 0
     integer               :: nspecies_tot_ = 0
     integer, allocatable  :: range_nspecies_(:)
     real(r8), allocatable :: bin_centers_(:) ! radii at bin center (nm)
     real(r8), allocatable :: bin_bounds_(:,:)! radii at bin bounds (nm)
     integer, allocatable  :: range_bounds_(:,:) ! index of bins at range bounds
     integer, allocatable  :: bins2ranges_(:) ! range index for each bin
     real(r8), allocatable :: particle_volume_(:) ! volume of a single particle in a bin (center radius)
     type(aerosol_species_properties), allocatable :: aer_spec_prop(:)

   contains
     procedure :: number_transported
     procedure :: get
     procedure :: amcube
     procedure :: density
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
     procedure :: nspecies_tot
     procedure :: range_nspecies
     procedure :: nranges
     procedure :: bin_centers
     procedure :: bin_bounds
     procedure :: particle_volume
     procedure :: spec_range_ndx
     procedure :: spec_bin_ndx
     procedure :: spec_nrange
     procedure :: spec_nbin
     procedure :: spec_tracernames
     procedure :: bins2ranges
     procedure :: nspecies_per_bin_rlist
     procedure :: alogsig_rlist
     procedure :: soluble
     procedure :: min_mass_mean_rad
     procedure :: bin_name
     procedure :: scav_diam
     procedure :: resuspension_resize
     procedure :: rebin_bulk_fluxes
     procedure :: hydrophilic
     procedure :: range_bounds
     procedure :: kappa
     procedure :: molecular_weight

     final :: destructor
  end type sectional_aerosol_properties

  interface sectional_aerosol_properties
     procedure :: constructor
  end interface sectional_aerosol_properties

  logical, parameter :: debug = .false.
  type(sectional_aerosol_properties), pointer :: prop_obj => null()


contains
  !------------------------------------------------------------------------------
  function constructor(nlfile) result(newobj)

    use mpi,               only: mpi_integer, mpi_real8, mpi_character, mpi_logical, MPI_SUCCESS
    use spmd_utils,        only: mstrid=>masterprocid, mpicom
    use string_utils,      only: int2str
    use namelist_utils,    only: find_group_name
    use infnan,            only: nan, assignment(=)

    type(sectional_aerosol_properties), pointer :: newobj

    character(len=*),optional, intent(in) :: nlfile
    integer                      :: ncnst_tot=0
    ! TODO: ncnst_tot = tracers ??
    integer,allocatable          :: nspecies(:) ! nspecies per bin (given by base_object)
    real(r8),allocatable         :: alogsig(:) ! given by base obj
    real(r8),allocatable         :: f1(:) ! given by base obj: Abdul-Razzak 1998 eq 28
    real(r8),allocatable         :: f2(:) ! given by base obj: Abdul-Razzak 1998 eq 29
    real(r8),allocatable         :: bin_centers(:)
    real(r8),allocatable         :: bin_bounds(:,:)
    integer,allocatable          :: range_bounds(:,:)
    integer,allocatable          :: bins2ranges(:)
    integer                      :: bin_ndx(100)
    integer                      :: ierr, unitn, pos, lower, upper
    integer                      :: ind, ibin, irange, ispec

    ! namelist variables
    integer               :: oslo_sectional_nbins
    integer               :: oslo_sectional_nranges
    integer               :: oslo_sectional_nspecies_tot
    integer               :: oslo_sectional_nspecies(500) ! range_nspecies TODO: make allocatable!!
    integer, parameter    :: strlen=50

    character(len=strlen) :: oslo_sectional_bin_centers(500)
    character(len=strlen) :: oslo_sectional_bin_bounds(500)
    character(len=strlen) :: oslo_sectional_range_bounds(500)

    ! namelist aerosol species variables
    type(aerosol_species_properties), allocatable :: oslo_sectional_species_properties(:)
    character(len=10)                     :: oslo_sectional_aerosol_name
    character(len=10)                     :: oslo_sectional_aerosol_range
    real(r8)                              :: oslo_sectional_aerosol_density
    real(r8)                              :: oslo_sectional_aerosol_weight
    real(r8)                              :: oslo_sectional_aerosol_kappa
    logical                               :: oslo_sectional_aerosol_mixed

    character(len=aero_name_len) :: spectype

    character(len=*), parameter :: subname = 'constructor'

!==================================================================================================
! Read Namelists
!==================================================================================================

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
                                            oslo_sectional_aerosol_density, &
                                            oslo_sectional_aerosol_weight, &
                                            oslo_sectional_aerosol_mixed, &
                                            oslo_sectional_aerosol_kappa
    if ( present(nlfile) ) then

    if ( associated(prop_obj) ) then
        call endrun(subname//':: ERROR sectional_aerosol_properties has already been initialized')
    end if
    ! initialize variables
    oslo_sectional_nspecies_tot = 0
    oslo_sectional_nspecies = 0
    oslo_sectional_nbins = 0
    oslo_sectional_nranges = 0
    oslo_sectional_bin_centers = ''
    oslo_sectional_bin_bounds = ''
    oslo_sectional_range_bounds = ''

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

!==================================================================================================
! Broadcast oslo_sectional_properties
!==================================================================================================
    call MPI_Bcast(oslo_sectional_nspecies_tot, 1, mpi_integer, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_nspecies_tot'")
    end if

    call MPI_Bcast(oslo_sectional_nspecies, size(oslo_sectional_nspecies), mpi_integer, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_nspecies'")
    end if

    call MPI_Bcast(oslo_sectional_nbins, 1, mpi_integer, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_nbins'")
    end if

    call MPI_Bcast(oslo_sectional_nranges, 1, mpi_integer, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_nranges'")
    end if

    call MPI_Bcast(oslo_sectional_bin_bounds, strlen*size(oslo_sectional_bin_bounds), mpi_character, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_bin_bounds'")
    end if

    call MPI_Bcast(oslo_sectional_bin_centers, strlen*size(oslo_sectional_bin_centers), mpi_character, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_bin_centers")
    end if

    call MPI_Bcast(oslo_sectional_range_bounds, strlen*size(oslo_sectional_range_bounds), mpi_character, mstrid, mpicom, ierr)
    if ( ierr /= MPI_SUCCESS ) then
        call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_range_bounds")
    end if

!==================================================================================================
! Read and broadcast aerosol species properties namelist in sequence
!==================================================================================================
 !    allocate array for species objects
    allocate(oslo_sectional_species_properties(oslo_sectional_nspecies_tot), stat=ierr)
    if(ierr/=0) then
        if (masterproc) close(unitn)
        call endrun(subname// ": ERROR "//int2str(ierr)//" allocating oslo_sectional_species_properties")
    end if

    do ispec=1,oslo_sectional_nspecies_tot

        ! namelist variables
        oslo_sectional_aerosol_name = ''
        oslo_sectional_aerosol_range = ''
        oslo_sectional_aerosol_density = 0.0_r8
        oslo_sectional_aerosol_weight = 0.0_r8
        oslo_sectional_aerosol_kappa = 0.0_r8
        oslo_sectional_aerosol_mixed = .false.
        lower = 0
        upper = 0

        ! read namelist
        if (masterproc) then
            call find_group_name(unitn, 'oslo_sectional_properties_aerosol_nl', ierr)
            if (ierr /= 0) then
                close(unitn)
                call endrun(subname//":: ERROR could not find group 'oslo_sectional_properties_aerosol_nl' for species number " &
                //int2str(ispec)//'of '//int2str(oslo_sectional_nspecies_tot) )
            end if

            read(unitn, oslo_sectional_properties_aerosol_nl, iostat=ierr)
            if (ierr /= 0) then
                close(unitn)
                call endrun(subname // ":: ERROR reading 'oslo_sectional_properties_aerosol_nl' for species number " &
                //int2str(ispec)//'of '//int2str(oslo_sectional_nspecies_tot) )
            end if

        end if

        ! broadcast
        call MPI_Bcast(oslo_sectional_aerosol_name, len(oslo_sectional_aerosol_name), mpi_character, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_aerosol_name'")
        end if

        call MPI_Bcast(oslo_sectional_aerosol_range, len(oslo_sectional_aerosol_range), mpi_character, mstrid, mpicom, ierr)
        if ( ierr /= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_aerosol_range'")
        end if

        call MPI_Bcast(oslo_sectional_aerosol_density, 1, mpi_real8, mstrid, mpicom, ierr)
        if ( ierr/= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_aerosol_density'")
        end if

        call MPI_Bcast(oslo_sectional_aerosol_weight, 1, mpi_real8, mstrid, mpicom, ierr)
        if ( ierr/= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_aerosol_weight'")
        end if

        call MPI_Bcast(oslo_sectional_aerosol_kappa, 1, mpi_real8, mstrid, mpicom, ierr)
        if ( ierr/= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_aerosol_kappa'")
        end if

        call MPI_Bcast(oslo_sectional_aerosol_mixed, 1, mpi_logical, mstrid, mpicom, ierr)
        if ( ierr/= MPI_SUCCESS ) then
            call endrun(subname//": Error "//int2str(ierr)//" broadcasting 'oslo_sectional_aerosol_mixed'")
        end if

        ! initialize properties object variables
        oslo_sectional_species_properties(ispec)%specname = ''
        oslo_sectional_species_properties(ispec)%density = 0.0_r8
        oslo_sectional_species_properties(ispec)%molecular_weight = 0.0_r8
        oslo_sectional_species_properties(ispec)%kappa = 0.0_r8
        oslo_sectional_species_properties(ispec)%mixed = .false.
        oslo_sectional_species_properties(ispec)%nbin = 0
        oslo_sectional_species_properties(ispec)%nrange = 0

        pos = index(oslo_sectional_aerosol_range, ':')

        if ( pos == 0 ) then
           call endrun(subname//":: ERROR invalid format for 'oslo_sectional_aerosol_range'")
        end if

        read(oslo_sectional_aerosol_range(1:pos-1), *) lower
        read(oslo_sectional_aerosol_range(pos+1:), *) upper

        oslo_sectional_species_properties(ispec)%nrange = upper - lower + 1
        oslo_sectional_species_properties(ispec)%specname = oslo_sectional_aerosol_name
        oslo_sectional_species_properties(ispec)%density = oslo_sectional_aerosol_density
        oslo_sectional_species_properties(ispec)%molecular_weight = oslo_sectional_aerosol_weight
        oslo_sectional_species_properties(ispec)%kappa = oslo_sectional_aerosol_kappa
        oslo_sectional_species_properties(ispec)%mixed = oslo_sectional_aerosol_mixed

        ! allocate and initialize range_ndx, and tracernames
        allocate(oslo_sectional_species_properties(ispec)%range_ndx(oslo_sectional_species_properties(ispec)%nrange), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating oslo_sectional_species_properties range_ndx")
        end if
        oslo_sectional_species_properties(ispec)%range_ndx = 0

        do ind = 1, (upper - lower + 1)
            oslo_sectional_species_properties(ispec)%range_ndx(ind) = lower + ind - 1
        end do


    end do

    if (masterproc) then
        close(unitn)
    end if

!==================================================================================================
! Allocate local arrays
!==================================================================================================

    ! allocate with size nbins
    allocate( nspecies(oslo_sectional_nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating nspecies")
            return
        end if

    allocate( bins2ranges(oslo_sectional_nbins), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": Error "//int2str(ierr)//" allocating bins2ranges'")
            return
        end if

    allocate( bin_centers(oslo_sectional_nbins), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating bin_centers'")
           return
       end if
    allocate( bin_bounds(oslo_sectional_nbins, 2), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating bin_bounds'")
           return
       end if

    allocate( range_bounds(oslo_sectional_nranges, 2), stat=ierr)
       if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating range_bounds'")
           return
       end if

    ! TODO: find out what to do with these three
    allocate( alogsig(oslo_sectional_nbins), stat=ierr)
        if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating alogsig'")
           return
       end if
    allocate( f1(oslo_sectional_nbins), stat=ierr)
           if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating f1'")
           return
       end if
    allocate( f2(oslo_sectional_nbins), stat=ierr)
           if( ierr /= 0 ) then
           call endrun(subname// ": Error "//int2str(ierr)//" allocating f2'")
           return
       end if
!==================================================================================================
! parse bin info
!==================================================================================================
    ! TODO: figure out what to do with these
    ! alogsig -> alogsig(m) = log(sigmag(m))
    ! f1 -> f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
    ! f2 -> f2(m) = 1._r8 + 0.25_r8*alogsig(m)
    alogsig = nan
    f1 = nan
    f2 = nan

    ! TODO: check indexer_ var in aerosol_properties mod. same as the indices from chemical pp?
    do ibin=1,oslo_sectional_nbins
        read(oslo_sectional_bin_centers(ibin),'(F10.5)') bin_centers(ibin)
        pos = index(oslo_sectional_bin_bounds(ibin), ':')
        read(oslo_sectional_bin_bounds(ibin)(1:pos-1), '(F10.5)') bin_bounds(ibin,1)
        read(oslo_sectional_bin_bounds(ibin)(pos+1:), '(F10.5)') bin_bounds(ibin,2)
    end do

    ! parse range bounds
    do irange=1,oslo_sectional_nranges
        pos = index(oslo_sectional_range_bounds(irange), ':')
        read(oslo_sectional_range_bounds(irange)(1:pos-1), '(I3)') range_bounds(irange,1)
        read(oslo_sectional_range_bounds(irange)(pos+1:), '(I3)') range_bounds(irange,2)
    end do

    ncnst_tot = oslo_sectional_nbins + sum(oslo_sectional_nspecies(:oslo_sectional_nranges)) ! nbins + sum(nspecies)

    ! initialize bins2ranges array
    do ibin=1,oslo_sectional_nbins
        do irange=1,oslo_sectional_nranges
            if (ibin >= range_bounds(irange,1) .and. ibin <= range_bounds(irange,2)) then
                bins2ranges(ibin) = irange
                nspecies(ibin) = oslo_sectional_nspecies(irange)
                exit
            end if
        end do
    end do

    ! fill species objects with bin info
    do ispec = 1, oslo_sectional_nspecies_tot
        ind = 0
        oslo_sectional_species_properties(ispec)%nbin = 0
        do ibin = 1, oslo_sectional_nbins
            if (bins2ranges(ibin) >= oslo_sectional_species_properties(ispec)%range_ndx(1) .and. &
                bins2ranges(ibin) <= maxval(oslo_sectional_species_properties(ispec)%range_ndx)) then
                ! nr of bins containing each species (e.g. dust_nbin in dust_model.F90)
                oslo_sectional_species_properties(ispec)%nbin = oslo_sectional_species_properties(ispec)%nbin + 1
                ind = ind+1
                ! indices of bins containing each species
                bin_ndx(ind) = ibin
            end if
        end do
! TODO: I moved the 6 lines below here, after species%nbin is initialized to make the bin_ndx allocatable but feels a bit messy..
        allocate(oslo_sectional_species_properties(ispec)%bin_ndx(oslo_sectional_species_properties(ispec)%nbin), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating oslo_sectional_species_properties bin_ndx")
        end if
!        oslo_sectional_species_properties(ispec)%bin_ndx = 0 ! TODO: does this need to be initialized?
        oslo_sectional_species_properties(ispec)%bin_ndx = bin_ndx(:oslo_sectional_species_properties(ispec)%nbin)

        allocate(oslo_sectional_species_properties(ispec)%tracernames(oslo_sectional_species_properties(ispec)%nrange), stat=ierr)
        if(ierr/=0) then
            call endrun(subname// ": ERROR "//int2str(ierr)//" allocating oslo_sectional_species_properties tracernames")
        end if
        oslo_sectional_species_properties(ispec)%tracernames = ''

        do irange = 1, oslo_sectional_species_properties(ispec)%nrange
            ! tracernames for each species - potentially check if consistent with cnst_get_ind?
            oslo_sectional_species_properties(ispec)%tracernames(irange) = &
            trim(oslo_sectional_species_properties(ispec)%specname)//'_R'//&
            trim(int2str(oslo_sectional_species_properties(ispec)%range_ndx(irange)))
        end do
    end do

!==================================================================================================
! Allocate and initialize newobj
!==================================================================================================

    allocate(newobj,stat=ierr)
    if ( ierr/=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%bins2ranges_(oslo_sectional_nbins), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%range_nspecies_(oslo_sectional_nranges), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%bin_centers_(oslo_sectional_nbins), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%bin_bounds_(oslo_sectional_nbins, 2), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%range_bounds_(oslo_sectional_nranges, 2), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%aer_spec_prop(oslo_sectional_nspecies_tot), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    allocate(newobj%particle_volume_(oslo_sectional_nbins), stat=ierr)
    if( ierr /=0 ) then
        nullify(newobj)
        return
    end if

    newobj%bins2ranges_ = bins2ranges
    newobj%nranges_ = oslo_sectional_nranges
    newobj%nspecies_tot_ = oslo_sectional_nspecies_tot
    newobj%range_nspecies_ = oslo_sectional_nspecies(:oslo_sectional_nranges)
    newobj%bin_centers_ = bin_centers(:oslo_sectional_nbins)
    newobj%bin_bounds_ = bin_bounds(:oslo_sectional_nbins, :)
    newobj%range_bounds_ = range_bounds(:oslo_sectional_nranges, :)
    newobj%aer_spec_prop = oslo_sectional_species_properties(:oslo_sectional_nspecies_tot)
! TODO: allocate aer_spec_props%bin_ndx and range_ndx?
    newobj%particle_volume_ = 4/3*pi*newobj%bin_centers_**3

    ! deallocate local variables
    if (allocated(bin_centers)) deallocate(bin_centers)
    if (allocated(bin_bounds)) deallocate(bin_bounds)
    if (allocated(range_bounds)) deallocate(range_bounds)
    if (allocated(bins2ranges)) deallocate(bins2ranges)
    if (allocated(oslo_sectional_species_properties)) deallocate(oslo_sectional_species_properties)

    call newobj%initialize(oslo_sectional_nbins, ncnst_tot, nspecies, nspecies, alogsig, f1, f2, ierr)
!==================================================================================================
! Report
!==================================================================================================
!    call MPI_Barrier(mpicom, ierr)
   if (masterproc) then
        write(iulog,*) 'sectional aerosol properties: '
        write(iulog,*) 'nbins = ', newobj%nbins()
        write(iulog,*) 'nranges = ', newobj%nranges_
        write(iulog,*) 'ncnst_tot = ', newobj%ncnst_tot()
        write(iulog,*) 'nspecies_tot = ', newobj%nspecies_tot_
        write(iulog,*) 'nspecies = ', newobj%nspecies()

        do irange=1,oslo_sectional_nranges,5
            write(iulog,*) 'range_nspecies = ', newobj%range_nspecies_(irange:min(irange+4, oslo_sectional_nranges))
        end do

        do ibin=1,oslo_sectional_nbins,5
            write(iulog,*) 'bin_centers = ', newobj%bin_centers_(ibin:min(ibin+4, oslo_sectional_nbins))
        end do

        do ibin=1,oslo_sectional_nbins,5
            write(iulog,*) 'bins2ranges = ',newobj%bins2ranges_(ibin:min(ibin+4, oslo_sectional_nbins))
        end do

        ! TODO: figure out what to do with these
        ! alogsig -> alogsig(m) = log(sigmag(m))
        ! f1 -> f1(m) = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
        ! f2 -> f2(m) = 1._r8 + 0.25_r8*alogsig(m)

        do ibin=1,oslo_sectional_nbins
            write(iulog,*) 'bin_bounds = ', newobj%bin_bounds_(ibin,1), &
                                     ' : ', newobj%bin_bounds_(ibin,2)
        end do
        do irange=1,oslo_sectional_nranges !TODO FIX format
            write(iulog,*) 'range_bounds = ', newobj%range_bounds_(irange,1), &
                                       ' : ', newobj%range_bounds_(irange,2)
        end do

        do ind = 1,oslo_sectional_nspecies_tot
            ! TODO (low priority): fix the format :D
            write(iulog ,*) 'sectional aerosol species properties: '
            write(iulog ,*) 'species name = ', newobj%aer_spec_prop(ind)%specname
            write(iulog ,*) 'range indices = ', newobj%aer_spec_prop(ind)%range_ndx(1), ' : ', newobj%aer_spec_prop(ind)%range_ndx(newobj%aer_spec_prop(ind)%nrange) ! TODO: is there a nicer way?
            write(iulog ,*) 'density = ', newobj%aer_spec_prop(ind)%density
            write(iulog ,*) 'molecular_weight = ', newobj%aer_spec_prop(ind)%molecular_weight
            write(iulog ,*) 'kappa = ', newobj%aer_spec_prop(ind)%kappa
            write(iulog ,*) 'mixed = ', newobj%aer_spec_prop(ind)%mixed
            write(iulog ,*) 'nrange = ', newobj%aer_spec_prop(ind)%nrange
            write(iulog ,*) 'nbin = ', newobj%aer_spec_prop(ind)%nbin
            write(iulog ,*) 'bin_indices = ', newobj%aer_spec_prop(ind)%bin_ndx(1), ' : ', maxval(newobj%aer_spec_prop(ind)%bin_ndx) ! TODO: is there a nicer way?
            write(iulog ,*) 'tracer names = ', newobj%aer_spec_prop(ind)%tracernames(:newobj%aer_spec_prop(ind)%nrange)
       end do
    end if
        prop_obj => newobj
    else
        if ( .not. associated(prop_obj) ) then
            call endrun("Internal Error: sectional_aerosol_properties has not been initialized")
        end if
        newobj => prop_obj

    end if

    ! deallocate local variables
    if (allocated(nspecies)) deallocate(nspecies)
    if (allocated(alogsig)) deallocate(alogsig)
    if (allocated(f1)) deallocate(f1)
    if (allocated(f2)) deallocate(f2)

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


    if (present(list_ndx)) then
        call endrun(subname//' list_ndx is not yet implemented')
    end if

    if (present(density)) then
        call endrun(subname//' density is not yet implemented')
    end if

    if (present(hygro)) then
        call endrun(subname//' hygro is not yet implemented')
    end if

    if (present(spectype)) then
        call endrun(subname//' spectype is not yet implemented')
    end if

    if (present(specname)) then
        specname = self%aer_spec_prop(species_ndx)%specname
    end if

    if (present(specmorph)) then
        call endrun(subname//' specmorph is not yet implemented')
    end if

    if (present(refindex_sw)) then
        call endrun(subname//' refindex_sw is not yet implemented')
    end if

    if (present(refindex_lw)) then
        call endrun(subname//' refindex_lw is not yet implemented')
    end if

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
  ! returns density for a species
  !------------------------------------------------------------------------------
  real(r8) function density(self, species_ndx)

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx

    density = self%aer_spec_prop(species_ndx)%density

  end function density

  !------------------------------------------------------------------------------
  ! returns kappa for a species (hygroscopicity)
  !------------------------------------------------------------------------------
  real(r8) function kappa(self, species_ndx)

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx

    kappa = self%aer_spec_prop(species_ndx)%kappa

  end function kappa

  !------------------------------------------------------------------------------
  ! returns molecular_weight for a species
  !------------------------------------------------------------------------------
  real(r8) function molecular_weight(self, species_ndx)

    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx

    molecular_weight = self%aer_spec_prop(species_ndx)%molecular_weight

  end function molecular_weight

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
  ! returns the total number of species objects
  !------------------------------------------------------------------------------
  function nspecies_tot(self)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer :: res
    character(len=*), parameter :: subname = 'nspecies_tot'

    res = self%nspecies_tot_

  end function nspecies_tot

  !------------------------------------------------------------------------------
  ! returns the total number of species in each range
  !------------------------------------------------------------------------------
  function range_nspecies(self, nranges)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: nranges
    integer :: res(nranges)
    character(len=*), parameter :: subname = 'range_nspecies'

    res = self%range_nspecies_(:min(nranges,size(self%range_nspecies_)))

  end function range_nspecies

  !------------------------------------------------------------------------------
  ! returns the total number of ranges
  !------------------------------------------------------------------------------
  function nranges(self)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer :: res
    character(len=*), parameter :: subname = 'nranges'

    res = self%nranges_

  end function nranges

  !------------------------------------------------------------------------------
  ! returns bin centers
  !------------------------------------------------------------------------------
  function bin_centers(self, nbins) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: nbins
    real(r8) :: res(nbins)
    character(len=*), parameter :: subname = 'bin_centers'

    res = self%bin_centers_(:min(nbins,size(self%bin_centers_)))

  end function bin_centers

  !------------------------------------------------------------------------------
  ! returns bin bounds
  !------------------------------------------------------------------------------
  function bin_bounds(self, nbins) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: nbins
    real(r8) :: res(nbins,2)
    character(len=*), parameter :: subname = 'bin_bounds'

    res = self%bin_bounds_(:min(nbins,size(self%bin_bounds_)),:)

  end function bin_bounds

  !------------------------------------------------------------------------------
  ! returns range bounds
  !------------------------------------------------------------------------------
  function range_bounds(self, nranges) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: nranges
    integer :: res(nranges,2)
    character(len=*), parameter :: subname = 'range_bounds'

    res = self%range_bounds_(:min(nranges,size(self%range_bounds_)),:)

  end function range_bounds

  !------------------------------------------------------------------------------
  ! returns volume of a particle in a bin
  !------------------------------------------------------------------------------
  function particle_volume(self, ibin) result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: ibin
    real(r8) :: res
    character(len=*), parameter :: subname = 'particle_volume'

    res = self%particle_volume_(ibin)

  end function particle_volume
  !------------------------------------------------------------------------------
  ! returns the upper or lower range idx TODO: change to array of idices?
  !------------------------------------------------------------------------------
  function spec_range_ndx(self, species_ndx, nranges)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx, nranges
    integer :: res(nranges)
    character(len=*), parameter :: subname = 'spec_range_ndx'

    res = self%aer_spec_prop(species_ndx)%range_ndx

  end function spec_range_ndx

  !------------------------------------------------------------------------------
  ! returns the upper or lower range idx TODO: change to array of idices?
  !------------------------------------------------------------------------------
  function spec_bin_ndx(self, species_ndx, nbins)  result(res)
    ! TODO: needed?
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx, nbins
    integer :: ibin
    integer :: res(nbins)
    character(len=*), parameter :: subname = 'spec_bin_ndx'

    res = self%aer_spec_prop(species_ndx)%bin_ndx

  end function spec_bin_ndx

  !------------------------------------------------------------------------------
  ! returns the upper or lower range idx TODO: change to array of idices?
  !------------------------------------------------------------------------------
  function spec_nrange(self, species_ndx)  result(res)
    ! TODO: needed?
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx
    integer :: res
    character(len=*), parameter :: subname = 'spec_nrange'

    res = self%aer_spec_prop(species_ndx)%nrange

  end function spec_nrange

  !------------------------------------------------------------------------------
  ! returns number of bins containing species
  !------------------------------------------------------------------------------
  function spec_nbin(self, species_ndx)  result(res)
    ! TODO: needed?
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx
    integer :: res
    character(len=*), parameter :: subname = 'spec_nbin'

    res = self%aer_spec_prop(species_ndx)%nbin

  end function spec_nbin

  !------------------------------------------------------------------------------
  ! returns number of bins containing species
  !------------------------------------------------------------------------------
  function spec_tracernames(self, species_ndx, nranges)  result(res)
    ! TODO: needed?
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in) :: species_ndx, nranges
    character(len=10) :: res(nranges)
    character(len=*), parameter :: subname = 'spec_tracernames'

    res = self%aer_spec_prop(species_ndx)%tracernames

  end function spec_tracernames

  !------------------------------------------------------------------------------
  ! returns the total number of species objects
  !------------------------------------------------------------------------------

  function bins2ranges(self, nbins)  result(res)
    class(sectional_aerosol_properties), intent(in) :: self
    integer, intent(in)         :: nbins
    integer                     :: res(nbins)
    character(len=*), parameter :: subname = 'bins2ranges'

    res = self%bins2ranges_(:min(nbins,size(self%bins2ranges_)))

  end function bins2ranges
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
