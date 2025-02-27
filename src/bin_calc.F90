module bin_calc
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use cam_abortutils,  only: endrun
   use cam_logfile,     only: iulog

   implicit none
   private
   public :: bin_readnl

   contains
      subroutine bin_readnl(nlfile, nbin, nrange, r, r_bnds, organic, dust, sulfate, &
         blackcarbon, blackcarbon_insoluble, seasalt, nitrate, ranges, range_bnds, organic_range, &
         dust_range, sulfate_range, blackcarbon_range, blackcarbon_insoluble_range, seasalt_range, nitrate_range)
         use mpi,             only: mpi_character, mpi_real8, mpi_integer, MPI_SUCCESS
         use spmd_utils,      only: mstrid=>masterproc
         use namelist_utils,  only: find_group_name
         use ifnan,           only: nan
         use bin_calc_utils,  only: calc_radii, check_species_bnds, find_species_bnds, find_range_bnds

         ! Namelist variables
         character(len=*), intent(in) :: nlfile
         integer, intent(out) :: nbin, nrange
         logical, intent(out) :: organic, dust, sulfate, blackcarbon, blackcarbon_insoluble, &
                              seasalt, nitrate, ranges
         real(r8) :: radius_1, radius_N
         real(r8) :: range_bounds_list(100)
         real(r8), intent(out) :: organic_range(2), dust_range(2), sulfate_range(2)
         real(r8), intent(out) :: blackcarbon_range(2), blackcarbon_insoluble_range(2)
         real(r8), intent(out) :: seasalt_range(2), nitrate_range(2)

         ! Derived from namelist variables
         real(r8), allocatable, intent(out) :: range_bnds(:), r(:), r_bnds(:)

         ! Local variables
         integer :: i, ierr, unitn
         character(len=*), parameter :: subname = 'bin_readnl'

         ! Define namelist

         ! Find way to make species more general?

         namelist /oslo_sectional_nl/ nbin, nrange, organic, dust, sulfate, &
                                    blackcarbon, blackcarbon_insoluble, seasalt, &
                                    nitrate, ranges, radius_1, radius_N, range_bounds_list, &
                                    organic_range, dust_range, sulfate_range, blackcarbon_range, &
                                    blackcarbon_insoluble_range, seasalt_range, nitrate_range

         !-----------------------------------------------------------------------
         ! Initialize all namelist variables
         !-----------------------------------------------------------------------
         nbin = nan
         organic = .false.
         dust = .false.
         sulfate = .false.
         blackcarbon = .false.
         blackcarbon_insoluble = .false.
         seasalt = .false.
         nitrate = .false.
         ranges = .false.
         radius_1 = nan
         radius_N = nan
         range_bounds_list = nan
         organic_range = nan
         dust_range = nan
         sulfate_range = nan
         blackcarbon_range = nan
         blackcarbon_insoluble_range = nan
         seasalt_range = nan
         nitrate_range = nan

         !-----------------------------------------------------------------------
         ! Read namelist from file
         !-----------------------------------------------------------------------

         if (masterproc) then
            open(newunit=unitn, file=trim(nlfile), status='old')
            call find_group_name(unitn, 'oslo_sectional_nl', ierr)
            if (ierr == 0) then
               read(unitn, oslo_sectional_nl, iostat=ierr)
               if (ierr /= 0) then
                  call endrun(subname//':: ERROR reading namelist')
               end if
            end if
            close(unitn)
         end if

         ! allocate range_bnds
         allocate(range_bnds(nrange+1))
         range_bnds = range_bounds_list(1:nrange+1)
         !-----------------------------------------------------------------------
         ! Calculate bins
         !-----------------------------------------------------------------------
         call calc_radii(nbin, radius_1, radius_N, r, r_bnds)
         !-----------------------------------------------------------------------
         ! Check species bounds
         !-----------------------------------------------------------------------
         check_species_bnds(organic_range, range_bnds, nrange)
         check_species_bnds(dust_range, range_bnds, nrange)
         check_species_bnds(sulfate_range, range_bnds, nrange)
         check_species_bnds(blackcarbon_range, range_bnds, nrange)
         check_species_bnds(blackcarbon_insoluble_range, range_bnds, nrange)
         check_species_bnds(seasalt_range, range_bnds, nrange)
         check_species_bnds(nitrate_range, range_bnds, nrange)
         !-----------------------------------------------------------------------
         ! Find bin bounds closest to range bounds
         !-----------------------------------------------------------------------
         call find_range_bnds(r_bnds, nbin, nrange, range_bnds, ranges)
         !-----------------------------------------------------------------------
         ! Find species bounds closest to range bounds
         !-----------------------------------------------------------------------
         call find_species_bnds(range_bnds, organic_range, nrange)
         call find_species_bnds(range_bnds, dust_range, nrange)
         call find_species_bnds(range_bnds, sulfate_range, nrange)
         call find_species_bnds(range_bnds, blackcarbon_range, nrange)
         call find_species_bnds(range_bnds, blackcarbon_insoluble_range, nrange)
         call find_species_bnds(range_bnds, seasalt_range, nrange)
         call find_species_bnds(range_bnds, nitrate_range, nrange)
         !-----------------------------------------------------------------------
         ! Insert something smart to make the species be in the correct bins/ranges
         ! solubility
         ! range indices
         !-----------------------------------------------------------------------
         !-----------------------------------------------------------------------
         ! Broadcast values to all MPI tasks
         !-----------------------------------------------------------------------
         ! nbin, nrange, r, r_bnds, range_bnds
         ! species_ranges
         call MPI_Bcast(nbin, 1, mpi_integer, mstrid, mpicom, ierr)
         if (ierr /= MPI_SUCCESS) then
            call endrun(subname//":: ERROR "//int2str(ierr)" broadcasting 'nbin'")
         end if
         call MPI_Bcast(nrange, 1, mpi_integer, mstrid, mpicom, ierr)
         if (ierr /= MPI_SUCCESS) then
            call endrun(subname//":: ERROR "//int2str(ierr)" broadcasting 'nrange'")
         end if
         call MPI_Bcast(r, nbin, mpi_real8, mstrid, mpicom, ierr)
         if (ierr /= MPI_SUCCESS) then
            call endrun(subname//":: ERROR "//int2str(ierr)" broadcasting 'r'")
         end if
         call MPI_Bcast(r_bnds, nbin+1, mpi_real8, mstrid, mpicom, ierr)
         if (ierr /= MPI_SUCCESS) then
            call endrun(subname//":: ERROR "//int2str(ierr)" broadcasting 'r_bnds'")
         end if
         call MPI_Bcast(nrange, 1, mpi_integer, mstrid, mpicom, ierr)
         if (ierr /= MPI_SUCCESS) then
            call endrun(subname//":: ERROR "//int2str(ierr)" broadcasting 'nrange'")
         end if
         call MPI_Bcast(range_bnds, nrange+1, mpi_real8, mstrid, mpicom, ierr)
         if (ierr /= MPI_SUCCESS) then
            call endrun(subname//":: ERROR "//int2str(ierr)" broadcasting 'range_bnds'")
         end if

         !-----------------------------------------------------------------------
         ! Report the settings
         !-----------------------------------------------------------------------
         if (masterproc) then
            write(iulog ,*) 'nbin: ', nbin
            write(iulog ,*) 'nrange: ', nrange
            write(iulog ,*) 'ranges: ', ranges
            write(iulog ,*) 'organic: ', organic
            write(iulog ,*) 'dust: ', dust
            write(iulog ,*) 'sulfate: ', sulfate
            write(iulog ,*) 'blackcarbon: ', blackcarbon
            write(iulog ,*) 'blackcarbon_insoluble: ', blackcarbon_insoluble
            write(iulog ,*) 'seasalt: ', seasalt
            write(iulog ,*) 'nitrate: ', nitrate
            write(iulog ,*) 'organic_range: ', organic_range
            write(iulog ,*) 'dust_range: ', dust_range
            write(iulog ,*) 'sulfate_range: ', sulfate_range
            write(iulog ,*) 'blackcarbon_range: ', blackcarbon_range
            write(iulog ,*) 'blackcarbon_insoluble_range: ', blackcarbon_insoluble_range
            write(iulog ,*) 'seasalt_range: ', seasalt_range
            write(iulog ,*) 'nitrate_range: ', nitrate_range
            do i = 1, nbin, 5
               write(iulog ,*) 'r: ', r(i:min(i+4,nbin))
            end do
            do i = 1, nbin+1, 5
               write(iulog ,*) 'r_bnds: ', r_bnds(i:min(i+4,nbin+1))
            end do
            write(iulog ,*) 'range_bnds: ', range_bnds
         end if


         end subroutine bin_readnl

   end module bin_calc