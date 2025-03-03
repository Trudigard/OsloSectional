module readnl
   implicit none
   private

   public sec_aer_readnl

contains
   subroutine sec_aer_readnl(nlfile)

      use mpi,            only: mpi_character, mpi_logical, mpi_integer, mpi_real
      use mpi,            only: mpi_bcast, MPI_SUCCESS
      use spmd_utils      only: masterproc, mstrid=>masterprocid, mpicom
      use string_utils,   only: int2str
      use namelist_utlis, only: find_group_name
      use cam_abortutils, only: endrun
      use cam_abortutils, only: iulog

      ! path to namelist file
      character(len=*), intent(in) :: nlfile

      ! namelist variables
      integer :: nbin, nrange, nspec
      !logical :: OA_active, DU_active, SO4_active, BC_active
      !logical :: BC_is_active, SS_active, NO3_active
      character(len=20) :: species_name_list(500), bin_name_list(500)
      real :: bin_list(500), bin_bnd_list(500), range_bnd_list(500)
      character(len=20), allocatable :: spec_names(:), bin_names(:)
      real, allocatable :: bins(:), bin_bnds(:), range_bnds(:)

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'sec_aer_readnl'

      namelist /oslo_sectional_nl/ nbin, nrange, nspec!, OA_active, DU_active, SO4_active
      !namelist /oslo_sectional_nl/ BC_active, BC_is_active, SS_active, NO3_active
      namelist /oslo_sectional_nl/ bin_list, bin_bnd_list, range_bnd_list, species_name_list, bin_name_list

      !-----------------------------------------------------------------------
      ! Initialize all namelist variables
      !-----------------------------------------------------------------------
      nbin = 0 ! nan
      nrange = 0 ! nan
      nspec = 0 ! nan
      !OA_active = .false.
      !DU_active = .false.
      !SO4_active = .false.
      !BC_active = .false.
      !BC_is_active = .false.
      !SS_active = .false.
      !NO3_active = .false.

      !-----------------------------------------------------------------------
      ! Read namelist from file
      !-----------------------------------------------------------------------

      if (masterproc) then
         open(newunit=unitn, file=trim(nlfile), status='old')
         call find_group_name(unitn, 'oslo_sectional_nl', ierr)
         if (ierr == 0) then
            read(unitn, oslo_sectional_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname//'ERROR reading namelist')
            end if
         end if
         close(unitn)
      end if

   !allocate bins, bin_bnds, range_bnds, spec_array
   allocate(bins(nbin))
   allocate(bin_bnds(nbin+1))
   allocate(range_bnds(nrange+1))
   allocate(bin_names(nbin))
   allocate(spec_names(nspec))
   bins = bin_list(1:nbin)
   bin_bnds = bin_bnd_list(1:nbin+1)
   range_bnds = range_bnd_list(1:nrange+1)
   bin_names = bin_name_list(1:nbin)
   spec_names = species_name_list(1:nspec)

   !-----------------------------------------------------------------------
   ! Broadcast values to all MPI tasks
   !-----------------------------------------------------------------------

   call mpi_bcast(nbin, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
     call endrun(subname//"ERROR "//int2str(ierr)//                     &
          " broadcasting 'nbin'")
   call mpi_bcast(nrange, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
      call endrun(subname//"ERROR "//int2str(ierr)//                     &
           " broadcasting 'nrange'")
   call mpi_bcast(nspec, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= MPI_SUCCESS) then
      call endrun(subname//"ERROR "//int2str(ierr)//                     &
           " broadcasting 'nsped'")
   call mpi_bcast(bins, nbin, mpi_real, mstrid, mpicom, ierr)
      if (ierr /= MPI_SUCCESS) then
      call endrun(subname//"ERROR "//int2str(ierr)//                     &
           " broadcasting 'bins'")
   call mpi_bcast(bin_bnds, nbin+1, mpi_real, mstrid, mpicom, ierr)
      if (ierr /= MPI_SUCCESS) then
         call endrun(subname//"ERROR "//int2str(ierr)//                     &
           " broadcasting 'bin_bnds'")
   call mpi_bcast(range_bnds, nrange+1, mpi_real, mstrid, mpicom, ierr)
      if (ierr /= MPI_SUCCESS) then
         call endrun(subname//"ERROR "//int2str(ierr)//                     &
           " broadcasting 'range_bnds'")
   call mpi_bcast(bin_names, nbin, mpi_character, mstrid, mpicom, ierr)
      if (ierr /= MPI_SUCCESS) then
     call endrun(subname//"ERROR "//int2str(ierr)//                     &
          " broadcasting 'bin_names'")
   call mpi_bcast(spec_names, nspec, mpi_character, mstrid, mpicom, ierr)
      if (ierr /= MPI_SUCCESS) then
     call endrun(subname//"ERROR "//int2str(ierr)//                     &
          " broadcasting 'spec_names'")

   !-----------------------------------------------------------------------
   ! Report the settings
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'nbin = ', nbin
      write(iulog,*) 'nrange = ', nrange
      write(iulog,*) 'nspec = ', nspec
      write(iulog,*) 'bins = ', bins
      write(iulog,*) 'bin_bnds = ', bin_bnds
      write(iulog,*) 'range_bnds = ', range_bnds
      write(iulog,*) 'spec_names = ', spec_names
      write(iulog,*) 'bin_names = ', bin_names
   end if

   end subroutine sec_aer_readnl
end module readnl
