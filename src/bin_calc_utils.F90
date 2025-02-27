module bin_calc_utils
   use shr_kind_mod, only : r8 => shr_kind_r8
   use cam_abortutils, only : endrun
   use cam_logfile, only : iulog

   implicit none
   private

   public :: calc_radii
   public :: check_species_bnds
   public :: find_species_bnds
   public :: find_range_bnds

contains
   subroutine calc_radii(nbin, radius_1, radius_N, r, r_bnds)
      use physconst, only : pi

      integer, intent(in) :: nbin
      real(r8), intent(in) :: radius_1, radius_N
      real(r8) :: v_rat, v(nbin)
      real(r8), allocatable, intent(out) :: r_bnds(:), r(:)
      real(r8) :: v_bnds(nbin+1)

      ! Local variables
      integer :: i
      character(len=*), parameter :: subname = 'calc_radii'

      ! Allocate arrays for radii and radius bounds
      allocate(r(nbin))
      allocate(r_bnds(nbin+1))

      ! Calculate volume and radius for each bin
      v_rat = (radius_N / radius_1) ** (3._r8 / real(nbin-1, r8)) ! calculate volume ratio
      v(1) = (3._r8/4._r8) * pi * (radius_1) ** 3 ! calculate smallest volume
      do i = 2, nbin
         v(i) = v(1) * v_rat ** real(i-1, r8) ! calculate all volumes for bins
      end do

      r = (v*4._r8/3._r8/pi)**(1._r8/3._r8) ! calculate radii for bins

      v_bnds(1:nbin) = (2._r8*v) / (1._r8+V_rat) ! lower radius bounds
      v_bnds(nbin+1) = V_rat*v_bnds(nbin)        ! upper radius bound

      r_bnds = (v_bnds*4._r8/3._r8/pi)**(1._r8/3._r8) ! calculate radus bounds from volume

   end subroutine calc_radii

   function check_species_bnds(species_bnds, range_bnds, nrange)
      integer, intent(in) :: nrange
      real(r8), intent(in) :: species_bnds(2), range_bnds(nrange+1)

      character(len=*), parameter :: subname = 'check_species_bnds'

      ! Check if species bounds are equal to any range bounds
      ! Use before adjusting any range bounds

      if (.not. any(species_bnds(1) == range_bnds) .or. .not. any(species_bnds(2) == range_bnds)) then
         call endrun('ERROR: Species range bound not equal to range bounds')
      end if

   end function check_species_bnds

   subroutine find_range_bnds(r_bnds, nbin, nrange, range_bnds, ranges)
      ! Find range bounds from bin bounds
      integer, intent(in) :: nbin, nrange
      real(r8), intent(in) :: r_bnds(nbin+1)
      real(r8), intent(inout) :: range_bnds(nrange+1)
      logical, intent(in) :: ranges

      integer :: i, j
      character(len=*), parameter :: subname = 'find_range_bnds'

      if (ranges) then
         range_bnds(1) = r_bnds(1)
         range_bnds(nrange+1) = r_bnds(nbin+1)
         do i = 2, nrange
            do j = 1, nbin
               if (r_bnds(j) < range_bnds(i) .and. r_bnds(j+1) > range_bnds(i)) then
                  if (abs(r_bnds(j)-range_bnds(i)) < abs(r_bnds(j+1)-range_bnds(i))) then
                     range_bnds(i) = r_bnds(j)
                  else
                     range_bnds(i) = r_bnds(j+1)
                  end if
               end if
            end do
         end do
      else
         range_bnds = r_bnds
      end if

   end subroutine find_range_bnds

   subroutine find_species_bnds(range_bnds, spec_bnds, nrange)
      ! Set species bounds to the closest range bound
      integer, intent(in) :: nrange
      real(r8), intent(in) :: range_bnds(nrange+1)
      real(r8), intent(inout) :: spec_bnds(2)
      integer :: idx0, idx1
      real(r8) :: min_diff0, min_diff1

      integer :: i
      character(len=*), parameter :: subname = 'find_species_bnds'

      ! Initialize the minimum differences with a large number
      min_diff0 = abs(range_bnds(1) - spec_bnds(1))
      min_diff1 = abs(range_bnds(1) - spec_bnds(2))
      idx0 = 1
      idx1 = 1

      ! find lower bound for species
      do i = 1, nrange
          if (abs(range_bnds(i) - spec_bnds(1)) < min_diff0) then
              min_diff0 = abs(range_bnds(i) - spec_bnds(1))
              idx0 = i
          end if
      end do

      ! find upper bound for species
      do i = 1, nrange+1
          if (abs(range_bnds(i) - spec_bnds(2)) < min_diff1) then
              min_diff1 = abs(range_bnds(i) - spec_bnds(2))
              idx1 = i
          end if
      end do

      ! Update spec_bnds with the closest range bounds
      spec_bnds(1) = range_bnds(idx0)
      spec_bnds(2) = range_bnds(idx1)

   end subroutine find_species_bnds

end module bin_calc_utils