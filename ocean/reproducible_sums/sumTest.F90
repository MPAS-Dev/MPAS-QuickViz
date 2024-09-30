!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! \file sumTest.F
!
! \brief  This is a test driver for reproducible sum algorithms. It
!         tests global sum for multiple data types and array dimensions.
!         Compile and run with, e.g. 
!            mpif90 mpas_global_sum_mod.F sumTest.F90 -fno-range-check -ffree-form
!            mpirun -n 36 a.out
!
!***********************************************************************

program sumTest

   use mpas_global_sum_mod

   implicit none

   include "mpif.h"


   !--------------------------------------------------------------------
   ! Local variables
   !--------------------------------------------------------------------

   integer, parameter :: & 
      R8 = selected_real_kind(12)

   ! Set local array sizes (global size is simply the local size*nRanks)
   ! for arrays up to five dimensions. Choose smaller sizes for higher
   ! dimensions to prevent too large a footprint.

   integer, parameter :: & 
      nx = 100,         &
      ny = 100,         &
      nz = 100,         &
      nl =  10,         &
      nm =   2,         &
      nn =   2

   integer, dimension(12) :: &
      indxRange = (/ 3, 90, 5, 95, 10, 90, 2, 9, 1, 2, 1, 2 /) 

   integer, parameter :: &
      nFields = 3,       &! test multiple field interface (keep=3)
      nIters  = 100       ! use a number of iterations for timing

   integer :: &
      i,j,k,l,m,n,r,nfld,iter, &! loop iterators
      i2,j2,k2,l2,m2,          &! loop iterators
      i3,j3,k3,l3,m3,          &! loop iterators
      ierr,                    &! local error code
      drange, rrange,          &! range of exponents
      irange, lrange,          &! range of integer kinds
      ipower,                  &! exponent to use
      myRank, nRanks,          &! MPI rank, number of ranks
      masterTask                ! MPI rank to use as master
      
   real :: &
      t1, t2, dt          ! timer variables

   real, dimension(:), allocatable :: &
      randoms          ! an array of random numbers

   ! arrays and reference vals for testing various types and dims
   double precision :: dref, dref2, dscalar, dsum1, dsum2, dsum3, dsum4
   double precision, dimension(:), allocatable :: dsum1n, dsum2n, dsum3n, dsum4n
   double precision, dimension(:), allocatable :: darray1d, dmask1d, darray0n
   double precision, dimension(:,:), allocatable :: darray2d, dmask2d, &
                                                    darray1n, dmask1n
   double precision, dimension(:,:,:), allocatable :: darray3d, dmask3d, &
                                                      darray2n, dmask2n
   double precision, dimension(:,:,:,:), allocatable :: darray4d, dmask4d, &
                                                        darray3n, dmask3n
   double precision, dimension(:,:,:,:,:), allocatable :: darray5d, dmask5d, &
                                                          darray4n, dmask4n
   double precision, dimension(:,:,:,:,:,:), allocatable :: darray6d, dmask6d, &
                                                            darray5n, dmask5n

   real :: rref, rref2, rscalar, rsum1, rsum2, rsum3, rsum4
   real, dimension(:), allocatable :: rsum1n, rsum2n, rsum3n, rsum4n
   real, dimension(:), allocatable :: rarray1d, rmask1d, rarray0n
   real, dimension(:,:), allocatable :: rarray2d, rmask2d, &
                                        rarray1n, rmask1n
   real, dimension(:,:,:), allocatable :: rarray3d, rmask3d, &
                                          rarray2n, rmask2n
   real, dimension(:,:,:,:), allocatable :: rarray4d, rmask4d, &
                                            rarray3n, rmask3n
   real, dimension(:,:,:,:,:), allocatable :: rarray5d, rmask5d, &
                                              rarray4n, rmask4n
   real, dimension(:,:,:,:,:,:), allocatable :: rarray6d, rmask6d, &
                                                rarray5n, rmask5n

   integer :: iref, iref2, iscalar, isum1, isum2, isum3, isum4, maxint
   integer, dimension(:), allocatable :: isum1n, isum2n, isum3n, isum4n
   integer, dimension(:), allocatable :: iarray1d, imask1d, iarray0n
   integer, dimension(:,:), allocatable :: iarray2d, imask2d, &
                                           iarray1n, imask1n
   integer, dimension(:,:,:), allocatable :: iarray3d, imask3d, &
                                             iarray2n, imask2n
   integer, dimension(:,:,:,:), allocatable :: iarray4d, imask4d, &
                                               iarray3n, imask3n
   integer, dimension(:,:,:,:,:), allocatable :: iarray5d, imask5d, &
                                                 iarray4n, imask4n
   integer, dimension(:,:,:,:,:,:), allocatable :: iarray6d, imask6d, &
                                                    iarray5n, imask5n

   integer (KIND=8) :: lref, lref2, lscalar, lsum1, lsum2, lsum3, lsum4
   integer (KIND=8), dimension(:), allocatable :: lsum1n, lsum2n, lsum3n, lsum4n
   integer (KIND=8), dimension(:), allocatable :: larray1d, lmask1d, larray0n
   integer (KIND=8), dimension(:,:), allocatable :: larray2d, lmask2d, &
                                                    larray1n, lmask1n
   integer (KIND=8), dimension(:,:,:), allocatable :: larray3d, lmask3d, &
                                                      larray2n, lmask2n
   integer (KIND=8), dimension(:,:,:,:), allocatable :: larray4d, lmask4d, &
                                                        larray3n, lmask3n
   integer (KIND=8), dimension(:,:,:,:,:), allocatable :: larray5d, lmask5d, &
                                                          larray4n, lmask4n
   integer (KIND=8), dimension(:,:,:,:,:,:), allocatable :: larray6d, lmask6d, &
                                                            larray5n, lmask5n

   ! End preamble
   !-------------
   ! Begin code

   ! Initialize MPI and retrieve some info
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myRank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)
   masterTask = 0

   ! Initialize random number generator
   call random_seed()

   ! Compute a max value for i4 to avoid overflow
   maxint = (huge(iref)-1)/nRanks

   !--------------------------------------------------------------------
   ! First a quick sanity check for all by doing a simple sum
   ! Test all single-field interfaces, including index range and
   ! sum-product (mask)
   !--------------------------------------------------------------------

   !*** Scalars
   if (myRank == masterTask) print *,'Testing scalar sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   dscalar = real(myRank,R8)
   rscalar = real(myRank)
   iscalar = myRank
   lscalar = myRank
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   do r=0,nRanks-1
      dref = dref + real(r,R8)
      rref = rref + real(r)
      iref = iref + r
      lref = lref + r
   end do
   dsum1 = mpas_global_sum(dscalar, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rscalar, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iscalar, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(lscalar, MPI_COMM_WORLD)
   if (dsum1 == dref) then
      print *, 'Global sum scalar double: PASS'
   else
      print *, 'Global sum scalar double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum scalar real: PASS'
   else
      print *, 'Global sum scalar real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum scalar integer: PASS'
   else
      print *, 'Global sum scalar integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum scalar integer(8): PASS'
   else
      print *, 'Global sum scalar integer(8): FAIL', lsum1, lref
   endif

   !*** 1-d arrays
   if (myRank == masterTask) print *,'Testing 1d sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray1d(nx), dmask1d(nx), rarray1d(nx), rmask1d(nx), &
            iarray1d(nx), imask1d(nx), larray1d(nx), lmask1d(nx))
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   dref2 = 0.d0
   rref2 = 0.0
   iref2 = 0
   lref2 = 0
   do i=1,nx
      lscalar = i*(myRank+1)
      darray1d(i) = real(lscalar,R8)
      rarray1d(i) = real(lscalar)
      larray1d(i) = lscalar
      iarray1d(i) = min(lscalar,maxint)
      do r=0,nRanks-1
         dref = dref + real(i*(r+1),R8)
         lref = lref + i*(r+1)
      end do
      if (i >= indxRange(1) .and. i <= indxRange(2)) then
         dmask1d(i) = 1.d0
         rmask1d(i) = 1.0
         lmask1d(i) = 1
         imask1d(i) = 1
         do r=0,nRanks-1
            dref2 = dref2 + real(i*(r+1),R8)
            lref2 = lref2 + i*(r+1)
         end do
      else
         dmask1d(i) = 0.d0
         rmask1d(i) = 0.0
         lmask1d(i) = 0
         imask1d(i) = 0
      endif
   end do
   rref  = dref
   rref2 = dref2
   iref  = lref
   iref2 = lref2
   dsum1 = mpas_global_sum(darray1d, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rarray1d, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iarray1d, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(larray1d, MPI_COMM_WORLD)
   dsum2 = mpas_global_sum(darray1d, MPI_COMM_WORLD, indxRange)
   rsum2 = mpas_global_sum(rarray1d, MPI_COMM_WORLD, indxRange)
   isum2 = mpas_global_sum(iarray1d, MPI_COMM_WORLD, indxRange)
   lsum2 = mpas_global_sum(larray1d, MPI_COMM_WORLD, indxRange)
   dsum3 = mpas_global_sum(darray1d, dmask1d, MPI_COMM_WORLD)
   rsum3 = mpas_global_sum(rarray1d, rmask1d, MPI_COMM_WORLD)
   isum3 = mpas_global_sum(iarray1d, imask1d, MPI_COMM_WORLD)
   lsum3 = mpas_global_sum(larray1d, lmask1d, MPI_COMM_WORLD)
   dsum4 = mpas_global_sum(darray1d, dmask1d, MPI_COMM_WORLD, indxRange)
   rsum4 = mpas_global_sum(rarray1d, rmask1d, MPI_COMM_WORLD, indxRange)
   isum4 = mpas_global_sum(iarray1d, imask1d, MPI_COMM_WORLD, indxRange)
   lsum4 = mpas_global_sum(larray1d, lmask1d, MPI_COMM_WORLD, indxRange)
   if (dsum1 == dref) then
      print *, 'Global sum array1d double: PASS'
   else
      print *, 'Global sum array1d double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum array1d real: PASS'
   else
      print *, 'Global sum array1d real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum array1d integer: PASS'
   else
      print *, 'Global sum array1d integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum array1d integer(8): PASS'
   else
      print *, 'Global sum array1d integer(8): FAIL', lsum1, lref
   endif
   if (dsum2 == dref2) then
      print *, 'Global sum irange array1d double: PASS'
   else
      print *, 'Global sum irange array1d double: FAIL', dsum2, dref2
   endif
   if (rsum2 == rref2) then
      print *, 'Global sum irange array1d real: PASS'
   else
      print *, 'Global sum irange array1d real: FAIL', rsum2, rref2
   endif
   if (isum2 == iref2) then
      print *, 'Global sum irange array1d integer: PASS'
   else
      print *, 'Global sum irange array1d integer: FAIL', isum2, iref2
   endif
   if (lsum2 == lref2) then
      print *, 'Global sum irange array1d integer(8): PASS'
   else
      print *, 'Global sum irange array1d integer(8): FAIL', lsum2, lref2
   endif
   if (dsum3 == dref2) then
      print *, 'Global sum prod array1d double: PASS'
   else
      print *, 'Global sum prod array1d double: FAIL', dsum3, dref2
   endif
   if (rsum3 == rref2) then
      print *, 'Global sum prod array1d real: PASS'
   else
      print *, 'Global sum prod array1d real: FAIL', rsum3, rref2
   endif
   if (isum3 == iref2) then
      print *, 'Global sum prod array1d integer: PASS'
   else
      print *, 'Global sum prod array1d integer: FAIL', isum3, iref2
   endif
   if (lsum3 == lref2) then
      print *, 'Global sum prod array1d integer(8): PASS'
   else
      print *, 'Global sum prod array1d integer(8): FAIL', lsum3, lref2
   endif
   if (dsum4 == dref2) then
      print *, 'Global sum prod irange array1d double: PASS'
   else
      print *, 'Global sum prod irange array1d double: FAIL', dsum4, dref2
   endif
   if (rsum4 == rref2) then
      print *, 'Global sum prod irange array1d real: PASS'
   else
      print *, 'Global sum prod irange array1d real: FAIL', rsum4, rref2
   endif
   if (isum4 == iref2) then
      print *, 'Global sum prod irange array1d integer: PASS'
   else
      print *, 'Global sum prod irange array1d integer: FAIL', isum4, iref2
   endif
   if (lsum4 == lref2) then
      print *, 'Global sum prod irange array1d integer(8): PASS'
   else
      print *, 'Global sum prod irange array1d integer(8): FAIL', lsum4, lref2
   endif
   deallocate(darray1d, dmask1d, rarray1d, rmask1d, &
              iarray1d, imask1d, larray1d, lmask1d)

   !*** 2-d arrays
   if (myRank == masterTask) print *,'Testing 2d sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray2d(nx,ny), dmask2d(nx,ny), rarray2d(nx,ny), rmask2d(nx,ny), &
            iarray2d(nx,ny), imask2d(nx,ny), larray2d(nx,ny), lmask2d(nx,ny))
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   dref2 = 0.d0
   rref2 = 0.0
   iref2 = 0
   lref2 = 0
   do j=1,ny
   do i=1,nx
      lscalar = (i+j)*(myRank+1)
      darray2d(i,j) = real(lscalar,R8)
      rarray2d(i,j) = real(lscalar)
      larray2d(i,j) = lscalar
      iarray2d(i,j) = min(lscalar,maxint)
      do r=0,nRanks-1
         dref = dref + real((i+j)*(r+1),R8)
         lref = lref + (i+j)*(r+1)
      end do
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4)) then
         dmask2d(i,j) = 1.d0
         rmask2d(i,j) = 1.0
         lmask2d(i,j) = 1
         imask2d(i,j) = 1
         do r=0,nRanks-1
            dref2 = dref2 + real((i+j)*(r+1),R8)
            lref2 = lref2 + (i+j)*(r+1)
         end do
      else
         dmask2d(i,j) = 0.d0
         rmask2d(i,j) = 0.0
         lmask2d(i,j) = 0
         imask2d(i,j) = 0
      endif
   end do
   end do
   rref  = dref
   rref2 = dref2
   iref  = lref
   iref2 = lref2
   dsum1 = mpas_global_sum(darray2d, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rarray2d, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iarray2d, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(larray2d, MPI_COMM_WORLD)
   dsum2 = mpas_global_sum(darray2d, MPI_COMM_WORLD, indxRange)
   rsum2 = mpas_global_sum(rarray2d, MPI_COMM_WORLD, indxRange)
   isum2 = mpas_global_sum(iarray2d, MPI_COMM_WORLD, indxRange)
   lsum2 = mpas_global_sum(larray2d, MPI_COMM_WORLD, indxRange)
   dsum3 = mpas_global_sum(darray2d, dmask2d, MPI_COMM_WORLD)
   rsum3 = mpas_global_sum(rarray2d, rmask2d, MPI_COMM_WORLD)
   isum3 = mpas_global_sum(iarray2d, imask2d, MPI_COMM_WORLD)
   lsum3 = mpas_global_sum(larray2d, lmask2d, MPI_COMM_WORLD)
   dsum4 = mpas_global_sum(darray2d, dmask2d, MPI_COMM_WORLD, indxRange)
   rsum4 = mpas_global_sum(rarray2d, rmask2d, MPI_COMM_WORLD, indxRange)
   isum4 = mpas_global_sum(iarray2d, imask2d, MPI_COMM_WORLD, indxRange)
   lsum4 = mpas_global_sum(larray2d, lmask2d, MPI_COMM_WORLD, indxRange)
   if (dsum1 == dref) then
      print *, 'Global sum array2d double: PASS'
   else
      print *, 'Global sum array2d double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum array2d real: PASS'
   else
      print *, 'Global sum array2d real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum array2d integer: PASS'
   else
      print *, 'Global sum array2d integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum array2d integer(8): PASS'
   else
      print *, 'Global sum array2d integer(8): FAIL', lsum1, lref
   endif
   if (dsum2 == dref2) then
      print *, 'Global sum irange array2d double: PASS'
   else
      print *, 'Global sum irange array2d double: FAIL', dsum2, dref2
   endif
   if (rsum2 == rref2) then
      print *, 'Global sum irange array2d real: PASS'
   else
      print *, 'Global sum irange array2d real: FAIL', rsum2, rref2
   endif
   if (isum2 == iref2) then
      print *, 'Global sum irange array2d integer: PASS'
   else
      print *, 'Global sum irange array2d integer: FAIL', isum2, iref2
   endif
   if (lsum2 == lref2) then
      print *, 'Global sum irange array2d integer(8): PASS'
   else
      print *, 'Global sum irange array2d integer(8): FAIL', lsum2, lref2
   endif
   if (dsum3 == dref2) then
      print *, 'Global sum prod array2d double: PASS'
   else
      print *, 'Global sum prod array2d double: FAIL', dsum3, dref2
   endif
   if (rsum3 == rref2) then
      print *, 'Global sum prod array2d real: PASS'
   else
      print *, 'Global sum prod array2d real: FAIL', rsum3, rref2
   endif
   if (isum3 == iref2) then
      print *, 'Global sum prod array2d integer: PASS'
   else
      print *, 'Global sum prod array2d integer: FAIL', isum3, iref2
   endif
   if (lsum3 == lref2) then
      print *, 'Global sum prod array2d integer(8): PASS'
   else
      print *, 'Global sum prod array2d integer(8): FAIL', lsum3, lref2
   endif
   if (dsum4 == dref2) then
      print *, 'Global sum prod irange array2d double: PASS'
   else
      print *, 'Global sum prod irange array2d double: FAIL', dsum4, dref2
   endif
   if (rsum4 == rref2) then
      print *, 'Global sum prod irange array2d real: PASS'
   else
      print *, 'Global sum prod irange array2d real: FAIL', rsum4, rref2
   endif
   if (isum4 == iref2) then
      print *, 'Global sum prod irange array2d integer: PASS'
   else
      print *, 'Global sum prod irange array2d integer: FAIL', isum4, iref2
   endif
   if (lsum4 == lref2) then
      print *, 'Global sum prod irange array2d integer(8): PASS'
   else
      print *, 'Global sum prod irange array2d integer(8): FAIL', lsum4, lref2
   endif
   deallocate(darray2d, dmask2d, rarray2d, rmask2d, &
              iarray2d, imask2d, larray2d, lmask2d)


   !*** 3-d arrays
   if (myRank == masterTask) print *,'Testing 3d sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray3d(nx,ny,nz), dmask3d(nx,ny,nz), &
            rarray3d(nx,ny,nz), rmask3d(nx,ny,nz), &
            iarray3d(nx,ny,nz), imask3d(nx,ny,nz), &
            larray3d(nx,ny,nz), lmask3d(nx,ny,nz))
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   dref2 = 0.d0
   rref2 = 0.0
   iref2 = 0
   lref2 = 0
   do k=1,nz
   do j=1,ny
   do i=1,nx
      lscalar = (i+j+k)*(myRank+1)
      darray3d(i,j,k) = real(lscalar,R8)
      rarray3d(i,j,k) = real(lscalar)
      larray3d(i,j,k) = lscalar
      iarray3d(i,j,k) = min(lscalar,maxint)
      do r=0,nRanks-1
         dref = dref + real((i+j+k)*(r+1),R8)
         lref = lref + (i+j+k)*(r+1)
      end do
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6)) then
         dmask3d(i,j,k) = 1.d0
         rmask3d(i,j,k) = 1.0
         lmask3d(i,j,k) = 1
         imask3d(i,j,k) = 1
         do r=0,nRanks-1
            dref2 = dref2 + real((i+j+k)*(r+1),R8)
            lref2 = lref2 + (i+j+k)*(r+1)
         end do
      else
         dmask3d(i,j,k) = 0.d0
         rmask3d(i,j,k) = 0.0
         lmask3d(i,j,k) = 0
         imask3d(i,j,k) = 0
      endif
   end do
   end do
   end do
   rref  = dref
   rref2 = dref2
   iref  = lref
   iref2 = lref2
   dsum1 = mpas_global_sum(darray3d, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rarray3d, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iarray3d, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(larray3d, MPI_COMM_WORLD)
   dsum2 = mpas_global_sum(darray3d, MPI_COMM_WORLD, indxRange)
   rsum2 = mpas_global_sum(rarray3d, MPI_COMM_WORLD, indxRange)
   isum2 = mpas_global_sum(iarray3d, MPI_COMM_WORLD, indxRange)
   lsum2 = mpas_global_sum(larray3d, MPI_COMM_WORLD, indxRange)
   dsum3 = mpas_global_sum(darray3d, dmask3d, MPI_COMM_WORLD)
   rsum3 = mpas_global_sum(rarray3d, rmask3d, MPI_COMM_WORLD)
   isum3 = mpas_global_sum(iarray3d, imask3d, MPI_COMM_WORLD)
   lsum3 = mpas_global_sum(larray3d, lmask3d, MPI_COMM_WORLD)
   dsum4 = mpas_global_sum(darray3d, dmask3d, MPI_COMM_WORLD, indxRange)
   rsum4 = mpas_global_sum(rarray3d, rmask3d, MPI_COMM_WORLD, indxRange)
   isum4 = mpas_global_sum(iarray3d, imask3d, MPI_COMM_WORLD, indxRange)
   lsum4 = mpas_global_sum(larray3d, lmask3d, MPI_COMM_WORLD, indxRange)
   if (dsum1 == dref) then
      print *, 'Global sum array3d double: PASS'
   else
      print *, 'Global sum array3d double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum array3d real: PASS'
   else
      print *, 'Global sum array3d real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum array3d integer: PASS'
   else
      print *, 'Global sum array3d integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum array3d integer(8): PASS'
   else
      print *, 'Global sum array3d integer(8): FAIL', lsum1, lref
   endif
   if (dsum2 == dref2) then
      print *, 'Global sum irange array3d double: PASS'
   else
      print *, 'Global sum irange array3d double: FAIL', dsum2, dref2
   endif
   if (rsum2 == rref2) then
      print *, 'Global sum irange array3d real: PASS'
   else
      print *, 'Global sum irange array3d real: FAIL', rsum2, rref2
   endif
   if (isum2 == iref2) then
      print *, 'Global sum irange array3d integer: PASS'
   else
      print *, 'Global sum irange array3d integer: FAIL', isum2, iref2
   endif
   if (lsum2 == lref2) then
      print *, 'Global sum irange array3d integer(8): PASS'
   else
      print *, 'Global sum irange array3d integer(8): FAIL', lsum2, lref2
   endif
   if (dsum3 == dref2) then
      print *, 'Global sum prod array3d double: PASS'
   else
      print *, 'Global sum prod array3d double: FAIL', dsum3, dref2
   endif
   if (rsum3 == rref2) then
      print *, 'Global sum prod array3d real: PASS'
   else
      print *, 'Global sum prod array3d real: FAIL', rsum3, rref2
   endif
   if (isum3 == iref2) then
      print *, 'Global sum prod array3d integer: PASS'
   else
      print *, 'Global sum prod array3d integer: FAIL', isum3, iref2
   endif
   if (lsum3 == lref2) then
      print *, 'Global sum prod array3d integer(8): PASS'
   else
      print *, 'Global sum prod array3d integer(8): FAIL', lsum3, lref2
   endif
   if (dsum4 == dref2) then
      print *, 'Global sum prod irange array3d double: PASS'
   else
      print *, 'Global sum prod irange array3d double: FAIL', dsum4, dref2
   endif
   if (rsum4 == rref2) then
      print *, 'Global sum prod irange array3d real: PASS'
   else
      print *, 'Global sum prod irange array3d real: FAIL', rsum4, rref2
   endif
   if (isum4 == iref2) then
      print *, 'Global sum prod irange array3d integer: PASS'
   else
      print *, 'Global sum prod irange array3d integer: FAIL', isum4, iref2
   endif
   if (lsum4 == lref2) then
      print *, 'Global sum prod irange array3d integer(8): PASS'
   else
      print *, 'Global sum prod irange array3d integer(8): FAIL', lsum4, lref2
   endif
   deallocate(darray3d, dmask3d, rarray3d, rmask3d, &
              iarray3d, imask3d, larray3d, lmask3d)

   !*** 4-d arrays
   if (myRank == masterTask) print *,'Testing 4d sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray4d(nx,ny,nz,nl), dmask4d(nx,ny,nz,nl), &
            rarray4d(nx,ny,nz,nl), rmask4d(nx,ny,nz,nl), &
            iarray4d(nx,ny,nz,nl), imask4d(nx,ny,nz,nl), &
            larray4d(nx,ny,nz,nl), lmask4d(nx,ny,nz,nl))
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   dref2 = 0.d0
   rref2 = 0.0
   iref2 = 0
   lref2 = 0
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      lscalar = (i+j+k+l)*(myRank+1)
      darray4d(i,j,k,l) = real(lscalar,R8)
      rarray4d(i,j,k,l) = real(lscalar)
      larray4d(i,j,k,l) = lscalar
      iarray4d(i,j,k,l) = min(lscalar,maxint)
      do r=0,nRanks-1
         dref = dref + real((i+j+k+l)*(r+1),R8)
         lref = lref + (i+j+k+l)*(r+1)
      end do
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6) .and. &
          l >= indxRange(7) .and. l <= indxRange(8)) then
         dmask4d(i,j,k,l) = 1.d0
         rmask4d(i,j,k,l) = 1.0
         lmask4d(i,j,k,l) = 1
         imask4d(i,j,k,l) = 1
         do r=0,nRanks-1
            dref2 = dref2 + real((i+j+k+l)*(r+1),R8)
            lref2 = lref2 + (i+j+k+l)*(r+1)
         end do
      else
         dmask4d(i,j,k,l) = 0.d0
         rmask4d(i,j,k,l) = 0.0
         lmask4d(i,j,k,l) = 0
         imask4d(i,j,k,l) = 0
      endif
   end do
   end do
   end do
   end do
   rref  = dref
   rref2 = dref2
   iref  = lref
   iref2 = lref2
   dsum1 = mpas_global_sum(darray4d, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rarray4d, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iarray4d, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(larray4d, MPI_COMM_WORLD)
   dsum2 = mpas_global_sum(darray4d, MPI_COMM_WORLD, indxRange)
   rsum2 = mpas_global_sum(rarray4d, MPI_COMM_WORLD, indxRange)
   isum2 = mpas_global_sum(iarray4d, MPI_COMM_WORLD, indxRange)
   lsum2 = mpas_global_sum(larray4d, MPI_COMM_WORLD, indxRange)
   dsum3 = mpas_global_sum(darray4d, dmask4d, MPI_COMM_WORLD)
   rsum3 = mpas_global_sum(rarray4d, rmask4d, MPI_COMM_WORLD)
   isum3 = mpas_global_sum(iarray4d, imask4d, MPI_COMM_WORLD)
   lsum3 = mpas_global_sum(larray4d, lmask4d, MPI_COMM_WORLD)
   dsum4 = mpas_global_sum(darray4d, dmask4d, MPI_COMM_WORLD, indxRange)
   rsum4 = mpas_global_sum(rarray4d, rmask4d, MPI_COMM_WORLD, indxRange)
   isum4 = mpas_global_sum(iarray4d, imask4d, MPI_COMM_WORLD, indxRange)
   lsum4 = mpas_global_sum(larray4d, lmask4d, MPI_COMM_WORLD, indxRange)
   if (dsum1 == dref) then
      print *, 'Global sum array4d double: PASS'
   else
      print *, 'Global sum array4d double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum array4d real: PASS'
   else
      print *, 'Global sum array4d real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum array4d integer: PASS'
   else
      print *, 'Global sum array4d integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum array4d integer(8): PASS'
   else
      print *, 'Global sum array4d integer(8): FAIL', lsum1, lref
   endif
   if (dsum2 == dref2) then
      print *, 'Global sum irange array4d double: PASS'
   else
      print *, 'Global sum irange array4d double: FAIL', dsum2, dref2
   endif
   if (rsum2 == rref2) then
      print *, 'Global sum irange array4d real: PASS'
   else
      print *, 'Global sum irange array4d real: FAIL', rsum2, rref2
   endif
   if (isum2 == iref2) then
      print *, 'Global sum irange array4d integer: PASS'
   else
      print *, 'Global sum irange array4d integer: FAIL', isum2, iref2
   endif
   if (lsum2 == lref2) then
      print *, 'Global sum irange array4d integer(8): PASS'
   else
      print *, 'Global sum irange array4d integer(8): FAIL', lsum2, lref2
   endif
   if (dsum3 == dref2) then
      print *, 'Global sum prod array4d double: PASS'
   else
      print *, 'Global sum prod array4d double: FAIL', dsum3, dref2
   endif
   if (rsum3 == rref2) then
      print *, 'Global sum prod array4d real: PASS'
   else
      print *, 'Global sum prod array4d real: FAIL', rsum3, rref2
   endif
   if (isum3 == iref2) then
      print *, 'Global sum prod array4d integer: PASS'
   else
      print *, 'Global sum prod array4d integer: FAIL', isum3, iref2
   endif
   if (lsum3 == lref2) then
      print *, 'Global sum prod array4d integer(8): PASS'
   else
      print *, 'Global sum prod array4d integer(8): FAIL', lsum3, lref2
   endif
   if (dsum4 == dref2) then
      print *, 'Global sum prod irange array4d double: PASS'
   else
      print *, 'Global sum prod irange array4d double: FAIL', dsum4, dref2
   endif
   if (rsum4 == rref2) then
      print *, 'Global sum prod irange array4d real: PASS'
   else
      print *, 'Global sum prod irange array4d real: FAIL', rsum4, rref2
   endif
   if (isum4 == iref2) then
      print *, 'Global sum prod irange array4d integer: PASS'
   else
      print *, 'Global sum prod irange array4d integer: FAIL', isum4, iref2
   endif
   if (lsum4 == lref2) then
      print *, 'Global sum prod irange array4d integer(8): PASS'
   else
      print *, 'Global sum prod irange array4d integer(8): FAIL', lsum4, lref2
   endif
   deallocate(darray4d, dmask4d, rarray4d, rmask4d, &
              iarray4d, imask4d, larray4d, lmask4d)


   !*** 5-d arrays
   if (myRank == masterTask) print *,'Testing 5d sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray5d(nx,ny,nz,nl,nm), dmask5d(nx,ny,nz,nl,nm), &
            rarray5d(nx,ny,nz,nl,nm), rmask5d(nx,ny,nz,nl,nm), &
            iarray5d(nx,ny,nz,nl,nm), imask5d(nx,ny,nz,nl,nm), &
            larray5d(nx,ny,nz,nl,nm), lmask5d(nx,ny,nz,nl,nm))
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   dref2 = 0.d0
   rref2 = 0.0
   iref2 = 0
   lref2 = 0
   do m=1,nm
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      lscalar = (i+j+k+l+m)*(myRank+1)
      darray5d(i,j,k,l,m) = real(lscalar,R8)
      rarray5d(i,j,k,l,m) = real(lscalar)
      larray5d(i,j,k,l,m) = lscalar
      iarray5d(i,j,k,l,m) = min(lscalar,maxint)
      do r=0,nRanks-1
         dref = dref + real((i+j+k+l+m)*(r+1),R8)
         lref = lref + (i+j+k+l+m)*(r+1)
      end do
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6) .and. &
          l >= indxRange(7) .and. l <= indxRange(8) .and. &
          m >= indxRange(9) .and. m <= indxRange(10)) then
         dmask5d(i,j,k,l,m) = 1.d0
         rmask5d(i,j,k,l,m) = 1.0
         lmask5d(i,j,k,l,m) = 1
         imask5d(i,j,k,l,m) = 1
         do r=0,nRanks-1
            dref2 = dref2 + real((i+j+k+l+m)*(r+1),R8)
            lref2 = lref2 + (i+j+k+l+m)*(r+1)
         end do
      else
         dmask5d(i,j,k,l,m) = 0.d0
         rmask5d(i,j,k,l,m) = 0.0
         lmask5d(i,j,k,l,m) = 0
         imask5d(i,j,k,l,m) = 0
      endif
   end do
   end do
   end do
   end do
   end do
   rref  = dref
   rref2 = dref2
   iref  = lref
   iref2 = lref2
   dsum1 = mpas_global_sum(darray5d, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rarray5d, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iarray5d, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(larray5d, MPI_COMM_WORLD)
   dsum2 = mpas_global_sum(darray5d, MPI_COMM_WORLD, indxRange)
   rsum2 = mpas_global_sum(rarray5d, MPI_COMM_WORLD, indxRange)
   isum2 = mpas_global_sum(iarray5d, MPI_COMM_WORLD, indxRange)
   lsum2 = mpas_global_sum(larray5d, MPI_COMM_WORLD, indxRange)
   dsum3 = mpas_global_sum(darray5d, dmask5d, MPI_COMM_WORLD)
   rsum3 = mpas_global_sum(rarray5d, rmask5d, MPI_COMM_WORLD)
   isum3 = mpas_global_sum(iarray5d, imask5d, MPI_COMM_WORLD)
   lsum3 = mpas_global_sum(larray5d, lmask5d, MPI_COMM_WORLD)
   dsum4 = mpas_global_sum(darray5d, dmask5d, MPI_COMM_WORLD, indxRange)
   rsum4 = mpas_global_sum(rarray5d, rmask5d, MPI_COMM_WORLD, indxRange)
   isum4 = mpas_global_sum(iarray5d, imask5d, MPI_COMM_WORLD, indxRange)
   lsum4 = mpas_global_sum(larray5d, lmask5d, MPI_COMM_WORLD, indxRange)
   if (dsum1 == dref) then
      print *, 'Global sum array5d double: PASS'
   else
      print *, 'Global sum array5d double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum array5d real: PASS'
   else
      print *, 'Global sum array5d real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum array5d integer: PASS'
   else
      print *, 'Global sum array5d integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum array5d integer(8): PASS'
   else
      print *, 'Global sum array5d integer(8): FAIL', lsum1, lref
   endif
   if (dsum2 == dref2) then
      print *, 'Global sum irange array5d double: PASS'
   else
      print *, 'Global sum irange array5d double: FAIL', dsum2, dref2
   endif
   if (rsum2 == rref2) then
      print *, 'Global sum irange array5d real: PASS'
   else
      print *, 'Global sum irange array5d real: FAIL', rsum2, rref2
   endif
   if (isum2 == iref2) then
      print *, 'Global sum irange array5d integer: PASS'
   else
      print *, 'Global sum irange array5d integer: FAIL', isum2, iref2
   endif
   if (lsum2 == lref2) then
      print *, 'Global sum irange array5d integer(8): PASS'
   else
      print *, 'Global sum irange array5d integer(8): FAIL', lsum2, lref2
   endif
   if (dsum3 == dref2) then
      print *, 'Global sum prod array5d double: PASS'
   else
      print *, 'Global sum prod array5d double: FAIL', dsum3, dref2
   endif
   if (rsum3 == rref2) then
      print *, 'Global sum prod array5d real: PASS'
   else
      print *, 'Global sum prod array5d real: FAIL', rsum3, rref2
   endif
   if (isum3 == iref2) then
      print *, 'Global sum prod array5d integer: PASS'
   else
      print *, 'Global sum prod array5d integer: FAIL', isum3, iref2
   endif
   if (lsum3 == lref2) then
      print *, 'Global sum prod array5d integer(8): PASS'
   else
      print *, 'Global sum prod array5d integer(8): FAIL', lsum3, lref2
   endif
   if (dsum4 == dref2) then
      print *, 'Global sum prod irange array5d double: PASS'
   else
      print *, 'Global sum prod irange array5d double: FAIL', dsum4, dref2
   endif
   if (rsum4 == rref2) then
      print *, 'Global sum prod irange array5d real: PASS'
   else
      print *, 'Global sum prod irange array5d real: FAIL', rsum4, rref2
   endif
   if (isum4 == iref2) then
      print *, 'Global sum prod irange array5d integer: PASS'
   else
      print *, 'Global sum prod irange array5d integer: FAIL', isum4, iref2
   endif
   if (lsum4 == lref2) then
      print *, 'Global sum prod irange array5d integer(8): PASS'
   else
      print *, 'Global sum prod irange array5d integer(8): FAIL', lsum4, lref2
   endif
   deallocate(darray5d, dmask5d, rarray5d, rmask5d, &
              iarray5d, imask5d, larray5d, lmask5d)

   !*** 6-d arrays
   if (myRank == masterTask) print *,'Testing 6d sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray6d(nx,ny,nz,nl,nm,nn), dmask6d(nx,ny,nz,nl,nm,nn), &
            rarray6d(nx,ny,nz,nl,nm,nn), rmask6d(nx,ny,nz,nl,nm,nn), &
            iarray6d(nx,ny,nz,nl,nm,nn), imask6d(nx,ny,nz,nl,nm,nn), &
            larray6d(nx,ny,nz,nl,nm,nn), lmask6d(nx,ny,nz,nl,nm,nn))
   dref = 0.d0
   rref = 0.0
   iref = 0
   lref = 0
   dref2 = 0.d0
   rref2 = 0.0
   iref2 = 0
   lref2 = 0
   do n=1,nn
   do m=1,nm
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      lscalar = (i+j+k+l+m+n)*(myRank+1)
      darray6d(i,j,k,l,m,n) = real(lscalar,R8)
      rarray6d(i,j,k,l,m,n) = real(lscalar)
      larray6d(i,j,k,l,m,n) = lscalar
      iarray6d(i,j,k,l,m,n) = min(lscalar,maxint)
      do r=0,nRanks-1
         dref = dref + real((i+j+k+l+m+n)*(r+1),R8)
         lref = lref + (i+j+k+l+m+n)*(r+1)
      end do
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6) .and. &
          l >= indxRange(7) .and. l <= indxRange(8) .and. &
          m >= indxRange(9) .and. m <= indxRange(10) .and. &
          n >= indxRange(11) .and. n <= indxRange(12)) then
         dmask6d(i,j,k,l,m,n) = 1.d0
         rmask6d(i,j,k,l,m,n) = 1.0
         lmask6d(i,j,k,l,m,n) = 1
         imask6d(i,j,k,l,m,n) = 1
         do r=0,nRanks-1
            dref2 = dref2 + real((i+j+k+l+m+n)*(r+1),R8)
            lref2 = lref2 + (i+j+k+l+m+n)*(r+1)
         end do
      else
         dmask6d(i,j,k,l,m,n) = 0.d0
         rmask6d(i,j,k,l,m,n) = 0.0
         lmask6d(i,j,k,l,m,n) = 0
         imask6d(i,j,k,l,m,n) = 0
      endif
   end do
   end do
   end do
   end do
   end do
   end do
   rref  = dref
   rref2 = dref2
   iref  = lref
   iref2 = lref2
   dsum1 = mpas_global_sum(darray6d, MPI_COMM_WORLD)
   rsum1 = mpas_global_sum(rarray6d, MPI_COMM_WORLD)
   isum1 = mpas_global_sum(iarray6d, MPI_COMM_WORLD)
   lsum1 = mpas_global_sum(larray6d, MPI_COMM_WORLD)
   dsum2 = mpas_global_sum(darray6d, MPI_COMM_WORLD, indxRange)
   rsum2 = mpas_global_sum(rarray6d, MPI_COMM_WORLD, indxRange)
   isum2 = mpas_global_sum(iarray6d, MPI_COMM_WORLD, indxRange)
   lsum2 = mpas_global_sum(larray6d, MPI_COMM_WORLD, indxRange)
   dsum3 = mpas_global_sum(darray6d, dmask6d, MPI_COMM_WORLD)
   rsum3 = mpas_global_sum(rarray6d, rmask6d, MPI_COMM_WORLD)
   isum3 = mpas_global_sum(iarray6d, imask6d, MPI_COMM_WORLD)
   lsum3 = mpas_global_sum(larray6d, lmask6d, MPI_COMM_WORLD)
   dsum4 = mpas_global_sum(darray6d, dmask6d, MPI_COMM_WORLD, indxRange)
   rsum4 = mpas_global_sum(rarray6d, rmask6d, MPI_COMM_WORLD, indxRange)
   isum4 = mpas_global_sum(iarray6d, imask6d, MPI_COMM_WORLD, indxRange)
   lsum4 = mpas_global_sum(larray6d, lmask6d, MPI_COMM_WORLD, indxRange)
   if (dsum1 == dref) then
      print *, 'Global sum array6d double: PASS'
   else
      print *, 'Global sum array6d double: FAIL', dsum1, dref
   endif
   if (rsum1 == rref) then
      print *, 'Global sum array6d real: PASS'
   else
      print *, 'Global sum array6d real: FAIL', rsum1, rref
   endif
   if (isum1 == iref) then
      print *, 'Global sum array6d integer: PASS'
   else
      print *, 'Global sum array6d integer: FAIL', isum1, iref
   endif
   if (lsum1 == lref) then
      print *, 'Global sum array6d integer(8): PASS'
   else
      print *, 'Global sum array6d integer(8): FAIL', lsum1, lref
   endif
   if (dsum2 == dref2) then
      print *, 'Global sum irange array6d double: PASS'
   else
      print *, 'Global sum irange array6d double: FAIL', dsum2, dref2
   endif
   if (rsum2 == rref2) then
      print *, 'Global sum irange array6d real: PASS'
   else
      print *, 'Global sum irange array6d real: FAIL', rsum2, rref2
   endif
   if (isum2 == iref2) then
      print *, 'Global sum irange array6d integer: PASS'
   else
      print *, 'Global sum irange array6d integer: FAIL', isum2, iref2
   endif
   if (lsum2 == lref2) then
      print *, 'Global sum irange array6d integer(8): PASS'
   else
      print *, 'Global sum irange array6d integer(8): FAIL', lsum2, lref2
   endif
   if (dsum3 == dref2) then
      print *, 'Global sum prod array6d double: PASS'
   else
      print *, 'Global sum prod array6d double: FAIL', dsum3, dref2
   endif
   if (rsum3 == rref2) then
      print *, 'Global sum prod array6d real: PASS'
   else
      print *, 'Global sum prod array6d real: FAIL', rsum3, rref2
   endif
   if (isum3 == iref2) then
      print *, 'Global sum prod array6d integer: PASS'
   else
      print *, 'Global sum prod array6d integer: FAIL', isum3, iref2
   endif
   if (lsum3 == lref2) then
      print *, 'Global sum prod array6d integer(8): PASS'
   else
      print *, 'Global sum prod array6d integer(8): FAIL', lsum3, lref2
   endif
   if (dsum4 == dref2) then
      print *, 'Global sum prod irange array6d double: PASS'
   else
      print *, 'Global sum prod irange array6d double: FAIL', dsum4, dref2
   endif
   if (rsum4 == rref2) then
      print *, 'Global sum prod irange array6d real: PASS'
   else
      print *, 'Global sum prod irange array6d real: FAIL', rsum4, rref2
   endif
   if (isum4 == iref2) then
      print *, 'Global sum prod irange array6d integer: PASS'
   else
      print *, 'Global sum prod irange array6d integer: FAIL', isum4, iref2
   endif
   if (lsum4 == lref2) then
      print *, 'Global sum prod irange array6d integer(8): PASS'
   else
      print *, 'Global sum prod irange array6d integer(8): FAIL', lsum4, lref2
   endif
   deallocate(darray6d, dmask6d, rarray6d, rmask6d, &
              iarray6d, imask6d, larray6d, lmask6d)

   if (myRank == masterTask) print *,'Testing basic sums complete '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   !--------------------------------------------------------------------
   ! Test multi-field interfaces and check for reproducibility by
   ! using a different ordering for each of three fields.
   ! Fill the arrays with random values that span the full range of a
   ! datatype to induce roundoff errors during the sums.
   !--------------------------------------------------------------------

   !*** Scalars
   if (myRank == masterTask) print *,'Testing scalar multi-field sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(dsum1n(nFields), dsum2n(nFields), dsum3n(nFields), dsum4n(nFields), &
            rsum1n(nFields), rsum2n(nFields), rsum3n(nFields), rsum4n(nFields), &
            isum1n(nFields), isum2n(nFields), isum3n(nFields), isum4n(nFields), &
            lsum1n(nFields), lsum2n(nFields), lsum3n(nFields), lsum4n(nFields))
   allocate(darray0n(nFields), rarray0n(nFields), & 
            iarray0n(nFields), larray0n(nFields)) 
   allocate (randoms(nRanks))
   drange = 15
   rrange = 7
   irange = huge(iref)/100
   lrange = huge(lref)/100

   ! create a random number for each MPI rank and broadcast so we
   ! can check reproducibility by using the same numbers in different
   ! orders
   call random_number(randoms)
   call MPI_Bcast(randoms, nRanks, MPI_REAL, masterTask, MPI_COMM_WORLD, ierr)

   rscalar = randoms(myRank+1)
   ipower = int(rscalar*drange) - 5
   darray0n(1) = real(rscalar,R8)*(10.d0**ipower)
   ipower = int(rscalar*rrange) - 5
   rarray0n(1) = rscalar*(10.0**ipower)
   iarray0n(1) = nint(rscalar*irange)
   larray0n(1) = nint(rscalar*lrange)

   ! Second field values are in inverse order
   rscalar = randoms(nRanks-myRank)
   ipower = int(rscalar*drange) - 5
   darray0n(2) = real(rscalar,R8)*(10.d0**ipower)
   ipower = int(rscalar*rrange) - 5
   rarray0n(2) = rscalar*(10.0**ipower)
   iarray0n(2) = nint(rscalar*irange)
   larray0n(2) = nint(rscalar*lrange)

   ! Third field has values another order
   iref = mod(nRanks/2+myRank ,nRanks) + 1
   rscalar = randoms(iref)
   ipower = int(rscalar*drange) - 5
   darray0n(3) = real(rscalar,R8)*(10.d0**ipower)
   ipower = int(rscalar*rrange) - 5
   rarray0n(3) = rscalar*(10.0**ipower)
   iarray0n(3) = nint(rscalar*irange)
   larray0n(3) = nint(rscalar*lrange)

   dsum1n(:) = 0.d0
   rsum1n(:) = 0.0
   isum1n(:) = 0
   lsum1n(:) = 0
   do r=0,nRanks-1
      iref = r+1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      dsum1n(1) = dsum1n(1) + real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rsum1n(1) = rsum1n(1) + rscalar*(10.0**ipower)
      isum1n(1) = isum1n(1) + nint(rscalar*irange)
      lsum1n(1) = lsum1n(1) + nint(rscalar*lrange)

      iref = nRanks - r
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      dsum1n(2) = dsum1n(2) + real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rsum1n(2) = rsum1n(2) + rscalar*(10.0**ipower)
      isum1n(2) = isum1n(2) + nint(rscalar*irange)
      lsum1n(2) = lsum1n(2) + nint(rscalar*lrange)

      iref = mod(nRanks/2+r ,nRanks) + 1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      dsum1n(3) = dsum1n(3) + real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rsum1n(3) = rsum1n(3) + rscalar*(10.0**ipower)
      isum1n(3) = isum1n(3) + nint(rscalar*irange)
      lsum1n(3) = lsum1n(3) + nint(rscalar*lrange)

   end do
   ! Compute a reference sum using the single field calls from above
   dscalar = darray0n(1)
   rscalar = rarray0n(1)
   iscalar = iarray0n(1)
   lscalar = larray0n(1)
   dref = mpas_global_sum(dscalar, MPI_COMM_WORLD)
   rref = mpas_global_sum(rscalar, MPI_COMM_WORLD)
   iref = mpas_global_sum(iscalar, MPI_COMM_WORLD)
   lref = mpas_global_sum(lscalar, MPI_COMM_WORLD)

   ! Check that the reference sums for floats and doubles
   ! are dependent on order of summation and then ints/longs are not
   if (lsum1n(1) /= lref .or. lsum1n(2) /= lref .or. &
       lsum1n(3) /= lref) then
      print *, 'Multi-field setup scalar integer(8): FAIL', lref, lsum1n
   endif
   if (isum1n(1) /= iref .or. isum1n(2) /= iref .or. &
       isum1n(3) /= iref) then
      print *, 'Multi-field setup scalar integer: FAIL', iref, isum1n
   endif
   if (dsum1n(1) == dsum1n(2) .and. dsum1n(1) == dsum1n(3) .and. &
       dsum1n(2) == dsum1n(3)) then
       print *, 'Multi-field double scalar not sensitive to roundoff : WARN', dref, dsum1n
       print *, 'Multi-field double scalar data:', myRank, darray0n
   endif
   if (rsum1n(1) == rsum1n(2) .and. rsum1n(1) == rsum1n(3) .and. &
       rsum1n(2) == rsum1n(3)) then
       print *, 'Multi-field real scalar not sensitive to roundoff : WARN', rref, rsum1n
       print *, 'Multi-field real scalar data:', myRank, rarray0n
   endif

   ! Now do actual multi-field and reproducibility tests
   dsum1n = mpas_global_sum_nfld(darray0n, MPI_COMM_WORLD)
   rsum1n = mpas_global_sum_nfld(rarray0n, MPI_COMM_WORLD)
   isum1n = mpas_global_sum_nfld(iarray0n, MPI_COMM_WORLD)
   lsum1n = mpas_global_sum_nfld(larray0n, MPI_COMM_WORLD)
   if (dsum1n(1) == dref .and. dsum1n(2) == dref .and. &
       dsum1n(3) == dref) then
      print *, 'Reproducibility and nfld scalar double: PASS'
   else
      print *, 'Reproducibility and nfld scalar double: FAIL', dref, dsum1n
   endif
   if (rsum1n(1) == rref .and. rsum1n(2) == rref .and. &
       rsum1n(3) == rref) then
      print *, 'Reproducibility and nfld scalar real: PASS'
   else
      print *, 'Reproducibility and nfld scalar real: FAIL', rref, rsum1n
   endif
   if (isum1n(1) == iref .and. isum1n(2) == iref .and. &
       isum1n(3) == iref) then
      print *, 'Reproducibility and nfld scalar integer: PASS'
   else
      print *, 'Reproducibility and nfld scalar integer: FAIL', iref, isum1n
   endif
   if (lsum1n(1) == lref .and. lsum1n(2) == lref .and. &
       lsum1n(3) == lref) then
      print *, 'Reproducibility and nfld scalar integer(8): PASS'
   else
      print *, 'Reproducibility and nfld scalar integer(8): FAIL', lref, lsum1n
   endif
   deallocate(darray0n, rarray0n, iarray0n, larray0n)
   deallocate(randoms)

   !*** 1-d arrays
   if (myRank == masterTask) print *,'Testing 1d multi-field sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray1n(nx,nFields), dmask1n(nx,nFields), &
            rarray1n(nx,nFields), rmask1n(nx,nFields), &
            iarray1n(nx,nFields), imask1n(nx,nFields), &
            larray1n(nx,nFields), lmask1n(nx,nFields))
   allocate(randoms(nx))
   call random_number(randoms)

   iref = 0
   do i=1,nx
      iref = iref + 1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      darray1n(i,1) = real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rarray1n(i,1) = rscalar*(10.0**ipower)
      iarray1n(i,1) = nint(rscalar*irange)
      larray1n(i,1) = nint(rscalar*lrange)
   end do
   do i=1,nx
      i2 = nx+1-i
      i3 = mod(nx/2+i,nx) + 1
      darray1n(i,2) = darray1n(i2,1)
      darray1n(i,3) = darray1n(i3,1)
      rarray1n(i,2) = rarray1n(i2,1)
      rarray1n(i,3) = rarray1n(i3,1)
      iarray1n(i,2) = iarray1n(i2,1)
      iarray1n(i,3) = iarray1n(i3,1)
      larray1n(i,2) = larray1n(i2,1)
      larray1n(i,3) = larray1n(i3,1)
   end do

   do nfld=1,nFields
   do i=1,nx
      if (i >= indxRange(1) .and. i <= indxRange(2)) then
         dmask1n(i,nfld) = 1.d0
         rmask1n(i,nfld) = 1.0
         lmask1n(i,nfld) = 1
         imask1n(i,nfld) = 1
      else
         dmask1n(i,nfld) = 0.d0
         rmask1n(i,nfld) = 0.0
         lmask1n(i,nfld) = 0
         imask1n(i,nfld) = 0
      endif
   end do
   end do
   ! Use the single-field interface to compute reference sums
   dref = mpas_global_sum(darray1n(:,1), MPI_COMM_WORLD)
   rref = mpas_global_sum(rarray1n(:,1), MPI_COMM_WORLD)
   iref = mpas_global_sum(iarray1n(:,1), MPI_COMM_WORLD)
   lref = mpas_global_sum(larray1n(:,1), MPI_COMM_WORLD)

   ! Now compute the various multi-field sums
   dsum1n = mpas_global_sum_nfld(darray1n, MPI_COMM_WORLD)
   rsum1n = mpas_global_sum_nfld(rarray1n, MPI_COMM_WORLD)
   isum1n = mpas_global_sum_nfld(iarray1n, MPI_COMM_WORLD)
   lsum1n = mpas_global_sum_nfld(larray1n, MPI_COMM_WORLD)
   dsum2n = mpas_global_sum_nfld(darray1n, MPI_COMM_WORLD, indxRange)
   rsum2n = mpas_global_sum_nfld(rarray1n, MPI_COMM_WORLD, indxRange)
   isum2n = mpas_global_sum_nfld(iarray1n, MPI_COMM_WORLD, indxRange)
   lsum2n = mpas_global_sum_nfld(larray1n, MPI_COMM_WORLD, indxRange)
   dsum3n = mpas_global_sum_nfld(darray1n, dmask1n, MPI_COMM_WORLD)
   rsum3n = mpas_global_sum_nfld(rarray1n, rmask1n, MPI_COMM_WORLD)
   isum3n = mpas_global_sum_nfld(iarray1n, imask1n, MPI_COMM_WORLD)
   lsum3n = mpas_global_sum_nfld(larray1n, lmask1n, MPI_COMM_WORLD)
   dsum4n = mpas_global_sum_nfld(darray1n, dmask1n, MPI_COMM_WORLD, indxRange)
   rsum4n = mpas_global_sum_nfld(rarray1n, rmask1n, MPI_COMM_WORLD, indxRange)
   isum4n = mpas_global_sum_nfld(iarray1n, imask1n, MPI_COMM_WORLD, indxRange)
   lsum4n = mpas_global_sum_nfld(larray1n, lmask1n, MPI_COMM_WORLD, indxRange)

   if (dsum1n(1) == dref .and. dsum1n(2) == dref .and. &
       dsum1n(3) == dref) then
      print *, 'Multi-field and reproducibility 1d double: PASS'
   else
      print *, 'Multi-field and reproducibility 1d double: FAIL', dref, dsum1n
   endif
   if (rsum1n(1) == rref .and. rsum1n(2) == rref .and. &
       rsum1n(3) == rref) then
      print *, 'Multi-field and reproducibility 1d real: PASS'
   else
      print *, 'Multi-field and reproducibility 1d real: FAIL', rref, rsum1n
   endif
   if (isum1n(1) == iref .and. isum1n(2) == iref .and. &
       isum1n(3) == iref) then
      print *, 'Multi-field and reproducibility 1d int: PASS'
   else
      print *, 'Multi-field and reproducibility 1d int: FAIL', iref, isum1n
   endif
   if (lsum1n(1) == lref .and. lsum1n(2) == lref .and. &
       lsum1n(3) == lref) then
      print *, 'Multi-field and reproducibility 1d int(8): PASS'
   else
      print *, 'Multi-field and reproducibility 1d int(8): FAIL', lref, lsum1n
   endif
   isum1 = 0
   isum2 = 0
   isum3 = 0
   isum4 = 0
   do nfld=1,nFields
      if (dsum2n(nfld) /= dsum3n(nfld)) isum1 = isum1+1
      if (dsum2n(nfld) /= dsum4n(nfld)) isum1 = isum1+1
      if (rsum2n(nfld) /= rsum3n(nfld)) isum2 = isum2+1
      if (rsum2n(nfld) /= rsum4n(nfld)) isum2 = isum2+1
      if (isum2n(nfld) /= isum3n(nfld)) isum3 = isum3+1
      if (isum2n(nfld) /= isum4n(nfld)) isum3 = isum3+1
      if (lsum2n(nfld) /= lsum3n(nfld)) isum4 = isum4+1
      if (lsum2n(nfld) /= lsum4n(nfld)) isum4 = isum4+1
   end do
   if (isum1 > 0) then
      print *, 'Multifield range/mask 1d double: FAIL', dsum2n,dsum3n,dsum4n
   else
      print *, 'Multifield range/mask 1d double: PASS'
   endif
   if (isum2 > 0) then
      print *, 'Multifield range/mask 1d real: FAIL', rsum2n,rsum3n,rsum4n
   else
      print *, 'Multifield range/mask 1d real: PASS'
   endif
   if (isum3 > 0) then
      print *, 'Multifield range/mask 1d int: FAIL', isum2n,isum3n,isum4n
   else
      print *, 'Multifield range/mask 1d int: PASS'
   endif
   if (isum4 > 0) then
      print *, 'Multifield range/mask 1d int(8): FAIL', lsum2n,lsum3n,lsum4n
   else
      print *, 'Multifield range/mask 1d int(8): PASS'
   endif
   deallocate(darray1n, dmask1n, rarray1n, rmask1n, &
              iarray1n, imask1n, larray1n, lmask1n)
   deallocate(randoms)

   !*** 2-d arrays
   if (myRank == masterTask) print *,'Testing 2d multi-field sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray2n(nx,ny,nFields), dmask2n(nx,ny,nFields), &
            rarray2n(nx,ny,nFields), rmask2n(nx,ny,nFields), &
            iarray2n(nx,ny,nFields), imask2n(nx,ny,nFields), &
            larray2n(nx,ny,nFields), lmask2n(nx,ny,nFields))
   allocate(randoms(nx*ny))
   call random_number(randoms)

   iref = 0
   do j=1,ny
   do i=1,nx
      iref = iref + 1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      darray2n(i,j,1) = real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rarray2n(i,j,1) = rscalar*(10.0**ipower)
      iarray2n(i,j,1) = nint(rscalar*irange)
      larray2n(i,j,1) = nint(rscalar*lrange)
   end do
   end do
   do j=1,ny
      j2 = ny+1-j
      j3 = mod(ny/2+j,ny) + 1
      do i=1,nx
         i2 = nx+1-i
         i3 = mod(nx/2+i,nx) + 1
         darray2n(i,j,2) = darray2n(i2,j2,1)
         darray2n(i,j,3) = darray2n(i3,j3,1)
         rarray2n(i,j,2) = rarray2n(i2,j2,1)
         rarray2n(i,j,3) = rarray2n(i3,j3,1)
         iarray2n(i,j,2) = iarray2n(i2,j2,1)
         iarray2n(i,j,3) = iarray2n(i3,j3,1)
         larray2n(i,j,2) = larray2n(i2,j2,1)
         larray2n(i,j,3) = larray2n(i3,j3,1)
      end do
   end do

   do nfld=1,nFields
   do j=1,ny
   do i=1,nx
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4)) then
         dmask2n(i,j,nfld) = 1.d0
         rmask2n(i,j,nfld) = 1.0
         lmask2n(i,j,nfld) = 1
         imask2n(i,j,nfld) = 1
      else
         dmask2n(i,j,nfld) = 0.d0
         rmask2n(i,j,nfld) = 0.0
         lmask2n(i,j,nfld) = 0
         imask2n(i,j,nfld) = 0
      endif
   end do
   end do
   end do
   ! Use the single-field interface to compute reference sums
   dref = mpas_global_sum(darray2n(:,:,1), MPI_COMM_WORLD)
   rref = mpas_global_sum(rarray2n(:,:,1), MPI_COMM_WORLD)
   iref = mpas_global_sum(iarray2n(:,:,1), MPI_COMM_WORLD)
   lref = mpas_global_sum(larray2n(:,:,1), MPI_COMM_WORLD)

   ! Now compute the various multi-field sums
   dsum1n = mpas_global_sum_nfld(darray2n, MPI_COMM_WORLD)
   rsum1n = mpas_global_sum_nfld(rarray2n, MPI_COMM_WORLD)
   isum1n = mpas_global_sum_nfld(iarray2n, MPI_COMM_WORLD)
   lsum1n = mpas_global_sum_nfld(larray2n, MPI_COMM_WORLD)
   dsum2n = mpas_global_sum_nfld(darray2n, MPI_COMM_WORLD, indxRange)
   rsum2n = mpas_global_sum_nfld(rarray2n, MPI_COMM_WORLD, indxRange)
   isum2n = mpas_global_sum_nfld(iarray2n, MPI_COMM_WORLD, indxRange)
   lsum2n = mpas_global_sum_nfld(larray2n, MPI_COMM_WORLD, indxRange)
   dsum3n = mpas_global_sum_nfld(darray2n, dmask2n, MPI_COMM_WORLD)
   rsum3n = mpas_global_sum_nfld(rarray2n, rmask2n, MPI_COMM_WORLD)
   isum3n = mpas_global_sum_nfld(iarray2n, imask2n, MPI_COMM_WORLD)
   lsum3n = mpas_global_sum_nfld(larray2n, lmask2n, MPI_COMM_WORLD)
   dsum4n = mpas_global_sum_nfld(darray2n, dmask2n, MPI_COMM_WORLD, indxRange)
   rsum4n = mpas_global_sum_nfld(rarray2n, rmask2n, MPI_COMM_WORLD, indxRange)
   isum4n = mpas_global_sum_nfld(iarray2n, imask2n, MPI_COMM_WORLD, indxRange)
   lsum4n = mpas_global_sum_nfld(larray2n, lmask2n, MPI_COMM_WORLD, indxRange)

   if (dsum1n(1) == dref .and. dsum1n(2) == dref .and. &
       dsum1n(3) == dref) then
      print *, 'Multi-field and reproducibility 2d double: PASS'
   else
      print *, 'Multi-field and reproducibility 2d double: FAIL', dref, dsum1n
   endif
   if (rsum1n(1) == rref .and. rsum1n(2) == rref .and. &
       rsum1n(3) == rref) then
      print *, 'Multi-field and reproducibility 2d real: PASS'
   else
      print *, 'Multi-field and reproducibility 2d real: FAIL', rref, rsum1n
   endif
   if (isum1n(1) == iref .and. isum1n(2) == iref .and. &
       isum1n(3) == iref) then
      print *, 'Multi-field and reproducibility 2d int: PASS'
   else
      print *, 'Multi-field and reproducibility 2d int: FAIL', iref, isum1n
   endif
   if (lsum1n(1) == lref .and. lsum1n(2) == lref .and. &
       lsum1n(3) == lref) then
      print *, 'Multi-field and reproducibility 2d int(8): PASS'
   else
      print *, 'Multi-field and reproducibility 2d int(8): FAIL', lref, lsum1n
   endif
   isum1 = 0
   isum2 = 0
   isum3 = 0
   isum4 = 0
   do nfld=1,nFields
      if (dsum2n(nfld) /= dsum3n(nfld)) isum1 = isum1+1
      if (dsum2n(nfld) /= dsum4n(nfld)) isum1 = isum1+1
      if (rsum2n(nfld) /= rsum3n(nfld)) isum2 = isum2+1
      if (rsum2n(nfld) /= rsum4n(nfld)) isum2 = isum2+1
      if (isum2n(nfld) /= isum3n(nfld)) isum3 = isum3+1
      if (isum2n(nfld) /= isum4n(nfld)) isum3 = isum3+1
      if (lsum2n(nfld) /= lsum3n(nfld)) isum4 = isum4+1
      if (lsum2n(nfld) /= lsum4n(nfld)) isum4 = isum4+1
   end do
   if (isum1 > 0) then
      print *, 'Multifield range/mask 2d double: FAIL', dsum2n,dsum3n,dsum4n
   else
      print *, 'Multifield range/mask 2d double: PASS'
   endif
   if (isum2 > 0) then
      print *, 'Multifield range/mask 2d real: FAIL', rsum2n,rsum3n,rsum4n
   else
      print *, 'Multifield range/mask 2d real: PASS'
   endif
   if (isum3 > 0) then
      print *, 'Multifield range/mask 2d int: FAIL', isum2n,isum3n,isum4n
   else
      print *, 'Multifield range/mask 2d int: PASS'
   endif
   if (isum4 > 0) then
      print *, 'Multifield range/mask 2d int(8): FAIL', lsum2n,lsum3n,lsum4n
   else
      print *, 'Multifield range/mask 2d int(8): PASS'
   endif
   deallocate(darray2n, dmask2n, rarray2n, rmask2n, &
              iarray2n, imask2n, larray2n, lmask2n)
   deallocate(randoms)

   !*** 3-d arrays
   if (myRank == masterTask) print *,'Testing 3d multi-field sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray3n(nx,ny,nz,nFields), dmask3n(nx,ny,nz,nFields), &
            rarray3n(nx,ny,nz,nFields), rmask3n(nx,ny,nz,nFields), &
            iarray3n(nx,ny,nz,nFields), imask3n(nx,ny,nz,nFields), &
            larray3n(nx,ny,nz,nFields), lmask3n(nx,ny,nz,nFields))
   allocate(randoms(nx*ny*nz))
   call random_number(randoms)

   iref = 0
   do k=1,nz
   do j=1,ny
   do i=1,nx
      iref = iref + 1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      darray3n(i,j,k,1) = real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rarray3n(i,j,k,1) = rscalar*(10.0**ipower)
      iarray3n(i,j,k,1) = nint(rscalar*irange)
      larray3n(i,j,k,1) = nint(rscalar*lrange)
   end do
   end do
   end do
   do k=1,nz
      k2 = nz+1-k
      k3 = mod(nz/2+k,nz) + 1
      do j=1,ny
         j2 = ny+1-j
         j3 = mod(ny/2+j,ny) + 1
         do i=1,nx
            i2 = nx+1-i
            i3 = mod(nx/2+i,nx) + 1
            darray3n(i,j,k,2) = darray3n(i2,j2,k2,1)
            darray3n(i,j,k,3) = darray3n(i3,j3,k3,1)
            rarray3n(i,j,k,2) = rarray3n(i2,j2,k2,1)
            rarray3n(i,j,k,3) = rarray3n(i3,j3,k3,1)
            iarray3n(i,j,k,2) = iarray3n(i2,j2,k2,1)
            iarray3n(i,j,k,3) = iarray3n(i3,j3,k3,1)
            larray3n(i,j,k,2) = larray3n(i2,j2,k2,1)
            larray3n(i,j,k,3) = larray3n(i3,j3,k3,1)
         end do
      end do
   end do

   do nfld=1,nFields
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6)) then
         dmask3n(i,j,k,nfld) = 1.d0
         rmask3n(i,j,k,nfld) = 1.0
         lmask3n(i,j,k,nfld) = 1
         imask3n(i,j,k,nfld) = 1
      else
         dmask3n(i,j,k,nfld) = 0.d0
         rmask3n(i,j,k,nfld) = 0.0
         lmask3n(i,j,k,nfld) = 0
         imask3n(i,j,k,nfld) = 0
      endif
   end do
   end do
   end do
   end do
   ! Use the single-field interface to compute reference sums
   dref = mpas_global_sum(darray3n(:,:,:,1), MPI_COMM_WORLD)
   rref = mpas_global_sum(rarray3n(:,:,:,1), MPI_COMM_WORLD)
   iref = mpas_global_sum(iarray3n(:,:,:,1), MPI_COMM_WORLD)
   lref = mpas_global_sum(larray3n(:,:,:,1), MPI_COMM_WORLD)

   ! Now compute the various multi-field sums
   dsum1n = mpas_global_sum_nfld(darray3n, MPI_COMM_WORLD)
   rsum1n = mpas_global_sum_nfld(rarray3n, MPI_COMM_WORLD)
   isum1n = mpas_global_sum_nfld(iarray3n, MPI_COMM_WORLD)
   lsum1n = mpas_global_sum_nfld(larray3n, MPI_COMM_WORLD)
   dsum2n = mpas_global_sum_nfld(darray3n, MPI_COMM_WORLD, indxRange)
   rsum2n = mpas_global_sum_nfld(rarray3n, MPI_COMM_WORLD, indxRange)
   isum2n = mpas_global_sum_nfld(iarray3n, MPI_COMM_WORLD, indxRange)
   lsum2n = mpas_global_sum_nfld(larray3n, MPI_COMM_WORLD, indxRange)
   dsum3n = mpas_global_sum_nfld(darray3n, dmask3n, MPI_COMM_WORLD)
   rsum3n = mpas_global_sum_nfld(rarray3n, rmask3n, MPI_COMM_WORLD)
   isum3n = mpas_global_sum_nfld(iarray3n, imask3n, MPI_COMM_WORLD)
   lsum3n = mpas_global_sum_nfld(larray3n, lmask3n, MPI_COMM_WORLD)
   dsum4n = mpas_global_sum_nfld(darray3n, dmask3n, MPI_COMM_WORLD, indxRange)
   rsum4n = mpas_global_sum_nfld(rarray3n, rmask3n, MPI_COMM_WORLD, indxRange)
   isum4n = mpas_global_sum_nfld(iarray3n, imask3n, MPI_COMM_WORLD, indxRange)
   lsum4n = mpas_global_sum_nfld(larray3n, lmask3n, MPI_COMM_WORLD, indxRange)

   if (dsum1n(1) == dref .and. dsum1n(2) == dref .and. &
       dsum1n(3) == dref) then
      print *, 'Multi-field and reproducibility 3d double: PASS'
   else
      print *, 'Multi-field and reproducibility 3d double: FAIL', dref, dsum1n
   endif
   if (rsum1n(1) == rref .and. rsum1n(2) == rref .and. &
       rsum1n(3) == rref) then
      print *, 'Multi-field and reproducibility 3d real: PASS'
   else
      print *, 'Multi-field and reproducibility 3d real: FAIL', rref, rsum1n
   endif
   if (isum1n(1) == iref .and. isum1n(2) == iref .and. &
       isum1n(3) == iref) then
      print *, 'Multi-field and reproducibility 3d int: PASS'
   else
      print *, 'Multi-field and reproducibility 3d int: FAIL', iref, isum1n
   endif
   if (lsum1n(1) == lref .and. lsum1n(2) == lref .and. &
       lsum1n(3) == lref) then
      print *, 'Multi-field and reproducibility 3d int(8): PASS'
   else
      print *, 'Multi-field and reproducibility 3d int(8): FAIL', lref, lsum1n
   endif
   isum1 = 0
   isum2 = 0
   isum3 = 0
   isum4 = 0
   do nfld=1,nFields
      if (dsum2n(nfld) /= dsum3n(nfld)) isum1 = isum1+1
      if (dsum2n(nfld) /= dsum4n(nfld)) isum1 = isum1+1
      if (rsum2n(nfld) /= rsum3n(nfld)) isum2 = isum2+1
      if (rsum2n(nfld) /= rsum4n(nfld)) isum2 = isum2+1
      if (isum2n(nfld) /= isum3n(nfld)) isum3 = isum3+1
      if (isum2n(nfld) /= isum4n(nfld)) isum3 = isum3+1
      if (lsum2n(nfld) /= lsum3n(nfld)) isum4 = isum4+1
      if (lsum2n(nfld) /= lsum4n(nfld)) isum4 = isum4+1
   end do
   if (isum1 > 0) then
      print *, 'Multifield range/mask 3d double: FAIL', dsum2n,dsum3n,dsum4n
   else
      print *, 'Multifield range/mask 3d double: PASS'
   endif
   if (isum2 > 0) then
      print *, 'Multifield range/mask 3d real: FAIL', rsum2n,rsum3n,rsum4n
   else
      print *, 'Multifield range/mask 3d real: PASS'
   endif
   if (isum3 > 0) then
      print *, 'Multifield range/mask 3d int: FAIL', isum2n,isum3n,isum4n
   else
      print *, 'Multifield range/mask 3d int: PASS'
   endif
   if (isum4 > 0) then
      print *, 'Multifield range/mask 3d int(8): FAIL', lsum2n,lsum3n,lsum4n
   else
      print *, 'Multifield range/mask 3d int(8): PASS'
   endif
   deallocate(darray3n, dmask3n, rarray3n, rmask3n, &
              iarray3n, imask3n, larray3n, lmask3n)
   deallocate(randoms)

   !*** 4-d arrays
   if (myRank == masterTask) print *,'Testing 4d multi-field sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray4n(nx,ny,nz,nl,nFields), &
            rarray4n(nx,ny,nz,nl,nFields), &
            iarray4n(nx,ny,nz,nl,nFields), &
            larray4n(nx,ny,nz,nl,nFields), &
             dmask4n(nx,ny,nz,nl,nFields), &
             rmask4n(nx,ny,nz,nl,nFields), &
             imask4n(nx,ny,nz,nl,nFields), &
             lmask4n(nx,ny,nz,nl,nFields))
   allocate(randoms(nx*ny*nz*nl))
   call random_number(randoms)

   iref = 0
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      iref = iref + 1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      darray4n(i,j,k,l,1) = real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rarray4n(i,j,k,l,1) = rscalar*(10.0**ipower)
      iarray4n(i,j,k,l,1) = nint(rscalar*irange)
      larray4n(i,j,k,l,1) = nint(rscalar*lrange)
   end do
   end do
   end do
   end do
   do l=1,nl
      l2 = nl+1-l
      l3 = mod(nl/2+l,nl) + 1
      do k=1,nz
         k2 = nz+1-k
         k3 = mod(nz/2+k,nz) + 1
         do j=1,ny
            j2 = ny+1-j
            j3 = mod(ny/2+j,ny) + 1
            do i=1,nx
               i2 = nx+1-i
               i3 = mod(nx/2+i,nx) + 1
               darray4n(i,j,k,l,2) = darray4n(i2,j2,k2,l2,1)
               darray4n(i,j,k,l,3) = darray4n(i3,j3,k3,l3,1)
               rarray4n(i,j,k,l,2) = rarray4n(i2,j2,k2,l2,1)
               rarray4n(i,j,k,l,3) = rarray4n(i3,j3,k3,l3,1)
               iarray4n(i,j,k,l,2) = iarray4n(i2,j2,k2,l2,1)
               iarray4n(i,j,k,l,3) = iarray4n(i3,j3,k3,l3,1)
               larray4n(i,j,k,l,2) = larray4n(i2,j2,k2,l2,1)
               larray4n(i,j,k,l,3) = larray4n(i3,j3,k3,l3,1)
            end do
         end do
      end do
   end do

   do nfld=1,nFields
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6) .and. &
          l >= indxRange(7) .and. l <= indxRange(8)) then
         dmask4n(i,j,k,l,nfld) = 1.d0
         rmask4n(i,j,k,l,nfld) = 1.0
         lmask4n(i,j,k,l,nfld) = 1
         imask4n(i,j,k,l,nfld) = 1
      else
         dmask4n(i,j,k,l,nfld) = 0.d0
         rmask4n(i,j,k,l,nfld) = 0.0
         lmask4n(i,j,k,l,nfld) = 0
         imask4n(i,j,k,l,nfld) = 0
      endif
   end do
   end do
   end do
   end do
   end do
   ! Use the single-field interface to compute reference sums
   dref = mpas_global_sum(darray4n(:,:,:,:,1), MPI_COMM_WORLD)
   rref = mpas_global_sum(rarray4n(:,:,:,:,1), MPI_COMM_WORLD)
   iref = mpas_global_sum(iarray4n(:,:,:,:,1), MPI_COMM_WORLD)
   lref = mpas_global_sum(larray4n(:,:,:,:,1), MPI_COMM_WORLD)

   ! Now compute the various multi-field sums
   dsum1n = mpas_global_sum_nfld(darray4n, MPI_COMM_WORLD)
   rsum1n = mpas_global_sum_nfld(rarray4n, MPI_COMM_WORLD)
   isum1n = mpas_global_sum_nfld(iarray4n, MPI_COMM_WORLD)
   lsum1n = mpas_global_sum_nfld(larray4n, MPI_COMM_WORLD)
   dsum2n = mpas_global_sum_nfld(darray4n, MPI_COMM_WORLD, indxRange)
   rsum2n = mpas_global_sum_nfld(rarray4n, MPI_COMM_WORLD, indxRange)
   isum2n = mpas_global_sum_nfld(iarray4n, MPI_COMM_WORLD, indxRange)
   lsum2n = mpas_global_sum_nfld(larray4n, MPI_COMM_WORLD, indxRange)
   dsum3n = mpas_global_sum_nfld(darray4n, dmask4n, MPI_COMM_WORLD)
   rsum3n = mpas_global_sum_nfld(rarray4n, rmask4n, MPI_COMM_WORLD)
   isum3n = mpas_global_sum_nfld(iarray4n, imask4n, MPI_COMM_WORLD)
   lsum3n = mpas_global_sum_nfld(larray4n, lmask4n, MPI_COMM_WORLD)
   dsum4n = mpas_global_sum_nfld(darray4n, dmask4n, MPI_COMM_WORLD, indxRange)
   rsum4n = mpas_global_sum_nfld(rarray4n, rmask4n, MPI_COMM_WORLD, indxRange)
   isum4n = mpas_global_sum_nfld(iarray4n, imask4n, MPI_COMM_WORLD, indxRange)
   lsum4n = mpas_global_sum_nfld(larray4n, lmask4n, MPI_COMM_WORLD, indxRange)

   if (dsum1n(1) == dref .and. dsum1n(2) == dref .and. &
       dsum1n(3) == dref) then
      print *, 'Multi-field and reproducibility 4d double: PASS'
   else
      print *, 'Multi-field and reproducibility 4d double: FAIL', dref, dsum1n
   endif
   if (rsum1n(1) == rref .and. rsum1n(2) == rref .and. &
       rsum1n(3) == rref) then
      print *, 'Multi-field and reproducibility 4d real: PASS'
   else
      print *, 'Multi-field and reproducibility 4d real: FAIL', rref, rsum1n
   endif
   if (isum1n(1) == iref .and. isum1n(2) == iref .and. &
       isum1n(3) == iref) then
      print *, 'Multi-field and reproducibility 4d int: PASS'
   else
      print *, 'Multi-field and reproducibility 4d int: FAIL', iref, isum1n
   endif
   if (lsum1n(1) == lref .and. lsum1n(2) == lref .and. &
       lsum1n(3) == lref) then
      print *, 'Multi-field and reproducibility 4d int(8): PASS'
   else
      print *, 'Multi-field and reproducibility 4d int(8): FAIL', lref, lsum1n
   endif
   isum1 = 0
   isum2 = 0
   isum3 = 0
   isum4 = 0
   do nfld=1,nFields
      if (dsum2n(nfld) /= dsum3n(nfld)) isum1 = isum1+1
      if (dsum2n(nfld) /= dsum4n(nfld)) isum1 = isum1+1
      if (rsum2n(nfld) /= rsum3n(nfld)) isum2 = isum2+1
      if (rsum2n(nfld) /= rsum4n(nfld)) isum2 = isum2+1
      if (isum2n(nfld) /= isum3n(nfld)) isum3 = isum3+1
      if (isum2n(nfld) /= isum4n(nfld)) isum3 = isum3+1
      if (lsum2n(nfld) /= lsum3n(nfld)) isum4 = isum4+1
      if (lsum2n(nfld) /= lsum4n(nfld)) isum4 = isum4+1
   end do
   if (isum1 > 0) then
      print *, 'Multifield range/mask 4d double: FAIL', dsum2n,dsum3n,dsum4n
   else
      print *, 'Multifield range/mask 4d double: PASS'
   endif
   if (isum2 > 0) then
      print *, 'Multifield range/mask 4d real: FAIL', rsum2n,rsum3n,rsum4n
   else
      print *, 'Multifield range/mask 4d real: PASS'
   endif
   if (isum3 > 0) then
      print *, 'Multifield range/mask 4d int: FAIL', isum2n,isum3n,isum4n
   else
      print *, 'Multifield range/mask 4d int: PASS'
   endif
   if (isum4 > 0) then
      print *, 'Multifield range/mask 4d int(8): FAIL', lsum2n,lsum3n,lsum4n
   else
      print *, 'Multifield range/mask 4d int(8): PASS'
   endif
   deallocate(darray4n, dmask4n, rarray4n, rmask4n, &
              iarray4n, imask4n, larray4n, lmask4n)
   deallocate(randoms)

   !*** 5-d arrays
   if (myRank == masterTask) print *,'Testing 5d multi-field sums '
   flush(6)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   allocate(darray5n(nx,ny,nz,nl,nm,nFields), &
            rarray5n(nx,ny,nz,nl,nm,nFields), &
            iarray5n(nx,ny,nz,nl,nm,nFields), &
            larray5n(nx,ny,nz,nl,nm,nFields), &
             dmask5n(nx,ny,nz,nl,nm,nFields), &
             rmask5n(nx,ny,nz,nl,nm,nFields), &
             imask5n(nx,ny,nz,nl,nm,nFields), &
             lmask5n(nx,ny,nz,nl,nm,nFields))
   allocate(randoms(nx*ny*nz*nl*nm))
   call random_number(randoms)

   iref = 0
   do m=1,nm
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      iref = iref + 1
      rscalar = randoms(iref)
      ipower = int(rscalar*drange) - 5
      darray5n(i,j,k,l,m,1) = real(rscalar,R8)*(10.d0**ipower)
      ipower = int(rscalar*rrange) - 5
      rarray5n(i,j,k,l,m,1) = rscalar*(10.0**ipower)
      iarray5n(i,j,k,l,m,1) = nint(rscalar*irange)
      larray5n(i,j,k,l,m,1) = nint(rscalar*lrange)
   end do
   end do
   end do
   end do
   end do
   do m=1,nm
      m2 = nm+1-m
      m3 = mod(nm/2+m,nm) + 1
      do l=1,nl
         l2 = nl+1-l
         l3 = mod(nl/2+l,nl) + 1
         do k=1,nz
            k2 = nz+1-k
            k3 = mod(nz/2+k,nz) + 1
            do j=1,ny
               j2 = ny+1-j
               j3 = mod(ny/2+j,ny) + 1
               do i=1,nx
                  i2 = nx+1-i
                  i3 = mod(nx/2+i,nx) + 1
                  darray5n(i,j,k,l,m,2) = darray5n(i2,j2,k2,l2,m2,1)
                  darray5n(i,j,k,l,m,3) = darray5n(i3,j3,k3,l3,m3,1)
                  rarray5n(i,j,k,l,m,2) = rarray5n(i2,j2,k2,l2,m2,1)
                  rarray5n(i,j,k,l,m,3) = rarray5n(i3,j3,k3,l3,m3,1)
                  iarray5n(i,j,k,l,m,2) = iarray5n(i2,j2,k2,l2,m2,1)
                  iarray5n(i,j,k,l,m,3) = iarray5n(i3,j3,k3,l3,m3,1)
                  larray5n(i,j,k,l,m,2) = larray5n(i2,j2,k2,l2,m2,1)
                  larray5n(i,j,k,l,m,3) = larray5n(i3,j3,k3,l3,m3,1)
               end do
            end do
         end do
      end do
   end do

   do nfld=1,nFields
   do m=1,nm
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (i >= indxRange(1) .and. i <= indxRange(2) .and. &
          j >= indxRange(3) .and. j <= indxRange(4) .and. &
          k >= indxRange(5) .and. k <= indxRange(6) .and. &
          l >= indxRange(7) .and. l <= indxRange(8) .and. &
          m >= indxRange(9) .and. m <= indxRange(10)) then
         dmask5n(i,j,k,l,m,nfld) = 1.d0
         rmask5n(i,j,k,l,m,nfld) = 1.0
         lmask5n(i,j,k,l,m,nfld) = 1
         imask5n(i,j,k,l,m,nfld) = 1
      else
         dmask5n(i,j,k,l,m,nfld) = 0.d0
         rmask5n(i,j,k,l,m,nfld) = 0.0
         lmask5n(i,j,k,l,m,nfld) = 0
         imask5n(i,j,k,l,m,nfld) = 0
      endif
   end do
   end do
   end do
   end do
   end do
   end do
   ! Use the single-field interface to compute reference sums
   dref = mpas_global_sum(darray5n(:,:,:,:,:,1), MPI_COMM_WORLD)
   rref = mpas_global_sum(rarray5n(:,:,:,:,:,1), MPI_COMM_WORLD)
   iref = mpas_global_sum(iarray5n(:,:,:,:,:,1), MPI_COMM_WORLD)
   lref = mpas_global_sum(larray5n(:,:,:,:,:,1), MPI_COMM_WORLD)

   ! Now compute the various multi-field sums
   dsum1n = mpas_global_sum_nfld(darray5n, MPI_COMM_WORLD)
   rsum1n = mpas_global_sum_nfld(rarray5n, MPI_COMM_WORLD)
   isum1n = mpas_global_sum_nfld(iarray5n, MPI_COMM_WORLD)
   lsum1n = mpas_global_sum_nfld(larray5n, MPI_COMM_WORLD)
   dsum2n = mpas_global_sum_nfld(darray5n, MPI_COMM_WORLD, indxRange)
   rsum2n = mpas_global_sum_nfld(rarray5n, MPI_COMM_WORLD, indxRange)
   isum2n = mpas_global_sum_nfld(iarray5n, MPI_COMM_WORLD, indxRange)
   lsum2n = mpas_global_sum_nfld(larray5n, MPI_COMM_WORLD, indxRange)
   dsum3n = mpas_global_sum_nfld(darray5n, dmask5n, MPI_COMM_WORLD)
   rsum3n = mpas_global_sum_nfld(rarray5n, rmask5n, MPI_COMM_WORLD)
   isum3n = mpas_global_sum_nfld(iarray5n, imask5n, MPI_COMM_WORLD)
   lsum3n = mpas_global_sum_nfld(larray5n, lmask5n, MPI_COMM_WORLD)
   dsum4n = mpas_global_sum_nfld(darray5n, dmask5n, MPI_COMM_WORLD, indxRange)
   rsum4n = mpas_global_sum_nfld(rarray5n, rmask5n, MPI_COMM_WORLD, indxRange)
   isum4n = mpas_global_sum_nfld(iarray5n, imask5n, MPI_COMM_WORLD, indxRange)
   lsum4n = mpas_global_sum_nfld(larray5n, lmask5n, MPI_COMM_WORLD, indxRange)

   if (dsum1n(1) == dref .and. dsum1n(2) == dref .and. &
       dsum1n(3) == dref) then
      print *, 'Multi-field and reproducibility 5d double: PASS'
   else
      print *, 'Multi-field and reproducibility 5d double: FAIL', dref, dsum1n
   endif
   if (rsum1n(1) == rref .and. rsum1n(2) == rref .and. &
       rsum1n(3) == rref) then
      print *, 'Multi-field and reproducibility 5d real: PASS'
   else
      print *, 'Multi-field and reproducibility 5d real: FAIL', rref, rsum1n
   endif
   if (isum1n(1) == iref .and. isum1n(2) == iref .and. &
       isum1n(3) == iref) then
      print *, 'Multi-field and reproducibility 5d int: PASS'
   else
      print *, 'Multi-field and reproducibility 5d int: FAIL', iref, isum1n
   endif
   if (lsum1n(1) == lref .and. lsum1n(2) == lref .and. &
       lsum1n(3) == lref) then
      print *, 'Multi-field and reproducibility 5d int(8): PASS'
   else
      print *, 'Multi-field and reproducibility 5d int(8): FAIL', lref, lsum1n
   endif
   isum1 = 0
   isum2 = 0
   isum3 = 0
   isum4 = 0
   do nfld=1,nFields
      if (dsum2n(nfld) /= dsum3n(nfld)) isum1 = isum1+1
      if (dsum2n(nfld) /= dsum4n(nfld)) isum1 = isum1+1
      if (rsum2n(nfld) /= rsum3n(nfld)) isum2 = isum2+1
      if (rsum2n(nfld) /= rsum4n(nfld)) isum2 = isum2+1
      if (isum2n(nfld) /= isum3n(nfld)) isum3 = isum3+1
      if (isum2n(nfld) /= isum4n(nfld)) isum3 = isum3+1
      if (lsum2n(nfld) /= lsum3n(nfld)) isum4 = isum4+1
      if (lsum2n(nfld) /= lsum4n(nfld)) isum4 = isum4+1
   end do
   if (isum1 > 0) then
      print *, 'Multifield range/mask 5d double: FAIL', dsum2n,dsum3n,dsum4n
   else
      print *, 'Multifield range/mask 5d double: PASS'
   endif
   if (isum2 > 0) then
      print *, 'Multifield range/mask 5d real: FAIL', rsum2n,rsum3n,rsum4n
   else
      print *, 'Multifield range/mask 5d real: PASS'
   endif
   if (isum3 > 0) then
      print *, 'Multifield range/mask 5d int: FAIL', isum2n,isum3n,isum4n
   else
      print *, 'Multifield range/mask 5d int: PASS'
   endif
   if (isum4 > 0) then
      print *, 'Multifield range/mask 5d int(8): FAIL', lsum2n,lsum3n,lsum4n
   else
      print *, 'Multifield range/mask 5d int(8): PASS'
   endif
   deallocate(darray5n, dmask5n, rarray5n, rmask5n, &
              iarray5n, imask5n, larray5n, lmask5n)
   deallocate(randoms)

   !--------------------------------------------------------------------
   ! clean up some arrays
   deallocate(dsum1n, dsum2n, dsum3n, dsum4n, &
              rsum1n, rsum2n, rsum3n, rsum4n, &
              isum1n, isum2n, isum3n, isum4n, &
              lsum1n, lsum2n, lsum3n, lsum4n)

   !--------------------------------------------------------------------
   ! All done, exit

   call MPI_Finalize(ierr)

!***********************************************************************

end program sumTest

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
