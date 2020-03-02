c-----------------------------------------------------------------------
      subroutine userqtl_scig
      return
      end
c-----------------------------------------------------------------------
      subroutine drivep
      return
      end
c-----------------------------------------------------------------------
      subroutine offline_init(comm_out)

      include 'LVAR'
      include 'IO'

      integer comm_out
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /rdump/ ntdump

      real kwave2
      logical ifemati

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest

      ! set word size for REAL
      wdsize = sizeof(rtest)
      ! set word size for INTEGER
      isize = sizeof(itest)
      ! set word size for INTEGER*8
      isize8 = sizeof(itest8)
      ! set word size for LOGICAL
      lsize = sizeof(ltest)
      ! set word size for CHARACTER
      csize = sizeof(ctest)

      call setupcomm
      nekcomm  = intracomm
      comm_out = nekcomm
      call iniproc
      igeom = 2
      call genwz           ! Compute GLL points, weights, etc.

      return
      end
c-----------------------------------------------------------------------