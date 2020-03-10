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
      subroutine reset_times

      include 'TIMES'

      time_ip=0.
      time_read=0.
      time_dump=0.
      time_sort=0.
      time_bip=0.
      time_aip=0.
      time_cip=0.
      time_cipc=0.
      time_bop=0.
      time_aop=0.
      time_cop=0.
      time_cdea=0.
      time_cconv=0.
      time_qop2=0.
      time_qop3=0.
      time_disk=0.
      time_io=0.
      time_off=0.

      return
      end
c-----------------------------------------------------------------------
      subroutine final

      include 'TIMES'

      common /nekmpi/ mid

      ! local dot-products
      if (mid.eq.0) write (6,10) 'time_bip:',time_bip,time_bip/time_off
      if (mid.eq.0) write (6,10) 'time_aip:',time_aip,time_aip/time_off
      if (mid.eq.0) write (6,10) 'time_cip:',time_cip,time_cip/time_off

      ! local operator application
      if (mid.eq.0) write (6,10) 'time_bop:',time_bop,time_bop/time_off
      if (mid.eq.0) write (6,10) 'time_aop:',time_aop,time_aop/time_off
      if (mid.eq.0) write (6,10) 'time_cop:',time_cdea+time_cconv,
     $   (time_cdea+time_cconv)/time_off
      if (mid.eq.0) write (6,10) 'time_dea:',time_cdea,
     $   time_cdea/time_off
      if (mid.eq.0) write (6,10) 'time_cnv:',time_cconv,
     $   time_cconv/time_off

      ! inter-rank communication
      if (mid.eq.0) write (6,10) 'time_gop:',time_gop,time_gop/time_off
      if (mid.eq.0) write (6,10) 'time_cst:',time_cst,time_cst/time_off

      ! disk related tasks
      if (mid.eq.0) write (6,10) 'time_dsk:',time_disk,
     $   time_disk/time_off
      if (mid.eq.0) write (6,10) 'time_dmp:',time_dump,
     $   time_dump/time_off

      ! hot
      if (mid.eq.0) write (6,10) 'time_mxm:',time_mxm,time_mxm/time_off
      if (mid.eq.0) write (6,10) 'time_eig:',time_eig,time_eig/time_off
      if (mid.eq.0) write (6,10) 'time_srt:',time_sort,
     $   time_sort/time_off

      ! total offline-mode
      if (mid.eq.0) write (6,10) 'time_off:',time_off,time_off/time_off

   10 format (a9,1p2e11.3)

      return
      end
c-----------------------------------------------------------------------