c-----------------------------------------------------------------------
      subroutine addfid(fname,fid)

      character*1 fname(132)
      integer fid

      character*8  eight,fmt,s8
      save         eight
      data         eight / "????????" /

      do ipass=1,2      ! 2nd pass, in case 1 file/directory
         do k=8,1,-1
            i1 = indx1(fname,eight,k)
            if (i1.ne.0) then ! found k??? string
               write(fmt,1) k,k
               write(s8,fmt) fid
               call chcopy(fname(i1),s8,k)
               goto 10
            endif
         enddo
   10    continue
      enddo

    1 format('(i',i1,'.',i1,')')
      
      return
      end
c-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)

      character*1 a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine blank(a,n)

      character*1 a(1)
      character*1 blnk
      save        blnk
      data        blnk /' '/

      do i=1,n
         a(i)=blnk
		enddo

      return
      end
c-----------------------------------------------------------------------
		logical function if_byte_swap_test(a,ierr)

		real a(1)

	   return
	   end
c-----------------------------------------------------------------------
      subroutine mfi_parse_hdr(hdr,ierr)

      character*132 hdr

      if (indx2(hdr,132,'#std',4).eq.1) then
          call parse_std_hdr(hdr)
      else
         ierr = 1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr(hdr)

      include 'LVAR'
      include 'IO'

      character*132 hdr
      character*4 dummy
      logical if_press_mesh

      p0thr = -1
      if_press_mesh = .false.

      read(hdr,*,iostat=ierr) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94
     $         ,  p0thr, if_press_mesh

      if (ierr.gt.0) then ! try again without pressure format flag
        read(hdr,*,iostat=ierr) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94
     $         ,  p0thr
      endif

      if (ierr.gt.0) then ! try again without mean pressure
        read(hdr,*,err=99) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94
      endif

c     set if_full_pres flag
      if_full_pres = .false.
      if (.not.ifsplit) if_full_pres = if_press_mesh

c      ifgtim  = .true.  ! always get time
      ifgetxr = .false.
      ifgetur = .false.
      ifgetpr = .false.
      ifgettr = .false.
      do k=1,ldimt-1
         ifgtpsr(k) = .false.
      enddo

      NPSR = 0
      do i=1,10 
         if (rdcode1(i).eq.'X') ifgetxr = .true.
         if (rdcode1(i).eq.'U') ifgetur = .true.
         if (rdcode1(i).eq.'P') ifgetpr = .true.
         if (rdcode1(i).eq.'T') ifgettr = .true.
         if (rdcode1(i).eq.'S') then
            read(rdcode1(i+1),'(I1)') NPS1
            read(rdcode1(i+2),'(I1)') NPS0
            NPSR = 10*NPS1+NPS0
            NPS  = NPSR
            if(NPSR.gt.ldimt-1) NPS=ldimt-1
            do k=1,NPS
               ifgtpsr(k) = .true.
            enddo
            ! nothing will follow
            GOTO 50
         endif
      enddo

  50  if (NPS.lt.NPSR) then
         if (nid.eq.0) then 
           write(*,'(A,/,A)') 
     &      'WARNING: restart file has a NSPCAL > LDIMT',
     &      'read only part of the fld-data!'
         endif
      endif

      if (NPS.lt.NPSCAL) then
         if (nid.eq.0) then 
           write(*,'(A,/,A)') 
     &      'WARNING: NPSCAL read from restart file differs from ',
     &      'currently used NPSCAL!'
         endif
      endif

      p0th = 1 
      if (p0thr.gt.0) p0th = p0thr

      return

   99 continue   !  If we got here, then the May 2008 variant of std hdr
                 !  failed and we may have an older input file.

      call parse_std_hdr_2006(hdr,rdcode)  ! try the original header format

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr_2006(hdr,rlcode)

      include 'LVAR'
		include 'IO'

      character*132 hdr
      character*1 rlcode(20)

c                4  7  10  13   23    33    53    62     68     74
      read(hdr,1) wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         , ifiler,nfiler
     $         , (rlcode(k),k=1,20)                   ! 74+20=94
    1 format(4x,i2,3i3,2i10,e20.13,i9,2i6,20a1)

      if (nid.eq.0) write(6,*) 'WARNING: reading depreacted header!'

      if (nelr.gt.lelr) then
        write(6,*)nid,nelr,lelr,'parse_std_hdr06: inc. lelr in RESTART'
        call exitt
      endif

c     Assign read conditions, according to rdcode
c     NOTE: In the old hdr format: what you see in file is what you get.
c      ifgtim  = .true.  ! always get time
      ifgetxr = .false.
      ifgetur = .false.
      ifgetpr = .false.
      ifgettr = .false.
      do k=1,npscal
         ifgtpsr(k) = .false.
      enddo

      if (rlcode(1).eq.'X') ifgetxr = .true.
      if (rlcode(2).eq.'U') ifgetur = .true.
      if (rlcode(3).eq.'P') ifgetpr = .true.
      if (rlcode(4).eq.'T') ifgettr = .true.
      do k=1,npscal
         if (rlcode(4+k).ne.' ') ifgtpsr(k) = .true.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cfill(a,b,n)
      DIMENSION  A(1)
C
      DO 100 I = 1, N
 100     A(I) = B
      return
      END
c-----------------------------------------------------------------------
      subroutine chsign(a,n)
      REAL A(1)
C
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
      return
      END
C
c-----------------------------------------------------------------------
      function glsum (x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      subroutine i8copy(a,b,n)
      INTEGER*8 A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      function iglsum(a,n)
      integer a(1),tsum
      integer tmp(1),work(1)
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call igop(tmp,work,'+  ',1)
      iglsum=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*132 S1,S2
C
      N1=132-L2+1
      INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine mxm(a,n1,b,n2,c,n3)
c
c     Compute matrix-matrix product C = A*B
c     for contiguously packed matrices A,B, and C.
c
      real a(n1,n2),b(n2,n3),c(n1,n3)
c
      call mxmf2(a,n1,b,n2,c,n3)

      return
      end
c-----------------------------------------------------------------------
      subroutine opzero(ux,uy,uz)

      real ux(1),uy(1),uz(1)
c

      call exitti('called unsupported opzero$',1)
c     n = lx1*ly1*lz1*nelfld(ifield)
      call rzero(ux,n)
      call rzero(uy,n)
c     if (if3d) call rzero(uz,n)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine izero(a,n)
      integer a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rzero(a,n)
      real a(1)
      do i=1,n
         a(i)=0.0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rescale_x(x,x0,x1)

      real x(1)

      write (6,*) 'called rescale_x, unsupported'
      return
c     n = lx1*ly1*lz1*nelt
      xmin = glmin(x,n)
      xmax = glmax(x,n)

      if (xmax.le.xmin) return

      scale = (x1-x0)/(xmax-xmin)
      do i=1,n
         x(i) = x0 + scale*(x(i)-xmin)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function indx2(s1,l1,s2,l2)
      character*132 s1,s2

      n1=l1-l2+1
      indx2=0
      if (n1.lt.1) return

      do i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx2=i
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function log2(k)
      RK=(K)
      RLOG=LOG10(RK)
      RLOG2=LOG10(2.0)
      RLOG=RLOG/RLOG2+0.5
      LOG2=INT(RLOG)
      return
      END
c-----------------------------------------------------------------------
      subroutine drivep
      return
      end
c-----------------------------------------------------------------------
      subroutine flush_io
      return
      end
c-----------------------------------------------------------------------
      subroutine userqtl_scig
      return
      end
c-----------------------------------------------------------------------
      subroutine ione(a,n)
      INTEGER  A(1)
      DO 100 I = 1, N
 100     A(I ) = 1
      return
      END
c-----------------------------------------------------------------------
      subroutine isort(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(1),ind(1)
      integer aa
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine sort(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real a(1),aa
      integer ind(1)
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      function ltruncr(string,l)

      character*1 string(l)
      character*1   blnk
      data blnk/' '/

      do i=1,l
         l1=i-1
         if (string(i).eq.blnk) goto 200
      enddo
      l1=0

  200 continue
      ltruncr=l1

      return
      end
c-----------------------------------------------------------------------
      subroutine offline_init(comm_out)
C
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
      write (6,*) 'wp 1'
      nekcomm  = intracomm
      write (6,*) 'wp 2'
      comm_out = nekcomm
      write (6,*) 'wp 3'
      call iniproc
c     igeom = 2
c     call genwz           ! Compute GLL points, weights, etc.

      return
      end
c-----------------------------------------------------------------------