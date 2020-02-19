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
		logical function if_byte_swap_test(a,ierr) ! TODO: implement

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
      function iglmin(a,n)
      integer a(1),tmin
      integer tmp(1),work(1)
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call igop(tmp,work,'m  ',1)
      iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      function iglmax(a,n)
      integer a(1),tmax
      integer tmp(1),work(1)
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call igop(tmp,work,'M  ',1)
      iglmax=tmp(1)
      return
      end
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
      nekcomm  = intracomm
      comm_out = nekcomm
      call iniproc
      igeom = 2
      call genwz           ! Compute GLL points, weights, etc.

      return
      end
c-----------------------------------------------------------------------
      subroutine genwz
C-----------------------------------------------------------------
C
C     GENERATE
C
C            - DERIVATIVE OPERATORS
C            - INTERPOLATION OPERATORS
C            - WEIGHTS
C            - COLLOCATION POINTS
C
C     ASSOCIATED WITH THE
C
C            - GAUSS-LOBATTO LEGENDRE MESH (SUFFIX M1/M2/M3)
C            - GAUSS LEGENDRE         MESH (SUFFIX M2)
C            - GAUSS-LOBATTO JACOBI   MESH (SUFFIX M1/M2/M3)
C
C-----------------------------------------------------------------
C
      INCLUDE 'LVAR'
      INCLUDE 'INTEG'

      logical ifsplit,ifaxis
      REAL TMP(LY1,LY1),TMPT(LY1,LY1)

      ifsplit=lx1.eq.lx2
      ifaxis=.false.
C
      IF (ldim.EQ.2) THEN
C
C***  Two-dimensional case  **********************
C
C
C     Gauss-Lobatto Legendre mesh (suffix M1)
C     Generate collocation points and weights
C
      CALL ZWGLL (ZGM1(1,1),WXM1,lx1)
      CALL ZWGLL (ZGM1(1,2),WYM1,ly1)
      ZGM1(lz1,3) = 0.
      WZM1(lz1)   = 1.
      DO 100 IY=1,ly1
      DO 100 IX=1,lx1
      W3M1(IX,IY,1)=WXM1(IX)*WYM1(IY)
  100 CONTINUE
C
C     Compute derivative matrices
C
      CALL DGLL (DXM1,DXTM1,ZGM1(1,1),lx1,lx1)
      CALL DGLL (DYM1,DYTM1,ZGM1(1,2),ly1,ly1)
      CALL RZERO (DZM1 ,lz1*lz1)
      CALL RZERO (DZTM1,lz1*lz1)
C
C     Gauss Legendre mesh (suffix M2)
C     Generate collocation points and weights
C
      IF(IFSPLIT)THEN
         CALL ZWGLL (ZGM2(1,1),WXM2,lx2)
         CALL ZWGLL (ZGM2(1,2),WYM2,ly2)
      ELSE
         CALL ZWGL  (ZGM2(1,1),WXM2,lx2)
         CALL ZWGL  (ZGM2(1,2),WYM2,ly2)
      ENDIF
      ZGM2(lz2,3) = 0.
      WZM2(lz2)   = 1.
      DO 200 IY=1,ly2
      DO 200 IX=1,lx2
      W3M2(IX,IY,1)=WXM2(IX)*WYM2(IY)
  200 CONTINUE
C
C     Gauss-Lobatto Legendre mesh (suffix M3).
C     Generate collocation points and weights.
C
      CALL ZWGLL (ZGM3(1,1),WXM3,lx3)
      CALL ZWGLL (ZGM3(1,2),WYM3,ly3)
      ZGM3(lz3,3) = 0.
      WZM3(lz3)   = 1.
      DO 300 IY=1,ly3
      DO 300 IX=1,lx3
      W3M3(IX,IY,1)=WXM3(IX)*WYM3(IY)
  300 CONTINUE
C
C     Compute derivative matrices
C
      CALL DGLL (DXM3,DXTM3,ZGM3(1,1),lx3,lx3)
      CALL DGLL (DYM3,DYTM3,ZGM3(1,2),ly3,ly3)
      CALL RZERO (DZM3 ,lz3*lz3)
      CALL RZERO (DZTM3,lz3*lz3)
C
C     Generate interpolation operators for the staggered mesh
C
      CALL IGLLM (IXM12,IXTM12,ZGM1(1,1),ZGM2(1,1),lx1,lx2,lx1,lx2)
      CALL IGLLM (IYM12,IYTM12,ZGM1(1,2),ZGM2(1,2),ly1,ly2,ly1,ly2)
      IZM12 (lz2,lz1) = 1.
      IZTM12(lz1,lz2) = 1.
C
C     NOTE: The splitting scheme has only one mesh!!!!!
C
      IF (IFSPLIT) THEN
         CALL IGLLM (IXM21,IXTM21,ZGM1(1,1),ZGM2(1,1),lx1,lx2,lx1,lx2)
         CALL IGLLM (IYM21,IYTM21,ZGM1(1,2),ZGM2(1,2),ly1,ly2,ly1,ly2)
      ELSE
         CALL IGLM  (IXM21,IXTM21,ZGM2(1,1),ZGM1(1,1),lx2,lx1,lx2,lx1)
         CALL IGLM  (IYM21,IYTM21,ZGM2(1,2),ZGM1(1,2),ly2,ly1,ly2,ly1)
      ENDIF
      IZM21 (lz1,lz2) = 1.
      IZTM21(lz2,lz1) = 1.
C
C     Compute derivative operators for the staggered mesh
C
      IF(IFSPLIT)THEN
         CALL COPY (DXM12, DXM1, lx1*lx2)
         CALL COPY (DXTM12,DXTM1,lx1*lx2)
         CALL COPY (DYM12, DYM1, ly1*ly2)
         CALL COPY (DYTM12,DYTM1,ly1*ly2)
         CALL COPY (DZM12, DZM1, lz1*lz2)
         CALL COPY (DZTM12,DZTM1,lz1*lz2)
      ELSE
         CALL DGLLGL (DXM12,DXTM12,ZGM1(1,1),ZGM2(1,1),IXM12,
     $                                       lx1,lx2,lx1,lx2)
         CALL DGLLGL (DYM12,DYTM12,ZGM1(1,2),ZGM2(1,2),IYM12,
     $                                       ly1,ly2,ly1,ly2)
         DZM12 (lz2,lz1) = 0.
         DZTM12(lz2,lz1) = 0.
      ENDIF
C
C     Compute interpolation operators for the geometry mesh M3.
C
      CALL IGLLM (IXM13,IXTM13,ZGM1(1,1),ZGM3(1,1),lx1,lx3,lx1,lx3)
      CALL IGLLM (IYM13,IYTM13,ZGM1(1,2),ZGM3(1,2),ly1,ly3,ly1,ly3)
      CALL IGLLM (IXM31,IXTM31,ZGM3(1,1),ZGM1(1,1),lx3,lx1,lx3,lx1)
      CALL IGLLM (IYM31,IYTM31,ZGM3(1,2),ZGM1(1,2),ly3,ly1,ly3,ly1)
      IZM13 (lz3,lz1) = 1.
      IZTM13(lz1,lz3) = 1.
      IZM31 (lz1,lz3) = 1.
      IZTM31(lz3,lz1) = 1.
C
C
      IF (IFAXIS) THEN
C
C     Special treatment for the axisymmetric case
C     Generate additional points, weights, derivative operators and
C     interpolation operators required for elements close to the axis.
C
C
C     Gauss-Lobatto Jacobi mesh (suffix M1).
C     Generate collocation points and weights (alpha=0, beta=1).
C
      ALPHA = 0.
      BETA  = 1.
      CALL ZWGLJ (ZAM1,WAM1,ly1,ALPHA,BETA)
      DO 400 IY=1,ly1
      DO 400 IX=1,lx1
         W2AM1(IX,IY)=WXM1(IX)*WAM1(IY)
         W2CM1(IX,IY)=WXM1(IX)*WYM1(IY)
  400 CONTINUE
C
C     Compute derivative matrices
C
      CALL COPY (DCM1,DYM1,ly1*ly1)
      CALL COPY (DCTM1,DYTM1,ly1*ly1)
      CALL DGLJ (DAM1,DATM1,ZAM1,ly1,ly1,ALPHA,BETA)
C
C     Gauss Jacobi mesh (suffix M2)
C     Generate collocation points and weights
C
      IF(IFSPLIT)THEN
         CALL ZWGLJ (ZAM2,WAM2,ly2,ALPHA,BETA)
      ELSE
         CALL ZWGJ  (ZAM2,WAM2,ly2,ALPHA,BETA)
      ENDIF
      DO 500 IY=1,ly2
      DO 500 IX=1,lx2
         W2CM2(IX,IY)=WXM2(IX)*WYM2(IY)
         W2AM2(IX,IY)=WXM2(IX)*WAM2(IY)
  500 CONTINUE
C
C     Gauss-Lobatto Jacobi mesh (suffix M3).
C     Generate collocation points and weights.
C
      CALL ZWGLJ (ZAM3,WAM3,ly3,ALPHA,BETA)
      DO 600 IY=1,ly3
      DO 600 IX=1,lx3
         W2CM3(IX,IY)=WXM3(IX)*WYM3(IY)
         W2AM3(IX,IY)=WXM3(IX)*WAM3(IY)
  600 CONTINUE
C
C     Compute derivative matrices
C
      CALL COPY (DCM3,DYM3,ly3*ly3)
      CALL COPY (DCTM3,DYTM3,ly3*ly3)
      CALL DGLJ (DAM3,DATM3,ZAM3,ly3,ly3,ALPHA,BETA)
C
C     Generate interpolation operators for the staggered mesh
C
      CALL COPY  (ICM12,IYM12,ly2*ly1)
      CALL COPY  (ICTM12,IYTM12,ly1*ly2)
      CALL IGLJM (IAM12,IATM12,ZAM1,ZAM2,ly1,ly2,ly1,ly2,ALPHA,BETA)
      CALL COPY  (ICM21,IYM21,ly1*ly2)
      CALL COPY  (ICTM21,IYTM21,ly2*ly1)
      IF (IFSPLIT) THEN
      CALL IGLJM (IAM21,IATM21,ZAM2,ZAM1,ly1,ly2,ly1,ly2,ALPHA,BETA)
      ELSE
      CALL IGJM  (IAM21,IATM21,ZAM2,ZAM1,ly2,ly1,ly2,ly1,ALPHA,BETA)
      ENDIF
C
C     Compute derivative operators for the staggered mesh
C
      CALL COPY  (DCM12,DYM12,ly2*ly1)
      CALL COPY  (DCTM12,DYTM12,ly1*ly2)
      IF(IFSPLIT)THEN
         CALL COPY (DAM12, DAM1, ly1*ly2)
         CALL COPY (DATM12,DATM1,ly1*ly2)
      ELSE
         CALL DGLJGJ (DAM12,DATM12,ZAM1,ZAM2,IAM12,
     $                             ly1,ly2,ly1,ly2,ALPHA,BETA)
      ENDIF
C
C     Compute interpolation operators for the geometry mesh M3.
C
      CALL COPY  (ICM13,IYM13,ly3*ly1)
      CALL COPY  (ICTM13,IYTM13,ly1*ly3)
      CALL IGLJM (IAM13,IATM13,ZAM1,ZAM3,ly1,ly3,ly1,ly3,ALPHA,BETA)
      CALL COPY  (ICM31,IYM31,ly1*ly3)
      CALL COPY  (ICTM31,IYTM31,ly3*ly1)
      CALL IGLJM (IAM31,IATM31,ZAM3,ZAM1,ly3,ly1,ly3,ly1,ALPHA,BETA)
C
C     Compute interpolation operators between Gauss-Lobatto Jacobi
C     and Gauss-Lobatto Legendre (to be used in PREPOST).
C
      CALL IGLJM(IAJL1,IATJL1,ZAM1,ZGM1(1,2),ly1,ly1,ly1,ly1,ALPHA,BETA)
      IF (IFSPLIT) THEN
      CALL IGLJM(IAJL2,IATJL2,ZAM2,ZGM2(1,2),ly2,ly2,ly2,ly2,ALPHA,BETA)
      ELSE
      CALL IGJM (IAJL2,IATJL2,ZAM2,ZGM2(1,2),ly2,ly2,ly2,ly2,ALPHA,BETA)
      ENDIF

      CALL INVMT(IAJL1 ,IALJ1 ,TMP ,ly1)
      CALL INVMT(IATJL1,IATLJ1,TMPT,ly1)
      CALL MXM (IATJL1,ly1,IATLJ1,ly1,TMPT,ly1)
      CALL MXM (IAJL1 ,ly1,IALJ1 ,ly1,TMP ,ly1)

C
C     Compute interpolation operators between Gauss-Lobatto Legendre
C     and Gauss-Lobatto Jacobi (to be used in subr. genxyz IN postpre).
C
c
c     This call is not right, and these arrays are not used. 3/27/02. pff
c     CALL IGLLM(IALJ3,IATLJ3,ZGM3(1,2),ZAM3,ly3,ly3,ly3,ly3,ALPHA,BETA)
      CALL IGLJM(IALJ3,IATLJ3,ZGM3(1,2),ZAM3,ly3,ly3,ly3,ly3,ALPHA,BETA)
C
      ENDIF
C
C
      ELSE
C
C***  Three-dimensional case ************************************
C
C
C     Gauss-Lobatto Legendre mesh (suffix M1)
C     Generate collocation points and weights
C
      CALL ZWGLL (ZGM1(1,1),WXM1,lx1)
      CALL ZWGLL (ZGM1(1,2),WYM1,ly1)
      CALL ZWGLL (ZGM1(1,3),WZM1,lz1)
      DO 700 IZ=1,lz1
      DO 700 IY=1,ly1
      DO 700 IX=1,lx1
      W3M1(IX,IY,IZ)=WXM1(IX)*WYM1(IY)*WZM1(IZ)
  700 CONTINUE
C
C     Compute derivative matrices
C
      CALL DGLL (DXM1,DXTM1,ZGM1(1,1),lx1,lx1)
      CALL DGLL (DYM1,DYTM1,ZGM1(1,2),ly1,ly1)
      CALL DGLL (DZM1,DZTM1,ZGM1(1,3),lz1,lz1)
C
C     Gauss Legendre mesh (suffix M2)
C     Generate collocation points and weights
C
      IF(IFSPLIT)THEN
         CALL ZWGLL (ZGM2(1,1),WXM2,lx2)
         CALL ZWGLL (ZGM2(1,2),WYM2,ly2)
         CALL ZWGLL (ZGM2(1,3),WZM2,lz2)
      ELSE
         CALL ZWGL  (ZGM2(1,1),WXM2,lx2)
         CALL ZWGL  (ZGM2(1,2),WYM2,ly2)
         CALL ZWGL  (ZGM2(1,3),WZM2,lz2)
      ENDIF
      DO 800 IZ=1,lz2
      DO 800 IY=1,ly2
      DO 800 IX=1,lx2
      W3M2(IX,IY,IZ)=WXM2(IX)*WYM2(IY)*WZM2(IZ)
  800 CONTINUE
C
C     Gauss-Loabtto Legendre mesh (suffix M3).
C     Generate collocation points and weights.
C
      CALL ZWGLL (ZGM3(1,1),WXM3,lx3)
      CALL ZWGLL (ZGM3(1,2),WYM3,ly3)
      CALL ZWGLL (ZGM3(1,3),WZM3,lz3)
      DO 900 IZ=1,lz3
      DO 900 IY=1,ly3
      DO 900 IX=1,lx3
      W3M3(IX,IY,IZ)=WXM3(IX)*WYM3(IY)*WZM3(IZ)
  900 CONTINUE
C
C     Compute derivative matrices
C
      CALL DGLL (DXM3,DXTM3,ZGM3(1,1),lx3,lx3)
      CALL DGLL (DYM3,DYTM3,ZGM3(1,2),ly3,ly3)
      CALL DGLL (DZM3,DZTM3,ZGM3(1,3),lz3,lz3)
C
C     Generate interpolation operators for the staggered mesh
C
      CALL IGLLM (IXM12,IXTM12,ZGM1(1,1),ZGM2(1,1),lx1,lx2,lx1,lx2)
      CALL IGLLM (IYM12,IYTM12,ZGM1(1,2),ZGM2(1,2),ly1,ly2,ly1,ly2)
      CALL IGLLM (IZM12,IZTM12,ZGM1(1,3),ZGM2(1,3),lz1,lz2,lz1,lz2)
C
C     NOTE: The splitting scheme has only one mesh!!!!!
C
      IF (IFSPLIT) THEN
         CALL IGLLM (IXM21,IXTM21,ZGM1(1,1),ZGM2(1,1),lx1,lx2,lx1,lx2)
         CALL IGLLM (IYM21,IYTM21,ZGM1(1,2),ZGM2(1,2),ly1,ly2,ly1,ly2)
         CALL IGLLM (IZM21,IZTM21,ZGM1(1,3),ZGM2(1,3),lz1,lz2,lz1,lz2)
      ELSE
         CALL IGLM  (IXM21,IXTM21,ZGM2(1,1),ZGM1(1,1),lx2,lx1,lx2,lx1)
         CALL IGLM  (IYM21,IYTM21,ZGM2(1,2),ZGM1(1,2),ly2,ly1,ly2,ly1)
         CALL IGLM  (IZM21,IZTM21,ZGM2(1,3),ZGM1(1,3),lz2,lz1,lz2,lz1)
      ENDIF
C
C     Compute derivative operators for the staggered mesh
C
      IF(IFSPLIT)THEN
         CALL COPY (DXM12, DXM1, lx1*lx2)
         CALL COPY (DXTM12,DXTM1,lx1*lx2)
         CALL COPY (DYM12, DYM1, ly1*ly2)
         CALL COPY (DYTM12,DYTM1,ly1*ly2)
         CALL COPY (DZM12, DZM1, lz1*lz2)
         CALL COPY (DZTM12,DZTM1,lz1*lz2)
      ELSE
         CALL DGLLGL (DXM12,DXTM12,ZGM1(1,1),ZGM2(1,1),IXM12,
     $                                       lx1,lx2,lx1,lx2)
         CALL DGLLGL (DYM12,DYTM12,ZGM1(1,2),ZGM2(1,2),IYM12,
     $                                       ly1,ly2,ly1,ly2)
         CALL DGLLGL (DZM12,DZTM12,ZGM1(1,3),ZGM2(1,3),IZM12,
     $                                       lz1,lz2,lz1,lz2)
      ENDIF
c
c     Compute interpolation operators for the geometry mesh M3.
c
      call igllm (ixm13,ixtm13,zgm1(1,1),zgm3(1,1),lx1,lx3,lx1,lx3)
      call igllm (iym13,iytm13,zgm1(1,2),zgm3(1,2),ly1,ly3,ly1,ly3)
      call igllm (izm13,iztm13,zgm1(1,3),zgm3(1,3),lz1,lz3,lz1,lz3)
      call igllm (ixm31,ixtm31,zgm3(1,1),zgm1(1,1),lx3,lx1,lx3,lx1)
      call igllm (iym31,iytm31,zgm3(1,2),zgm1(1,2),ly3,ly1,ly3,ly1)
      call igllm (izm31,iztm31,zgm3(1,3),zgm1(1,3),lz3,lz1,lz3,lz1)
c
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine geodat1(nel)
C-----------------------------------------------------------------------
C
C     Routine to generate elemental geometric matrices on mesh 1
C     (Gauss-Legendre Lobatto mesh).
C
C-----------------------------------------------------------------------
      INCLUDE 'LVAR'
c     INCLUDE 'TSTEP'
      INCLUDE 'INTEG'
C
C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3 
C           share the same array structure in Scratch Common /SCRNS/.
C
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
     $ ,             XTM1(LX1,LY1,LZ1,LELT)
     $ ,             YTM1(LX1,LY1,LZ1,LELT)
     $ ,             ZRM1(LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ ZSM1(LX1,LY1,LZ1,LELT)
     $ ,             ZTM1(LX1,LY1,LZ1,LELT)
     $ ,             WJ   (LX1,LY1,LZ1,LELT)

      logical ifaxis

      ifaxis=.false.

      NXYZ1 = lxyz
      NTOT1 = lxyz*NEL

      IF (.NOT.IFAXIS) THEN
         CALL INVERS2 (WJ,JACM1,NTOT1)
      ELSE
         DO 500 IEL=1,NEL
           IF (IFRZER(IEL)) THEN
              DO 510 J=1,ly1
              DO 510 I=1,lx1
                IF (J.GT.1) THEN
                   WJ(I,J,1,IEL) = YM1(I,J,1,IEL)/
     $                            (JACM1(I,J,1,IEL)*(1.+ZAM1(J)))
                ELSE
                   WJ(I,J,1,IEL) = YSM1(I,J,1,IEL)/JACM1(I,J,1,IEL)
                ENDIF
 510          CONTINUE
           ELSE
              CALL INVCOL3 (WJ(1,1,1,IEL),YM1(1,1,1,IEL),
     $                      JACM1(1,1,1,IEL),NXYZ1)
           ENDIF
 500     CONTINUE
      ENDIF
C
C     Compute geometric factors for integrated del-squared operator.
C
      IF (ldim.EQ.2) THEN
         CALL VDOT2 (G1M1,RXM1,RYM1,RXM1,RYM1,NTOT1)
         CALL VDOT2 (G2M1,SXM1,SYM1,SXM1,SYM1,NTOT1)
         CALL VDOT2 (G4M1,RXM1,RYM1,SXM1,SYM1,NTOT1)
         CALL COL2  (G1M1,WJ,NTOT1)
         CALL COL2  (G2M1,WJ,NTOT1)
         CALL COL2  (G4M1,WJ,NTOT1)
         CALL RZERO (G3M1,NTOT1)
         CALL RZERO (G5M1,NTOT1)
         CALL RZERO (G6M1,NTOT1)
      ELSE
         CALL VDOT3 (G1M1,RXM1,RYM1,RZM1,RXM1,RYM1,RZM1,NTOT1)
         CALL VDOT3 (G2M1,SXM1,SYM1,SZM1,SXM1,SYM1,SZM1,NTOT1)
         CALL VDOT3 (G3M1,TXM1,TYM1,TZM1,TXM1,TYM1,TZM1,NTOT1)
         CALL VDOT3 (G4M1,RXM1,RYM1,RZM1,SXM1,SYM1,SZM1,NTOT1)
         CALL VDOT3 (G5M1,RXM1,RYM1,RZM1,TXM1,TYM1,TZM1,NTOT1)
         CALL VDOT3 (G6M1,SXM1,SYM1,SZM1,TXM1,TYM1,TZM1,NTOT1)
         CALL COL2  (G1M1,WJ,NTOT1)
         CALL COL2  (G2M1,WJ,NTOT1)
         CALL COL2  (G3M1,WJ,NTOT1)
         CALL COL2  (G4M1,WJ,NTOT1)
         CALL COL2  (G5M1,WJ,NTOT1)
         CALL COL2  (G6M1,WJ,NTOT1)
      ENDIF
C
C     Multiply the geometric factors GiM1,i=1,5 with the
C     weights on mesh M1.
C
      DO 580 IEL=1,NEL
         IF (IFAXIS) CALL SETAXW1 ( IFRZER(IEL) )
            CALL COL2 (G1M1(1,1,1,IEL),W3M1,NXYZ1)
            CALL COL2 (G2M1(1,1,1,IEL),W3M1,NXYZ1)
            CALL COL2 (G4M1(1,1,1,IEL),W3M1,NXYZ1)
         IF (ldim.EQ.3) THEN
            CALL COL2 (G3M1(1,1,1,IEL),W3M1,NXYZ1)
            CALL COL2 (G5M1(1,1,1,IEL),W3M1,NXYZ1)
            CALL COL2 (G6M1(1,1,1,IEL),W3M1,NXYZ1)
         ENDIF
  580 CONTINUE
C
C     Compute the mass matrix on mesh M1.
C
      DO 700 IEL=1,NEL
         IF (IFAXIS) CALL SETAXW1 ( IFRZER(IEL) )
            CALL COL3 (BM1  (1,1,1,IEL),JACM1(1,1,1,IEL),W3M1,NXYZ1)
         IF (IFAXIS) THEN 
             CALL COL3(BAXM1(1,1,1,IEL),JACM1(1,1,1,IEL),W3M1,NXYZ1)
          IF (IFRZER(IEL)) THEN
            DO 600 J=1,ly1
            IF (J.GT.1) THEN
               DO 610 I=1,lx1
                  BM1(I,J,1,IEL) = BM1(I,J,1,IEL)*YM1(I,J,1,IEL)
     $                                           /(1.+ZAM1(J))
                  BAXM1(I,J,1,IEL)=BAXM1(I,J,1,IEL)/(1.+ZAM1(J))
 610           CONTINUE
            ELSE
               DO 620 I=1,lx1
                  BM1(I,J,1,IEL) = BM1(I,J,1,IEL)*YSM1(I,J,1,IEL)
                  BAXM1(I,J,1,IEL)=BAXM1(I,J,1,IEL)
 620           CONTINUE
            ENDIF
 600        CONTINUE
          ELSE
            CALL COL2 (BM1(1,1,1,IEL),YM1(1,1,1,IEL),NXYZ1)
          ENDIF
         ENDIF
C
 700  CONTINUE

      IF(IFAXIS) THEN
        DO IEL=1,NEL
          IF(IFRZER(IEL)) THEN
            DO J=1,ly1
            DO I=1,lx1
              IF(J.EQ.1) THEN
                 YINVM1(I,J,1,IEL)=1.0D0/YSM1(I,J,1,IEL)
              ELSE
                 YINVM1(I,J,1,IEL)=1.0D0/YM1 (I,J,1,IEL)
              ENDIF
            ENDDO 
            ENDDO 
          ELSE
            CALL INVERS2(YINVM1(1,1,1,IEL),YM1(1,1,1,IEL),NXYZ1)
          ENDIF
        ENDDO
      ELSE
        CALL CFILL(YINVM1,1.0D0,NXYZ1*NEL)
      ENDIF
C
C     Compute normals, tangents, and areas on elemental surfaces
C
      CALL SETAREA(nel)
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine col2(a,b,n)

      real a(1),b(1)

      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine invers2(a,b,n)

      real a(1),b(1)
 
      do i=1,n
         a(i)=1./b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine vdot2 (dot,u1,u2,v1,v2,n)
C
C     Compute a Cartesian vector dot product. 2-d version
C
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1)
      DIMENSION V1(1),V2(1)
C
C
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) 
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine setarea(nel)
c
c     Compute surface data: areas, normals and tangents
c
      include 'LVAR'
      include 'INTEG'

      nsrf=6*lx1*lz1*nel

      call rzero  (area,nsrf)
      call rzero3 (unx,uny,unz,nsrf)
      call rzero3 (t1x,t1y,t1z,nsrf)      
      call rzero3 (t2x,t2y,t2z,nsrf)      

      if (ldim.eq.2) then
         call area2(nel)
      else
         call area3(nel)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine area2(nel)
C--------------------------------------------------------------------
C
C     Compute areas, normals and tangents (2D and Axisymmetric geom.)
C
C--------------------------------------------------------------------
      INCLUDE 'LVAR'
      INCLUDE 'INTEG'
C
C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3 
C           share the same array structure in Scratch Common /SCRNS/.
C
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
      COMMON /CTMP0/ WGTR1(LX1,LELT)
     $ ,             WGTR2(LY1,LELT)
     $ ,             WGTR3(LX1,LELT)
     $ ,             WGTR4(LY1,LELT)
C
      CALL SETWGTR(WGTR1,WGTR2,WGTR3,WGTR4,nel)
C
C     "R"
C
      DO 100 IEL=1,NEL
      DO 100 IY=1,ly1
         XS2  = XSM1(lx1,IY,1,IEL)
         YS2  = YSM1(lx1,IY,1,IEL)
         XS4  = XSM1(  1,IY,1,IEL)
         YS4  = YSM1(  1,IY,1,IEL)
         SS2  = SQRT( XS2**2 + YS2**2 )
         SS4  = SQRT( XS4**2 + YS4**2 )
         T1X (IY,1,2,IEL) =  XS2 / SS2
         T1Y (IY,1,2,IEL) =  YS2 / SS2
         T1X (IY,1,4,IEL) = -XS4 / SS4
         T1Y (IY,1,4,IEL) = -YS4 / SS4
         UNX (IY,1,2,IEL) =  T1Y(IY,1,2,IEL)
         UNY (IY,1,2,IEL) = -T1X(IY,1,2,IEL)
         UNX (IY,1,4,IEL) =  T1Y(IY,1,4,IEL)
         UNY (IY,1,4,IEL) = -T1X(IY,1,4,IEL)
         AREA(IY,1,2,IEL) =  SS2 * WGTR2(IY,IEL)
         AREA(IY,1,4,IEL) =  SS4 * WGTR4(IY,IEL)
  100 CONTINUE
C
C     "S"
C
      DO 200 IEL=1,NEL
      DO 200 IX=1,lx1
         XR1  = XRM1(IX,  1,1,IEL)
         YR1  = YRM1(IX,  1,1,IEL)
         XR3  = XRM1(IX,ly1,1,IEL)
         YR3  = YRM1(IX,ly1,1,IEL)
         RR1  = SQRT( XR1**2 + YR1**2 )
         RR3  = SQRT( XR3**2 + YR3**2 )
         T1X (IX,1,1,IEL) =  XR1 / RR1
         T1Y (IX,1,1,IEL) =  YR1 / RR1
         T1X (IX,1,3,IEL) = -XR3 / RR3
         T1Y (IX,1,3,IEL) = -YR3 / RR3
         UNX (IX,1,1,IEL) =  T1Y(IX,1,1,IEL)
         UNY (IX,1,1,IEL) = -T1X(IX,1,1,IEL)
         UNX (IX,1,3,IEL) =  T1Y(IX,1,3,IEL)
         UNY (IX,1,3,IEL) = -T1X(IX,1,3,IEL)
         AREA(IX,1,1,IEL) =  RR1 * WGTR1(IX,IEL)
         AREA(IX,1,3,IEL) =  RR3 * WGTR3(IX,IEL)
  200 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine area3(nel)
C--------------------------------------------------------------------
C
C     Compute areas, normals and tangents (3D geom.)
C
C--------------------------------------------------------------------
      include 'LVAR'
      include 'INTEG'
C
C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3 
C           share the same array structure in Scratch Common /SCRNS/.
C
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
     $ ,             XTM1(LX1,LY1,LZ1,LELT)
     $ ,             YTM1(LX1,LY1,LZ1,LELT)
     $ ,             ZRM1(LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ ZSM1(LX1,LY1,LZ1,LELT)
     $ ,             ZTM1(LX1,LY1,LZ1,LELT)
     $ ,             A  (LX1,LY1,LZ1,LELT)
     $ ,             B  (LX1,LY1,LZ1,LELT)
      COMMON /CTMP0/ C  (LX1,LY1,LZ1,LELT)
     $ ,             DOT(LX1,LY1,LZ1,LELT)
C
      NXY1  = lx1*ly1
      NFACE = 2*ldim
      NTOT  = lx1*ly1*lz1*NEL
      NSRF  = 6*lx1*ly1*NEL
C
C        "R"
C
      CALL VCROSS(A,B,C,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,NTOT)
      CALL VDOT3 (DOT,A,B,C,A,B,C,NTOT)
C
      DO 100 IEL=1,NEL
      DO 100 IZ=1,lz1
      DO 100 IY=1,ly1
         WEIGHT = WYM1(IY)*WZM1(IZ)
         AREA(IY,IZ,2,IEL) = SQRT(DOT(lx1,IY,IZ,IEL))*WEIGHT
         AREA(IY,IZ,4,IEL) = SQRT(DOT(  1,IY,IZ,IEL))*WEIGHT
         UNX (IY,IZ,4,IEL) = -A(  1,IY,IZ,IEL)
         UNX (IY,IZ,2,IEL) =  A(lx1,IY,IZ,IEL)
         UNY (IY,IZ,4,IEL) = -B(  1,IY,IZ,IEL)
         UNY (IY,IZ,2,IEL) =  B(lx1,IY,IZ,IEL)
         UNZ (IY,IZ,4,IEL) = -C(  1,IY,IZ,IEL)
         UNZ (IY,IZ,2,IEL) =  C(lx1,IY,IZ,IEL)
  100 CONTINUE
C
C        "S"
C
      CALL VCROSS(A,B,C,XRM1,YRM1,ZRM1,XTM1,YTM1,ZTM1,NTOT)
      CALL VDOT3 (DOT,A,B,C,A,B,C,NTOT)
      DO 200 IEL=1,NEL
      DO 200 IZ=1,lz1
      DO 200 IX=1,lx1
         WEIGHT=WXM1(IX)*WZM1(IZ)
         AREA(IX,IZ,1,IEL) = SQRT(DOT(IX,  1,IZ,IEL))*WEIGHT
         AREA(IX,IZ,3,IEL) = SQRT(DOT(IX,ly1,IZ,IEL))*WEIGHT
         UNX (IX,IZ,1,IEL) =  A(IX,  1,IZ,IEL)
         UNX (IX,IZ,3,IEL) = -A(IX,ly1,IZ,IEL)
         UNY (IX,IZ,1,IEL) =  B(IX,  1,IZ,IEL)
         UNY (IX,IZ,3,IEL) = -B(IX,ly1,IZ,IEL)
         UNZ (IX,IZ,1,IEL) =  C(IX,  1,IZ,IEL)
         UNZ (IX,IZ,3,IEL) = -C(IX,ly1,IZ,IEL)
  200 CONTINUE
C
C        "T"
C
      CALL VCROSS(A,B,C,XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,NTOT)
      CALL VDOT3 (DOT,A,B,C,A,B,C,NTOT)
      DO 300 IEL=1,NEL
      DO 300 IX=1,lx1
      DO 300 IY=1,ly1
         WEIGHT=WXM1(IX)*WYM1(IY)
         AREA(IX,IY,5,IEL) = SQRT(DOT(IX,IY,  1,IEL))*WEIGHT
         AREA(IX,IY,6,IEL) = SQRT(DOT(IX,IY,lz1,IEL))*WEIGHT
         UNX (IX,IY,5,IEL) = -A(IX,IY,  1,IEL)
         UNX (IX,IY,6,IEL) =  A(IX,IY,lz1,IEL)
         UNY (IX,IY,5,IEL) = -B(IX,IY,  1,IEL)
         UNY (IX,IY,6,IEL) =  B(IX,IY,lz1,IEL)
         UNZ (IX,IY,5,IEL) = -C(IX,IY,  1,IEL)
         UNZ (IX,IY,6,IEL) =  C(IX,IY,lz1,IEL)
  300 CONTINUE
C
      CALL UNITVEC (UNX,UNY,UNZ,NSRF)
C
C     COMPUTE UNIT TANGENT T1
C
      DO 600 IEL=1,NEL
      DO 600 IFC=1,NFACE
      IF (IFC.EQ.1 .OR. IFC.EQ.6) THEN
         CALL FACEXV (T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL),
     $                T1Z(1,1,IFC,IEL),
     $                XRM1(1,1,1,IEL),YRM1(1,1,1,IEL),
     $                ZRM1(1,1,1,IEL),IFC,0)
      ELSEIF (IFC.EQ.2 .OR. IFC.EQ.5) THEN
         CALL FACEXV (T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL),
     $                T1Z(1,1,IFC,IEL),
     $                XSM1(1,1,1,IEL),YSM1(1,1,1,IEL),
     $                ZSM1(1,1,1,IEL),IFC,0)
      ELSE
         CALL FACEXV (T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL),
     $                T1Z(1,1,IFC,IEL),
     $                XTM1(1,1,1,IEL),YTM1(1,1,1,IEL),
     $                ZTM1(1,1,1,IEL),IFC,0)
      ENDIF
  600 CONTINUE
C
      CALL UNITVEC (T1X,T1Y,T1Z,NSRF)
C
C     COMPUTE UNIT TANGENT T2  ( T2 = Normal X T1 )
C
      DO 700 IEL=1,NEL
      DO 700 IFC=1,NFACE
         CALL VCROSS (T2X(1,1,IFC,IEL),T2Y(1,1,IFC,IEL),
     $                T2Z(1,1,IFC,IEL),
     $                UNX(1,1,IFC,IEL),UNY(1,1,IFC,IEL),
     $                UNZ(1,1,IFC,IEL),
     $                T1X(1,1,IFC,IEL),T1Y(1,1,IFC,IEL),
     $                T1Z(1,1,IFC,IEL),NXY1)
  700 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine vdot3 (dot,u1,u2,u3,v1,v2,v3,n)
C
C     Compute a Cartesian vector dot product. 3-d version
C
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
C
C
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)
C
C     Compute a Cartesian vector cross product.
C
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
      DIMENSION W1(1),W2(1),W3(1)
C
C
      DO 100 I=1,N
         U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
         U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
         U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine setwgtr (wgtr1,wgtr2,wgtr3,wgtr4,nel)
C
      INCLUDE 'LVAR'
      INCLUDE 'INTEG'

      logical ifaxis
C
C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3 
C           share the same array structure in Scratch Common /SCRNS/.
C
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
C
      DIMENSION WGTR1(LX1,1)
     $ ,        WGTR2(LY1,1)
     $ ,        WGTR3(LX1,1)
     $ ,        WGTR4(LY1,1)

      ifaxis=.false.
C
      IF (IFAXIS) THEN
         DO 100 IEL=1,NEL
            DO 120 IX=1,lx1
               WGTR1(IX,IEL) = YM1(IX,  1,1,IEL) * WXM1(IX)
               WGTR3(IX,IEL) = YM1(IX,ly1,1,IEL) * WXM1(IX)
  120       CONTINUE
            IF ( IFRZER(IEL) ) THEN
               IY = 1
               WGTR2(IY,IEL) = YSM1(lx1,IY,1,IEL) * WAM1(IY)
               WGTR4(IY,IEL) = YSM1(  1,IY,1,IEL) * WAM1(IY)
               DO 160 IY=2,ly1
                  DNR = 1. + ZAM1(IY)
                  WGTR2(IY,IEL) = YM1(lx1,IY,1,IEL) * WAM1(IY) / DNR
                  WGTR4(IY,IEL) = YM1(  1,IY,1,IEL) * WAM1(IY) / DNR
  160          CONTINUE
            ELSE
               DO 180 IY=1,ly1
                  WGTR2(IY,IEL) = YM1(lx1,IY,1,IEL) * WYM1(IY)
                  WGTR4(IY,IEL) = YM1(  1,IY,1,IEL) * WYM1(IY)
  180          CONTINUE
            ENDIF
  100    CONTINUE
      ELSE
         DO 200 IEL=1,NEL
            CALL COPY (WGTR1(1,IEL),WXM1,lx1)
            CALL COPY (WGTR2(1,IEL),WYM1,ly1)
            CALL COPY (WGTR3(1,IEL),WXM1,lx1)
            CALL COPY (WGTR4(1,IEL),WYM1,ly1)
  200    CONTINUE
      ENDIF
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE RZERO3 (A,B,C,N)
      DIMENSION A(1),B(1),C(1)
      DO 100 I=1,N
         A(I)=0.0
         B(I)=0.0
         C(I)=0.0
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UNITVEC (X,Y,Z,N)
      DIMENSION X(1),Y(1),Z(1)
      DO 100 I=1,N
      XLNGTH = SQRT( X(I)**2 + Y(I)**2 + Z(I)**2 )
      IF (XLNGTH.NE.0.0) THEN
         X(I) = X(I)/XLNGTH
         Y(I) = Y(I)/XLNGTH
         Z(I) = Z(I)/XLNGTH
      ENDIF
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine xyzrst (xrm1,yrm1,zrm1,xsm1,ysm1,zsm1,
     $                   XTM1,YTM1,ZTM1,IFAXIS,nel)
C-----------------------------------------------------------------------
C
C     Compute global-to-local derivatives on mesh 1.
C
C-----------------------------------------------------------------------
      INCLUDE 'LVAR'
      INCLUDE 'INTEG'
C
      DIMENSION XRM1(LX1,LY1,LZ1,1),YRM1(LX1,LY1,LZ1,1)
     $        , ZRM1(LX1,LY1,LZ1,1),XSM1(LX1,LY1,LZ1,1)
     $        , YSM1(LX1,LY1,LZ1,1),ZSM1(LX1,LY1,LZ1,1)
     $        , XTM1(LX1,LY1,LZ1,1),YTM1(LX1,LY1,LZ1,1)
     $        , ZTM1(LX1,LY1,LZ1,1)
      LOGICAL IFAXIS

      ifaxis=.false.
C
      NXY1=lx1*ly1
      NYZ1=ly1*lz1
C
      DO 100 IEL=1,NEL
C
c     IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )

C
      CALL MXM (DXM1,lx1,XM1(1,1,1,IEL),lx1,XRM1(1,1,1,IEL),NYZ1)
      CALL MXM (DXM1,lx1,YM1(1,1,1,IEL),lx1,YRM1(1,1,1,IEL),NYZ1)
      CALL MXM (DXM1,lx1,ZM1(1,1,1,IEL),lx1,ZRM1(1,1,1,IEL),NYZ1)
C
      DO 10 IZ=1,lz1
      CALL MXM (XM1(1,1,IZ,IEL),lx1,DYTM1,ly1,XSM1(1,1,IZ,IEL),ly1)
      CALL MXM (YM1(1,1,IZ,IEL),lx1,DYTM1,ly1,YSM1(1,1,IZ,IEL),ly1)
      CALL MXM (ZM1(1,1,IZ,IEL),lx1,DYTM1,ly1,ZSM1(1,1,IZ,IEL),ly1)
   10 CONTINUE
C
      IF (ldim.EQ.3) THEN
         CALL MXM (XM1(1,1,1,IEL),NXY1,DZTM1,lz1,XTM1(1,1,1,IEL),lz1)
         CALL MXM (YM1(1,1,1,IEL),NXY1,DZTM1,lz1,YTM1(1,1,1,IEL),lz1)
         CALL MXM (ZM1(1,1,1,IEL),NXY1,DZTM1,lz1,ZTM1(1,1,1,IEL),lz1)
      ELSE
         CALL RZERO (XTM1(1,1,1,IEL),NXY1)
         CALL RZERO (YTM1(1,1,1,IEL),NXY1)
         CALL RONE  (ZTM1(1,1,1,IEL),NXY1)
      ENDIF
C
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      function vlsc3(x,y,b,n)

c     local inner product, with weight

      real x(1),y(1),b(1)

      vlsc3=0.

      do i=1,n
         vlsc3=vlsc3+x(i)*y(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rone(a,n)
      DIMENSION A(1)
      DO 100 I = 1, N
 100     A(I ) = 1.0
      return
      END
c-----------------------------------------------------------------------
      subroutine glmapm1(nel)
C-----------------------------------------------------------------------
C
C     Routine to generate mapping data based on mesh 1
C     (Gauss-Legendre Lobatto meshes).
C
C         XRM1,  YRM1,  ZRM1   -   dx/dr, dy/dr, dz/dr
C         XSM1,  YSM1,  ZSM1   -   dx/ds, dy/ds, dz/ds
C         XTM1,  YTM1,  ZTM1   -   dx/dt, dy/dt, dz/dt
C         RXM1,  RYM1,  RZM1   -   dr/dx, dr/dy, dr/dz
C         SXM1,  SYM1,  SZM1   -   ds/dx, ds/dy, ds/dz
C         TXM1,  TYM1,  TZM1   -   dt/dx, dt/dy, dt/dz
C         JACM1                -   Jacobian
C
C-----------------------------------------------------------------------
      include 'LVAR'
      include 'INTEG'
C
C     Note: Subroutines GLMAPM1, GEODAT1, AREA2, SETWGTR and AREA3
C           share the same array structure in Scratch Common /SCRNS/.
C
      COMMON /SCRNS/ XRM1(LX1,LY1,LZ1,LELT)
     $ ,             YRM1(LX1,LY1,LZ1,LELT)
     $ ,             XSM1(LX1,LY1,LZ1,LELT)
     $ ,             YSM1(LX1,LY1,LZ1,LELT)
     $ ,             XTM1(LX1,LY1,LZ1,LELT)
     $ ,             YTM1(LX1,LY1,LZ1,LELT)
     $ ,             ZRM1(LX1,LY1,LZ1,LELT)
      COMMON /CTMP1/ ZSM1(LX1,LY1,LZ1,LELT)
     $ ,             ZTM1(LX1,LY1,LZ1,LELT)
C
      logical ifaxis
      ifaxis=.false.

      NXY1  = lx1*ly1
      NYZ1  = ly1*lz1
      NXYZ1 = lx1*ly1*lz1
      NTOT1 = NXYZ1*NEL
C
      CALL XYZRST (XRM1,YRM1,ZRM1,XSM1,YSM1,ZSM1,XTM1,YTM1,ZTM1,
     $             IFAXIS,nel)
C
      IF (ldim.EQ.2) THEN
         CALL RZERO   (JACM1,NTOT1)
         CALL ADDCOL3 (JACM1,XRM1,YSM1,NTOT1)
         CALL SUBCOL3 (JACM1,XSM1,YRM1,NTOT1)
         CALL COPY    (RXM1,YSM1,NTOT1)
         CALL COPY    (RYM1,XSM1,NTOT1)
         CALL CHSIGN  (RYM1,NTOT1)
         CALL COPY    (SXM1,YRM1,NTOT1)
         CALL CHSIGN  (SXM1,NTOT1)
         CALL COPY    (SYM1,XRM1,NTOT1)
         CALL RZERO   (RZM1,NTOT1)
         CALL RZERO   (SZM1,NTOT1)
         CALL RONE    (TZM1,NTOT1)
      ELSE
         CALL RZERO   (JACM1,NTOT1)
         CALL ADDCOL4 (JACM1,XRM1,YSM1,ZTM1,NTOT1)
         CALL ADDCOL4 (JACM1,XTM1,YRM1,ZSM1,NTOT1)
         CALL ADDCOL4 (JACM1,XSM1,YTM1,ZRM1,NTOT1)
         CALL SUBCOL4 (JACM1,XRM1,YTM1,ZSM1,NTOT1)
         CALL SUBCOL4 (JACM1,XSM1,YRM1,ZTM1,NTOT1)
         CALL SUBCOL4 (JACM1,XTM1,YSM1,ZRM1,NTOT1)
         CALL ASCOL5  (RXM1,YSM1,ZTM1,YTM1,ZSM1,NTOT1)
         CALL ASCOL5  (RYM1,XTM1,ZSM1,XSM1,ZTM1,NTOT1)
         CALL ASCOL5  (RZM1,XSM1,YTM1,XTM1,YSM1,NTOT1)
         CALL ASCOL5  (SXM1,YTM1,ZRM1,YRM1,ZTM1,NTOT1)
         CALL ASCOL5  (SYM1,XRM1,ZTM1,XTM1,ZRM1,NTOT1)
         CALL ASCOL5  (SZM1,XTM1,YRM1,XRM1,YTM1,NTOT1)
         CALL ASCOL5  (TXM1,YRM1,ZSM1,YSM1,ZRM1,NTOT1)
         CALL ASCOL5  (TYM1,XSM1,ZRM1,XRM1,ZSM1,NTOT1)
         CALL ASCOL5  (TZM1,XRM1,YSM1,XSM1,YRM1,NTOT1)
      ENDIF
C
      kerr = 0
      DO 500 ie=1,NEL
         CALL CHKJAC(JACM1(1,1,1,ie),NXYZ1,ie,xm1(1,1,1,ie),
     $ ym1(1,1,1,ie),zm1(1,1,1,ie),ldim,ierr)
         if (ierr.ne.0) kerr = kerr+1
  500 CONTINUE
      kerr = iglsum(kerr,1)
      if (kerr.gt.0) then
         if (nid.eq.0) write(6,*) 'Jac error 1, setting p66=4, ifxyo=t'
         call exitt
      endif

      call invers2(jacmi,jacm1,ntot1)

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)

      real a(1),b(1),c(1)

      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chkjac(jac,n,iel,X,Y,Z,ND,IERR)
c
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
C
C     Check the array JAC for a change in sign.
C
      REAL JAC(N),x(1),y(1),z(1)
c
      ierr = 1
      SIGN = JAC(1)
      DO 100 I=2,N
         IF (SIGN*JAC(I).LE.0.0) THEN
            WRITE(6,101) mid,I,iel
            write(6,*) jac(i-1),jac(i)
            if (ldim.eq.3) then
               write(6,7) mid,x(i-1),y(i-1),z(i-1)
               write(6,7) mid,x(i),y(i),z(i)
            else
               write(6,7) mid,x(i-1),y(i-1)
               write(6,7) mid,x(i),y(i)
            endif
    7       format(i5,' xyz:',1p3e14.5)
            return
         ENDIF
  100 CONTINUE
  101 FORMAT(//,i5,2x
     $ ,'ERROR:  Vanishing Jacobian near',i7,'th node of element'
     $ ,I10,'.')
c
c
      ierr = 0
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine subcol3(a,b,c,n)
      REAL A(1),B(1),C(1)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE FACEXV (A1,A2,A3,B1,B2,B3,IFACE1,IOP)
C
C     IOP = 0
C     Extract vector (A1,A2,A3) from (B1,B2,B3) on face IFACE1.
C
C     IOP = 1
C     Extract vector (B1,B2,B3) from (A1,A2,A3) on face IFACE1.
C
C     A1, A2, A3 have the (NX,NY,NFACE) data structure
C     B1, B2, B3 have the (NX,NY,NZ)    data structure
C     IFACE1 is in the preprocessor notation
C     IFACE  is the dssum notation.
C
      include 'LVAR'
      include 'INTEG'
C
      DIMENSION A1(LX1,LY1),A2(LX1,LY1),A3(LX1,LY1),
     $          B1(LX1,LY1,LZ1),B2(LX1,LY1,LZ1),B3(LX1,LY1,LZ1)
C
      CALL DSSET(lx1,ly1,lz1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
      I = 0
C
      IF (IOP.EQ.0) THEN
         DO 100 J2=JS2,JF2,JSKIP2
         DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A1(I,1) = B1(J1,J2,1)
            A2(I,1) = B2(J1,J2,1)
            A3(I,1) = B3(J1,J2,1)
  100    CONTINUE
      ELSE
         DO 150 J2=JS2,JF2,JSKIP2
         DO 150 J1=JS1,JF1,JSKIP1
            I = I+1
            B1(J1,J2,1) = A1(I,1)
            B2(J1,J2,1) = A2(I,1)
            B3(J1,J2,1) = A3(I,1)
  150    CONTINUE
      ENDIF
C
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine col3(a,b,c,n)

      real a(1),b(1),c(1)

      do i=1,n
         a(i)=b(i)*c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dsset(nx,ny,nz)
C
C     Set up arrays IXCN,ESKIP,SKPDAT,NEDG,NOFFST for new NX,NY,NZ
C
      include 'LVAR'
      include 'INTEG'

      INTEGER NXO,NYO,NZO
      SAVE    NXO,NYO,NZO
      DATA    NXO,NYO,NZO /3*0/
C
C     Check if element surface counters are already set from last call...
C
      IF (NXO.EQ.NX.AND.NYO.EQ.NY.AND.NZO.EQ.NZ) RETURN
C
C     else, proceed....
C
      NXO = NX
      NYO = NY
      NZO = NZ
C
C     Establish corner to elemental node number mappings
C
      IC=0
      DO 10 ICZ=0,1
      DO 10 ICY=0,1
      DO 10 ICX=0,1
C       Supress vectorization to
c        IF(ICX.EQ.0)DUMMY=0
c        IF(ICX.EQ.1)DUMMY=1
c        DUMMY2=DUMMY2+DUMMY
        IC=IC+1
        IXCN(IC)= 1 + (NX-1)*ICX + NX*(NY-1)*ICY + NX*NY*(NZ-1)*ICZ
   10   CONTINUE
C
C     Assign indices for direct stiffness summation of arbitrary faces.
C
C
C     Y-Z Planes (Faces 1 and 2)
C
      SKPDAT(1,1)=1
      SKPDAT(2,1)=NX*(NY-1)+1
      SKPDAT(3,1)=NX
      SKPDAT(4,1)=1
      SKPDAT(5,1)=NY*(NZ-1)+1
      SKPDAT(6,1)=NY
C
      SKPDAT(1,2)=1             + (NX-1)
      SKPDAT(2,2)=NX*(NY-1)+1   + (NX-1)
      SKPDAT(3,2)=NX
      SKPDAT(4,2)=1
      SKPDAT(5,2)=NY*(NZ-1)+1
      SKPDAT(6,2)=NY
C
C     X-Z Planes (Faces 3 and 4)
C
      SKPDAT(1,3)=1
      SKPDAT(2,3)=NX
      SKPDAT(3,3)=1
      SKPDAT(4,3)=1
      SKPDAT(5,3)=NY*(NZ-1)+1
      SKPDAT(6,3)=NY
C
      SKPDAT(1,4)=1           + NX*(NY-1)
      SKPDAT(2,4)=NX          + NX*(NY-1)
      SKPDAT(3,4)=1
      SKPDAT(4,4)=1
      SKPDAT(5,4)=NY*(NZ-1)+1
      SKPDAT(6,4)=NY
C
C     X-Y Planes (Faces 5 and 6)
C
      SKPDAT(1,5)=1
      SKPDAT(2,5)=NX
      SKPDAT(3,5)=1
      SKPDAT(4,5)=1
      SKPDAT(5,5)=NY
      SKPDAT(6,5)=1
C
      SKPDAT(1,6)=1           + NX*NY*(NZ-1)
      SKPDAT(2,6)=NX          + NX*NY*(NZ-1)
      SKPDAT(3,6)=1
      SKPDAT(4,6)=1
      SKPDAT(5,6)=NY
      SKPDAT(6,6)=1
C
C     Set up skip indices for each of the 12 edges
C
C         Note that NXY = NX*NY even for 2-D since
C         this branch does not apply to the 2D case anyway.
C
C     ESKIP(*,1) = start location
C     ESKIP(*,2) = end
C     ESKIP(*,3) = stride
C
      NXY=NX*NY
      ESKIP( 1,1) = IXCN(1) + 1
      ESKIP( 1,2) = IXCN(2) - 1
      ESKIP( 1,3) = 1
      ESKIP( 2,1) = IXCN(3) + 1
      ESKIP( 2,2) = IXCN(4) - 1
      ESKIP( 2,3) = 1
      ESKIP( 3,1) = IXCN(5) + 1
      ESKIP( 3,2) = IXCN(6) - 1
      ESKIP( 3,3) = 1
      ESKIP( 4,1) = IXCN(7) + 1
      ESKIP( 4,2) = IXCN(8) - 1
      ESKIP( 4,3) = 1
      ESKIP( 5,1) = IXCN(1) + NX
      ESKIP( 5,2) = IXCN(3) - NX
      ESKIP( 5,3) = NX
      ESKIP( 6,1) = IXCN(2) + NX
      ESKIP( 6,2) = IXCN(4) - NX
      ESKIP( 6,3) = NX
      ESKIP( 7,1) = IXCN(5) + NX
      ESKIP( 7,2) = IXCN(7) - NX
      ESKIP( 7,3) = NX
      ESKIP( 8,1) = IXCN(6) + NX
      ESKIP( 8,2) = IXCN(8) - NX
      ESKIP( 8,3) = NX
      ESKIP( 9,1) = IXCN(1) + NXY
      ESKIP( 9,2) = IXCN(5) - NXY
      ESKIP( 9,3) = NXY
      ESKIP(10,1) = IXCN(2) + NXY
      ESKIP(10,2) = IXCN(6) - NXY
      ESKIP(10,3) = NXY
      ESKIP(11,1) = IXCN(3) + NXY
      ESKIP(11,2) = IXCN(7) - NXY
      ESKIP(11,3) = NXY
      ESKIP(12,1) = IXCN(4) + NXY
      ESKIP(12,2) = IXCN(8) - NXY
      ESKIP(12,3) = NXY
C
C     Load reverse direction edge arrays for reverse mappings...
C
      DO 20 IED=1,12
      IEDM=-IED
      ESKIP(IEDM,1) =  ESKIP(IED,2)
      ESKIP(IEDM,2) =  ESKIP(IED,1)
      ESKIP(IEDM,3) = -ESKIP(IED,3)
   20 CONTINUE
C
C     Compute offset for global edge vector given current element
C     dimensions.....
C
C     NGSPED(ITE,ICMP) = number of global (ie, distinct) special edges
C                        of type ITE (1,2, or 3)  for field ICMP.
C
C                        ITE = 1 implies an "X" edge
C                        ITE = 2 implies an "Y" edge
C                        ITE = 3 implies an "Z" edge
C
C     Set up number of nodes along each of the 3 types of edges
C     (endpoints excluded).
C
      NEDG(1)=NX-2
      NEDG(2)=NY-2
      NEDG(3)=NZ-2
C
C
      RETURN
      END
c-----------------------------------------------------------------------
      real function vlsc2(x,y,n)

      real x(1),y(1)

      vlsc2=0.
      do i=1,n
         vlsc2=vlsc2+x(i)*y(i)
      enddo

      return
      end
c-----------------------------------------------------------------------