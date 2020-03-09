c-----------------------------------------------------------------------
      subroutine setgeom(xyz,nel)

      include 'LVAR'
      include 'INTEG'

      real xyz(1)

      n=lxyz*nel

      call copy(xm1,xyz,n)
      call copy(ym1,xyz(n+1),n)
      if (ldim.eq.3) call copy(zm1,xyz(n*2+1),n)

      call glmapm1(nel)
      call geodat1(nel)
      call setrx(nel)

      return
      end
c-----------------------------------------------------------------------
      subroutine geom_check(nel)

      include 'LVAR'
      include 'INTEG'

      write (6,*) 'inside geom_check',nel,lxyz

      do ie=1,nel
      do i=1,lxyz
c        write (6,*) i,ie,bm1(i,1,1,ie),'bm1'
      enddo
      enddo

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
      call initds
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
      subroutine intp_rstd(ju,u,mx,md,if3d,idir) ! GLL->GL interpolation

c     GLL interpolation from mx to md.

c     If idir ^= 0, then apply transpose operator  (md to mx)

      include 'LVAR'

      real    ju(1),u(1)
      logical if3d

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      parameter (ld=2*lxd)
      common /ctmp0/ w(ld**ldim,2)

      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')

      ldw = 2*(ld**ldim)

      call get_int_ptr (i,mx,md)
c
      if (idir.eq.0) then
         call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else
         call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md) ! GLL-->GL pointer

c     Get pointer to jgl() for interpolation pair (mx,md)

      include 'LVAR'

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c
      if (ip.eq.0) then
c
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
c
         call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgl_ptr (ip,mx,md)
c
c     Get pointer to GL-GL interpolation dgl() for pair (mx,md)
c
      include 'LVAR'
c
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c
      parameter (ld=2*lxd)
      common /jgrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pdg (ij)

      if (ip.eq.0) then

         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c
         call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rst(ur,us,ut,u,md,if3d) ! Gauss-->Gauss grad

      include 'LVAR'
      include 'INTEG'

      real    ur(1),us(1),ut(1),u(1)
      logical if3d

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      m0 = md-1
      call get_dgl_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,u,m0,1,dg(ip),dgt(ip))
      else
         call local_grad2(ur,us   ,u,m0,1,dg(ip),dgt(ip))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_dgl(dgl,dgt,mp,np,w)
c
c     Generate derivative from np GL points onto mp GL points
c
c        dgl  = interpolation matrix, mapping from velocity nodes to pressure
c        dgt  = transpose of interpolation matrix
c        w    = work array of size (3*np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
c
      real dgl(mp,np),dgt(np*mp),w(1)
c
c
      iz = 1
      id = iz + np
c
      call zwgl  (w(iz),dgt,np)  ! GL points
      call zwgl  (w(id),dgt,mp)  ! GL points
c
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
         enddo
      enddo
c
      call transpose(dgt,np,dgl,mp)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
      include 'SIZE'            ! need nid
      character*5  avar5,lvar5
      character*10 sub_name10

      if (n.gt.m) then
         write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i12,2(1x,a5),1x,a10)
         call exitti('lim_chk problem. $',n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,m,c)
c
c     This routine evaluates the derivative based on all points
c     in the stencils.  It is more memory efficient than "fd_weights"
c
c     This set of routines comes from the appendix of
c     A Practical Guide to Pseudospectral Methods, B. Fornberg
c     Cambridge Univ. Press, 1996.   (pff)
c
c     Input parameters:
c       xx -- point at wich the approximations are to be accurate
c       x  -- array of x-ordinates:   x(0:n)
c       n  -- polynomial degree of interpolant (# of points := n+1)
c       m  -- highest order of derivative to be approxxmated at xi
c
c     Output:
c       c  -- set of coefficients c(0:n,0:m).
c             c(j,k) is to be applied at x(j) when
c             the kth derivative is approxxmated by a
c             stencil extending over x(0),x(1),...x(n).
c
c
      real x(0:n),c(0:n,0:m)
c
      c1       = 1.
      c4       = x(0) - xx
c
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
c
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo
c     call outmat(c,n+1,m+1,'fdw',n)
      return
      end
c-----------------------------------------------------------------------
      subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)
C
C     -  Spectral interpolation from A to B via tensor products
C     -  scratch arrays: w(na*na*nb + nb*nb*na)
C
C     5/3/00  -- this routine replaces specmp in navier1.f, which
c                has a potential memory problem
C
C
      logical if3d
c
      real b(nb,nb,nb),a(na,na,na)
      real w(ldw)
c
      ltest = na*nb
      if (if3d) ltest = na*na*nb + nb*na*na
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif
      return
      end
c
c-----------------------------------------------------------------------
      subroutine gen_int(jgl,jgt,mp,np,w)
c
c     Generate interpolation from np GLL points to mp GL points
c
c        jgl  = interpolation matrix, mapping from velocity nodes to pressure
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
      real jgl(mp,np),jgt(np*mp),w(1)
c
      iz = 1
      id = iz + np
c
      call zwgll (w(iz),jgt,np)
      call zwgl  (w(id),jgt,mp)
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo
c
      call transpose(jgt,np,jgl,mp)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine initds
C
C          -- Direct Stiffness Initialization Routine --
C
C     Set up required data for packing data on faces of spectral cubes.
C
      INCLUDE 'LVAR'
      INCLUDE 'INTEG'
C
C     Nominal ordering for direct stiffness summation of faces
C
      J=0
      DO 5 IDIM=1,ldim
      DO 5 IFACE=1,2
        J=J+1
         NOMLIS(IFACE,IDIM)=J
    5 CONTINUE
C
C     Assign Ed's numbering scheme to PF's scheme.
C
      EFACE(1)=4
      EFACE(2)=2
      EFACE(3)=1
      EFACE(4)=3
      EFACE(5)=5
      EFACE(6)=6
C
C     Assign inverse of Ed's numbering scheme to PF's scheme.
C
      EFACE1(1)=3
      EFACE1(2)=2
      EFACE1(3)=4
      EFACE1(4)=1
      EFACE1(5)=5
      EFACE1(6)=6
C
C     Assign group designation to each face to determine ordering of indices.
C
      GROUP(1)=0
      GROUP(2)=1
      GROUP(3)=1
      GROUP(4)=0
      GROUP(5)=0
      GROUP(6)=1
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine gen_igll(jgl,jgt,mp,np,w)

c     Generate interpolation from np GLL points to mp GLL points
c
c     jgl   = interpolation matrix, mapping from np to mp
c     jgt   = transpose of interpolation matrix
c     w     = work array of size (np+mp)
c     np,mp = number of points on each grid

      real jgl(mp,np),jgt(np*mp),w(1)

      iz = 1
      id = iz + np

      call zwgll (w(iz),jgt,np)
      call zwgll (w(id),jgt,mp)

      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo

      call transpose(jgt,np,jgl,mp)

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_proj

      include 'LVAR'

      common /sei/ pmat(lx1**2,lx1),ptmat(lx1**2,lx1),wk(lx1**2,5)

      l=1
      do i=lx1-1,2,-1
         call gen_igll(wk(1,1),wk(1,2),i,lx1,wk(1,5))
         call gen_igll(wk(1,3),wk(1,4),lx1,i,wk(1,5))
         call mxm(wk(1,3),lx1,wk(1,1),i,pmat(1,l),lx1)
         call mxm(wk(1,2),lx1,wk(1,4),i,ptmat(1,l),lx1)
         l=l+1
      enddo

      return
      end
c-----------------------------------------------------------------------