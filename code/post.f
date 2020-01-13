c-----------------------------------------------------------------------
      subroutine drivep

      include 'POST'
      include 'MASS'
      include 'PARALLEL'

      character*127 flist

      nelp=512
      nsnap=ns

      ! assigned snapshot range

      isg0=1
      isg1=nsg

      call reade_init(flist,ns)

      inel=1
      ieg1=0

      ns=ls
      nsg=lsg

      call rzero(gram,ns*nsg)
      call rzero(aa,ns*nsg)
      call rzero(bb,ns*nsg)
      call rzero(cc,ns*nsg*nsg)

      do while (ieg1+1.le.nelgv)
         ieg0=ieg1+1
         ieg1=min(ieg1+inel+nelp-1,nelgv)
         nel=ieg1-ieg0+1
         n=lxyz*(ieg1-ieg0+1)
         write (6,*) 'ieg0,ieg1,nelgv',ieg0,ieg1,nelgv
         call reade_dummy(uu,ieg0,ieg1)
         call setmass(mass,wv1,ieg0,ieg1,lxyz)
         call setgg(gram,uu,mass,wvf1,wvf2,ns,nsg,n,ndim,lxyz,nel)
      enddo

      call dump_serial(gram,ns*ns,'ops/gramp ',nid)

      ! dssum to ensure gram is symmetric
      ! shift each row in gram
      ! eigendecomposition here or external process

      inel=1
      ieg1=0

      call rone(visc,n)

      do while (ieg1+1.le.nelgv)
         ieg0=ie1+1
         ieg1=min(inel+nelp-1,nelgv)
         call reade(uu,ieg0,ieg1,'U',flist)

c        call setvisc(visc,w,ieg0,ieg1,lxyz,nid) ! not required
         call setgeom(gfac,w,ieg0,ieg1,lxyz,3*(ndim-1),nid)
         call setmass(mass,w,ieg0,ieg1,lxyz,nid)
         call setconv(rxd) ! TODO: implement

         call setzz(zz,uu,evec,n,ns)
         call setaa(aa,zz,visc,gfac,wvf1,wvf2,ns,nsg,n,ndim)
         call setbb(bb,zz,mass,wvf1,wvf2,ns,nsg,n,ndim)
c        call setc(c,z,wvf1,wvf2,ns,nsg,n,ndim,ndim) ! TODO: implement

c        call setkk(uk,zz,mass,wvf1,wvf2,ns,nsg,n,ndim)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine reade_init(flist,ns)

      character*127 flist

      return
      end
c-----------------------------------------------------------------------
      subroutine reade_dummy(v,ieg0,ieg1)

      include 'POST'

      parameter (ll=lx1*ly1*lz1)

      real v(lxyz,ieg1-ieg0+1,ldim,ns)

      do is=1,ns
      do idim=1,ndim
      do ie=ieg0,ieg1
         call copy(v(1,ie-ieg0+1,idim,is),us0(1+(ie-1)*ll,idim,is),ll)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine reade(u,ieg0,ieg1,clist,flist)

      real u(1)
      logical ifread(0:4)
      character*127 clist
      character*127 fname

      call parse_clist(ifread,clist)
c     call reade_helper(u,ieg0,ieg1,ifread,fname)

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_clist(ifread,clist)

      logical ifread(0:4)
      character*128 clist

      ! no support for thermal for now

      do i=0,4
         ifread(i)=.false.
      enddo

c     if (clist matches 'U') then
         ifread(1)=.true.
c     else
c        call exitti('unsupported clist$',1)
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setzz(z,u,evec,w1,w2,n,ns,nsg)

      real z(n,ns), u(n,ns)
      real evec(ns,nsg)
      real w1(n),w2(n)

      do isg=1,nsg
         is=iglls(isg)
         call mxm(u,n,evec(1,isg),ns,w1,1)
         call gop(w1,w2,'+  ',n)
         if (is.gt.0 .and. is.le.ns) call copy(z(1,is),w1,n)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function ilgls(is)

      ilgls=0

      return
      end
c-----------------------------------------------------------------------
      function iglls(isg)

      iglls=0

      return
      end
c-----------------------------------------------------------------------
      subroutine setaa(a,z,w1,w2,ns,nsg,n,ndim)

      real a(ns,nsg),z(n,ndim,ns),w1(n,ndim,ns),w2(n,ndim,ns)

      call copy(w1,z,n*ndim*ns)
      call copy(w2,z,n*ndim*ns)

      do i=1,ns*ndim
c        call 'axhelm'(w2(1,i,1))
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            a(is,js+ns*ioff)=a(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setbb(b,z,mass,w1,w2,ns,nsg,n,ndim)

      real b(ns,nsg),z(n,ndim,ns),mass(n),w1(n,ndim,ns),w2(n,ndim,ns)

      call copy(w1,z,n*ndim*ns)
      call copy(w2,z,n*ndim*ns)

      do i=1,ns*ndim
         call col2(w2(1,i,1),mass,n)
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            b(is,js+ns*ioff)=b(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc(c,z,t,w1,w2,w3,w4,ns,nsg,n,mdim,ndim)

      real c(ns,nsg,nsg),z(n,mdim,ns),t(n,ndim,ns),
     $     w1(n,mdim,ns),w2(n,ndim,ns),w3(n,mdim,ns),w4(n,mdim,ns)

      call copy(w1,t,m*mdim*ns)
      call copy(w2,z,n*ndim*ns)
      call copy(w3,t,n*mdim*ns)

      do i=1,ns*ndim
c        call col2(w2(1,i,1),mass,n)
         ! apply grad on w1
      enddo

      do jof=0,nsg/ns-1 ! assume ns is the same across all processors
         do iof=0,nsg/ns-1
            do ks=1,ns
            do js=1,ns
c              call conv(w4,w1,w2) ! w4= mass * (w1 * w4)
               do is=1,ns
                  c(is,js+ns*iof,ks+ns*jof)=c(is,js+ns*iof,ks+ns*jof)
     $               +vlsc2(w3(1,1,is),w4(1,1,js),n*ndim)
               enddo
            enddo
            enddo
            
            call shift(w3,n*ndim*ns)
         enddo
         call shift(w2,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setgeom(gfac,w,ieg0,ieg1,lxyz,ng)!,nid)

      include 'SIZE' ! nid
      include 'GEOM'
      include 'PARALLEL'

      real gfac(lxyz,ieg1-ieg0+1,6),w(lxyz*(ieg1-ieg0+1)*6)

      ! TODO: set g1..6m1

      call rzero(gfac,lxyz*(ieg1-ieg0+1)*ng)

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
            call copy(gfac(1,ieg-ieg0+1,1),g1m1(1,1,1,ie),lxyz)
            call copy(gfac(1,ieg-ieg0+1,2),g2m1(1,1,1,ie),lxyz)
            call copy(gfac(1,ieg-ieg0+1,4),g4m1(1,1,1,ie),lxyz)
            if (ng.eq.6) then
               call copy(gfac(1,ieg-ieg0+1,3),g3m1(1,1,1,ie),lxyz)
               call copy(gfac(1,ieg-ieg0+1,5),g5m1(1,1,1,ie),lxyz)
               call copy(gfac(1,ieg-ieg0+1,6),g6m1(1,1,1,ie),lxyz)
            endif
         endif
      enddo

      call gop(gfac,w,'+  ',lxyz*(ieg1-ieg0+1)*ng)

      return
      end
c-----------------------------------------------------------------------
      subroutine setvisc(visc,w,ieg0,ieg1,lxyz)!,nid)

      include 'SIZE'
      include 'PARALLEL'

      real visc(lxyz,ieg1-ieg0+1), w(lxyz*(ieg1-ieg0+1))

      call rzero(visc,lxyz*(ieg1-ieg0+1))

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
c           call copy(visc(1,i-ieg0+1),h1(1,1,1,ie),lxyz)
         endif
      enddo

      call gop(visc,w,'+  ',lxyz*(ieg1-ieg0+1))

      return
      end
c-----------------------------------------------------------------------
      subroutine setmass(mass,w,ieg0,ieg1,lxyz)!,nid)

      include 'SIZE' ! nid
      include 'MASS'
      include 'PARALLEL'

      real mass(lxyz,ieg1-ieg0+1), w(lxyz*(ieg1-ieg0+1))

      ! TODO: set bm1

      call rzero(mass,lxyz*(ieg1-ieg0+1))

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
            call copy(mass(1,ieg-ieg0+1),bm1(1,1,1,ie),lxyz)
         endif
      enddo

      call gop(mass,w,'+  ',lxyz*(ieg1-ieg0+1))

      return
      end
c-----------------------------------------------------------------------
      subroutine setconv(rxd,w,ieg0,ieg1,lxyz)!,ndim,nid)

      include 'SIZE' ! ndim & nid
      include 'GEOM'
      include 'PARALLEL'

      real rxd(lxyz,ndim*ndim,ieg1-ieg0+1), w1(lxyz,9,(ieg1-ieg0+1)),
     $     w(lxyz*(ieg1-ieg0+1))

c     logical if3d

      call rzero(w,lxyz*(ieg1-ieg0+1))
c     if3d=ndim.eq.3

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
            if (ndim.eq.2) then
               call copy(w1(1,1,ieg-ieg0+1),rxm1(1,1,1,ie),lxyz)
               call copy(w1(1,2,ieg-ieg0+1),rym1(1,1,1,ie),lxyz)
               call copy(w1(1,3,ieg-ieg0+1),sxm1(1,1,1,ie),lxyz)
               call copy(w1(1,4,ieg-ieg0+1),sym1(1,1,1,ie),lxyz)
            else
               call copy(w1(1,1,ieg-ieg0+1),rxm1(1,1,1,ie),lxyz)
               call copy(w1(1,2,ieg-ieg0+1),rym1(1,1,1,ie),lxyz)
               call copy(w1(1,3,ieg-ieg0+1),rzm1(1,1,1,ie),lxyz)
               call copy(w1(1,4,ieg-ieg0+1),sxm1(1,1,1,ie),lxyz)
               call copy(w1(1,5,ieg-ieg0+1),sym1(1,1,1,ie),lxyz)
               call copy(w1(1,6,ieg-ieg0+1),rzm1(1,1,1,ie),lxyz)
               call copy(w1(1,7,ieg-ieg0+1),txm1(1,1,1,ie),lxyz)
               call copy(w1(1,8,ieg-ieg0+1),tym1(1,1,1,ie),lxyz)
               call copy(w1(1,9,ieg-ieg0+1),tzm1(1,1,1,ie),lxyz)
            endif
         endif
      enddo

      call gop(w1,w2,'+  ',lxyz*(ieg1-ieg0+1)*ndim*ndim)

      do ie=1,ieg1-ieg0+1
         call intp_rstd(rxd(1,1,ie),w1(1,1,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         call intp_rstd(rxd(1,2,ie),w1(1,2,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         call intp_rstd(rxd(1,3,ie),w1(1,3,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         call intp_rstd(rxd(1,4,ie),w1(1,4,ie),lx1,lxd,if3d,0) ! 0 --> fwd

         if (ndim.eq.3) then
            call intp_rstd(rxd(1,5,ie),w1(1,5,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,6,ie),w1(1,6,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,7,ie),w1(1,7,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,8,ie),w1(1,8,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,9,ie),w1(1,9,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setgg(gram,u,mass,w1,w2,ns,nsg,n,ndim)

      real gram(ns,nsg),u(n,ndim,ns),mass(n),w1(n,ndim,ns),w2(n,ndim,ns)

      write (6,*) 'ns,nsg,n,ndim',ns,nsg,n,ndim

      call copy(w1,u,n*ndim*ns)
      call copy(w2,u,n*ndim*ns)

      do i=1,ns*ndim
         call col2(w2(1,i,1),mass,n)
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            gram(is,js+ns*ioff)=gram(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine shift(u,n)

      real u(n)

      return
      end
c-----------------------------------------------------------------------
c     subroutine ax(au,u,helm1,helm2,imesh,isd)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
c     include 'SIZE'
c     include 'WZ'
c     include 'DXYZ'
c     include 'GEOM'
c     include 'MASS'
c     include 'INPUT'
c     include 'PARALLEL'
c     include 'CTIMER'
c
c     COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
c     COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
c     LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
c
c     REAL           AU    (LX1,LY1,LZ1,1)
c    $ ,             U     (LX1,LY1,LZ1,1)
c    $ ,             HELM1 (LX1,LY1,LZ1,1)
c    $ ,             HELM2 (LX1,LY1,LZ1,1)
c     COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
c    $ ,             DUDS  (LX1,LY1,LZ1)
c    $ ,             DUDT  (LX1,LY1,LZ1)
c    $ ,             TMP1  (LX1,LY1,LZ1)
c    $ ,             TMP2  (LX1,LY1,LZ1)
c    $ ,             TMP3  (LX1,LY1,LZ1)

c     REAL           TM1   (LX1,LY1,LZ1)
c     REAL           TM2   (LX1,LY1,LZ1)
c     REAL           TM3   (LX1,LY1,LZ1)
c     REAL           DUAX  (LX1)
c     REAL           YSM1  (LX1)
c     EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

c     integer e

c     naxhm = naxhm + 1
c     etime1 = dnekclock()

c     nel=nelt
c     if (imesh.eq.1) nel=nelv

c     NXY=lx1*ly1
c     NYZ=ly1*lz1
c     NXZ=lx1*lz1
c     NXYZ=lx1*ly1*lz1
c     NTOT=NXYZ*NEL

c     IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
c     CALL RZERO (AU,NTOT)

c     do 100 e=1,nel
c
c       if (ifaxis) call setaxdy ( ifrzer(e) )
C
c       IF (ldim.EQ.2) THEN
c
c       2-d case ...............
c
c          if (iffast(e)) then
c
c          Fast 2-d mode: constant properties and undeformed element
c
c          h1 = helm1(1,1,1,e)
c          call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
c          call mxm   (u(1,1,1,e),lx1,wddyt,ly1,tm2,ly1)
c          call col2  (tm1,g4m1(1,1,1,e),nxyz)
c          call col2  (tm2,g5m1(1,1,1,e),nxyz)
c          call add3  (au(1,1,1,e),tm1,tm2,nxyz)
c          call cmult (au(1,1,1,e),h1,nxyz)
c
c          else
c
c
c          call mxm  (dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
c          call mxm  (u(1,1,1,e),lx1,dytm1,ly1,duds,ly1)
c          call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
c          call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
c          if (ifdfrm(e)) then
c             call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
c             call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
c          endif
c          call col2 (tmp1,helm1(1,1,1,e),nxyz)
c          call col2 (tmp2,helm1(1,1,1,e),nxyz)
c          call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
c          call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
c          call add2 (au(1,1,1,e),tm1,nxyz)
c          call add2 (au(1,1,1,e),tm2,nxyz)

c       endif
c
c       else
c
c       3-d case ...............
c
c          if (iffast(e)) then
c
c          Fast 3-d mode: constant properties and undeformed element
c
c          h1 = helm1(1,1,1,e)
c          call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
c          do 5 iz=1,lz1
c          call mxm   (u(1,1,iz,e),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
c5         continue
c          call mxm   (u(1,1,1,e),nxy,wddzt,lz1,tm3,lz1)
c          call col2  (tm1,g4m1(1,1,1,e),nxyz)
c          call col2  (tm2,g5m1(1,1,1,e),nxyz)
c          call col2  (tm3,g6m1(1,1,1,e),nxyz)
c          call add3  (au(1,1,1,e),tm1,tm2,nxyz)
c          call add2  (au(1,1,1,e),tm3,nxyz)
c          call cmult (au(1,1,1,e),h1,nxyz)
c
c          else
c
c
c          call mxm(dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
c          do 10 iz=1,lz1
c             call mxm(u(1,1,iz,e),lx1,dytm1,ly1,duds(1,1,iz),ly1)
c  10      continue
c          call mxm     (u(1,1,1,e),nxy,dztm1,lz1,dudt,lz1)
c          call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
c          call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
c          call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
c          if (ifdfrm(e)) then
c             call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
c             call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
c             call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
c             call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
c             call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
c             call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
c          endif
c          call col2 (tmp1,helm1(1,1,1,e),nxyz)
c          call col2 (tmp2,helm1(1,1,1,e),nxyz)
c          call col2 (tmp3,helm1(1,1,1,e),nxyz)
c          call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
c          do 20 iz=1,lz1
c             call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
c  20      continue
c          call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
c          call add2 (au(1,1,1,e),tm1,nxyz)
c          call add2 (au(1,1,1,e),tm2,nxyz)
c          call add2 (au(1,1,1,e),tm3,nxyz)
c
c          endif
c
c       endif
c
c100  continue
c
c     if (ifh2) call addcol4 (au,helm2,bm1,u,ntot)
c
c     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
c
c     if (ifaxis.and.(isd.eq.2)) then
c        do 200 e=1,nel
c
c           if (ifrzer(e)) then
c              call mxm(u  (1,1,1,e),lx1,datm1,ly1,duax,1)
c              call mxm(ym1(1,1,1,e),lx1,datm1,ly1,ysm1,1)
c           endif
c
c           do 190 j=1,ly1
c           do 190 i=1,lx1
c               if (ym1(i,j,1,e).ne.0.) then
c                 if (ifrzer(e)) then
c                    term1 = 0.0
c                    if(j.ne.1) 
c    $             term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
c                    term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i)
c    $                       *jacm1(i,1,1,e)/ysm1(i)
c                 else
c                  term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
c                    term2 = 0.
c                 endif
c                 au(i,j,1,e) = au(i,j,1,e)
c    $                          + helm1(i,j,1,e)*(term1+term2)
c               endif
c 190       continue
c 200    continue
c     endif
c
c     taxhm=taxhm+(dnekclock()-etime1)
c     return
c     end
c-----------------------------------------------------------------------
c     subroutine setk(uk,u,z,mass,w1,w2,ns,nsg,n,ndim)
c
c     real uk(ns,nsg),u(n,ndim,ns),z(n,ndim,ns),
c    $     mass(n),w1(n,ndim,ns),w2(n,ndim,ns)
c
c     call copy(w1,z,n*ndim*ns)
c     call copy(w2,u,n*ndim*ns)
c
c     do i=1,ns*ndim
c        call col2(w2(1,i,1),mass,n)
c     enddo
c
c     do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
c        do js=1,ns
c        do is=1,ns
c           b(is,js+ns*ioff)=b(is,js+ns*ioff)
c    $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
c        enddo
c        enddo
c        
c        call shift(w1,n*ndim*ns)
c     enddo
c
c     return
c     end
c-----------------------------------------------------------------------
