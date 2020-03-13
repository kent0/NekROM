c-----------------------------------------------------------------------
      subroutine setops(a,b,c,wk,t,u,nel,nb,mdim,ndim,ifm1)

      real a(1),b(1),c(1),t(1),u(1)
      logical ifm1

      call bip(b,t,t,nel,nb,mdim)
      if (ifm1) then
         call aip(a,t,t,nel,nb,mdim)
         call cip(c,t,t,u,nel,nb,mdim,ndim)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bip(b,u,v,nel,nb,ndim)

      include 'LVAR'
      include 'INTEG'
      include 'TIMES'

      common /tinteg/ tmp(lxyz,lel,ldim)

      real b(nb,nb)
      real u(lxyz,nel,nb,ndim),v(lxyz,nel,nb,ndim)

      do l=1,ndim
      do k=1,nb
         tt=dnekclock()
         call col3(tmp,u(1,1,k,l),bm1,lxyz*nel)
         time_bop=time_bop+(dnekclock()-tt)
         tt=dnekclock()
         do j=1,nb
            b(j,k)=b(j,k)+vlsc2(tmp,v(1,1,j,l),lxyz*nel)
         enddo
         time_bip=time_bip+(dnekclock()-tt)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine aip(a,u,v,nel,nb,ndim)

      include 'LVAR'
      include 'INTEG'
      include 'TIMES'

      common /tinteg/ tmp(lxyz,lel,ldim)

      real a(nb,nb)
      real u(lxyz,nel,nb,ndim),v(lxyz,nel,nb,ndim)

      do l=1,ndim
      do k=1,nb
         tt=dnekclock()
         call aop(tmp,u(1,1,k,l),nel)
         time_aop=time_aop+(dnekclock()-tt)
         tt=dnekclock()
         do j=1,nb
            a(j,k)=a(j,k)+vlsc2(tmp,v(1,1,j,l),lxyz*nel)
         enddo
         time_aip=time_aip+(dnekclock()-tt)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cip(c,u,v,w,nel,nb,mdim,ndim)

      include 'LVAR'
      include 'INTEG'
      include 'TIMES'

      logical ifuf,ifcf

      common /tinteg/ tmp(lxyz*lel*ldim),dw(lxyzd,lel,ldim),
     $                ud(lxyzd*lel*lsg*ldim),cd(lxyz*lel*lsg*ldim),
     $                wk(lsg*lsg,2)


      real c(nb,nb,nb)
      real u(lxyz,nel,nb,mdim),v(lxyz,nel,nb,mdim),w(lxyz,nel,nb,ndim)

      ifuf=.false.
      ifcf=.true.

      tt=dnekclock()
      call intp_rstd_all2(ud,u,nel,nb,mdim)
      time_cdea=time_cdea+(dnekclock()-tt)

      do k=1,nb
         tt=dnekclock()
         call set_convect_new_part(dw,dw(1,1,2),dw(1,1,ndim),
     $      w(1,1,k,1),w(1,1,k,2),w(1,1,k,ndim),nel)
         time_cdea=time_cdea+(dnekclock()-tt)
         call rzero(wk,nb*nb)
         do j=1,nb
            do idim=1,mdim
               tt=dnekclock()
               call conv(tmp(1+(idim-1)*lxyz*nel),
     $            ud(1+(j-1)*lxyzd*nel+(idim-1)*lxyzd*nel*nb),.true.,
     $            dw(1,1,1),dw(1,1,2),dw(1,1,ndim),.true.,nel)
               time_cconv=time_cconv+(dnekclock()-tt)
            enddo
            tt=dnekclock()
            do idim=1,mdim
            do i=1,nb
               wk(i+(j-1)*nb,1)=wk(i+(j-1)*nb,1)+
     $            vlsc2(tmp(1+(idim-1)*lxyz*nel),v(1,1,i,idim),lxyz*nel)
            enddo
            enddo
            time_cip=time_cip+(dnekclock()-tt)
         enddo
         call gop(wk,wk(1,2),'+  ',nb*nb)
         tt=dnekclock()
         call updatec(c,wk,nb,k)
         time_cip=time_cip+(dnekclock()-tt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine updatec(c,wk,nb,k)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real c(nb,nb,1),wk(nb,nb)

      integer iloc
      save iloc

c     if (k.eq.1) iloc=1
c     jloc=1

c     do j=1,nb
c     if (mod(j-1,mp).eq.mid) then
c        jloc=(j-1)/mp
c        call add2(c(1,iloc+jloc,1),wk(1,j),nb)
c     endif
c     enddo
c     iloc=iloc+jloc+1

      if (mod(k-1,mp).eq.mid) then
         kloc=(k-1)/mp+1
         call add2(c(1,1,kloc),wk,nb*nb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine aop(au,u,mel)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
      include 'LVAR'
      include 'INTEG'

      common /fastax/ wddx(lx1,lx1),wddyt(ly1,ly1),wddzt(lz1,lz1)
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv

      real           au    (lx1,ly1,lz1,1)
     $ ,             u     (lx1,ly1,lz1,1)
      common /ctmp1/ dudr  (lx1,ly1,lz1)
     $ ,             duds  (lx1,ly1,lz1)
     $ ,             dudt  (lx1,ly1,lz1)
     $ ,             tmp1  (lx1,ly1,lz1)
     $ ,             tmp2  (lx1,ly1,lz1)
     $ ,             tmp3  (lx1,ly1,lz1)

      real           tm1   (lx1,ly1,lz1)
      real           tm2   (lx1,ly1,lz1)
      real           tm3   (lx1,ly1,lz1)
      real           duax  (lx1)
      real           ysm1  (lx1)
      equivalence    (dudr,tm1),(duds,tm2),(dudt,tm3)
      logical ifaxis

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nxy=lx1*ly1
      nyz=ly1*lz1
      nxz=lx1*lz1
      nxyz=lx1*ly1*lz1
      ntot=nxyz*mel

      call rzero (au,ntot)

      do ie=1,mel

        if (ifaxis) call setaxdy ( ifrzer(ie) )

        if (ldim.eq.2) then

c       2-d case ...............

c          if (iffast(ie)) then
           if (.false.) then

c          Fast 2-d mode: constant properties and undeformed element

           h1 = 1.
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           call mxm   (u(1,1,1,ie),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,g4m1(1,1,1,ie),nxyz)
           call col2  (tm2,g5m1(1,1,1,ie),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)

           else


           call mxm  (dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           call mxm  (u(1,1,1,ie),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,g1m1(1,1,1,ie),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,ie),nxyz)
           call addcol3 (tmp1,duds,g4m1(1,1,1,ie),nxyz)
           call addcol3 (tmp2,dudr,g4m1(1,1,1,ie),nxyz)
c          call col2 (tmp1,visc(1,1,1,ie),nxyz)
c          call col2 (tmp2,visc(1,1,1,ie),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,ie),tm1,nxyz)
           call add2 (au(1,1,1,ie),tm2,nxyz)
        endif

        else

c       3-d case ...............

           if (iffast(ie)) then

c          Fast 3-d mode: constant properties and undeformed element

           h1 = 1.
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,ie),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,ie),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,g4m1(1,1,1,ie),nxyz)
           call col2  (tm2,g5m1(1,1,1,ie),nxyz)
           call col2  (tm3,g6m1(1,1,1,ie),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call add2  (au(1,1,1,ie),tm3,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)

           else

           call mxm(dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,ie),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,ie),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,g1m1(1,1,1,ie),nxyz)
           call col3    (tmp2,duds,g2m1(1,1,1,ie),nxyz)
           call col3    (tmp3,dudt,g3m1(1,1,1,ie),nxyz)
           if (ifdfrm(ie)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,ie),nxyz)
              call addcol3 (tmp1,dudt,g5m1(1,1,1,ie),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,ie),nxyz)
              call addcol3 (tmp2,dudt,g6m1(1,1,1,ie),nxyz)
              call addcol3 (tmp3,dudr,g5m1(1,1,1,ie),nxyz)
              call addcol3 (tmp3,duds,g6m1(1,1,1,ie),nxyz)
           endif
c          call col2 (tmp1,visc(1,1,1,ie),nxyz)
c          call col2 (tmp2,visc(1,1,1,ie),nxyz)
c          call col2 (tmp3,visc(1,1,1,ie),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
           call add2 (au(1,1,1,ie),tm1,nxyz)
           call add2 (au(1,1,1,ie),tm2,nxyz)
           call add2 (au(1,1,1,ie),tm3,nxyz)

           endif

        endif

      enddo

c     If axisymmetric, add a diagonal term in the radial direction (ISD=2)

      if (ifaxis.and.(isd.eq.2)) then
         do ie=1,mel

            if (ifrzer(ie)) then
               call mxm(u  (1,1,1,ie),lx1,datm1,ly1,duax,1)
               call mxm(ym1(1,1,1,ie),lx1,datm1,ly1,ysm1,1)
            endif

            do 190 j=1,ly1
            do 190 i=1,lx1
                if (ym1(i,j,1,ie).ne.0.) then
                  if (ifrzer(ie)) then
                     term1 = 0.0
                     if(j.ne.1)
     $             term1 = bm1(i,j,1,ie)*u(i,j,1,ie)/ym1(i,j,1,ie)**2
                     term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i)
     $                       *jacm1(i,1,1,ie)/ysm1(i)
                  else
                   term1 = bm1(i,j,1,ie)*u(i,j,1,ie)/ym1(i,j,1,ie)**2
                     term2 = 0.
                  endif
c                 au(i,j,1,ie) = au(i,j,1,ie)
c    $                          + visc(i,j,1,ie)*(term1+term2)
                  au(i,j,1,ie) = au(i,j,1,ie)
     $                          + (term1+term2)
                endif
  190       continue
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_convect_new_part(cr,cs,ct,ux,uy,uz,mel)!,nb,k)
C
C     Put vxd,vyd,vzd into rst form on fine mesh
C
C     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
C
      include 'LVAR'
      include 'INTEG'

      logical if3d

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real cr(ltd,1),cs(ltd,1),ct(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)
c     real uxyz(lxy,mel,nb,ldim)

      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e
      if3d=ldim.eq.3

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      ic = 1    ! pointer to vector field C

      do e=1,mel

c        Map coarse velocity to fine mesh (C-->F)

         call intp_rstd(fx,ux(1,e),lx1,lxd,if3d,0) ! 0 --> forward
         call intp_rstd(fy,uy(1,e),lx1,lxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(fz,uz(1,e),lx1,lxd,if3d,0) ! 0 --> forward

c        call intp_rstd(fx,uxyz(1,e,k,1),lx1,lxd,if3d,0)
c        call intp_rstd(fy,uxyz(1,e,k,2),lx1,lxd,if3d,0)
c        if (if3d) call intp_rstd(fz,uxyz(1,e,k,ldim),lx1,lxd,if3d,0)

c        Convert convector F to r-s-t coordinates

         if (if3d) then

           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
              cs(i,e)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
              ct(i,e)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
           enddo

         else

           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
              cs(i,e)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
           enddo

         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine conv(bdu,u,ifuf,cx,cy,cz,ifcf,mel)

c     Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
c
      include 'LVAR'
      include 'INTEG'

      real bdu(1),u(1),cx(1),cy(1),cz(1)
      logical ifuf,ifcf,if3d            ! u and/or c already on fine mesh?

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /iconv/ icc
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e
      if3d=ldim.eq.3

      icc=icc+1
      nxyz1=lx1*ly1*lz1
      nxyzd=lxd*lyd*lzd

      nxyzu=nxyz1
      if (ifuf) nxyzu=nxyzd

      nxyzc = nxyz1
      if (ifcf) nxyzc=nxyzd

      iu=1 ! pointer to scalar field u
      ic=1 ! pointer to vector field C
      ib=1 ! pointer to scalar field Bdu

      do e=1,mel
         if (ifcf) then
            call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
            call copy(tr(1,2),cy(ic),nxyzd)
            if (if3d) call copy(tr(1,3),cz(ic),nxyzd)
         else  ! map coarse velocity to fine mesh (C-->F)
            call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0)
            call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0)
            if (if3d) call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0)

            if (if3d) then  ! Convert convector F to r-s-t coordinates
               do i=1,nxyzd
                  tr(i,1)=rx(i,1,e)*fx(i)+
     $                    rx(i,2,e)*fy(i)+
     $                    rx(i,3,e)*fz(i)
                  tr(i,2)=rx(i,4,e)*fx(i)+
     $                    rx(i,5,e)*fy(i)+
     $                    rx(i,6,e)*fz(i)
                  tr(i,3)=rx(i,7,e)*fx(i)+
     $                    rx(i,8,e)*fy(i)+
     $                    rx(i,9,e)*fz(i)
               enddo
            else
               do i=1,nxyzd
                  tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
                  tr(i,2)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
               enddo
           endif
         endif

         if (ifuf) then
            call grad_rst(ur,us,ut,u(iu),lxd,if3d)
         else
            call intp_rstd(uf,u(iu),lx1,lxd,if3d,0)
            call grad_rst(ur,us,ut,uf,lxd,if3d)
         endif

         if (if3d) then
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i)=tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
            enddo
         else
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i)=tr(i,1)*ur(i)+tr(i,2)*us(i)
            enddo
         endif

         call intp_rstd(bdu(ib),uf,lx1,lxd,if3d,1)

         ic=ic+nxyzc
         iu=iu+nxyzu
         ib=ib+nxyz1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setrx(nel)
c
c     Eulerian scheme, add convection term to forcing function
c     at current time step.
c
      include 'LVAR'
      include 'INTEG'

      common /dealias1/ zd(lxd),wd(lxd)
      logical if3d

      nxyz1=lx1*ly1*lz1
      nxyzd=lxd*lyd*lzd

      call zwgl(zd,wd,lxd)  ! zwgl -- NOT zwgll!

      if3d=ldim.eq.3

      if (if3d) then
         do ie=1,nel
            call intp_rstd(rx(1,1,ie),rxm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,2,ie),rym1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,3,ie),rzm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,4,ie),sxm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,5,ie),sym1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,6,ie),szm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,7,ie),txm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,8,ie),tym1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,9,ie),tzm1(1,1,1,ie),lx1,lxd,if3d,0)
            l=0
            do k=1,lzd
            do j=1,lyd
            do i=1,lxd
               l=l+1
               w=wd(i)*wd(j)*wd(k)
               do ii=1,9
                  rx(l,ii,ie)=w*rx(l,ii,ie)
               enddo
            enddo
            enddo
            enddo
         enddo
      else ! 2D
         do ie=1,nel
            call intp_rstd(rx(1,1,ie),rxm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,2,ie),rym1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,3,ie),sxm1(1,1,1,ie),lx1,lxd,if3d,0)
            call intp_rstd(rx(1,4,ie),sym1(1,1,1,ie),lx1,lxd,if3d,0)
            l=0
            do j=1,lyd
            do i=1,lxd
               l=l+1
               w=wd(i)*wd(j)
               do ii=1,4
                  rx(l,ii,ie)=w*rx(l,ii,ie)
               enddo
            enddo
            enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gengrams(ga,gb,gc,gm,gf,gt,buf,tmpf,gvec,ieg,indxr,
     $   nsg,nb,mp,neg,mel,ldim,lx1,lxyz,iftherm,ifbuoy,ifgf,ifm1)

      include 'TIMES'

      common /nekmpi/ mid

      integer ieg(1),indxr(1)
      real ga((nb+1)**2,1),gb((nb+1)**2,1),gm(1)
      real gf(nb**2,lx1-2,1)
      real gc(((nb-1)/mp+2)*(nb+1)**2,1)
      real gt(((nb-1)/mp+2)*(nb+1)**2,1)
      real buf(1),tmpf(1),gvec(nsg,nb)

      logical iftherm,ifdebug,ifgf,ifbuoy,ifm1

      if (mid.eq.0) write (6,*) 'starting gengrams'

      ifdebug=.false.

      ie=1

      call rzero(gb,nb*nb)
      if (iftherm) call rzero(gb(1,2),nb*nb)

      if (ifm1) then
         call rzero(ga,nb*nb)
         call rzero(gc,nb*nb*((nb-1)/mp+1))

         if (iftherm) then
            call rzero(ga(1,2),nb*nb)
            call rzero(gc(1,2),nb*nb*((nb-1)/mp+1))
         endif

         if (ifbuoy) call rzero(gm,nb*nb)
      endif

      tt=dnekclock()

      do while (ie.le.neg)
         ieg(1)=ie
         ieg(2)=min(ie+mel*mp-1,neg)
         rate=(dnekclock()-tt)/(ieg(1)-1)
         if (mid.eq.0) then
            if (ieg(1).eq.1) then
               write (6,1) ieg(1),ieg(2),neg
            else
               write (6,2) ieg(1),ieg(2),neg,rate,(neg-ieg(1)+1)*rate
            endif
         endif
         call loadsnaps(buf,ieg,indxr,nsg,iftherm)
         if (ieg(2).lt.0) then
            neg=ieg(1)+ieg(2)
            exit
         endif

         nel=ieg(4)

         call setgeom(buf,nel,ifm1)
         call transop(buf(nel*lxyz*ldim+1),buf(nel*lxyz*ldim*(nsg+1)+1),
     $      tmpf,gvec,lxyz*nel,ldim,nsg,nb,iftherm)

         write (6,*) 'nel=',nel,mel
         call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),
     $      buf(nel*lxyz*ldim+1),nel,nb,ldim,ldim,ifm1)
         if (ifgf.and.ifm1)
     $      call setgf(gf,buf(nel*lxyz*ldim+1),tmpf,nel,nb,ldim)
         if (iftherm) then
            call setops(ga(1,2),gb(1,2),gc(1,2),gt,
     $         buf(nel*lxyz*ldim+nel*lxyz*ldim*nsg+1),
     $         buf(nel*lxyz*ldim+1),nel,nb,1,ldim,ifm1)
            if (ifgf) call setgf(gf(1,1,2),
     $         buf(nel*lxyz*ldim+nel*lxyz*ldim*nsg+1),tmpf,nel,nb,1)
         endif
         if (ifbuoy.and.ifm1) call setgm(gm,buf,buf(nel*lxyz*ldim+1),
     $      buf(nel*lxyz*ldim+nel*lxyz*ldim*nsg+1),tmpf,nel,nb)
         ie=ieg(2)+1
      enddo

      if (mid.eq.0) write (6,*) 'summing Gramians'

      call gop(gb,gt,'+  ',nb*nb)
      if (iftherm) call gop(gb(1,2),gt,'+  ',nb*nb)

      if (ifm1) then
         call gop(ga,gt,'+  ',nb*nb)
         if (ifgf) then
            do i=1,lx1-2
               call gop(gf(1,i,1),gt,'+  ',nb*nb)
            enddo
         endif

         if (iftherm) then
            call gop(ga(1,2),gt,'+  ',nb*nb)
            if (ifgf) then
               do i=1,lx1-2
                  call gop(gf(1,i,2),gt,'+  ',nb*nb)
               enddo
            endif
         endif

         if (ifbuoy) call gop(gm,gt,'+  ',nb*nb)
      endif

      if (mid.eq.0) write (6,*) 'ending gengrams'

    1 format ('reading [',i8,',',i8,'] / ',i8,' ...')
    2 format ('reading [',i8,',',i8,'] / ',i8,
     $   ', Rate: ',1p1e12.5,', ETA: ',1p1e12.5,' ...')

      return
      end
c-----------------------------------------------------------------------
      subroutine qop1(g1,gt,gvec,nsg,nsg1)

      real g1(1),gt(1),gvec(1)

      call mxm(g1,1,gvec,nsg,gt,nsg1)
      call copy(g1,gt,nsg1)

      return
      end
c-----------------------------------------------------------------------
      subroutine qop2(ga,gt,gvec,gvect,nsg,nsg1)

      real ga(1),gt(1),gvec(1),gvect(1)

      call mxm(gvect,nsg1,ga,nsg,gt,nsg)
      call mxm(gt,nsg1,gvec,nsg,ga,nsg1)

      return
      end
c-----------------------------------------------------------------------
      subroutine qop3(gc,gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid,pfx)

      include 'TIMES'

      real gc(1),gt(1),gvec(1),gvect(1),gvecc(1)
      character*1 pfx
      character*1 fname(2)
      character*2 fname2

      call chcopy(fname,pfx,1)
      fname(2)='c'
      call chcopy(fname2,fname,2)
      call mxm(gvect,nsg1,gc,nsg,gt,nsg*nsc)

      do i=1,nsc
         call mxm(gt(1+(i-1)*(nsg1)*nsg),nsg1,
     $        gvec,nsg,gc(1+(i-1)*(nsg1)**2),nsg1)
      enddo

      if (mid.eq.0) open (unit=10,file=fname2)

      ncount=0
      do k=1,nsg1
         call mxm(gc,(nsg1)**2,gvecc(1+(k-1)*nsc),nsc,gt,1)
         call gop(gt,gt((nsg1)**2+1),'+  ',(nsg1)**2)
         tt=dnekclock()
         do j=0,nsg1-1
         do i=0,nsg1-1
            if (mid.eq.0) write (10,*)
     $         gt(1+i+j*(nsg1))
         enddo
         enddo
         time_dump=time_dump+(dnekclock()-tt)
      enddo

      if (mid.eq.0) close (unit=10)

      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine setg(gram,gub,gram0,nsg,ifavg0)

      real gram(nsg,nsg),gub(nsg,1),gram0(nsg)
      integer ilgls(nsg)
      logical ifavg0

      call copy(gram,gub,nsg*nsg)

      if (ifavg0) then
         call rzero(gram0,nsg)

         do i=1,nsg
            call add2(gram0,gram(1,i),nsg)
         enddo

         s=1./real(nsg)

         call cmult(gram0,s,nsg)
         gram00=s*vlsum(gram0,nsg)

         do i=1,nsg
            call sub2(gram(1,i),gram0,nsg)
         enddo

         call cadd(gram0,-gram00,nsg)

         do i=1,nsg
            s=-gram0(i)
            call cadd(gram(1,i),s,nsg)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setq(qu,qut,quc,elam,gram,nsg,nsc,mp,mid,ifavg0,nsg1)

      logical ifavg0
      real qu(nsg,nsg+1),qut(nsg+1,nsg),quc(1)
      real elam(nsg),gram(nsg,nsg)

      call genevec(qu(1,2),elam,gram,nsg,1)

      do i=1,nsg
         s=1./sqrt(abs(elam(i)))
         call cmult(qu(1,i+1),s,nsg)
      enddo

      if (ifavg0) then
         s=1./nsg
         call cfill(qu,s,nsg)
         do i=1,nsg
            s=-(1./nsg)*vlsum(qu(1,i+1),nsg)
            call cadd(qu(1,i+1),s,nsg)
         enddo
         call rzero(qu(1,nsg+1),nsg)
      else
         call rzero(qu,nsg)
      endif

      nsg1=nsg+1
      if (ifavg0) nsg1=nsg

      call settranspose(qut,qu,nsg,nsg1)

      do i=1,nsg
         if (mod(i-1,mp).eq.mid) nsc=(i-1)/mp+1
      enddo

      do j=1,nsg1
      do i=1,nsc
         quc(i+(j-1)*nsc)=qu(1+(i-1)*mp+mid,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd_all(uf,u,nel)

      include 'LVAR'

      parameter (lxyz1=lxyz)

      real uf(lxyzd,lelt), u(lxyz1,lelt)

      do i=1,nel
         call intp_rstd(uf(1,i),u(1,i),lx1,lxd,ldim.eq.3,0)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd_all2(ud,u,nel,nb,mdim)

      include 'LVAR'

      real ud(lxyzd,nel,nb,mdim),u(lxyz,nel,nb,mdim)

      do idim=1,mdim
      do ib=1,nb
         call intp_rstd_all(ud(1,1,ib,idim),u(1,1,ib,idim),nel)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setgf(gf,buf,tmpf,nel,nsg,ndim)

      include 'LVAR'

      common /sei/ pmat(lx1**2,lx1),ptmat(lx1**2,lx1),
     $             wk(lx1**ldim,5)

      real gf(nsg*nsg,lx1-2),buf(lxyz,nel,nsg,ndim),
     $     tmpf(nel*lxyz*ndim*nsg,2)

      do i=1,lx1-2
         do ie=1,nel*nsg*ndim
            call specmpn(tmpf(1+(ie-1)*lxyz,1),lx1,buf(1,ie,1,1),lx1,
     $         ptmat(1,i),pmat(1,i),ldim.eq.3,wk,5*lx1**ldim)
         enddo
         call copy(tmpf(1,2),tmpf,lxyz*nel*ndim*nsg)
         s=-2.
         call add2s2(tmpf(1,2),buf,s,lxyz*nel*nsg*ndim)
         call bip(gf(1,i),tmpf,tmpf(1,2),nel,nsg,ndim)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setgm(gm,xyz,uvw,t,gxyz,nel,nsg)

      include 'LVAR'

      real gxyz(lxyz*nel*nsg,1)
      real xyz(lxyz,nel,1),uvw(lxyz,nel,1),t(lxyz,nel)

      ! implement gx,gy,gz

      call rzero(gxyz,nel*lxyz*(ldim-1)*nsg)
      call copy(gxyz(1,ldim),t,nel*lxyz*nsg)
      call chsign(gxyz(1,ldim),nel*lxyz*nsg)

c     call bip(gm,uvw,gxyz,nel,nsg,ldim)
      call bip(gm,gxyz,uvw,nel,nsg,ldim)

      return
      end
c-----------------------------------------------------------------------
      subroutine transop(uvw,t,tmpf,gvec,lxyze,ldim,nsg,nb,iftherm)

      real uvw(lxyze,nsg,ldim),t(lxyze,nsg),tmpf(1),gvec(1)
      logical iftherm

      do idim=1,ldim
         call mxm(uvw(1,1,idim),lxyze,gvec,nsg,tmpf,nb)
         call copy(uvw(1,1+(idim-1)*nb,1),tmpf,lxyze*nb)
      enddo

      if (iftherm) then
         call mxm(t,lxyze,gvec,nsg,tmpf,nb)
         call copy(t,tmpf,lxyze*nb)
      endif

      return
      end
c-----------------------------------------------------------------------