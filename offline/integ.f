c-----------------------------------------------------------------------
      subroutine setops(a,b,c,u,nel,nb,ndim)

      real a(1),b(1),c(1),u(1)

      call bip(b,u,u,nel,nb,ndim)

      do j=1,nb
      do i=1,nb
c        write (6,*) i,j,b(i+(j-1)*nb),'b'
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bip(b,u,v,nel,nb,ndim)

      include 'LVAR'
      include 'INTEG'

      common /tinteg/ tmp(lxyz,lel,ldim)

      real b(nb,nb)
      real u(lxyz,nel,nb,ndim),v(lxyz,nel,nb,ndim)

      do l=1,ndim
      do k=1,nb
         call col3(tmp,u(1,1,k,l),bm1,lxyz*nel)
         do j=1,nb
            b(j,k)=b(j,k)+vlsc2(tmp,v(1,1,j,l),lxyz*nel)
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine aip(a,u,v,nel,nb,ndim)

      include 'LVAR'
      include 'INTEG'

      common /tinteg/ tmp(lxyz,lel,ldim)

      real a(nb,nb)
      real u(lxyz,nel,nb,ndim),v(lxyz,nel,nb,ndim)

      do l=1,ndim
      do k=1,nb
         call aop(tmp,u(1,1,k,l),nel)
         do j=1,nb
            a(j,k)=a(j,k)+vlsc2(tmp,v(1,1,j,l),lxyz*nel)
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cip(c,u,v,w,nb,n)

      real c(nb,nb,nb)
      real u(n,nb),v(n,nb),w(n,nb)

      do l=1,nb
      do k=1,nb
      do j=1,nb
      do i=1,n
         c(j,k,l)=c(j,k,l)+u(i,l)*v(i,j)*w(i,k)
      enddo
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
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------