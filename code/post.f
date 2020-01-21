c-----------------------------------------------------------------------
      subroutine drivep

      include 'POST'
      include 'MASS'
      include 'GEOM'
      include 'PARALLEL'

      character*127 flist

      nelp=512
c     nelp=511
c     nelp=11
c     nelp=16
c     nelp=8
c     nelp=4
c     nelp=2
      nelp=1
c     nelp=34
      nsnap=ns

      ! assigned snapshot range

      isg0=1
      isg1=nsg

      call reade_init(fnames,ns)

      inel=1
      ieg1=0
c     ieg1=33

      ns=ls
      nsg=lsg

      call rzero(gram,ns*nsg)
      call rzero(aa,ns*nsg)
      call rzero(bb,ns*nsg)
      call rzero(cc,ns*nsg*nsg)
      call rzero(cc2,ns*nsg*nsg)

      ng=(ndim-1)*3

      do while (ieg1+1.le.nelgv)
         ieg0=ieg1+1
         ieg1=min(ieg1+inel+nelp-1,nelgv)
         nel=ieg1-ieg0+1
         n=lxyz*(ieg1-ieg0+1)
c        write (6,*) ieg0,ieg1,nel,n,'info'
c        call reade_dummy(uu,ieg0,ieg1)
c        call reade2(uu,ieg0,ieg1)
         call reade(uu,ieg0,ieg1)

c        call setgeom(gfac,w9,ieg0,ieg1,lxyz,ng,nid)
c        call setvisc(visc,w,ieg0,ieg1,lxyz,nid)
c        call setmass(mass,wv1,ieg0,ieg1,lxyz)
c        call dump_serial(rx,lxd*lyd*lzd*4*nelt,'ops/rx ',nid)
c        call setrxp(rxp,rxpt,ieg0,ieg1)
c        if (ieg0.eq.1) then
c           call dump_serial(rxp,lxd*lyd*lzd*4*nelp,'ops/rxp1 ',nid)
c        else
c           call dump_serial(rxp,lxd*lyd*lzd*4*nelp,'ops/rxp2 ',nid)
c        endif

c        call setbb(bb,uu,mass,wvf1,wvf2,wvf12,ns,nsg,n,ndim)
c        call setaa(aa,uu,visc,gfac,wvf1,wvf2,wvf12,
c    $      ns,nsg,n,nel,ndim,ng)
c        call setcc(cc,uu,uu,rxp,wvf1,wvf2,wvf12,wvf12,
c    $      ns,nsg,n,ndim,ndim,nel)
c        call setcc_snap(cc2)
      enddo

      call dump_serial(bb,ns*ns,'ops/graml2 ',nid)
      call dump_serial(aa,ns*ns,'ops/gramh10 ',nid)
      call dump_serial(cc,ns*ns*ns,'ops/gramc ',nid)
      call dump_serial(cc2,ns*ns*ns,'ops/gramc2 ',nid)

      ! eigendecomposition here or external process

      mmm=ns*ns
      call read_serial(evecp,mmm,'ops/evecp ',ug,nid)
      call read_serial(evecpt,mmm,'ops/evecpt ',ug,nid)

      call mxm(bb,ns,evecp,ns,wevec,ns)
      call mxm(evecpt,ns,wevec,ns,bb,ns)
      call dump_serial(bb,ns*ns,'ops/bup_new ',nid)

      call mxm(aa,ns,evecp,ns,wevec,ns)
      call mxm(evecpt,ns,wevec,ns,aa,ns)
      call dump_serial(aa,ns*ns,'ops/aup_new ',nid)

      call mxm(cc,ns*ns,evecp,ns,wevecc,ns)
      do i=1,ns
         call mxm(wevecc(2+(i-1)*ns*ns),ns,evecp,ns,
     $      cc(1+(i-1)*ns*ns),ns)
      enddo
      call mxm(evecpt,ns,cc,ns,wevecc,ns*ns)
      call copy(cc,wevecc,ns*ns*ns)
      call dump_serial(cc,ns*ns*ns,'ops/cup_new ',nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine reade(v,ieg0,ieg1)

      include 'POST'
      include 'MASS'
      include 'GEOM'

      common /screrr/ err(lxyz,lelt,ldim)

      real v(lxyz,ieg1-ieg0+1,ldim,ns)

      n=lxyz*(ieg1-ieg0+1)

      call mfip_setup
      is=1
      call mfip_init(fnames(1+(is-1)*132))

      do is=1,1
         call mfip_read(v(1,1,1,is),v(1,1,2,is),v(1,1,ldim,is),
     $      ieg0,ieg1)
      enddo

      do ieg=ieg0,ieg1
         ie=ieg-ieg0+1
         u1l2=vlsc2(v(1,ie,1,1),v(1,ie,1,1),lxyz)
         u2l2=vlsc2(us0(1+(ieg-1)*lxyz,1,1),
     $              us0(1+(ieg-1)*lxyz,1,1),lxyz)
         call sub3(err,v(1,ie,1,1),us0(1+(ieg-1)*lxyz,1,1),lxyz)
         uel2=vlsc2(err,err,lxyz)

         v1l2=vlsc2(v(1,ie,2,1),v(1,ie,2,1),lxyz)
         v2l2=vlsc2(us0(1+(ieg-1)*lxyz,2,1),
     $              us0(1+(ieg-1)*lxyz,2,1),lxyz)
         call sub3(err,v(1,ie,2,1),us0(1+(ieg-1)*lxyz,2,1),lxyz)
         vel2=vlsc2(err,err,lxyz)

         write (6,*) ieg,u1l2,u2l2,uel2,'uerr'
         write (6,*) ieg,v1l2,v2l2,vel2,'verr'
      enddo

      call mfip_end

      return
      end
c-----------------------------------------------------------------------
      subroutine reade2(v,ieg0,ieg1)

      include 'POST'
      include 'MASS'
      include 'GEOM'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /screrr/ err(lxyz,lelt,ldim)

      real v(lxyz,ieg1-ieg0+1,ldim,ns)

      n=lxyz*(ieg1-ieg0+1)

c     call mfip_setup
      is=1
      call cread_head(fnames(1+(is-1)*132))
      call byte_open(fnames(1+(is-1)*132),ierr)
c     call byte_open('cyl0.f01000',ierr)

      do is=1,1
         call cread(v(1,1,1,is),v(1,1,2,is),v(1,1,ldim,is),
     $      ieg0,ieg1)
      enddo
c     call exitt0

      do ieg=ieg0,ieg1
         ie=ieg-ieg0+1
         u1l2=vlsc2(v(1,ie,1,1),v(1,ie,1,1),lxyz)
         u2l2=vlsc2(us0(1+(ieg-1)*lxyz,1,1),
     $              us0(1+(ieg-1)*lxyz,1,1),lxyz)
         call sub3(err,v(1,ie,1,1),us0(1+(ieg-1)*lxyz,1,1),lxyz)
         uel2=vlsc2(err,err,lxyz)

         v1l2=vlsc2(v(1,ie,2,1),v(1,ie,2,1),lxyz)
         v2l2=vlsc2(us0(1+(ieg-1)*lxyz,2,1),
     $              us0(1+(ieg-1)*lxyz,2,1),lxyz)
         call sub3(err,v(1,ie,2,1),us0(1+(ieg-1)*lxyz,2,1),lxyz)
         vel2=vlsc2(err,err,lxyz)

         write (6,*) ieg,u1l2,u2l2,uel2,'uerr'
         write (6,*) ieg,v1l2,v2l2,vel2,'verr'
      enddo

      call byte_close(ierr)

c     call cread_end

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

c     include 'SIZE'
c     include 'PARALLEL'

      real z(n,ns),u(n,ns)
      real evec(ns,nsg)
      real w1(n),w2(n)

c     do isg=1,nsg
c        is=iglls(isg)
c        is=isg
c        call mxm(u,n,evec(1,isg),ns,w1,1)
c        call gop(w1,w2,'+  ',n)
c        if (is.gt.0.and.is.le.ns) call copy(z(1,is),w1,n)

c     enddo

      do i=1,ns
         call mxm(u,n,evec(1,i),ns,z(1,i),1)
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
      subroutine setaa(a,z,visc,gfac,w1,w2,w3,ns,nsg,n,nel,ndim,ng)

      real a(ns,nsg),z(n,ndim,ns),w1(n,ndim,ns),w2(n,ndim,ns)
      real w3(n,ndim,ns,2)
      real visc(n),gfac(n,ng)

      call copy(w1,z,n*ndim*ns)
      call copy(w2,z,n*ndim*ns)

      imesh=1
      isd=1

      call rzero(w3,n)

      do i=1,ns*ndim
         call aop(w2(1,i,1),z(1,i,1),visc,gfac,imesh,isd,nel)
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            a(is,js+ns*ioff)=a(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,w3,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setbb(b,u,mass,w1,w2,w3,ns,nsg,n,ndim)

      real b(ns,nsg),u(n,ndim,ns),mass(n),w1(n,ndim,ns),w2(n,ndim,ns)
      real w3(n,ndim,ns,3)

      write (6,*) 'ns,nsg,n,ndim',ns,nsg,n,ndim

      call copy(w1,u,n*ndim*ns)
      call copy(w2,u,n*ndim*ns)

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
         
         call shift(w1,w3,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc_snap(cc)

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrcc/ t1(lt),t2(lt),t3(lt)

      real cc(ns,ns,ns)

      n=lx1*ly1*lz1*nelv

      call rone(ones,lx1*ly1*lz1*nelv)
      ifield=1

      do ks=1,ns
         do js=1,ns
            call convect_new(t1,us0(1,1,js),.false.,
     $         us0(1,1,ks),us0(1,2,ks),us0(1,mdim,ks),.false.)
            call convect_new(t2,us0(1,2,js),.false.,
     $         us0(1,1,ks),us0(1,2,ks),us0(1,mdim,ks),.false.)
            do is=1,ns
               cc(is,js,ks)=cc(is,js,ks)
     $            +vlsc2(t1,us0(1,1,is),n)+vlsc2(t2,us0(1,2,is),n)
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc(c,z,t,rxp,w1,w2,w3,w4,ns,nsg,n,mdim,ndim,nel)

      real rxp(1)
      real c(ns,nsg,nsg),z(n,ndim,ns),t(n,mdim,ns),
     $     w1(n,mdim,ns),w2(n,ndim,ns),w3(n,mdim,ns),w4(n,mdim,ns)

      do ks=1,ns
         do js=1,ns
c           call convect_new(w3,t(1,1,js),.false.,
c    $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.)
c           call convect_new(w3(1,2,1),t(1,2,js),.false.,
c    $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.)
            call conv(w3,t(1,1,js),.false.,
     $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.,rxp,nel)
            call conv(w3(1,2,1),t(1,2,js),.false.,
     $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.,rxp,nel)
            do is=1,ns
               c(is,js,ks)=c(is,js,ks)
     $            +vlsc2(w3,t(1,1,is),n*mdim)
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setgeom(gfac,w,ieg0,ieg1,lxyz,ng)!,nid)

      include 'SIZE' ! nid
      include 'GEOM'
      include 'PARALLEL'

      real gfac(lxyz,ieg1-ieg0+1,ng),w(lxyz*(ieg1-ieg0+1)*ng)

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

      call rone(visc,lxyz*(ieg1-ieg0+1))
      return

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
      subroutine setgg(gram,u,mass,w1,w2,w3,ns,nsg,n,ndim)

      real gram(ns,nsg),u(n,ndim,ns),mass(n),w1(n,ndim,ns),w2(n,ndim,ns)
      real w3(n,ndim,ns,3)

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
         
         call shift(w1,w3,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine shift_setup(iw,n)

      include 'SIZE'

      common /pcomm/ igsh(2)
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer*8 iw(n,2)

      if (mp.eq.1) return

      if (mod(mp,2).eq.1)
     $   call exitti('no support for odd number of mpi-ranks$',mp)

      call izero(igsh,2)
      do i=1,n
         i1=nid/2
         i2=mod(nid-1+mp,mp)/2
         irem=mod(nid,2).eq.0
         iw(i,1)=i+i1*n
         iw(i,2)=i+i2*n
      enddo

      do i=1,mp
         if (nid.eq.(i-1)) then
            do j=1,n
               write (6,*) nid,i,j,iw(j,1),iw(j,2)
            enddo
         endif
         call nekgsync
      enddo

      call fgslib_gs_setup(igsh,iw,n,nekcomm,mp)
      call fgslib_gs_setup(igsh(2),iw(1,2),n,nekcomm,mp)

      return
      end
c-----------------------------------------------------------------------
      subroutine shift(u,w,n)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /pcomm/ igsh(2)

      real u(n),w(n,2)

      write (6,*) 'mp=',mp

      if (mp.eq.1) return

      if (mod(mp,2).eq.1)
     $   call exitti('no support for odd number of mpi-ranks$',mp)

      if (mod(mid,2).eq.0) then
         call copy(w,u,n)
         call rzero(w(1,2),n)
      else
         call copy(w(1,2),u,n)
         call rzero(w,n)
      endif

      call fgslib_gs_op(igsh,w,1,1,0)
      call fgslib_gs_op(igsh(2),w(1,2),1,1,0)

      if (mod(mid,2).eq.0) then
         call copy(u,w(1,2),n)
      else
         call copy(u,w,n)
      endif

      return
      end
c-----------------------------------------------------------------------
c     subroutine ax(au,u,helm1,helm2,imesh,isd)
      subroutine aop(au,u,visc,gfac,imesh,isd,mel)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
c     include 'PARALLEL'
c     include 'CTIMER'

      common /fastax/ wddx(lx1,lx1),wddyt(ly1,ly1),wddzt(lz1,lz1)
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv

      real           au    (lx1,ly1,lz1,1)
     $ ,             u     (lx1,ly1,lz1,1)
     $ ,             visc (lx1,ly1,lz1,1)
     $ ,             helm2 (lx1,ly1,lz1,1)
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

      real gfac(lx1,ly1,lz1,mel,ng)

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nxy=lx1*ly1
      nyz=ly1*lz1
      nxz=lx1*lz1
      nxyz=lx1*ly1*lz1
      ntot=nxyz*mel

c     if (.not.ifsolv) call setfast(visc,helm2,imesh)
      call rzero (au,ntot)

      do ie=1,mel
 
        if (ifaxis) call setaxdy ( ifrzer(ie) )

        if (ldim.eq.2) then

c       2-d case ...............

c          if (iffast(ie)) then
           if (.false.) then

c          Fast 2-d mode: constant properties and undeformed element

           h1 = visc(1,1,1,ie)
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           call mxm   (u(1,1,1,ie),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,gfac(1,1,1,ie,4),nxyz)
           call col2  (tm2,gfac(1,1,1,ie,5),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)

           else


           call mxm  (dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           call mxm  (u(1,1,1,ie),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,gfac(1,1,1,ie,1),nxyz)
           call col3 (tmp2,duds,gfac(1,1,1,ie,2),nxyz)
c          if (ifdfrm(ie)) then
              call addcol3 (tmp1,duds,gfac(1,1,1,ie,4),nxyz)
              call addcol3 (tmp2,dudr,gfac(1,1,1,ie,4),nxyz)
c          endif
           call col2 (tmp1,visc(1,1,1,ie),nxyz)
           call col2 (tmp2,visc(1,1,1,ie),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,ie),tm1,nxyz)
           call add2 (au(1,1,1,ie),tm2,nxyz)
        endif

        else

c       3-d case ...............

           if (iffast(ie)) then

c          Fast 3-d mode: constant properties and undeformed element

           h1 = visc(1,1,1,ie)
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,ie),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,ie),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,gfac(1,1,1,ie,4),nxyz)
           call col2  (tm2,gfac(1,1,1,ie,5),nxyz)
           call col2  (tm3,gfac(1,1,1,ie,6),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call add2  (au(1,1,1,ie),tm3,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)

           else

           call mxm(dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,ie),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,ie),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,gfac(1,1,1,ie,1),nxyz)
           call col3    (tmp2,duds,gfac(1,1,1,ie,2),nxyz)
           call col3    (tmp3,dudt,gfac(1,1,1,ie,3),nxyz)
           if (ifdfrm(ie)) then
              call addcol3 (tmp1,duds,gfac(1,1,1,ie,4),nxyz)
              call addcol3 (tmp1,dudt,gfac(1,1,1,ie,5),nxyz)
              call addcol3 (tmp2,dudr,gfac(1,1,1,ie,4),nxyz)
              call addcol3 (tmp2,dudt,gfac(1,1,1,ie,6),nxyz)
              call addcol3 (tmp3,dudr,gfac(1,1,1,ie,5),nxyz)
              call addcol3 (tmp3,duds,gfac(1,1,1,ie,6),nxyz)
           endif
           call col2 (tmp1,visc(1,1,1,ie),nxyz)
           call col2 (tmp2,visc(1,1,1,ie),nxyz)
           call col2 (tmp3,visc(1,1,1,ie),nxyz)
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
                  au(i,j,1,ie) = au(i,j,1,ie)
     $                          + visc(i,j,1,ie)*(term1+term2)
                endif
  190       continue
         enddo
      endif

      return
      end
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
      subroutine conv(bdu,u,ifuf,cx,cy,cz,ifcf,rxp,mel)

c     Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
c
      include 'SIZE'
      include 'TOTAL'

      real bdu(1),u(1),cx(1),cy(1),cz(1)
      real rxp(lxd*lyd*lzd,ldim*ldim,mel)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e

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
                  tr(i,1)=rxp(i,1,e)*fx(i)+
     $                    rxp(i,2,e)*fy(i)+
     $                    rxp(i,3,e)*fz(i)
                  tr(i,2)=rxp(i,4,e)*fx(i)+
     $                    rxp(i,5,e)*fy(i)+
     $                    rxp(i,6,e)*fz(i)
                  tr(i,3)=rxp(i,7,e)*fx(i)+
     $                    rxp(i,8,e)*fy(i)+
     $                    rxp(i,9,e)*fz(i)
               enddo
            else
               do i=1,nxyzd
                  tr(i,1)=rxp(i,1,e)*fx(i)+rxp(i,2,e)*fy(i)
                  tr(i,2)=rxp(i,3,e)*fx(i)+rxp(i,4,e)*fy(i)
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
      subroutine setrxp(rxp,tmp,ieg0,ieg1)
c
c     Eulerian scheme, add convection term to forcing function
c     at current time step.
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP' ! for istep
      include 'PARALLEL'

      common /dealias1/ zd(lxd),wd(lxd)
      real rxp(lxd*lyd*lzd,ldim*ldim,ieg1-ieg0+1)
      real tmp(lxd*lyd*lzd,ldim*ldim,ieg1-ieg0+1)
      integer e

      integer ilstep
      save    ilstep
      data    ilstep /-1/

      write (6,*) ieg0,ieg1,'rxp'

c     if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
c     if (ifgeom.and.ilstep.eq.istep)  return  ! already computed
      ilstep=istep

      nxyz1=lx1*ly1*lz1
      nxyzd=lxd*lyd*lzd

      call zwgl(zd,wd,lxd)  ! zwgl -- NOT zwgll!
      call rzero(rxp,lxd*lyd*lzd*ldim*ldim*(ieg1-ieg0+1))

      if (if3d) then
         do ieg=ieg0,ieg1
            ie=ieg-ieg0+1
            if (gllnid(ieg).eq.nid) then
               e=gllel(ieg)

c              Interpolate z+ and z- into fine mesh, translate to r-s-t coords

               call intp_rstd(rxp(1,1,ie),rxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,2,ie),rym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,3,ie),rzm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,4,ie),sxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,5,ie),sym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,6,ie),szm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,7,ie),txm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,8,ie),tym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,9,ie),tzm1(1,1,1,e),lx1,lxd,if3d,0)

               l=0
               do k=1,lzd
               do j=1,lyd
               do i=1,lxd
                  l=l+1
                  w=wd(i)*wd(j)*wd(k)
                  do ii=1,9
                     rxp(l,ii,ie)=w*rxp(l,ii,ie)
                  enddo
               enddo
               enddo
               enddo
            endif
         enddo

      else ! 2D
         do ieg=ieg0,ieg1
            ie=ieg-ieg0+1
            if (gllnid(ieg).eq.nid) then
               e=gllel(ieg)

c              Interpolate z+ and z- into fine mesh, translate to r-s-t coords

               call intp_rstd(rxp(1,1,ie),rxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,2,ie),rym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,3,ie),sxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,4,ie),sym1(1,1,1,e),lx1,lxd,if3d,0)

               l=0
               do j=1,lyd
               do i=1,lxd
                  l=l+1
                  w=wd(i)*wd(j)
                  do ii=1,4
                     rxp(l,ii,ie)=w*rxp(l,ii,ie)
                  enddo
               enddo
               enddo
            endif
         enddo
      endif

      call gop(rxp,tmp,'+  ',lxd*lyd*lzd*ldim*ldim*(ieg1-ieg0+1))

      return
      end
c-----------------------------------------------------------------------
      subroutine mfip_setup
c----------------------------------------------------------------------
c
c     (1) Open restart file(s)
c     (2) Check previous spatial discretization 
c     (3) Map (K1,N1) => (K2,N2) if necessary
c
c     nfiles > 1 has several implications:
c
c     i.   For std. run, data is taken from last file in list, unless
c          explicitly specified in argument list of filename
c
c     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
C          subsequent files are for B-field or perturbation fields
c
c
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'RESTART'

      common /inelr/ nelrr

      parameter (lxr=lx1+6)
      parameter (lyr=ly1+6)
      parameter (lzr=lz1+6)
      parameter (lxyzr=lxr*lyr*lzr)
      parameter (lxyzt=lx1*ly1*lz1*lelt)
      parameter (lpsc9=ldimt+9)

      common /scrcg/ pm1(lx1*ly1*lz1,lelv)
      COMMON /SCRNS/ SDUMP(LXYZT,7)
      integer mesg(40)

c     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 9.
      COMMON /CTMP1/ TDUMP(LXYZR,LPSC9)
      real*4         tdump

      REAL SDMP2(LXYZT,LDIMT)

c     cdump comes in via PARALLEL (->TOTAL)

      character*30 excoder
      character*1  excoder1(30)
      equivalence (excoder,excoder1)

      character*132 fname
      character*1  fname1(132)
      equivalence (fname1,fname)

      integer       hnami (30)
      character*132 hname
      character*1   hname1(132)
      equivalence  (hname,hname1)
      equivalence  (hname,hnami )

      character*132 header

c     Local logical flags to determine whether to copy data or not.
      logical ifok,iffmat
      integer iposx,iposz,iposu,iposw,iposp,ipost,ipsps(ldimt1)

      logical ifbytsw, if_byte_swap_test
      real*4   bytetest

      ifok=.false.
      ifbytsw = .false.

      nfiles=1

      if (nio.eq.0) write(6,*) 'Reading checkpoint data '

c use new reader (only binary support)
      p67 = abs(param(67))
      if (p67.eq.6.0) then
         ifile=1
         call sioflag(ndumps,fname,initc(ifile))
      else
         call exitti('non-binary reading not supported$',nint(p67))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mfip_init(fname_in)
c
c     (1) Open restart file(s)
c     (2) Check previous spatial discretization 
c     (3) Map (K1,N1) => (K2,N2) if necessary
c
c     nfiles > 1 has several implications:
c
c     i.   For std. run, data is taken from last file in list, unless
c          explicitly specified in argument list of filename
c
c     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
c          subsequent files are for B-field or perturbation fields
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      character*132  fname_in

      character*132  fname
      character*1    fnam1(132)
      equivalence   (fnam1,fname)

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)

      ! add path
      call blank(fname,132)
      lenp = ltrunc(path,132)
      lenf = ltrunc(fname_in,132)
      call chcopy(fnam1(1),path,lenp)
      call chcopy(fnam1(lenp+1),fname_in,lenf)

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping 

      return
      end
c-----------------------------------------------------------------------
      subroutine cread_head(fname_in)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      character*132  fname_in

      character*132  fname
      character*1    fnam1(132)
      equivalence   (fnam1,fname)

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)

      ! add path
      call blank(fname,132)
      lenp = ltrunc(path,132)
      lenf = ltrunc(fname_in,132)
      call chcopy(fnam1(1),path,lenp)
      call chcopy(fnam1(lenp+1),fname_in,lenf)

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping 

      return
      end
c-----------------------------------------------------------------------
      subroutine mfip_read(ux,uy,uz,ieg0,ieg1)
c
c     (1) Open restart file(s)
c     (2) Check previous spatial discretization 
c     (3) Map (K1,N1) => (K2,N2) if necessary
c
c     nfiles > 1 has several implications:
c
c     i.   For std. run, data is taken from last file in list, unless
c          explicitly specified in argument list of filename
c
c     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
c          subsequent files are for B-field or perturbation fields
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)
      common /scrread/ i1(lelt),i2(lelt),i3(lelt),i4(lelt)

      real ux(lx1,ly1,lz1,ieg1-ieg0+1)
      real uy(lx1,ly1,lz1,ieg1-ieg0+1)
      real uz(lx1,ly1,lz1,ieg1-ieg0+1)

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      call rzero(wk,7*lx1*ly1*lz1*nelt*ldim)

      iofldsr=0
      if (ifgetxr) iofldsr=ldim

      nelr=ieg1-ieg0+1

      icount=1
      ie=1

      do while (icount.le.nelr)
         if (er(ie).ge.ieg0.and.er(ie).le.ieg1) then
            i1(icount)=ie
            icount=icount+1
         endif
         ie=ie+1
      enddo

      call isort(i1,i2,nelr)

      ngroup=1

      ic=1
      ng=1
      i3(1)=i1(1)

      call ione(i4,nelr)

      do ie=2,nelr
         if (i1(ie)-i1(ie-1).gt.1) then ! if jump
            ng=ng+1
            i3(ng)=i1(ie)
         else
            i4(ng)=i4(ng)+1
         endif
      enddo

      neltmp=nelr

      iloc=1

      do ig=1,ng
         ieg=i3(ig)
         nelr=i4(ig)
         nwk=nelr*ldim*nxyzr8
         offs = offs0 + iofldsr*stride + ldim*strideB + 
     $      ldim*(ieg-1)*nxyzr8*wdsizr

         write (6,*) ieg,nelr,nwk,'info'
         call byte_set_view(offs,ifh_mbyte)
         call mfi_getw(wk(iloc),nwk,.false.)
         iloc=iloc+nwk
      enddo

      nelr=neltmp

      lxyz=lx1*ly1*lz1
      do ie=1,nelr
         write (6,*) ieg0+(ie-1),i1(ie),i2(ie),'info1'
         call copy(ux(1,1,1,i2(ie)),wk(1+(ie-1)*lxyz*ldim),lxyz)
         call copy(uy(1,1,1,i2(ie)),wk(1+(ie-1)*lxyz*ldim+lxyz),lxyz)
         if (ldim.eq.3) call
     $      copy(uz(1,1,1,i2(ie)),wk(1+(ie-1)*lxyz*ldim+2*lxyz),lxyz)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cread(ux,uy,uz,ieg0,ieg1)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      parameter (lxyz=lx1*ly1*lz1)

      common /scrns/ wk(lwk)

      real ux(lxyz*lelt,ieg1-ieg0+1)
      real uy(lxyz*lelt,ieg1-ieg0+1)
      real uz(lxyz*lelt,ieg1-ieg0+1)

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      iofldsr=0
      if (ifgetxr) iofldsr=ldim

      nelr=ieg1-ieg0+1
      offs = offs0 + iofldsr*stride + ldim*strideB + 
     $   ldim*(ieg0-1)*nxyzr8*wdsizr

      call byte_seek(offs/4,ierr)
      call byte_read(wk,nelr*ldim*nxyzr8*wdsizr/4,ierr)

      do ie=1,nelr
         call copy(ux(1,ie),wk(1+(ie-1)*lxyz*ldim),lxyz)
         call copy(uy(1,ie),wk(1+(ie-1)*lxyz*ldim+lxyz),lxyz)
         if (ldim.eq.3)
     $      call copy(uz(1,ie),wk(1+(ie-1)*lxyz+2*lxyz),lxyz)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mfip_end

      include 'SIZE'
      include 'INPUT'
      include 'RESTART'

      if (ifmpiio) then
         if (nid.eq.pid0r) call byte_close_mpi(ifh_mbyte,ierr)
      else
         if (nid.eq.pid0r) call byte_close(ierr)
      endif

      call err_chk(ierr,'Error closing restart file, in mfi.$')

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_getv_p(u,v,w,wk,lwk,iskip)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,1),v(lx1*ly1*lz1,1),w(lx1*ly1*lz1,1)
      logical iskip

      real*4 wk(lwk) ! message buffer
      parameter(lrbs=50*lx1*ly1*lz1*lelt)
      common /vrthov/ w2(lrbs) ! read buffer
      real*4 w2

      integer e,ei,eg,msg_id(lelt)
      integer*8 i8tmp

      call nekgsync() ! clear outstanding message queues.

      nxyzr  = ldim*nxr*nyr*nzr
c     dnxyzr = nxyzr
c     len    = nxyzr*wdsizr             ! message length in bytes
      if (wdsizr.eq.8) nxyzr = 2*nxyzr

      ! check message buffer
c     num_recv  = len
c     num_avail = lwk*wdsize
c     call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv a')

c     ! setup read buffer
c     i8tmp = int(nxyzr,8)*int(nelr,8)
c     nread = i8tmp/int(lrbs,8)
c     if (mod(i8tmp,int(lrbs,8)).ne.0) nread = nread + 1
c     if(ifmpiio) nread = iglmax(nread,1) ! needed because of collective read
c     nelrr = nelr/nread
c     call bcast(nelrr,4)
c     call lim_chk(nxyzr*nelrr,lrbs,'     ','     ','mfi_getv b')

      ierr = 0

      call byte_read_mpi(wk,nxyzr*nelr,-1,ifh_mbyte,ierr)

      nxyzr = nxr*nyr*nzr
      nxyzw = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzw = 2*nxyzw

      l = 1
      do e=1,nelr
c        ei = er(e) 
         call copy  (u(1,e),wk(l        ),nxyzr)
         call copy  (v(1,e),wk(l+  nxyzw),nxyzr)
         l = l+ldim*nxyzw
      enddo

 100  call err_chk(ierr,'Error reading restart data, in getv.$')
      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_getw(wk,lwk,iskip)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      logical iskip

      real*4 wk(lwk) ! message buffer
      parameter(lrbs=50*lx1*ly1*lz1*lelt)
      common /vrthov/ w2(lrbs) ! read buffer
      real*4 w2

      integer e,ei,eg,msg_id(lelt)
      integer*8 i8tmp

      call nekgsync() ! clear outstanding message queues.

      nxyzr  = ldim*nxr*nyr*nzr
c     dnxyzr = nxyzr
c     len    = nxyzr*wdsizr             ! message length in bytes
      if (wdsizr.eq.8) nxyzr = 2*nxyzr

      ! check message buffer
c     num_recv  = len
c     num_avail = lwk*wdsize
c     call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv a')

c     ! setup read buffer
c     i8tmp = int(nxyzr,8)*int(nelr,8)
c     nread = i8tmp/int(lrbs,8)
c     if (mod(i8tmp,int(lrbs,8)).ne.0) nread = nread + 1
c     if(ifmpiio) nread = iglmax(nread,1) ! needed because of collective read
c     nelrr = nelr/nread
c     call bcast(nelrr,4)
c     call lim_chk(nxyzr*nelrr,lrbs,'     ','     ','mfi_getv b')

      ierr = 0

      call byte_read_mpi(wk,nxyzr*nelr,-1,ifh_mbyte,ierr)

 100  call err_chk(ierr,'Error reading restart data, in getv.$')
      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_getv_p2(u,v,w,wk,lwk,iskip)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,1),v(lx1*ly1*lz1,1),w(lx1*ly1*lz1,1)
      logical iskip

      real*4 wk(lwk) ! message buffer
      parameter(lrbs=20*lx1*ly1*lz1*lelt)
      common /vrthov/ w2(lrbs) ! read buffer
      real*4 w2

      integer e,ei,eg,msg_id(lelt)
      integer*8 i8tmp

      call nekgsync() ! clear outstanding message queues.

      nxyzr  = ldim*nxr*nyr*nzr
      dnxyzr = nxyzr
      len    = nxyzr*wdsizr             ! message length in bytes
      if (wdsizr.eq.8) nxyzr = 2*nxyzr

      ! check message buffer
      num_recv  = len
      num_avail = lwk*wdsize
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv a')

      ! setup read buffer
      if(nid.eq.pid0r) then
c         dtmp  = dnxyzr*nelr
c         nread = dtmp/lrbs
c         if(mod(dtmp,1.0*lrbs).ne.0) nread = nread + 1
         i8tmp = int(nxyzr,8)*int(nelr,8)
         nread = i8tmp/int(lrbs,8)
         if (mod(i8tmp,int(lrbs,8)).ne.0) nread = nread + 1
         if(ifmpiio) nread = iglmax(nread,1) ! needed because of collective read
         nelrr = nelr/nread
      endif
      call bcast(nelrr,4)
      call lim_chk(nxyzr*nelrr,lrbs,'     ','     ','mfi_getv b')

      ! pre-post recieves (one mesg per element)
      ! this assumes we never pre post more messages than supported
      if (np.gt.1) then
         l = 1
         do e=1,nelt
            msg_id(e) = irecv(e,wk(l),len)
            l = l+nxyzr
         enddo
      endif

      ierr = 0
      if (nid.eq.pid0r .and. np.gt.1) then ! only i/o nodes
         k = 0
         do i = 1,nread
            if(i.eq.nread) then ! clean-up 
              nelrr = nelr - (nread-1)*nelrr
              if(nelrr.lt.0) nelrr = 0
            endif

            if(ierr.eq.0) then
              if(ifmpiio) then 
                call byte_read_mpi(w2,nxyzr*nelrr,-1,ifh_mbyte,ierr)
              else
                call byte_read (w2,nxyzr*nelrr,ierr)
              endif
            endif

            ! redistribute data based on the current el-proc map
            l = 1
            do e = k+1,k+nelrr
               jnid = gllnid(er(e))                ! where is er(e) now?
               jeln = gllel(er(e))
               if(ierr.ne.0) call rzero(w2(l),len)
               call csend(jeln,w2(l),len,jnid,0)  ! blocking send
               l = l+nxyzr
            enddo
            k  = k + nelrr
         enddo
      elseif (np.eq.1) then
         if(ifmpiio) then 
           call byte_read_mpi(wk,nxyzr*nelr,-1,ifh_mbyte,ierr)
         else
           call byte_read(wk,nxyzr*nelr,ierr)
         endif
      endif

      if (iskip) then
         call nekgsync() ! clear outstanding message queues.
         goto 100     ! don't assign the data we just read
      endif

      nxyzr = nxr*nyr*nzr
      nxyzv = ldim*nxr*nyr*nzr
      nxyzw = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzw = 2*nxyzw

      l = 1
      do e=1,nelr
         if (np.gt.1) then
            call msgwait(msg_id(e))
            ei = e
         else if(np.eq.1) then
            ei = er(e) 
         endif
         if (if_byte_sw) then
            if(wdsizr.eq.8) then
               call byte_reverse8(wk(l),nxyzv*2,ierr)
            else
               call byte_reverse(wk(l),nxyzv,ierr)
            endif
         endif
         if (nxr.eq.lx1.and.nyr.eq.ly1.and.nzr.eq.lz1) then
            if (wdsizr.eq.4) then         ! COPY
               call copy4r(u(1,ei),wk(l        ),nxyzr)
               call copy4r(v(1,ei),wk(l+  nxyzw),nxyzr)
               if (if3d) 
     $         call copy4r(w(1,ei),wk(l+2*nxyzw),nxyzr)
            else
               call copy  (u(1,ei),wk(l        ),nxyzr)
               call copy  (v(1,ei),wk(l+  nxyzw),nxyzr)
               if (if3d) 
     $         call copy  (w(1,ei),wk(l+2*nxyzw),nxyzr)
            endif
         else                             ! INTERPOLATE
            if (wdsizr.eq.4) then
               call mapab4r(u(1,ei),wk(l        ),nxr,1)
               call mapab4r(v(1,ei),wk(l+  nxyzw),nxr,1)
               if (if3d) 
     $         call mapab4r(w(1,ei),wk(l+2*nxyzw),nxr,1)
            else
               call mapab  (u(1,ei),wk(l        ),nxr,1)
               call mapab  (v(1,ei),wk(l+  nxyzw),nxr,1)
               if (if3d) 
     $         call mapab  (w(1,ei),wk(l+2*nxyzw),nxr,1)
            endif
         endif
         l = l+ldim*nxyzw
      enddo

 100  call err_chk(ierr,'Error reading restart data, in getv.$')
      return
      end
c-----------------------------------------------------------------------
      subroutine reade_init(fnames,ns)

      include 'SIZE'

      character*1 fnames(132*ns)
      character*132 fname

      open (unit=2,file='file.list')

      do i=1,ns
         call blank(fnames(1+(i-1)*132),132)
         call blank(fname,132)
         fname=':'
         read (2,1) fname
         if (indx1(fname,':',1).ne.0) goto 3
         call chcopy(fnames(1+(i-1)*132),fname,132)
         itmp=i
      enddo

      close (unit=2)

    3 continue
    1 format(a132)

      return
      end
c-----------------------------------------------------------------------
      subroutine reader_test

      include 'SIZE'
      include 'TOTAL'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*132 fname

      call blank(fname,132)
      fname='r0.f00001'

      ! init communicator
      call mpi_comm_split(nekcomm,mid,mid,mycomm,ierr)
      if (ierr.ne.0) call exitti('error splitting comm$',ierr)

      ! open
      call my_byte_open_mpi(mycomm,fname,mpi_fh,.true.,ierr)

      ! read header
      call byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)

      ! read elements sequentially

      do i=1,nel
         call byte_set_view(ioff_in,mpi_fh)
         call byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)
      enddo

      ! close
      call byte_close_mpi(mpi_fh,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_byte_open_mpi(nekcomm,fnamei,mpi_fh,ifro,ierr)

      include 'mpif.h'

      character fnamei*(*)
      logical ifro

      character*132 fname
      character*1   fname1(132)
      equivalence  (fname1,fname)

      l = ltrunc(fnamei,len(fnamei))
      if(l+1.gt.len(fname))
     $ call exitti('invalid string length$',l)

      call chcopy(fname1     ,fnamei ,l)
      call chcopy(fname1(l+1),char(0),1)

      imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
      if(ifro) then
        imode = MPI_MODE_RDONLY
      endif

#ifndef NOMPIIO
      call MPI_file_open(nekcomm,fname,imode,
     &                   MPI_INFO_NULL,mpi_fh,ierr)
#else
      call exitti('MPI_file_open unsupported!$',0)
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bopen(sourcefld)
c
c     generic field file reader
c     reads sourcefld and interpolates all avaiable fields
c     onto current mesh
c
c     memory requirement:
c     nelgs*nxs**ldim < np*(4*lelt*lx1**ldim)
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'

      character sourcefld*(*)

      common /scrcg/  pm1(lx1*ly1*lz1,lelv)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      character*1   hdr(iHeaderSize)

      integer*8 dtmp8

      logical if_byte_swap_test
      real*4 bytetest

      etime_t = dnekclock_sync()
      if(nio.eq.0) write(6,*) 'call gfldr ',trim(sourcefld)

      ! open source field file
      ierr = 0
      if(nid.eq.0) then
        open (90,file=sourcefld,status='old',err=100)
        close(90)
        goto 101
 100    ierr = 1
 101  endif
      call err_chk(ierr,' Cannot open source fld file!$')
      call byte_open_mpi(sourcefld,fldh_gfldr,.true.,ierr)

      ! read and parse header
      call byte_read_mpi(hdr,iHeaderSize/4,0,fldh_gfldr,ierr)
      call byte_read_mpi(bytetest,1,0,fldh_gfldr,ierr)

      call mfi_parse_hdr(hdr,ierr)
      call err_chk(ierr,' Invalid endian tag!$')

      return
      end

c-----------------------------------------------------------------------
      subroutine bread

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'

      nelgs   = nelgr
      nxs     = nxr
      nys     = nyr
      nzs     = nzr
      if(nzs.gt.1) then
        ldims = 3
      else
        ldims = 2
      endif
      if (ifgtim) time = timer

      ! distribute elements across all ranks
      nels = nelgs/np
      do i = 0,mod(nelgs,np)-1
         if(i.eq.nid) nels = nels + 1
      enddo
      nxyzs      = nxs*nys*nzs
      dtmp8      = nels
      ntots_b    = dtmp8*nxyzs*wdsizr
      rankoff_b  = igl_running_sum(nels) - dtmp8
      rankoff_b  = rankoff_b*nxyzs*wdsizr
      dtmp8      = nelgs
      nSizeFld_b = dtmp8*nxyzs*wdsizr
      noff0_b    = iHeaderSize + iSize + iSize*dtmp8

      ! do some checks
      if(ldims.ne.ldim)
     $ call exitti('ldim of source does not match target!$',0)
      if(ntots_b/wdsize .gt. ltots) then
        dtmp8 = nelgs
        lelt_req = dtmp8*nxs*nys*nzs / (np*ltots/lelt)
        lelt_req = lelt_req + 1
        if(nio.eq.0) write(6,*)
     $   'ABORT: buffer too small, increase lelt > ', lelt_req
        call exitt
      endif

      ifldpos = 0
      if(ifgetxr) then
        ! read source mesh coordinates
        call gfldr_getxyz(xm1s,ym1s,zm1s)
        ifldpos = ldim
      else
        call exitti('source does not contain a mesh!$',0)
      endif

      if(if_full_pres) then
        call exitti('no support for if_full_pres!$',0)
      endif

      ! initialize interpolation tool using source mesh
      nxf   = 2*nxs
      nyf   = 2*nys
      nzf   = 2*nzs
      nhash = nels*nxs*nys*nzs
      nmax  = 128

      ! read source fields and interpolate
c     if(ifgetur) then
c       if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading vel'
c       ntot = nx1*ny1*nz1*nelv
c       call gfldr_getfld(vx,vy,vz,ntot,ldim,ifldpos+1)
c       ifldpos = ifldpos + ldim
c     endif

      call byte_close_mpi(fldh_gfldr,ierr)
      etime_t = dnekclock_sync() - etime_t
      call fgslib_findpts_free(inth_gfldr)
      if(nio.eq.0) write(6,'(A,1(1g9.2),A)')
     &                   ' done :: gfldr  ', etime_t, ' sec'

      return
      end
c-----------------------------------------------------------------------
