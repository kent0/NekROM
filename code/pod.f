c-----------------------------------------------------------------------
      subroutine setbases

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real u0(lt,3)
      common /scrk2/ a1(lt),a2(lt),a3(lt)

      if (nio.eq.0) write (6,*) 'inside setbases'

      n=lx1*ly1*lz1*nelt

      if (ifread) then
         call loadbases(ub,vb,wb,nb)
      else
         n=lx1*ly1*lz1*nelt

         do i=1,nb
            call rzero(ub(1,i),n)
            call rzero(vb(1,i),n)
            if (ldim.eq.3) call rzero(wb(1,i),n)
         enddo

         ! ub, vb, wb, are the modes
         do j=1,ns
         do i=1,nb
            call opadds(ub(1,i),vb(1,i),wb(1,i),
     $         us0(1,1,j),us0(1,2,j),us0(1,ldim,j),evec(j,i,1),n,2)
         enddo
         enddo

         call vnorm(ub,vb,wb)

         if (ifpod(2)) then
            do i=1,nb
               call rzero(tb(1,i),n)
               do j=1,ns
                  call add2s2(tb(1,i),ts0(1,j),evec(j,i,2),n)
               enddo
            enddo
            call snorm(tb)
         endif

         if (ifcintp) then
            do i=1,nb
               call rzero(cxb(1,i),n)
               call rzero(cyb(1,i),n)
               if (ldim.eq.3) call rzero(czb(1,i),n)
            enddo
            do j=1,ns
            do i=1,nb
               call opadds(cxb(1,i),cyb(1,i),czb(1,i),
     $            cs0(1,1,j),cs0(1,2,j),cs0(1,ldim,j),evec(j,i,0),n,2)
            enddo
            enddo
            call vnorm(cxb,cyb,czb)
         endif
      endif

      if (ifcdrag) then
         do i=0,nb
            call lap2d(a1,ub(1,i))
            call lap2d(a2,vb(1,i))
            if (ldim.eq.3) call lap2d(a3,wb(1,i))
            call opcmult(a1,a2,a3,param(2))
            call cint(fd1(1,i),ub(1,i),vb(1,i),wb(1,i))
            call cint(fd3(1,i),a1,a2,a3)
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exiting setbases'

      return
      end
c-----------------------------------------------------------------------
      subroutine ps2k(ck,ux,uub)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ck(0:nb,ls),ux(lt,ls),uub(lt,0:nb)

      nio=-1
      do i=1,ns
         call ps2b(ck(0,i),ux(1,i),uub)
      enddo
      nio=nid

      return
      end
c-----------------------------------------------------------------------
      subroutine pv2k(ck,usnap,uub,vvb,wwb)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ck(0:nb,ls),usnap(lt,ldim,ls),
     $     uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      nio=-1
      do i=1,ns
         call pv2b(ck(0,i),usnap(1,1,i),usnap(1,2,i),usnap(1,ldim,i),
     $        uub,vvb,wwb)
      enddo
      nio=nid

      return
      end
c-----------------------------------------------------------------------
      subroutine ps2b(coef,tt,sb)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrp/ t1(lt),t2(lt),t3(lt)

      real coef(0:nb),tt(lt),sb(lt,0:nb)

      if (nio.eq.0) write (6,*) 'inside ps2b'

      n=lx1*ly1*lz1*nelt

      call sub3(t1,tt,sb,n)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         ww=sip(sb(1,i),sb(1,i))
         vv=sip(sb(1,i),t1)
         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww,ips
      enddo

      if (nio.eq.0) write (6,*) 'exiting ps2b'

    1 format(' h10coef',1p3e16.8,1x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine pv2b(coef,ux,uy,uz,uub,vvb,wwb)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      common /scrp/ t1(lt),t2(lt),t3(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside pv2b'

      n=lx1*ly1*lz1*nelt

      call opsub3(t1,t2,t3,ux,uy,uz,uub,vvb,wwb)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         ww=vip(uub(1,i),vvb(1,i),wwb(1,i),uub(1,i),vvb(1,i),wwb(1,i))
         vv=vip(uub(1,i),vvb(1,i),wwb(1,i),t1,t2,t3)
         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww,ips
      enddo

      if (nio.eq.0) write (6,*) 'exiting pv2b'

    1 format(' h10coef',1p3e16.8,1x,a3)

      return
      end
c-----------------------------------------------------------------------
      function sip(t1,t2)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      if (ips.eq.'L2 ') then
         sip=wl2sip(t1,t2)
      else if (ips.eq.'H10') then
         sip=h10sip(t1,t2)
      else if (ips.eq.'HLM') then
         sip=hlmsip(t1,t2)
      else
         call exitti('did not provide supported inner product space$')
      endif

      return
      end
c-----------------------------------------------------------------------
      function vip(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      if (ips.eq.'L2 ') then
         vip=wl2vip(t1,t2,t3,t4,t5,t6)
      else if (ips.eq.'H10') then
         vip=h10vip(t1,t2,t3,t4,t5,t6)
      else if (ips.eq.'HLM') then
         vip=hlmvip(t1,t2,t3,t4,t5,t6)
      else
         call exitti('did not provide supported inner product space$')
      endif

      return
      end
c-----------------------------------------------------------------------
      function h10sip_vd(t1,t2)

      include 'SIZE'
      include 'SOLN'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      common /scrip/ t3(lt)

      call axhelm(t3,t1,vdiff,zeros,1,1)
      h10sip_vd=glsc2(t3,t2,lx1*ly1*lz1*nelt)

      return
      end
c-----------------------------------------------------------------------
      function h10sip(t1,t2)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      common /scrip/ t3(lt)

      call axhelm(t3,t1,ones,zeros,1,1)
      h10sip=glsc2(t3,t2,lx1*ly1*lz1*nelt)

      return
      end
c-----------------------------------------------------------------------
      function h10vip(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      common /scrip/ t7(lt),t8(lt),t9(lt)

      n=lx1*ly1*lz1*nelt

      call axhelm(t7,t1,ones,zeros,1,1)
      h10vip=glsc2(t7,t4,n)

      call axhelm(t8,t2,ones,zeros,1,2)
      h10vip=h10vip+glsc2(t8,t5,n)

      if (ldim.eq.3) then
         call axhelm(t9,t3,ones,zeros,1,3)
         h10vip=h10vip+glsc2(t9,t6,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      function hlmsip(t1,t2)

      include 'SIZE'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      hlmsip=h10sip(t1,t2)/ad_re+wl2sip(t1,t2)*ad_beta(1,3)/ad_dt

      return
      end
c-----------------------------------------------------------------------
      function hlmvip(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      hlmvip=h10vip(t1,t2,t3,t4,t5,t6)/ad_re
     $      +wl2vip(t1,t2,t3,t4,t5,t6)*ad_beta(1,3)/ad_dt

      return
      end
c-----------------------------------------------------------------------
      function wl2sip(t1,t2)

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt)

      n=lx1*ly1*lz1*nelt

      wl2sip = glsc3(t1,t2,bm1,n)

      return
      end
c-----------------------------------------------------------------------
      function wl2vip(t1,t2,t3,t4,t5,t6)

      include 'SIZE'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      wl2vip=op_glsc2_wt(t1,t2,t3,t4,t5,t6,bm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine setgram

      include 'SIZE'
      include 'MOR'
      include 'SOLN'

      if (.not.ifread) then
         jfield=ifield
         ifield=1
         if (ifpod(1)) call gengram(ug(1,1,1),us0,ns,ldim)
         ifield=2
         if (ifpod(2)) call gengram(ug(1,1,2),ts0,ns,1)
         ifield=jfield
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine setevec

      include 'SIZE'
      include 'MOR'

      if (.not.ifread) then
         do i=0,ldimt1
            if (ifpod(i)) call
     $         genevec(evec(1,1,i),eval(1,i),ug(1,1,i),i)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hlmgg(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengram HLM'

      n=lx1*ly1*lz1*nelt

      s2=ad_beta(1,3)/ad_dt

      if (ifield.eq.1) then
         s1=1./ad_re
      else
         s1=1./ad_pe
      endif

      do j=1,ms
         isd=1
         if (ifaxis) isd=2
         call axhelm(uu,s(1,1,j),ones,zeros,1,isd)
         if (mdim.ge.2) then
            isd=2
            if (ifaxis) isd=1
            call axhelm(vv,s(1,2,j),ones,zeros,1,isd)
         endif
         if (mdim.eq.3) call axhelm(ww,s(1,3,j),ones,zeros,1,3)
         do i=1,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
            gram(i,j)=s1*glsc2(uu,s(1,1,i),n)
     $               +s2*glsc3(s(1,1,i),s(1,1,j),bm1,n)
            if (mdim.ge.2) then
               gram(i,j)=gram(i,j)+s1*glsc2(vv,s(1,2,i),n)
     $                            +s2*glsc3(s(1,2,i),s(1,2,j),bm1,n)
            endif
            if (mdim.eq.3) then
               gram(i,j)=gram(i,j)+s1*glsc2(ww,s(1,3,i),n)
     $                            +s2*glsc3(s(1,3,i),s(1,3,j),bm1,n)
            endif
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j),'HLM'
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengram'
    1 format (' gram',i5,1p1e16.6,2x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine h10gg(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengram H10'

      n=lx1*ly1*lz1*nelt

      do j=1,ms
         isd=1
         if (ifaxis) isd=2
         call axhelm(uu,s(1,1,j),ones,zeros,1,isd)
         if (mdim.ge.2) then
            isd=2
            if (ifaxis) isd=1
            call axhelm(vv,s(1,2,j),ones,zeros,1,isd)
         endif
         if (mdim.eq.3) call axhelm(ww,s(1,3,j),ones,zeros,1,3)
         do i=1,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
            gram(i,j)=glsc2(uu,s(1,1,i),n)
            if (mdim.ge.2) then
               gram(i,j)=gram(i,j)+glsc2(vv,s(1,2,i),n)
            endif
            if (mdim.eq.3) then
               gram(i,j)=gram(i,j)+glsc2(ww,s(1,3,i),n)
            endif
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j),'H10'
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengram H10'
    1 format (' gram',i5,1p1e16.6,2x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine wl2gg(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgg/ uu(lt),vv(lt),ww(lt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengram'

      n=lx1*ly1*lz1*nelt

      do j=1,ms
         do i=1,ms ! Form the Gramian, U=U_K^T A U_K using H^1_0 Norm
            if (mdim.eq.1) then
               gram(i,j)=sip(s(1,1,i),s(1,1,j))
            else
               gram(i,j)=vip(s(1,1,i),s(1,2,i),s(1,ldim,i),
     $                       s(1,1,j),s(1,2,j),s(1,ldim,j))
            endif
         enddo
         if (nio.eq.0) write(6,1) j,gram(1,j),'L2 '
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengram'
    1 format (' gram',i5,1p1e16.6,2x,a3)

      return
      end
c-----------------------------------------------------------------------
      subroutine gengraml2(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrgram/ uw(lt),vw(lt),ww(lt)
      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (nio.eq.0) write (6,*) 'inside gengraml2'

      n=lx1*ly1*lz1*nelv

      do j=1,ms ! Form the Gramian, U=U_K^T A U_K using L2 Norm
      do i=1,ms
         gram(i,j)=glsc3(s(1,1,i),s(1,1,j),bm1,n)
         if (mdim.ge.2)
     $      gram(i,j)=gram(i,j)+glsc3(s(1,2,i),s(1,2,j),bm1,n)
         if (mdim.ge.3)
     $      gram(i,j)=gram(i,j)+glsc3(s(1,3,i),s(1,3,j),bm1,n)
      enddo
         if (nio.eq.0) write (6,1) j,gram(1,j)
      enddo

      if (nio.eq.0) write (6,*) 'exiting gengraml2'

    1 format (' gram',i5,' ',1p1e16.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine gengram(gram,s,ms,mdim)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real gram(ms,ms)
      real s(lt,mdim,ms)

      if (ips.eq.'L2 ') then
         call gengraml2(gram,s,ms,mdim)
      else if (ips.eq.'H10') then
         call h10gg(gram,s,ms,mdim)
      else
         call hlmgg(gram,s,ms,mdim)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine genevec(vec,val,gram,ifld)

      !!! does not work if ns.lt.ls !!!

      include 'SIZE'
      include 'TSTEP'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      integer icalld
      save    icalld
      data    icalld /0/
      
      common /scrgvec/ eye(ls,ls),gc(ls,ls),wk(ls,ls),eigv(ls,ls)

      real gram(ls,ls),vec(ls,nb),val(ls)

      if (nio.eq.0) write (6,*) 'inside genevec'

      if (icalld.eq.0) then
         icalld=1
      endif

      call rzero(eye,ls*ls)
      do j=1,ls
         eye(j,j) = 1.
      enddo

      call copy(gc,gram,ls*ls)

      call generalev(gc,eye,val,ls,wk)
      call copy(eigv,gc,ls*ls)

      do l = 1,nb
         call copy(vec(1,l),eigv(1,ls-l+1),ls) ! reverse order of eigv
      enddo

      do i=1,ns
         if (nio.eq.0) write (6,'(i5,1p1e16.6,3x,a,i1)')
     $      i,val(ns-i+1),'eval',ifld
      enddo

      if (nio.eq.0) write (6,*) 'exiting genevec'

      return
      end
c-----------------------------------------------------------------------
      subroutine vnorm(uub,vvb,wwb)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      jfield=ifield
      ifield=1
      nio=-1
      do i=1,nb
         p=vip(uub(1,i),vvb(1,i),wwb(1,i),uub(1,i),vvb(1,i),wwb(1,i))
         s=1./sqrt(p)
         call opcmult(uub(1,i),vvb(1,i),wwb(1,i),s)
      enddo
      nio=nid
      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine snorm(ssb)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ssb(lt,0:nb)

      nio=-1
      do i=1,nb
         p=sip(ssb(1,i),ssb(1,i))
         s=1./sqrt(p)
         call cmult(ssb(1,i),s,lx1*ly1*lz1*nelt)
      enddo
      nio=nid

      return
      end
c-----------------------------------------------------------------------
      subroutine h10pv2b(coef,ux,uy,uz,uub,vvb,wwb)

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      common /scrp/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside h10vpv2b'

      n=lx1*ly1*lz1*nelt

      call opsub3(t1,t2,t3,ux,uy,uz,uub,vvb,wwb)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      do i=1,nb
         call axhelm(t4,uub(1,i),ones,zeros,1,1)
         call axhelm(t5,vvb(1,i),ones,zeros,1,1)

         ww = glsc2(t4,uub(1,i),n)+glsc2(t5,vvb(1,i),n)
         vv = glsc2(t4,t1,n)+glsc2(t5,t2,n)

         if (ldim.eq.3) then
            call axhelm(t6,wwb(1,i),ones,zeros,1,1)
            ww = ww + glsc2(t6,wwb(1,i),n)
            vv = vv + glsc2(t6,t3,n)
         endif

         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww
      enddo

      if (nio.eq.0) write (6,*) 'exiting h10pv2b'

    1 format(' h10coef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
      subroutine hlmpv2b(coef,ux,uy,uz,uub,vvb,wwb)

      include 'SIZE'
      include 'MASS'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ux(lt),uy(lt),uz(lt),uub(lt,0:nb),vvb(lt,0:nb),wwb(lt,0:nb)

      common /scrp/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt)

      real coef(0:nb)

      if (nio.eq.0) write (6,*) 'inside hlmvpv2b'

      n=lx1*ly1*lz1*nelt

      call opsub3(t1,t2,t3,ux,uy,uz,uub,vvb,wwb)

      coef(0) = 1.
      if (nio.eq.0) write (6,1) coef(0),coef(0),1.

      s1=1./ad_re
      s2=ad_beta(1,3)/ad_dt

      do i=1,nb
         call axhelm(t4,uub(1,i),ones,zeros,1,1)
         call axhelm(t5,vvb(1,i),ones,zeros,1,1)

         ww=s1*(glsc2(t4,uub(1,i),n)+glsc2(t5,vvb(1,i),n))
         vv=s1*(glsc2(t4,t1,n)+glsc2(t5,t2,n))

         vv=vv+s2*(glsc3(t1,uub(1,i),bm1,n)+glsc3(t2,vvb(1,i),bm1,n))
         ww=ww+s2*(glsc3(uub(1,i),uub(1,i),bm1,n)
     $            +glsc3(vvb(1,i),vvb(1,i),bm1,n))

         if (ldim.eq.3) then
            call axhelm(t6,wwb(1,i),ones,zeros,1,1)
            ww=ww+s1*glsc2(t6,wwb(1,i),n)
            vv=vv+s1*glsc2(t6,t3,n)
            vv=vv+s2*glsc3(t3,wwb(1,i),bm1,n)
            ww=ww+s2*glsc3(wwb(1,i),wwb(1,i),bm1,n)
         endif

         coef(i) = vv/ww
         if (nio.eq.0) write (6,1) coef(i),vv,ww
      enddo

      if (nio.eq.0) write (6,*) 'exiting hlmpv2b'

    1 format(' hlmcoef',1p3e16.8)

      return
      end
c-----------------------------------------------------------------------
