c-----------------------------------------------------------------------
      subroutine rsnapsm(v,s,ieg0,ieg1)

      include 'POST'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real v(lxyz*(ieg1-ieg0+1),ldim,ns)
      real s(lxyz*(ieg1-ieg0+1),ns)

      n=lxyz*(ieg1-ieg0+1)

      melt=(lelt-3)*lxyz

      ifcread=.true.
      ttime=dnekclock()

      nio=-1
      do is=1,ms(mid+1)
         js=ilgls(is)
         call rfldm_open(fnames(1+(js-1)*132),ieg0,ieg1,ifcread)
         call rfldm_read(v(1,1,is),v(1,2,is),v(1,ldim,is),s(1,is),
     $      ieg0,ieg1,ifcread,iftherm)
         call rfldm_close(ifcread)
      enddo
      nio=nid

      read_time=read_time+dnekclock()-ttime

      return
      end
c-----------------------------------------------------------------------
      subroutine rsnapsd(v,ieg0,ieg1)

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
      subroutine shift_setup(igsh,ncomm,iw,n,ns)

      include 'SIZE'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myshift/ m

      integer*8 iw(n,2)
      integer igsh(2)

      if (ns.eq.1) return

      if (mod(ns,2).eq.1)
     $   call exitti('no support for odd number of mpi-ranks$',ns)

      m=n

      if (mid.le.(ns-1)) then
         call izero(igsh,2)
         do i=1,n
            i1=nid/2
            i2=mod(nid-1+ns,ns)/2
            irem=mod(nid,2).eq.0
            iw(i,1)=i+i1*n
            iw(i,2)=i+i2*n
         enddo
         call fgslib_gs_setup(igsh,iw,n,ncomm,mp)
         call fgslib_gs_setup(igsh(2),iw(1,2),n,ncomm,mp)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine shift(igsh,u,w,n)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myshift/ m
      common /ptime/ aa_time,bb_time,cc_time,zz_time,read_time,comm_time

      integer igsh(2)
      real u(n),w(m,2)

      if (mp.eq.1) return
      if (n.gt.m) call exitti('n > m in shift$',n-m)

      if (mod(mp,2).eq.1)
     $   call exitti('no support for odd number of mpi-ranks$',mp)

      if (mod(mid,2).ne.0) then
         call copy(w,u,n)
         call rzero(w(1,2),n)
      else
         call copy(w(1,2),u,n)
         call rzero(w,n)
      endif

      ttime=dnekclock()
      call fgslib_gs_op(igsh,w,1,1,0)
      call fgslib_gs_op(igsh(2),w(1,2),1,1,0)
      comm_time=comm_time+dnekclock()-ttime

      if (mod(mid,2).ne.0) then
         call copy(u,w(1,2),n)
      else
         call copy(u,w,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rfldm_setup

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
      subroutine rfldm_open(fname_in,ieg0,ieg1,ifcread)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      character*132  fname_in
      logical ifcread

      character*132  fname
      character*1    fnam1(132)
      equivalence   (fnam1,fname)

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)

c     if (nio.eq.0) write (6,1) nid,fname_in
      ! add path
      call blank(fname,132)
      lenp = ltrunc(path,132)
      lenf = ltrunc(fname_in,132)
      call chcopy(fnam1(1),path,lenp)
      call chcopy(fnam1(lenp+1),fname_in,lenf)

      nio=-1
      call my_mfi_prepare(fname,ieg0,ieg1,ifcread) ! determine reader nodes +
                                                   ! read hdr + element mapping 
      nio=nid

    1 format (i8,2x,a60,' fname_in')

      return
      end
c-----------------------------------------------------------------------
      subroutine rfldm_read(ux,uy,uz,ut,ieg0,ieg1,ifcread,iftherm)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)
      common /scrread/ i1(lelt),i2(lelt),i3(lelt),i4(lelt)

      logical ifcread,iftherm

      real ux(lx1*ly1*lz1,ieg1-ieg0+1)
      real uy(lx1*ly1*lz1,ieg1-ieg0+1)
      real uz(lx1*ly1*lz1,ieg1-ieg0+1)
      real ut(lx1*ly1*lz1,ieg1-ieg0+1)

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      call rzero(wk,7*lx1*ly1*lz1*nelt)

      iofldsr=0
      if (ifgetxr) iofldsr=ldim

      nelr=ieg1-ieg0+1

      icount=1

      call isort(i1,i2,nelr)

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
         offs = offs0 + iofldsr*stride + !ldim*strideB +
     $      ldim*(ieg-1)*nxyzr8*wdsizr
         if (ifcread) then
            call byte_seek(offs/4,ierr)
         else
            call byte_set_view(offs,ifh_mbyte)
         endif
         call mfi_getw(wk(iloc),nwk,ifcread)

         iloc=iloc+nwk
      enddo

      nelr=neltmp
      nwk=nelr*ldim*nxyzr8

      if (if_byte_sw) then
      if (wdsizr.eq.8) then
         call byte_reverse8(wk,nwk*2,ierr)
      else
         call byte_reverse(wk,nwk,ierr)
      endif
      endif

      lxyz=lx1*ly1*lz1
      do ie=1,nelr
         call copy(ux(1,i2(ie)),wk(1+(ie-1)*lxyz*ldim),lxyz)
         call copy(uy(1,i2(ie)),wk(1+(ie-1)*lxyz*ldim+lxyz),lxyz)
         if (ldim.eq.3) call
     $      copy(uz(1,i2(ie)),wk(1+(ie-1)*lxyz*ldim+2*lxyz),lxyz)
      enddo

      if (iftherm) then
         iofldsr=iofldsr+ldim
         if (ifgetpr) iofldsr=iofldsr+1
         iloc=1
         do ig=1,ng
            ieg=i3(ig)
            nelr=i4(ig)
            nwk=nelr*nxyzr8
            offs = offs0 + iofldsr*stride + !ldim*strideB +
     $         (ieg-1)*nxyzr8*wdsizr
            if (ifcread) then
               call byte_seek(offs/4,ierr)
            else
               call byte_set_view(offs,ifh_mbyte)
            endif
            call mfi_getw(wk(iloc),nwk,ifcread)
            iloc=iloc+nwk
         enddo

         nelr=neltmp
         nwk=nelr*nxyzr8

         if (if_byte_sw) then
         if (wdsizr.eq.8) then
            call byte_reverse8(wk,nwk*2,ierr)
         else
            call byte_reverse(wk,nwk,ierr)
         endif
         endif

         do ie=1,nelr
            call copy(ut(1,i2(ie)),wk(1+(ie-1)*lxyz),lxyz)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rfldm_close(ifcread)

      include 'SIZE'
      include 'INPUT'
      include 'RESTART'

      logical ifcread

      if (.not.ifcread) then
         if (nid.eq.pid0r) call byte_close_mpi(ifh_mbyte,ierr)
      else
         if (nid.eq.pid0r) call byte_close(ierr)
      endif

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
      subroutine mfi_getw(wk,n,ifcread)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real wk(1)
      logical ifcread

      nxyzr  = n
      if (wdsizr.eq.8) nxyzr = 2*nxyzr

      ierr = 0
      if (ifcread) then
         call byte_read(wk,nxyzr,ierr)
         ierr=ierr*10*(nid+1)
      else
         call byte_read_mpi(wk,nxyzr,-1,ifh_mbyte,ierr)
      endif

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
      subroutine rflist(fnames,ns)

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
      subroutine setup_mycomm

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
      subroutine my_mfi_prepare(hname,ieg0,ieg1,ifcread)  ! determine which nodes are readers
      character*132 hname

      include 'SIZE'
      include 'PARALLEL'
      include 'RESTART'
      include 'INPUT'

      common /ipparallel/ nps,lenpb,melt,itmp(lx1*ly1*lz1*lelt)
      common /scrread/ i1(lelt),i2(lelt),i3(lelt),i4(lelt)

      integer stride
      character*132 hdr, hname_
      logical if_byte_swap_test,ifcread
      real*4 bytetest

      integer*8 offs0,offs

      ierr = 0

      call chcopy(hname_,hname,132)
      call addfid(hname_,0)
      call byte_open(hname_,ierr)
      if(ierr.ne.0) goto 101

      call blank     (hdr,iHeaderSize)
      call byte_read (hdr,iHeaderSize/4,ierr)
      if(ierr.ne.0) goto 101

      call byte_read (bytetest,1,ierr)
      if(ierr.ne.0) goto 101

      if_byte_sw = if_byte_swap_test(bytetest,ierr) ! determine endianess
      if(ierr.ne.0) goto 101

      call byte_close(ierr)
      if(ierr.ne.0) goto 101

 101  continue

      call mfi_parse_hdr(hdr,ierr)

      ifmpiio=.not.ifcread

      if (ifmpiio) then
         if (nelt.gt.lelr) then
            write(6,*) 'ERROR: increase lelr in SIZE!', lelr, nelt
            call exitt
         endif
      else
         if (nelr.gt.lelr) then
            write(6,*) 'ERROR: increase lelr in SIZE!', lelr, nelr
            call exitt
         endif
      endif

      if (.not.ifmpiio) then

        stride = np / nfiler
        if (stride.lt.1) then
           write(6,*) nfiler,np,'  TOO MANY FILES, mfi_prepare'
           call exitt
        endif

        pid0r = nid
        pid1r = nid + stride
        fid0r = nid / stride

        call blank(hdr,iHeaderSize)

        call addfid(hname,fid0r)
        call byte_open(hname,ierr)
        if(ierr.ne.0) goto 102

        call byte_read(hdr,iHeaderSize/4,ierr)  

        if(ierr.ne.0) goto 102

        call byte_read(bytetest,1,ierr) 
        if(ierr.ne.0) goto 102

        call mfi_parse_hdr(hdr,ierr) ! replace hdr with correct one 

        ic=1
        ie=1

        do while (ic.le.nelgv)
           mel=min(melt,nelgv)
           call byte_read(er,mel,ierr)! get element mapping
           if (if_byte_sw) call byte_reverse(er,mel,ierr)
           jc=1
           do while (ie.le.mel)
              if (er(ie).ge.ieg0.and.er(ie).le.ieg1) then
                 i1(er(ie)-ieg0+1)=ie
                 jc=jc+1
              endif
              ie=ie+1
           enddo
           ic=ic+mel
        enddo
      else

        pid0r  = nid
        pid1r  = nid
        offs0  = iHeaderSize + 4
        nfiler = np
 
        ! number of elements to read 
        nelr = nelgr/np
        do i = 0,mod(nelgr,np)-1
           if(i.eq.nid) nelr = nelr + 1
        enddo
        nelBr = igl_running_sum(nelr) - nelr 
        offs = offs0 + nelBr*isize

        call addfid(hname,fid0r)
        if(nio.eq.0) write(6,*) '      FILE:',hname
        call byte_open_mpi(hname,ifh_mbyte,.true.,ierr)

        if(ierr.ne.0) goto 102
        call byte_set_view(offs,ifh_mbyte)
        call byte_read_mpi(er,nelr,-1,ifh_mbyte,ierr)
        if(ierr.ne.0) goto 102
        if(if_byte_sw) call byte_reverse(er,nelr,ierr)

      endif

 102  continue
c     call err_chk(ierr,'Error reading header/element map.$')

      ifmpiio = .false.
      if (nfiler.eq.1 .and. abs(param(67)).eq.6) ifmpiio = .true.
#ifdef NOMPIIO
      ifmpiio = .false.
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_parallel(a,n,fname,nid)

      real a(n)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*128 fname
      character*128 fntrunc

      call blank(fntrunc,128)

      len=ltruncr(fname,128)
      call chcopy(fntrunc,fname,len)

      do id=0,mp-1
         if (id.eq.nid) call dump_parallel_helper(a,n,fntrunc,id.eq.0)
         call nekgsync
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_parallel_helper(a,n,fname,if0)

      real a(n)
      character*128 fname
      logical if0

      if (if0) then
         open (unit=12,file=fname)
      else
         open (unit=12,file=fname,access='append',status='old')
      endif

      do i=1,n
         write (12,1) a(i)
      enddo

      close (unit=12)
    1 format(1pe24.16)

      return
      end
c-----------------------------------------------------------------------
      subroutine shift_test

      include 'POST'
      include 'MASS'
      include 'GEOM'
      include 'PARALLEL'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*127 flist

      nelp=47

      nsnap=ns

      ns=ls
      call rflist(fnames,ns)

      iftherm=.false.

      call ilgls_setup(ilgls,msr,ms,ns,np,nid)
      call iglls_setup(iglls,itmp,ms,ns,np,nid)

      call ilgls_setup(ilglsp,msrp,msp,nb+1,min(np,nb+1),nid)

      call shift_setup(igsh,nekcomm,itmp,3,np)

      do i=1,3
         iw1(i)=i+mid*3
      enddo
      call shift(igsh,iw1,iw2,3)

      do id=0,mp-1
         if (id.eq.mid) then
            do i=1,3
               write (6,*) id,i,iw1(i),'shift'
            enddo
         endif
         call nekgsync
      enddo

      return
      end
c-----------------------------------------------------------------------
