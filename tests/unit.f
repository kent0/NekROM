c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      include 'rom.f'
c-----------------------------------------------------------------------
      subroutine grammian_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      logical iflag
      real vv(ls,ls)

      param(33) = 1
      if (iflag) param(33) = 0
      param(34) = 1
      param(35) = 2

      call rom_setup

      iexit=0

      call copy(vv,uu,ls*ls)
      call gengram

      s1=0.
      s2=0.
      s3=0.

      do j=1,ls
      do i=1,ls
         s1=s1+(uu(i,j)-uu(j,i))**2
         s2=s2+(uu(i,j)-vv(i,j))**2
         s3=s3+vv(i,j)**2
         if (nio.eq.0) write (6,*) 'gram',i,j,uu(i,j),vv(i,j)
      enddo
      enddo

      esym=sqrt(s1/s3)
      edif=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'esym',esym,s1,s3
      if (nio.eq.0) write (6,*) 'edif',edif,s2,s3

      if (esym.gt.1e-14) iexit=iexit+1
      if (edif.gt.1e-14) iexit=iexit+2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine initial_condition_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      logical iflag
      real u0(0:nb)

      param(33) = 1
      if (iflag) param(33) = 0
      param(34) = 1
      param(35) = 2

      call rom_setup

      call copy(u0,u,nb+1)
      call setu

      s1=0.
      s2=0.

      do i=0,nb
         write (6,*) u0(i),u(i,1),'ic'
      enddo

      do i=1,nb
         s1=s1+(u0(i)-u(i,1))**2
         s2=s2+u0(i)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.5e-15) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine a0_unit(iflag)

      include 'SIZE'
      include 'INPUT'
      include 'MOR'

      logical iflag

      real aa(0:nb,0:nb)

      param(33) = 1
      if (iflag) param(33) = 0
      param(34) = 1
      param(35) = 2

      call rom_setup

      iexit=0

      call copy(aa,a0,(nb+1)**2)
      call seta

      s1=0.
      s2=0.
      s3=0.

      do j=0,nb
      do i=0,nb
         s1=s1+(aa(i,j)-a0(i,j))**2
         s2=s2+(a0(i,j)-a0(j,i))**2
         s3=s3+aa(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      if (edif.gt.1.e-16) iexit=iexit+1
      if (esym.gt.2.4e-15) iexit=iexit+2

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+a0(i,j)**2
         s2=s2+a0(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2
      if (.not.ifl2.and.edia.gt.5.8e-15) iexit=iexit+4

      s1=0.

      do i=1,nb
         s1=s1+(a0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2
      if (.not.ifl2.and.euni.gt.4.8e-15) iexit=iexit+8

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine b0_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      logical iflag

      real bb(0:nb,0:nb)

      param(33) = 1
      if (iflag) param(33) = 0
      param(34) = 1
      param(35) = 2

      call rom_setup

      iexit=0

      call copy(bb,b0,(nb+1)**2)
      call setb

      s1=0.
      s2=0.
      s3=0.

      do j=0,nb
      do i=0,nb
         s1=s1+(bb(i,j)-b0(i,j))**2
         s2=s2+(b0(i,j)-b0(j,i))**2
         s3=s3+bb(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s3)
      esym=sqrt(s2/s3)

      if (edif.gt.1.e-16) iexit=iexit+1
      if (esym.gt.1.e-16) iexit=iexit+2

      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3
      if (nio.eq.0) write (6,*) 'esym',esym,s2,s3

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,nb
         if (i.ne.j) s1=s1+b0(i,j)**2
         s2=s2+bb(i,j)**2
      enddo
      enddo

      edia=sqrt(s1/s2)

      if (ifl2.and.edia.gt.1.0e-14) iexit=iexit+4
      if (nio.eq.0) write (6,*) 'edia',edia,s1,s2

      s1=0.

      do i=1,nb
         s1=s1+(b0(i,i)-1.)**2
      enddo

      euni=sqrt(s1/s2)

      if (ifl2.and.euni.gt.1.0e-14) iexit=iexit+8
      if (nio.eq.0) write (6,*) 'euni',euni,s1,s2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine c0_unit(iflag)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MOR'

      logical iflag

      real cc(lcloc), cglob(nb,nb+1,nb+1)

      param(33) = 1
      if (iflag) param(33) = 0
      param(34) = 1
      param(35) = 2

      call rom_setup

      iexit=0

      call copy(cc,clocal,nb*(nb+1)**2)
      call setc

      s1=0.
      s2=0.
      s3=0.

      call rzero(cglob,nb*(nb+1)**2)

      do jc=1,nb*(nb+1)**2
         i=icloc(1,jc)
         j=icloc(2,jc)
         k=icloc(3,jc)

         cglob(i,j,k)=cglob(i,j,k)+clocal(jc)
         cglob(k,j,i)=cglob(k,j,i)-cc(jc)

         s1=s1+(cc(jc)-clocal(jc))**2
         s3=s3+cc(jc)**2
         write (6,*) 'cc',clocal(jc),cc(jc)
      enddo

      do k=1,nb
      do j=1,nb
      do i=1,nb
         s2=s2+(cglob(i,j,k)+cglob(k,j,i))**2
      enddo
      enddo
      enddo

      edif=sqrt(s1/s3)
      if (edif.gt.1.e-16) iexit=iexit+1
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s3

      eskew=sqrt(s2/s3)
c     if (eskew.gt.1.e-16) iexit=iexit+2
      if (nio.eq.0) write (6,*) 'eskew',eskew,s1,s2

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
