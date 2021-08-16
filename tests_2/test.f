c-----------------------------------------------------------------------
c     include 't.f'
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,eg

      udiff =0.
      utrans=0.
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,eg

      ffx = 0.0
      ffy = 0.0
      ffz = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,eg

      qvol   = 0.0
      source = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux=1-y*y
      uy   = 0
      uz   = 0
      temp = 0
      pa   = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      n=lx1*ly1*lz1*nelv

      s=.5/pi

      do i=1,lx1*ly1*lz1*nelt
         x=xm1(i,1,1,1)
         y=ym1(i,1,1,1)
         sobj(i,1,1,1)=mod(atan2(y,x)*s+1.,1.)
      enddo

      call set_obj

      return
      end
c----------------------------------------------------------------------
      subroutine usrdat3
      return
      end
c----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,eg)
c     NOTE ::: This is not guaranteed to be called by every process
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer e,eg

      e = gllel(eg)

      ux=1-y*y
      uy=0.0
      uz=0.0
      temp=0.0
      pa  =0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

c     call test
      call rom_update
      call exitt0

      return
      end
c-----------------------------------------------------------------------
      subroutine set_obj  ! define objects for surface integrals

      include 'SIZE'
      include 'TOTAL'

      integer e,f

c     Define new objects

      nobj = 1
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F' ! 'F'
         hcode( 2,ii) = 'F' ! 'F'
         hcode( 3,ii) = 'F' ! 'F'
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj

      if (maxobj.lt.nobj) write(6,*) 'increase maxobj in SIZE'
      if (maxobj.lt.nobj) call exitt

      do e=1,nelv
      do f=1,2*ndim
         if (cbc(f,e,1).eq.'W  ') then
            iobj = 0
            iobj=1              ! cylinder wall
            if (iobj.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               ieg = lglel(e)
               object(iobj,mem,1) = ieg
               object(iobj,mem,2) = f
    1          format(6i9,a4)
            endif
c
         endif
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
