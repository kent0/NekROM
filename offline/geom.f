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
      call setrxp(rx,nel)

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