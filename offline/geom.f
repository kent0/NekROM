c-----------------------------------------------------------------------
      subroutine setgeom(xyz,nel)

      include 'LVAR'
      include 'INTEG'

      real xyz(1)

      n=lxyz*nel

      call copy(xm1,xyz,n)
      call copy(ym1,xyz(n+1),n)
      if (ldim.eq.3) call copy(zm1,xyz(n*2+1),n)

      call geodat1

      return
      end
c-----------------------------------------------------------------------
