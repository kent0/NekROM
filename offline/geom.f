c-----------------------------------------------------------------------
		  subroutine setgeom

        include 'LVAR'
        include 'OFFLINE'
        include 'INTEG'

        n=lxyz*nel

        call copy(xm1,xupt,n)
        call copy(ym1,xupt(n+1),n)
        if (ldim.eq.3) call copy(zm1,xupt(n*2+1),n)

        call geodat1

		  return
		  end
c-----------------------------------------------------------------------
