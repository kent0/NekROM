c-----------------------------------------------------------------------
      subroutine mode_normalize(aa,m,n)

      real aa(m,n)
      real vlngth
      integer m,n

      do jj=1,n

         vlngth = vlsc2(aa(1,jj),aa(1,jj),m)
         vlngth = 1./sqrt(vlngth)
         do ii=1,m
            aa(ii,jj) = aa(ii,jj)*vlngth
         enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine rand_initial()

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real mode_unfold(nb*ltr,3)
      integer,parameter :: seed = 86456
  
      call srand(seed)
      
      do jj=1,3
         do ii=1,nb*ltr
            mode_unfold(ii,jj) = rand()
         enddo
         call mode_normalize(mode_unfold(1,jj),nb,ltr)
      enddo

      do jj=1,3
         call check_normalize(mode_unfold(1,jj),nb,ltr)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine check_normalize(aa,m,n)

      real aa(m,n)
      real vlngth
      integer m,n

      do jj=1,n

         vlngth = vlsc2(aa(1,jj),aa(1,jj),m)
         write(6,*)jj,vlngth,'norm of jj colum'
      enddo
      
      return
      end
