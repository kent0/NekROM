c-----------------------------------------------------------------------
      subroutine setb(b,u,v,nb,n)

      real b(nb,nb)
      real u(n,nb),v(n,nb)

      do k=1,nb
      do j=1,nb
      do i=1,n
         b(j,k)=b(j,k)+u(i,k)*v(i,j)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine seta(a,u,v,nb,n)

      real a(nb,nb)
      real u(n,nb),v(n,nb)

      do k=1,nb
      do j=1,nb
      do i=1,n
         a(j,k)=a(j,k)+u(i,k)*v(i,j)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setc(c,u,v,w,nb,n)

      real c(nb,nb,nb)
      real u(n,nb),v(n,nb),w(n,nb)

      do l=1,nb
      do k=1,nb
      do j=1,nb
      do i=1,n
         c(j,k,l)=c(j,k,l)+u(i,l)*v(i,j)*w(i,k)
      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------