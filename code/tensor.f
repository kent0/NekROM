c-----------------------------------------------------------------------
      subroutine CP_ALS(cl)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real fcm(0:(lub+1)*ltr-1,3)
      real fcmpm(ltr*ltr,3)
      real lsm(ltr*ltr,3),lsminv(ltr*ltr,3)
      real tmp(ltr*ltr),tmp_wrk(ltr)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real lsr((lub+1)*ltr)
      integer mode,maxit

      maxit = 1000

      call rand_initial(fcm)

      do mode=2,3
         call set_product_matrix(fcm(0,mode),fcmpm(1,mode),lub+1,ltr)
      enddo

      do ii=1,1!maxit
         do mode=1,3
            call mttkrp(lsr,cl,fcm,mode)
            call set_lsm(lsm,fcmpm,mode,ltr)
            call invmat(lsminv(1,mode),tmp,lsm(1,mode),tmp_wrk,ltr)
            if (mode.eq.1) then
               do jj=1,ltr
                  call mxm(lsr,lub+1,lsminv(1,mode),ltr,
     $            fcm(0+(lub+1)*(jj-1),mode),1)
               enddo
               do jj=0,(lub+1)*ltr-1
                  write(6,*)jj,fcm(jj,mode),'check'
               enddo
            elseif (mode.ne.1) then
               do jj=1,ltr
                  call mxm(lsr,lub+1,lsminv(1,mode),ltr,
     $            fcm(0+(lub+1)*(jj-1),mode),1)
               enddo
            endif
            call set_product_matrix(fcm(0,mode),fcmpm(1,mode),lub+1,ltr)
         enddo
      enddo

      

      return
      end
c-----------------------------------------------------------------------
      subroutine mttkrp(lsr,cl,fcm,mode)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real lsr((lub+1)*ltr)
      real fcm(0:(lub+1)*ltr-1,3)
      real cl(ic1:ic2,jc1:jc2,kc1:kc2)
      real cm(ic1:ic2,jc1:jc2,ltr)
      real cm2(ic1:ic2,ltr,kc1:kc2)
      integer mode,tr

      write(6,*)ic1,ic2,jc1,jc2,kc1,kc2,'index'

 
      if (mode.eq.1) then

         call rzero(lsr,(lub+1)*ltr)

         ! construct temporary mttkrp
         do tr=1,ltr
            call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $               fcm(kc1+(lub+1)*(tr-1),3),(kc2-kc1+1),cm(1,0,tr),1)
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,ltr
            do jj=jc1,jc2
               call add2s2(lsr(1+(lub+1)*(tr-1)),cm(1,jj,tr),
     $                     fcm(jj+(lub+1)*(tr-1),2),lub)
            enddo
         enddo

         write(6,*) 'First cp gradient'
         do ii=1,(lub+1)*ltr
            write(6,*)ii,lsr(ii)
         enddo

      elseif (mode.eq.2) then

         call rzero(lsr,(lub+1)*ltr)

         ! construct temporary mttkrp
         do tr=1,ltr
            call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
     $               fcm(kc1+(lub+1)*(tr-1),3),(kc2-kc1+1),cm(1,0,tr),1)
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,ltr
            do jj=jc1,jc2
               lsr((jj+1)+(lub+1)*(tr-1)) = vlsc2(cm(1,jj,tr),
     $         fcm(1+(lub+1)*(tr-1),1),lub)
            enddo
         enddo

         write(6,*) 'Second cp gradient'
         do ii=1,(lub+1)*ltr
            write(6,*)ii,lsr(ii)
         enddo

      elseif (mode.eq.3) then

         call rzero(lsr,(lub+1)*ltr)

         ! construct temporary mttkrp
         do kk=kc1,kc2
         do tr=1,ltr
            call mxm(cl(ic1,jc1,kk),(ic2-ic1+1),
     $               fcm(jc1+(lub+1)*(tr-1),2),(jc2-jc1+1),
     $               cm2(1,tr,kk),1) 
c           call mxm(cl,(ic2-ic1+1)*(jc2-jc1+1),
c    $               fcm(kc1+(lub+1)*(tr-1),3),(kc2-kc1+1),cm(1,0,tr),1)
         enddo
         enddo

         ! temporary mttkrp with factor matrix
         do tr=1,ltr
            do kk=kc1,kc2
               lsr((kk+1)+(lub+1)*(tr-1)) = vlsc2(cm2(1,tr,kk),
     $         fcm(1+(lub+1)*(tr-1),1),lub)
            enddo
         enddo

         write(6,*) 'Third cp gradient'
         do ii=1,(lub+1)*ltr
            write(6,*)ii,lsr(ii)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_lsm(lsm,fcmpm,mode,nn)

      real fcmpm(nn*nn,3)
      real lsm(nn*nn,3)
      integer mode,idx,nn

      ! Hadamard product

      if (mode.eq.1) then 
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,1),fcmpm(idx,2),fcmpm(idx,3),nn)
         enddo
c        do ii=1,nn*nn
c           write(6,*)ii,lsm(ii,1),fcmpm(ii,2),fcmpm(ii,3),'check'
c        enddo
      elseif (mode.eq.2) then
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,2),fcmpm(idx,1),fcmpm(idx,3),nn)
         enddo
      elseif (mode.eq.3) then
         do ii=1,nn
            idx = 1+(ii-1)*nn
            call col3(lsm(idx,3),fcmpm(idx,1),fcmpm(idx,2),nn)
         enddo
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_relerr(relerr,lsr,fcm,lsm,fcmpm,norm_c,mm,nn)

      real relerr,norm_c
      real inner_prod,norm_approx 
      real lsr(mm,nn),lsm(nn,nn)
      real fcm(mm,nn),fcmpm(nn,nn)
      real tmp1(mm,nn),tmp2(nn,nn)
      real tmp3(nn),tmp4(nn)
      integer mm,nn

      inner_prod=0.
      norm_approx=0.

      write(6,*)'check dimension',mm,nn
      write(6,*)'check lsr,fcm'
      do jj=1,nn
      do ii=1,mm
         write(6,*)ii,jj,lsr(ii,jj),fcm(ii,jj),'check'
      enddo
      enddo
      do ii=1,nn
         call col3(tmp1(1,ii),lsr(1,ii),fcm(1,ii),mm)
      enddo
      write(6,*)'check tmp1'
      do jj=1,nn
      do ii=1,mm
         write(6,*)ii,jj,tmp1(ii,jj),'check'
      enddo
      enddo

      inner_prod = -2*vlsc2(tmp1,tmp1,mm*nn)

      do ii=1,nn
         call col3(tmp2(1,ii),lsm(1,ii),fcmpm(1,ii),nn)
      enddo
      call rone(tmp3,nn)
      call mxm(tmp2,nn,tmp3,nn,tmp4,1)
      norm_approx = vlsc2(tmp4,tmp3,nn)

      write(6,*)'inner_prod',inner_prod
      write(6,*)'norm_c',norm_c
      write(6,*)'norm_approx',norm_approx
      relerr = sqrt(norm_c + inner_prod + norm_approx)/sqrt(norm_c)
      write(6,*)'relerr',relerr



      return
      end
c-----------------------------------------------------------------------
      subroutine compute_cp_weight(cp_weight,aa,m,n)

      real aa(m,n)
      real cp_weight(n)
      integer m,n

      do jj=1,n
         cp_weight(jj) = sqrt(vlsc2(aa(1,jj),aa(1,jj),m))
      enddo

      return
      end
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
      subroutine rand_initial(fcm)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real fcm((lub+1)*ltr,3)
      integer,parameter :: seed = 86456
  
      call srand(seed)
      
      do jj=2,3
         do ii=1,(lub+1)*ltr
            fcm(ii,jj) = rand()
         enddo
         call mode_normalize(fcm(1,jj),lub+1,ltr)

         call check_normalize(fcm(1,jj),lub+1,ltr)
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
         write(6,*)jj,vlngth,'jj'
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_product_matrix(aa,bb,m,n)

      real aa(m,n)
      real bb(n,n)

      do jj=1,n
         do ii=1,n
            bb(ii,jj) = vlsc2(aa(1,ii),aa(1,jj),m)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine read_cp_weight

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetmode/ wk(ltr+1)

      logical ifexist

      inquire (file='./ops/lambda',exist=ifexist)
      if (ifexist) call read_serial(cp_w,ntr,'./ops/lambda ',wk,nid)

      return
      end
c-----------------------------------------------------------------------
      subroutine read_cp_mode

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrsetmode/ wk1(lt)

      if (nio.eq.0) write (6,*) 'reading A1...'
      call read_mat_serial(cua,nb,ntr,'./ops/cua ',mb,ntr,wk1,nid)

      if (nio.eq.0) write (6,*) 'reading A2...'
      call read_mat_serial(cub,nb+1,ntr,'./ops/cub ',mb+1,ntr,wk1,nid)

      if (nio.eq.0) write (6,*) 'reading A3...'
      call read_mat_serial(cuc,nb+1,ntr,'./ops/cuc ',mb+1,ntr,wk1,nid)

      return
      end
