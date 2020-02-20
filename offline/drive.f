c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        nsg=lsg
        nsg=5

        ie=1

        neg=512
        mel=8
        ng=neg/mel

        do i=1,ng
           write (6,*) i,'offline'
           ieg(1)=ie
           ieg(2)=ie+mel-1
           nel=ieg(2)-ieg(1)+1

           call loadsnaps(buf,ieg,indxr,nsg)
           call setgeom(buf,nel)
           call setops(ga,gb,gc,gt,buf(nel*lxyz*ldim+1),nel,nsg,ldim)
           ie=ie+mel
        enddo

        call gop(ga,gt,'+  ',nsg*nsg)
        call gop(gb,gt,'+  ',nsg*nsg)
c       call gop(gc,gt,'+  ',nsg*nsg*nsg)

        do j=1,nsg
        do i=1,nsg
           write (6,*) i,j,gb(i+(j-1)*nsg),'gb'
        enddo
        enddo

        write (6,*) ' '

        do j=1,nsg
        do i=1,nsg
           write (6,*) i,j,ga(i+(j-1)*nsg),'ga'
        enddo
        enddo

        write (6,*) ' '

        do k=1,nsg
        do j=1,nsg
        do i=1,nsg
           write (6,*) i,j,k,gc(i+(j-1)*nsg+(k-1)*nsg*nsg),'gc'
        enddo
        enddo
        enddo

        return
        end
c-----------------------------------------------------------------------