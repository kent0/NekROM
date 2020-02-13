c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        call blank(fname,132)
        fname='/Users/kaneko/Developer/MOR/offline/r0.f00001 '
        call chcopy(fnames(1),fname,132)

        call blank(fname,132)
        fname='/Users/kaneko/Developer/MOR/offline/r0.f00002 '
        call chcopy(fnames(2),fname,132)

        call setindxr(indxr,ldimt)

        ns=1

        mel=16
        ng=32

        mel=512
        ng=1

        call rzero(gb,ns*ns)

        do ig=1,ng
           ieg(1)=mel*(ig-1)+1
           ieg(2)=mel*ig
           call rxupt(buf,ieg,indxr,fnames(1))
           call setb(gb,buf,buf,ns,mel*lxyz*ldim)
c          do k=1,mel
c          do j=1,ldim
c          do i=1,lxyz
c          write (6,*) ig,i,j,k,buf(i+(j-1)*lxyz+(k-1)*lxyz*ldim),
c    $             'wp buf'
c          enddo
c          enddo
c          enddo
           write (6,*) ig,gb(1),buf(1),'wp deb'
        enddo


        end
c-----------------------------------------------------------------------