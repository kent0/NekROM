c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        nsg=64
        neg=512
        nel=512
        nep=512

        call loadflist(fnames,nsg,nsl)

        call setindxr(indxr,ldimt)

        do is=1,3
           write (6,'(a132)') fnames(is)
        enddo

        ns=1

        call rzero(gb,ns*ns)

        iloc=1

        is0=1
        is1=2

        ng=neg/nep
        ng=1
        mel=1
        nel=mel

        nsg=lsg
        ieg(1)=1
        ieg(2)=1
        call loadsnaps(buf,ieg,indxr,nsg)

        do i=1,lxyz*4
            write (6,*) is,i,ibuf(i),ibuf8(i),'ibuf'
        enddo

        return
        end
c-----------------------------------------------------------------------