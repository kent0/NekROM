c-----------------------------------------------------------------------
        program offline

        include 'LVAR'
        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        call blank(fname,132)
        fname='/Users/kaneko/Developer/MOR/offline/r0.f00001 '

        nel=512

        ieg(1)=1
        ieg(2)=nel

        call setindxr(indxr,ldimt)
        write (6,*) 'buf(1)',buf(1)
        call rxupt(buf,ieg,indxr,fname)
        write (6,*) 'buf(1)',buf(1)
        call setb(gb,buf,buf,1,512*lxyz)

        write (6,*) 'gb(1)',gb(1)

        end
c-----------------------------------------------------------------------