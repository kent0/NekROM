c-----------------------------------------------------------------------
        program offline

        include 'OFFLINE'

        character*132 fname

        call offline_init(icomm)

        call blank(fname,132)
        fname='/Users/kaneko/Developer/MOR/Offline/r0.f00001 '

        ieg(1)=1
        ieg(2)=2

        call setifread(ifread,ldimt)
        write (6,*) 'upt(1)',upt(1)
        call rxupt(xyz,upt,ieg,ifread,fname)
        write (6,*) 'upt(1)',upt(1)

        end
c-----------------------------------------------------------------------