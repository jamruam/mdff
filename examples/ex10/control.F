&controltag
        calc          = 'md'
        cutshortrange = 7.0d0
        cutlongrange  = 7.0d0
        skindiff      = 1.0d0
        !lvnlist       = .false.
        lreduced      = .false.
        lbmhftd       = .false.
        lcoulomb      = .true.
        lstatic       = .true.
        longrange     = 'ewald'
        ltraj         =.true.
        itraj_start   = 0
        itraj_period  = 1
&end
&mdtag
        temp          = 0.0d0
&end
&fieldtag
        ctrunc = 'notrunc'
        qch=-2.0d0 4.0d0 3.0d0
        mass=15.99994d0 28.0855d0 10.811d0

        lautoES = .true.
        epsw = 1.e-10

&end
