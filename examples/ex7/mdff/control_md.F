&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        lreduced      = .false.
        lrestart      = .false.
        lbmlj         = .true.
        lvnlist       = .true.
&end
&mdtag
        integrator    = 'nve-vv'
        nequil        = 5000
        nequil_period = 10
        npas          = 50000
        nprint        = 100
        fprint        = 100
        dt            = 0.003
        temp          = 86.56d0
        ltraj         = .true.
        itraj_start   = 1
        itraj_period  = 10
&end
&fieldtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 119.8d0
        sigmalj(1,1)= 3.405d0
&end
