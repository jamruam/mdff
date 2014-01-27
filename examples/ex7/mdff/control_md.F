&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        lreduced      = .false.
        lrestart      = .false.
        lbmlj         = .true.
        lvnlist       = .true.
        ltraj         = .true.
        itraj_start   = 10000
        itraj_period  = 10
&end
&mdtag
        integrator    = 'nve-vv'
        nequil        = 10000 
        nequil_period = 10
        npas          = 20000
        nprint        = 100
        fprint        = 1
        dt            = 0.003
        temp          = 86.56d0
&end
&fieldtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 119.8d0
        sigmalj(1,1)= 3.405d0
        mass=1.0
&end
