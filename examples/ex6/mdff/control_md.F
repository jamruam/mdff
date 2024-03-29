&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        restart_data  = 'rnn'
        lnmlj         = .true.
        lvnlist       = .false.
        ltraj         = .true.
        itraj_start   = 0 
        itraj_period  = 5
&end
&mdtag
        npas          = 40000 
        nprint        = 100
        fprint        = 100

        integrator    = 'nvt-nhcn'
        dt            = 0.006
        temp          = 119.8d0
        nhc_n         = 2 
        nhc_yosh_order= 1 
        nhc_mults     = 1
        timesca_thermo= 2.0 
&end
&fieldtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 0.010323576 
        sigmalj(1,1)= 3.405d0
        mass=39.95d0
&end
