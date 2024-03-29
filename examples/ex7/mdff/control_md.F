&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        restart_data  = 'rnn'
        lreduced      = .false.
        lrestart      = .false.
        lnmlj         = .true.
        lvnlist       = .true.
        ltraj         = .true.
        itraj_start   = 20000
        itraj_period  = 10
&end
&mdtag
        npas          = 40000
        nequil        = 20000
        nequil_period = 10
        nprint        = 100
        fprint        = 1
        
        dt            = 0.001 
        temp          = 86.56d0 
&end
&fieldtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 0.010323576 
        sigmalj(1,1)= 3.405d0
        mass=39.95d0
&end
