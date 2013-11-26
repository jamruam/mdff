&controltag
        calc    = 'vib'
        lbmlj   = .true.
        lreduced= .false.
&end

&vibtag
        lwrite_vectff = .true.
        nconf        = 1
        resdos       = 1.0d0    
        omegamax     = 30.0d0
        ncell        = 4
&end

&mdtag
&end

&fieldtag
&end
