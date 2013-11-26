&controltag
        calc     = 'vib+dos'
        lbmlj    = .true.
        lreduced = .false.
&end

&vibtag
        nconf        = 1
        resdos       = 0.1d0    
        omegamax     = 30.0d0
        nkphon       = 50 ! (nkphon+1) * (nkphon+1) * (nkphon+1) generated in IBZKPTFF
        ncell        = 4
&end

&mdtag
&end

&fieldtag
&end
