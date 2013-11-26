&controltag
 calc      = 'vib+band'
 cutshortrange    = 3.0d0
 lbmlj     = .true.
 lreduced  = .false.
&end

&vibtag
  nconf        = 1
  resdos       = 2.0d0    
  omegamax     = 30.0d0
  nkphon       = 100
  path='X-G'
  ks = 0.0d0 0.0d0 1.0d0,
  kf = 0.0d0 0.0d0 0.0d0,
  ncell=4
 
&end

&mdtag
&end

&proptag 
&end
&fieldtag
&end
