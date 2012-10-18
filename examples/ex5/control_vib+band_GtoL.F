&controltag
 calc      = 'vib+band'
 cutoff    = 3.0d0
 lbmlj     = .true.
 lreduced  = .false.
&end

&vibtag
  ncvib        = 1
  resdos       = 2.0d0    
  omegamax     = 30.0d0
  nkphon       = 100
! GAMMA      ->     X
  ksx=0.0d0,    kfx=0.5d0,
  ksy=0.0d0,    kfy=0.5d0,
  ksz=0.0d0,    kfz=0.5d0
 
&end

&mdtag
&end

&configtag
lgenconf=.false.
ncell=4
&end
&proptag 
&end
&fieldtag
&end
