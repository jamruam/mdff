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
! GAMMA       ->       K 
  ksx=0.0d0,       kfx=1.0d0,
  ksy=0.0d0,       kfy=1.0d0, 
  ksz=0.0d0,       kfz=0.0d0
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
