&controltag
 calc='vib'
 lbmlj     = .true.
 lreduced  = .false.
&end

&vibtag
  lwrite_vectff = .true.
  ncvib        = 1
  resdos       = 1.0d0    
  omegamax     = 30.0d0
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
