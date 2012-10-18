&controltag
 calc='vib+dos'
 lbmlj     = .true.
 lreduced  = .false.
&end

&vibtag
  ncvib        = 1
  resdos       = 0.50d0    
  omegamax     = 30.0d0
  nkphon       = 40 ! (nkphon+1) * (nkphon+1) * (nkphon+1) generated in IBZKPTFF
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
