&controltag
 calc='opt'
 cutoff = 3.0d0
 lbmlj     = .true.
 lreduced  = .false.
&end

&opttag
ncopt        = 1
nmaxopt      = 1
!optalgo ='sastry' 
optalgo ='lbfgs'
!optalgo ='m1qn3'

&end

&mdtag
&end

&configtag
lgenconf=.false.
&end
&proptag 
&end
&fieldtag
&end
