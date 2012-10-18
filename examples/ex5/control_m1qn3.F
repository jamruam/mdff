&controltag
calc='opt'
cutoff =3.0d0
lbmlj     = .true.
lreduced  = .false.
&end

&opttag
ncopt        = 1
nmaxopt      = 1
!optalgo ='sastry' 
!optalgo ='lbfgs'
optalgo ='m1qn3'
epsrel_m1qn3 = 1.0e-5  ! default is 1.0e-6 some problems occurs in this example
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
