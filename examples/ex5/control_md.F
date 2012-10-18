&controltag
 calc      = 'md'
 cutoff    = 3.0d0
 lbmlj     = .true.
 lcoulomb  = .false.
 lreduced  = .false.
&end

&mdtag
npas=1
ltraj=.true.
itraj_start=0
itraj_period=1
&end

&configtag
lgenconf=.true.
lfcc=.true.
ntype=1
ncell=4
rho=1.0d0
&end
&proptag 
&end
&fieldtag
&end
