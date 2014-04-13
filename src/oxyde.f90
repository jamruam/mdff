MODULE oxyde

  USE constants,        ONLY :  dp 
! pour le moment il n'y a que 6 oxydes 
! mais l'ajout de nouveaux oxydes ne devraient pas poser de probleme.
! à modifier (dans le cas d'ajout d'un ou plusieurs oxydes) :
!        - le nombre d'oxyde = "noxyde"
!        - creer une variable reel pour etre lu en entré : real(kind=dp) :: sio2 , na2o , b2o3 ...
!        - ajouter le nouvel element à la liste des elements "ele"
!        - ajouter le degree d'oxydation "numoxyd"
!        - ajouter les infos pour les oxydes aux listes lox, ele_ox et nel_ox
!        - attention lox, ele_ox et nel_ox sont definies dans le meme ordre

! note:
! il y aurai certainement un moyen de reduire le nombre d'info 
! necessaire... juste à partir de la liste des elements et du numero
! d'oxydation par exemple.

  integer , PARAMETER :: noxyde = 6
  integer , PARAMETER :: nelem  = noxyde + 1 ! + oxygen
  real(kind=dp) :: g_to_am = 1.660538782_dp ! grammes -> unité de masse atomique

  character(len=4)    :: lox     ( noxyde )       ! label de l'oxyde ex: SiO2, Na2O
  character(len=2)    :: ele_ox  ( 2 , noxyde )   ! element de l'oxyde ex:Si , Na . le second est toujours l'oxygene
  integer             :: nel_ox  ( 2 , noxyde )   ! nombre d'element dans
  real(kind=dp)       :: sio2 , na2o , b2o3 , cao , p2o5 , al2o3 ! valeur en entré

  character(len=2)    :: ele     ( nelem ) ! liste des elements chimique distincts
  integer             :: numoxyd ( nelem ) ! degree d'oxydation
  integer             :: valence ( nelem ) ! valence electrons ( as recommended in VASP )
  real(kind=dp)       :: massele ( nelem ) ! Standard atomic weight 

  data ele      / 'Si' , &
                  'Na' , &
                  ' B' , &
                  'Ca' , &
                  ' P' , &
                  'Al' , &
                  ' O' /

  data numoxyd  /   4  , &
                    1  , & 
                    3  , &
                    2  , & 
                    5  , & 
                    3  , & 
                   -2  /

  data valence  /   4  , &
                    7  , & 
                    3  , &
                    8  , & 
                    5  , & 
                    3  , & 
                    6  /

  data massele  / 28.0855_dp     , &
                  22.98976928_dp , &
                  10.811_dp      , &
                  40.078_dp      , & 
                  30.973762_dp   , & 
                  26.9815386_dp  , & 
                  15.9994_dp      /

  data lox    / 'SiO2' , 'Na2O' , 'B2O3' , 'CaO' , 'P2O5' , 'Al2O3' /

  data ele_ox / 'Si' , ' O' , &
                'Na' , ' O' , & 
                ' B' , ' O' , & 
                'Ca' , ' O' , & 
                ' P' , ' O' , & 
                'Al' , ' O' / 

  data nel_ox /  1 , 2 , & ! Si 1 O 2
                 2 , 1 , & ! Na 2 O 1
                 2 , 3 , & !  B 2 O 3 
                 1 , 1 , & ! Ca 1 O 1 
                 2 , 5 , & !  P 2 O 5
                 2 , 3   / ! Al 2 O 3


END MODULE oxyde
