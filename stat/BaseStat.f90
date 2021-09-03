

 Module BaseStatistics

   Implicit None
   Interface Mean
     module procedure MeanS,MeanV
   End Interface Mean
   Interface Variance
     module procedure VarianceS, VarianceV,VarianceM
   End Interface Variance

   Contains
     !!!
     !!! compute average of scalar array
     Function MeanS( x )
       
       Implicit None
       Real*8,dimension( : )  ::x
       Real*8 ::MeanS
       MeanS = 0d0
       MeanS = Sum( x )
       MeanS = MeanS / Size( x , 1 )
     End Function MeanS

     !!
     !! computes the averages of 
     !! a matrix over second index
     Function MeanV( x )
       
       Implicit None
       Real*8,dimension( : , : )  ::x
       Real*8,dimension( 1:Size( x , 1 ) ) ::MeanV

       MeanV  =  0d0
       MeanV  =  Sum( x , dim = 2 )
       MeanV  =  MeanV / Dble( Size( x , 2 ) )

     End Function MeanV

     !!
     !! computes the average of 
     !! a matrix
     Function MeanM( x )
       
       Implicit None
       Real*8,dimension( : , : )  ::x
       Real*8  ::MeanM

       MeanM = 0d0
       MeanM = Sum( x )
       MeanM = MeanM / Dble( Size( x , 2 ) *  Size( x , 1 ) )

     End Function MeanM

     !!
     !! computes the variance of
     !! of a vector x;
     !! the average value can be obtained if 
     !! a variable is supplied
     !!

     Subroutine VarianceS( x , var , mean )

       Implicit None
       Real*8,dimension( : ) ::x
       Real*8  ::var
       Real*8,optional ::mean

       Real*8  ::x2
       Real*8  ::xbar


       Integer ::i

       x2  =  0d0
       xbar =  0d0

       Do i = 1, Size( x  , 1 )
         xbar  =  xbar + x( i )
         x2    =  x2   + x( i ) * x( i )
       End Do

       xbar = xbar / Dble( Size( x , 1 ) )
       x2   = x2 / Dble( Size( x , 1 ) )

       var = x2 - xbar*xbar

       if ( present( mean ) ) mean = xbar

     End Subroutine VarianceS

     !!
     !! the variance of a matrix is computed 
     !! over the second index optionally the 
     !! average value for every column of the matrix is computed
     !!

     Subroutine VarianceV( x , var , mean )

       Implicit None
       Real*8,dimension( : , : ) ::x
       Real*8,dimension( : )     ::var
       Real*8,dimension( : ),optional ::mean

       Real*8,dimension( 1:Size( x , 1 ) )  ::x2
       Real*8,dimension( 1:Size( x , 1 ) ) ::xbar


       Integer ::i

       x2  =  0d0
       xbar =  0d0

       Do i = 1, Size( x  , 1 )
         xbar  =  xbar + x( : , i )
         x2    =  x2   + x( : , i ) * x( : , i )
       End Do

       xbar = xbar / Dble( Size( x , 1 ) )
       x2   = x2 / Dble( Size( x , 1 ) )

       var = x2 - xbar*xbar

       if ( present( mean ) ) mean = xbar

     End Subroutine VarianceV


     !!
     !! compute variance over a 
     !! matrix; mean value is optional
     !!

     Subroutine VarianceM( x , var , mean )

       Implicit None
       Real*8,dimension( : , : ) ::x
       Real*8  ::var
       Real*8,optional ::mean

       Real*8  ::x2
       Real*8  ::xbar


       Integer ::i
       Integer ::j

       xbar =  0d0
       x2   =  0d0

       Do j = 1, Size( x , 2 )
         Do i = 1, Size( x  , 1 )
           xbar  =  xbar + x( i , j )
           x2    =  x2   + x( i , j ) * x( i , j )
         End Do
       End Do

       xbar = xbar / Dble( Size( x , 1 ) ) / Dble( Size( x , 2 ) )
       x2   = x2 / Dble( Size( x , 1 ) ) / Dble( Size( x , 2 ) )

       var = x2 - xbar*xbar

       if ( present( mean ) ) mean = xbar

     End Subroutine VarianceM



     Subroutine StudentTDist( grid , DOF , Dist )

       Implicit None
       Real*8,dimension( : )  ::grid
       Integer ::DOF
       Real*8,dimension( 1:Size( Grid ) )  ::Dist

       Real*8 ::factor

       factor = DGamma( ( Dble( DOF ) + 1d0 ) / 2d0 ) / &
                  Dsqrt( Dble( DOF ) ) / DGamma( Dble( DOF ) / 2d0 )

       Dist = factor * ( 1d0 + grid*grid / Dble( DOF ) )**(-(Dble( DOF )+1d0)/2d0 )

     End Subroutine StudentTDist








 End Module BaseStatistics
