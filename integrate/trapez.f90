  Module trapez_int
   Interface
     Function trapezoid( dx , func )
       Implicit None
       Real*8                ::dx
       Real*8                ::trapezoid
       Real*8,dimension( : ) ::func
       Integer               ::i
     End Function trapezoid
   End Interface


  End Module trapez_int


  Function trapezoid( dx , func )
    Implicit None
    Real*8                ::dx
    Real*8                ::trapezoid
    Real*8,dimension( : ) ::func
    Integer               ::i

    trapezoid  =  0.5d0 * func( 1 )
    trapezoid  =  trapezoid + 0.5d0 * func( Size( func , 1 ) )

    Do i = 2 , Size( func , 1 ) - 1
       trapezoid  =  trapezoid + func( i )
    End Do

    trapezoid  =  trapezoid  *  dx
  
  End Function



  Module trapez_int_2d

     Implicit None

     Contains

     Function trapezoid_2d( func , dx , dy )

       use trapez_int, only:trapezoid
       Implicit None
       Real*8,dimension( : , : ) ::func
       Real*8 ::dx
       Real*8 ::dy
       Real*8 ::trapezoid_2d
       Real*8,dimension( 1 : Size( func , 2 ) ) ::temp
       Integer ::i

       Do i = 1 , Size( func , 2 )
         temp( i )  =  trapezoid( dx , func( : , i ) )
       End Do

       trapezoid_2d = trapezoid( dy , temp )
     End Function trapezoid_2d
  End Module trapez_int_2d


  Module trapez_int_3d

     Implicit None

     Contains

     Function trapezoid_3d( func , dx , dy , dz )

       use trapez_int, only:trapezoid
       use trapez_int_2d, only:trapezoid_2d
       Implicit None
       Real*8,dimension( : , : , : ) ::func
       Real*8 ::dx
       Real*8 ::dy
       Real*8 ::dz
       Real*8 ::trapezoid_3d
       Real*8,dimension( 1 : Size( func , 3 ) ) ::temp
       Integer ::i

       Do i = 1 , Size( func , 3 )
         temp( i )  =  trapezoid_2d( func( : , : , i ) , dx , dy )
       End Do

       trapezoid_3d = trapezoid( dz , temp )
     End Function trapezoid_3d
  End Module trapez_int_3d




