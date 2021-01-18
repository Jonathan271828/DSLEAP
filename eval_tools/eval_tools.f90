


  
  Module FlagReader

    Implicit None
    Interface ReadFlags
    module procedure flag_finder,value_finder
    End Interface ReadFlags
    Contains
    Subroutine flag_finder( tag , flag , fnameInput )
   
       Implicit None
       Character( len = * )          ::tag
       Logical                       ::flag
       Character( len = * ),optional ::fnameInput
       Character( len = : ),allocatable ::fname
       Character ( len = 70 ) ::inline
       Integer                ::ioerror
       Logical                ::check
       Integer                ::j

       fname  =  "eval_input"
       If ( present( fnameInput ) ) then
           fname = fnameInput
       Else
           Write( * , "(A)" ) 'Taking default file name'
           Write( * , "(A)" ) fname
       End If

   
       ioerror = 0
       flag = .False.
      open( unit = 25 , action = 'Read' , file = Trim( fname ) , status='old' )
       check = .False.
       RL: Do
         Read( 25 , "( A )" , iostat = ioerror ) inline
         If ( ioerror .ne. 0 ) Exit RL
         If ( Trim( inline( 1 : len( tag ) ) ) .eq. Trim( tag )  ) then
            Do j = 1 , len( inline )
              If ( Trim( inline( j : j ) ) .eq. "=" ) then
                  check = .True.
              End If
              If ( check ) then
                If ( Trim( inline( j : j ) ) .eq. "T" .or. Trim( inline( j : j ) ) .eq. "t" ) then
                  flag  = .True.
                  Exit RL
                End If
                If ( Trim( inline( j : j ) ) .eq. "F" .or. Trim( inline( j : j ) ) .eq. "f" ) then
                  flag  = .False.
                  Exit RL
                End If
              End If
            End Do
         End If
       End Do RL
      close( 25 )
   
    End Subroutine
  
 
    Subroutine value_finder( tag , val , check_found , fnameInput )
   
     Implicit None
     Character( len = * )   ::tag
     Character ( len = 70 ) ::inline
     Character( len = * ),optional ::fnameInput
     Character( len = : ),allocatable ::fname
     Integer                ::ioerror
     Logical                ::check
     Integer                ::j
     Real*8                 ::val
     Character( len = 10 )  ::string_value
     Logical                ::check_found
     Integer                ::kk
   
       fname  =  "eval_input"
       If ( present( fnameInput ) ) then
           fname = fnameInput
       Else
           Write( * , "(A)" ) 'Taking default file name'
           Write( * , "(A)" ) fname
       End If

   
      open( unit = 25 , action = 'Read' , file = Trim( fname ) , status='old' )
       check = .False.
       string_value  =  ""
       check_found = .False.
       RL: Do
         Read( 25 , "( A )" , iostat = ioerror ) inline
         If ( ioerror .ne. 0 ) Exit RL
         If ( Trim( inline( 1 : len( tag ) ) ) .eq. Trim( tag )  ) then
            Do j = 1 , len( inline )
              If ( Trim( inline( j : j ) ) .eq. "=" ) then
                  check = .True.
                  kk = j
              End If
              If ( check .and. kk .lt. j ) then
                 If ( Trim( inline( j : j ) ) .ne. " " ) then
                   string_value  =  Trim( string_value ) // inline( j : j )
                   check_found  =  .True.
                 End If
              End If
            End Do
            Exit RL
         End If
       End Do RL
      close( 25 )
   
      If ( check_found ) then
        Read( string_value , * ) val
      End If
   
   
    End Subroutine value_finder
  End Module FlagReader





  Module Binning

    Implicit None
    Interface get_bin
       module procedure bin,compute_binLB
    End Interface get_bin
    Contains

    Function bin( theta , dtheta , Ntheta )
      implicit none
      Real*8            ::theta
      Real*8            ::dtheta
      Integer           ::Ntheta
      Integer           ::bin

      bin = Ceiling( theta / dtheta )
      If ( bin .lt. 1 )       bin  =  1
      If ( bin .gt. Ntheta )  bin  =  Ntheta
    End Function bin

    Function compute_binLB( value , dx , N , lower ) Result( bin )

      Implicit None
      Real*8  ::value
      Real*8  ::dx
      Integer ::N
      Real*8  ::lower
      Integer ::bin


      bin = Ceiling( value / dx ) - Ceiling( lower / dx )
      If ( bin .lt. 1 ) then
        bin = 1
      ElseIf ( bin .gt. N ) then
        bin = N
      End If
    End Function compute_binLB

  End Module Binning



   Module tool_interfaces

     Interface
        Function get_norm( in_vector )
          Real*8,dimension( : )  ::in_vector
          Real*8                 ::get_norm
        End Function get_norm
     End Interface

     Interface
        Function get_periodic_images_and_link( data ) Result( images )
          Logical     ::link_list
          Real*8,dimension( : , : )  ::data
          Integer ::i
          Integer ::aa
          Integer ::bb
          Integer ::cc
          !Real*8,dimension( 1 : 27 * Size( data , 1 ) , 1 : Size( data , 2 ) + 1 )  ::images
          Real*8,dimension( 1 : Size( data , 1 ) + 1 , 1 : 27 * Size( data , 2 ) )  ::images
        End Function get_periodic_images_and_link
     End Interface

     Interface
        Function get_periodic_images_and_link_cart( data , lattice ) Result( images )
          Logical     ::link_list
          Real*8,dimension( : , : )  ::data
          Real*8,dimension( : , : )  ::lattice
          Integer ::i
          Integer ::aa
          Integer ::bb
          Integer ::cc
          !Real*8,dimension( 1 : 27 * Size( data , 1 ) , 1 : Size( data , 2 ) + 1 )  ::images
          Real*8,dimension( 1 : Size( data , 1 ) + 1 , 1 : 27 * Size( data , 2 ) )  ::images
        End Function get_periodic_images_and_link_cart
     End Interface

     Interface
        Function get_periodic_images( data ) Result( images )
          Logical     ::link_list
          Real*8,dimension( : , : )  ::data
          Integer ::i
          Integer ::aa
          Integer ::bb
          Integer ::cc
          Real*8,dimension( 1 : 3 ,  27 * Size( data , 2 ) )  ::images
        End Function get_periodic_images
     End Interface

     Interface
        Function get_nearest_image( vectorA , vectorB )  Result(  vectorC )
          Real*8,dimension( : )  ::vectorA
          Real*8,dimension( : )  ::vectorB
          Real*8,dimension( 1 : Size( vectorA , 1 ) ) ::vectorC
        End Function get_nearest_image
     End Interface

     Interface
        Function calc_spherical_theta( vector )
          Real*8,dimension( : )      ::vector
          Real*8                     ::calc_spherical_theta
        End Function calc_spherical_theta
     End Interface

     Interface
        Function calc_spherical_phi( vector )
          Real*8,dimension( : )      ::vector
          Real*8                     ::calc_spherical_phi
        End Function calc_spherical_phi
     End Interface

     Interface
       Function calc_pearson_terms( xnew , ynew ) Result( data )
         Real*8         ::xnew
         Real*8         ::ynew
         Real*8,dimension( 1 : 5 ) ::data
       End Function calc_pearson_terms
     End interface

     Interface
       Function mean_value_IIxII_array( in_array , N ) Result( out_array )
         Integer   ::N
         Real*8,dimension( : , : )  ::in_array
         Real*8,dimension( 1 : Size( in_array , 1 ) , 1 : Size( in_array , 2 ) ) &
                         ::out_array
       End Function mean_value_IIxII_array
     End Interface

     Interface
       Function calc_single_pearson_from_terms( in_data , N )
         Integer       ::N
         Real*8,dimension( 1 : 5 )  ::in_data
         Real*8                     ::calc_single_pearson_from_terms
       End Function calc_single_pearson_from_terms
     End Interface

     Interface
       Function periodic_dist( vector_in ) Result( vector_out )
         Real*8,dimension( : )  ::vector_in
         Real*8,dimension( 1 : Size( vector_in , 1 ) ) ::vector_out
       End Function periodic_dist
     End interface

     Interface
       Function bin( theta , dtheta , Ntheta )
         Real*8            ::theta
         Real*8            ::dtheta
         Integer           ::Ntheta
         Integer           ::bin
       End Function
     End Interface

     Interface
       Function shift_atom_to_prim_cell( atom ) Result( atom_in_cell )
         Implicit None
         Real*8,dimension( : )  ::atom
         Real*8,dimension( 1 : Size( atom , 1 ) )  ::atom_in_cell
       End Function shift_atom_to_prim_cell
     End Interface

     Interface
       Subroutine write_poscar_format( cp , latc , species , Natoms , coord_type , coords )
         Implicit None
         Real*8    ::cp
         Real*8,dimension( : , : )  ::latc
         Character( len = * )   ::species
         Integer,dimension( : ) ::Natoms
         Character( len = * )   ::coord_type
         Real*8,dimension( : , : ) ::coords
       End Subroutine write_poscar_format
     End Interface

     Interface
       Function get_N_ats_unit_cell( Ntot , unit_size )
         Implicit None
         Integer                 ::Ntot
         Real*8,dimension( : )   ::unit_size
         Integer                 ::get_N_ats_unit_cell
       End Function get_N_ats_unit_cell
     End interface

     Interface
       Function get_nearest_neighbour_no_periodicity( atoms , center , N , Nself ) Result( list )
         Implicit None
         Real*8,dimension( : , : )  ::atoms
         Real*8,dimension( : )      ::center
         Integer                    ::N
         Integer,dimension( 1 : N ) ::list
         Integer  ::Nself
       End Function get_nearest_neighbour_no_periodicity
     End Interface

     Interface
       Function get_cart_from_spher( vector ) Result( cart_vec )
         Implicit None
         Real*8,dimension( : )  ::vector
         Real*8,dimension( 1 : 3 ) ::cart_vec
       End Function get_cart_from_spher
     End interface

     Interface
       Function spherical_sin_correction( in_data , dimen , u_limit ) Result( out_data )
         Implicit None
         Integer    ::dimen
         Real*8,dimension( : , : )  ::in_data
         Real*8,dimension( 1 : Size( in_data , 1 ) , &
                           1 : Size( in_data , 2 ) )  ::out_data
         Real*8                     ::u_limit
       End function spherical_sin_correction
     End interface
     Interface
       Subroutine gen_file_number( Num , file_string )
         Implicit None
         Character( len = * ) ::file_string
         Integer  ::Num
       End Subroutine gen_file_number
     End Interface

     Interface
       Function histogram_2d( data , min1 , min2 , max1 , max2 , N1 , N2 ) Result( hist )
         Implicit None
         Real*8,dimension( : , : )  ::data
         Real*8     ::min1
         Real*8     ::min2
         Real*8     ::max1
         Real*8     ::max2
         Integer    ::N1
         Integer    ::N2
         Integer,dimension( 1 : N1 , 1 : N2 ) ::hist
       End Function histogram_2d
     End Interface

     Interface
       Function Cross_Product( a , b ) Result( c )
         Implicit None
         Real*8,dimension( : ) ::a
         Real*8,dimension( : ) ::b
         Real*8,dimension( 1 : Size( a , 1 ) ) ::c
       End Function Cross_Product
     End interface

     Interface
         Function get_2d_bin( theta , phi , dtheta , dphi , Ntheta , Nphi ) Result( bins )
           Implicit None
           Real*8   ::theta
           Real*8   ::phi
           Real*8   ::dtheta
           Real*8   ::dphi
           Integer  ::Ntheta
           Integer  ::Nphi
           Integer,dimension( 1 : 2 )  ::bins
         End Function get_2d_bin
     End Interface
     Interface
         Function give_closest_lower_power_of_2( N1 ) Result( N2 )
           Implicit None
           Integer ::N1
           Integer ::N2
         End Function give_closest_lower_power_of_2
     End Interface

     Interface
         Function get_nearest_image_cartesian( vectorA , vectorB , lattice ) Result( vectorC )
           Implicit None
           Real*8,dimension( : )  ::vectorA
           Real*8,dimension( : )  ::vectorB
           Real*8,dimension( : , : )  ::lattice
           Real*8,dimension( 1 : Size( vectorA , 1 ) ) ::vectorC
         End function get_nearest_image_cartesian
     End Interface

     Interface
         Function assign_closest_atoms_to_array( atomsA , atomsB ) Result( ClosestAtoms )
           Real*8,dimension( : , : )  ::atomsA
           Real*8,dimension( : , : )  ::atomsB
           Real*8,dimension( 1:3 , 1:Size( atomsA , 2 ) )  ::ClosestAtoms
           Real*8,dimension( 1:3 , 1:27*Size( atomsA , 2 ) )  ::images
           Real*8   ::Sdist  !Shortest dist
           Real*8   ::Adist  !Actual dist
           Integer ::i
           Integer ::j
         End Function assign_closest_atoms_to_array
     End Interface

     Interface
         Subroutine Drift_correction( atoms )
           Implicit None
           Real*8,dimension( : , : ) ::atoms
           Real*8,Save,dimension( 1:3 )  ::center = [0.5d0, 0.5d0, 0.5d0 ]
           Real*8,dimension( 1:3 ) :: COM
           Integer ::i
         End Subroutine Drift_correction
     End Interface
     Interface
       Subroutine transformation_matrix( latVa , latVb , tm )
         Implicit None
         Real*8,dimension( : , : )  ::latVa
         Real*8,dimension( : , : )  ::latVb
         Real*8,dimension( : , : )  ::tm
         Integer                    ::i
         Integer                    ::j
       End Subroutine transformation_matrix
     End Interface
     Interface
       Subroutine coordinate_transformation( A , x , xprime )
         implicit none
         Real*8,dimension( : , : )  ::A
         Real*8,dimension( : , : )  ::x
         Real*8,dimension( : , : )  ::xprime
         Integer                    ::i
       End Subroutine coordinate_transformation
     End Interface
     Interface
       Subroutine coordinate_transformationIP( A , x )
         implicit none
         Real*8,intent( in ),dimension( : , : )  ::A
         Real*8,intent( inout ),dimension( : , : )  ::x
         Integer                    ::i
       End Subroutine coordinate_transformationIP
     End Interface
     Interface
       Subroutine lattice_consts( a , b , c , alpha , beta , gamma , bravais )
         Implicit None
         Real*8                    ::a,b,c
         Real*8                    ::alpha,beta,gamma
         Real*8,dimension( 3 , 3 ) ::bravais
       End Subroutine lattice_consts
     End Interface
     Interface
       Function make_lattice( a , b, c ) Result( lattice )
         Implicit None
         Real*8   ::a
         Real*8   ::b
         Real*8   ::c
         Real*8,dimension( 1 : 3 , 1 : 3 )  ::lattice
       End Function
     End Interface
     Interface
       Function lattice_volume( Bravais ) Result( volume ) !j.L 20.07.2018
           Implicit None
           Real*8,dimension( : , : )  ::Bravais
           Real*8  ::volume
       End Function lattice_volume
     End interface

     Interface
       Function compute_angle( vecA , vecB ) Result( angle )
         Implicit None
         Real*8,dimension( : )  ::vecA
         Real*8,dimension( : )  ::vecB
         Real*8   ::angle
       End Function compute_angle
     End Interface

     Interface
       Function count_words( fname , word ) Result( Nwords )

         Implicit None
         Character( len = * )  ::fname
         Character( len = * )  ::word
         Integer   ::Nwords
       End Function
     End Interface

     Interface
       Function PeriodicImageCartS( Pos , lattice , mode ) Result( PPositions ) !J.L 23.01.2019
         Implicit None
         Real*8,dimension( : )   ::Pos
         Real*8,dimension( : , : )  ::lattice
         Integer ::mode
         Real*8,allocatable,dimension( : , : ) ::PPositions   !!Periodic positions
       End Function PeriodicImageCartS
     End Interface

     Interface
       Function ComputeQgrid( Basis , Nx , Ny , Nz ) Result( Grid )
         Implicit None
         Real*8,dimension( : , : )  ::Basis
         Integer ::Nx
         Integer ::Ny
         Integer ::Nz
         Real*8,allocatable,dimension( : , : )  ::Grid
       End Function ComputeQgrid
     End Interface

     Interface
       Function MatMulArray( Array , matrix ) Result( Prod )
         Implicit None
         Real*8,dimension( : , : )  ::Array
         Real*8,dimension( : , : )  ::matrix
         Real*8,allocatable,dimension( : , : )  ::Prod
       End Function MatMulArray
     End Interface


     Interface
       Function ComputeRotationMatrix( A , B ) Result( RotMat )
         Implicit None
         Real*8,dimension( : ) ::A
         Real*8,dimension( : ) ::B
         Real*8,dimension( 1:3 , 1:3 ) ::RotMat

       End Function ComputeRotationMatrix
     End Interface

     Interface
       Subroutine RotatePolarDist( Polar , dtheta , dphi , RotMat )
         Implicit None
         Real*8,dimension( : , : )  ::Polar
         Real*8 ::dtheta
         Real*8 ::dphi
         Real*8,dimension( : , : )  ::RotMat
       End Subroutine RotatePolarDist
     End Interface

     Interface
       Function MyFindLoc( x , val ) Result( indx )
          Implicit None
          Integer,dimension( : ) ::x
          Integer ::val
          Integer ::indx
       End Function MyFindLoc
     End Interface

     InterFace
       Subroutine WriteCHGCAR( lattice , atoms , fname , types , Ntypes , density3d )
          Implicit None
          Real*8,dimension( : , : )  ::lattice
          Real*8,dimension( : , : )  ::atoms
          Character( len=* ) ::fname
          Character( len=* ) ::types
          Integer,dimension( : ) ::Ntypes
          Real*8,dimension( : , : , : ) ::density3d
       End Subroutine WriteCHGCAR
     End Interface

     Interface 
       Function GetSelfCorrSteps( geostart , EndGeo , every ) Result( steps )
         Implicit None
         Integer ::geostart
         Integer ::endgeo
         Integer ::every
         Integer ::steps
         Integer ::start
         Integer ::delta
       End Function GetSelfCorrSteps
     End Interface

     Interface
       Subroutine ReplicateCell( atoms , Nx , Ny , Nz , Natoms )
         Implicit None
         Real*8,intent(inout),dimension( : , : )  ::atoms
         Integer,intent( in ) ::Nx
         Integer,intent( in ) ::Ny
         Integer,intent( in ) ::Nz
         Integer,intent( in ) ::Natoms
       End Subroutine ReplicateCell
     End Interface
   End Module tool_interfaces



   !!
   !! supply atoms in direct coordinates [0,1]
   !! routine will make Nx replicas of this cell
   !! in x direction, Ny in y direction, Nz
   !! in z direction 
   !! if Nx,Ny,Nz=1 original cell will be returned
   !! input array will be resized

   Subroutine ReplicateCell( atoms , Nx , Ny , Nz , Natoms )

     Implicit None
     Real*8,allocatable,intent(inout),dimension( : , : )  ::atoms
     Integer,intent( in ) ::Nx
     Integer,intent( in ) ::Ny
     Integer,intent( in ) ::Nz
     Integer,intent( in ) ::Natoms
     Real*8,allocatable,dimension( : , : )  ::TempAtoms
     Integer ::xx
     Integer ::yy
     Integer ::zz
     Integer ::i
     Integer ::col

     Allocate( TempAtoms( 1:Size( atoms , 1 ) , &
               1:Nx*Ny*Nz*Natoms ) )


     col = 1
     Do xx = 0 , Nx - 1
       Do yy = 0 , Ny - 1
         Do zz = 0 , Nz - 1
           Do i = 1, Natoms
              TempAtoms( 1 , col )  =  atoms( 1 , i ) + xx
              TempAtoms( 2 , col )  =  atoms( 2 , i ) + yy
              TempAtoms( 3 , col )  =  atoms( 3 , i ) + zz
              col = col + 1
           End Do
         End Do
       End Do
     End Do

     If ( Natoms .ne. Size( TempAtoms , 2 ) ) then
       Deallocate( atoms )
       Allocate( atoms( 1:Size( TempAtoms , 1 ) ,&
                      1:Size( TempAtoms , 2 ) ) )
     End If
     atoms = TempAtoms
     Deallocate( TempAtoms )
   End Subroutine ReplicateCell



   Function GetSelfCorrSteps( geostart , endgeo , every ) Result( steps )

     Implicit None
     Integer ::geostart
     Integer ::endgeo
     Integer ::every
     Integer ::steps
     Integer ::start
     Integer ::delta

     If ( geostart .eq. 0 ) then
       start = 1
     Else
       start = geostart
     End If
     
     delta = ( endgeo - start ) + 1
     
     
     If ( Modulo( start , every ) .eq. 0 .and. Modulo( delta , every ) .ne.0 ) then
         Steps = ( endgeo - start + 1 ) / every + 1
     Else
         Steps = ( endgeo - start + 1 ) / every
     End If
   End Function GetSelfCorrSteps




   Subroutine WriteCHGCAR( lattice , atoms , fname , types , Ntypes , density3d )


     Implicit None
     Real*8,dimension( : , : )  ::lattice
     Real*8,dimension( : , : )  ::atoms
     Character( len=* ) ::fname
     Character( len=* ) ::types
     Integer,dimension( : ) ::Ntypes
     Real*8,dimension( : , : , : ) ::density3d

     Integer ::i

     Integer ::Nx
     Integer ::Ny
     Integer ::Nz


     open( unit = 312 , action='write' , status='unknown' , file = Trim( fname ) )
       Write( 312 , "(A)" ) "3d Sitribution eval routine"
       Write( 312 , * ) 1.0
       Do i = 1, Size( lattice , 2 )
          Write( 312 , "(3F15.8)" ) lattice( : , i )
       End Do
       Write( 312 , "(A)" ) types
       Write( 312 , "(100I5)" ) Ntypes
       Write( 312 , "(A)" ) "Direct"
       Do i = 1, Size( atoms , 2 )
         Write( 312 , "(3F15.8)" ) atoms( : , i )
       End Do
       Write( 312 , * )
       Write( 312 , * ) Size(density3d , 1 ) , Size(density3d , 2 ) , Size( density3d , 3 )
       Write( 312 , '(1(1X,E17.11))') ( ( ( density3d( NX , NY , NZ ) , NX = 1 , Size( density3d , 1 ) ) ,&
                                           NY = 1 , Size( density3d , 2 ) ) ,&
                                           NZ = 1 , Size( density3d , 3 ) )
     close( 312 )
   End Subroutine WriteCHGCAR




   Function MyFindLoc( x , val ) Result( indx )

      Implicit None
      Integer,dimension( : ) ::x
      Integer ::val

      Integer ::indx
      Integer ::i

      indx = -Huge( indx )
      Do i = 1, Size( x , 1 )
        If ( x( i ) .eq. val ) then
           indx = i
           Return
        End If
      End Do

   End Function MyFindLoc



   Subroutine RotatePolarDist( Polar , dtheta , dphi , RotMat )

     use tool_interfaces , only:get_cart_from_spher,&
                             calc_spherical_theta,calc_spherical_phi,&
                             get_norm
     use Binning, only :get_bin
     Implicit None
     Real*8,dimension( : , : )  ::Polar
     Real*8 ::dtheta
     Real*8 ::dphi
     Real*8,dimension( : , : )  ::RotMat

     Real*8,allocatable,dimension( : , : )  ::Temp
     Integer ::N1
     Integer ::N2

     Integer ::i
     Integer ::j

     Real*8  ::theta
     Real*8  ::phi

     Real*8,dimension( 1:3 ) ::cart
     Real*8,dimension( 1:3 ) ::spher

     Integer ::bin1
     Integer ::bin2


     N1  =  Size( Polar , 1 )
     N2  =  Size( Polar , 2 )

     Allocate( Temp( 1:N1 , 1:N2 ) )

     Temp = Polar


     Do j = 1, N2
       Do i = 1 , N1
         phi    =  Dble( j ) * dphi - dphi / 2d0
         theta  =  Dble( i ) * dtheta - dtheta / 2d0
         spher( 1 ) = theta
         spher( 2 ) = phi
         spher( 3 ) = 1d0
         cart = get_cart_from_spher( spher )

         cart = MatMul( RotMat , cart )

         theta =  calc_spherical_theta( cart )
         phi   =  calc_spherical_phi( cart )
         !Write( 777 , "(2F8.5)" ) theta , phi
         bin1 = get_bin( theta , dtheta , N1 )
         bin2 = get_bin( phi , dphi , N2 )
         !Write( 778 , "(2I3)" ) bin1 , bin2
         Polar( bin1 , bin2 ) = Temp( i , j )
       End Do
     End Do



     Deallocate( Temp )



   End Subroutine RotatePolarDist




   !!
   !! supply vectors A and B the function returns
   !! a rotation matrix from vector A to vector B
   !! works only for 3d
   !! obtained rotation matrix can be used as x'=Matmul( RotMat , x )
   Function ComputeRotationMatrix( A , B ) Result( RotMat )

     use tool_interfaces, only :Cross_Product,get_norm
     Implicit None
     Real*8,dimension( : ) ::A
     Real*8,dimension( : ) ::B
     Real*8,dimension( 1:3 , 1:3 ) ::RotMat

     Real*8,dimension( 1:3 ) ::v
     Real*8 ::s
     Real*8 ::c


     Integer ::i

     RotMat = 0d0

     Do i = 1, 3
       RotMat( i , i )  = 1d0
     End Do


     v = Cross_Product( A , B )

     s = get_norm( v )
     s = s*s
     c = Dot_Product( A , B )



     v = v / get_norm( v )
     s = Dsqrt( s )

     RotMat( 1 , 1 ) = v(1)*v(1)*(1-c)+c
     RotMat( 1 , 2 ) = v(1)*v(2)*(1-c)-s*v(3)
     RotMat( 1 , 3 ) = v(1)*v(3)*(1-c)+s*v(2)

     RotMat( 2 , 1 ) = v(1)*v(2)*(1-c)+v(3)*s
     RotMat( 2 , 2 ) = v(2)*v(2)*(1-c)+c
     RotMat( 2 , 3 ) = v(2)*v(3)*(1-c)-s*v(1)

     RotMat( 3 , 1 ) = v( 1 )*v( 3 )*( 1 - c ) - v( 2 )*s
     RotMat( 3 , 2 ) = v(2)*v(3)*(1-c)+s*v(1)
     RotMat( 3 , 3 ) = v(3)*v(3)*(1-c)+c

   End Function ComputeRotationMatrix




   !!
   !!  multiplies an vector of vectors times
   !!  a matrix
   !!
   !!

   Function MatMulArray( Array , matrix ) Result( Prod )

     Implicit None
     Real*8,dimension( : , : )  ::Array
     Real*8,dimension( : , : )  ::matrix
     Real*8,allocatable,dimension( : , : )  ::Prod

     Integer ::i

     Allocate( Prod( 1:size( Array , 1 ) , 1:Size( Array , 2 ) ) )

     Do i = 1, Size( Array, 2 )
       Prod( : , i )  =  MatMul( matrix , Array( : , i ) )
     End Do

   End Function


   !!
   !! compute fourier grid for given reciprocal
   !! lattice
   !!
   Function ComputeQgrid( Basis , Nx , Ny , Nz ) Result( Grid )

     Implicit None
     Real*8,dimension( : , : )  ::Basis
     Integer ::Nx
     Integer ::Ny
     Integer ::Nz

     Real*8,allocatable,dimension( : , : )  ::Grid

     Integer ::kx
     Integer ::ky
     Integer ::kz

     Real*8 ::dx
     Real*8 ::dy
     Real*8 ::dz


     Integer ::counter

     Allocate( Grid( 1:3 , 1 : Nx * Ny * Nz ) )

     dx = 1d0 / Dble( Nx )
     dy = 1d0 / Dble( Ny )
     dz = 1d0 / Dble( Nz )


     counter = 1
     Do kx = 0 , Nx-1
       Do ky = 0 , Ny-1
         Do kz = 0 , Nz-1
           Grid( : , counter )  =  Dble( kx ) * Basis( : , 1 ) * dx
           Grid( : , counter )  =  Grid( : , counter ) + Dble( ky ) * Basis( : , 2 ) * dy
           Grid( : , counter )  =  Grid( : , counter ) + Dble( kz ) * Basis( : , 3 ) * dz
           counter  =  counter + 1
         End Do
       End Do
     End Do
   End Function ComputeQgrid






   !!
   !! supply position in direct coordinates
   !! and the lattice to which those direct coords
   !! belong
   !! function will return the periodic images
   !! depending on the mode the function
   !! mode = 1:
   !! function returns only periodic images along x,y,z
   !! mode= 2
   !! function also returns periodic images along
   !! x,y,z,xy,xz,yz
   !! mode = 3
   !! x,y,z,xy,xz,yz,xyz
   !! takinginto acount +-

   Function PeriodicImageCartS( Pos , lattice , mode ) Result( PPositions ) !J.L 23.01.2019

      Implicit None
      Real*8,dimension( : )   ::Pos
      Real*8,dimension( : , : )  ::lattice
      Integer ::mode
      Real*8,allocatable,dimension( : , : ) ::PPositions   !!Periodic positions
      Real*8,dimension( 1:3 )  ::TempPos

      Integer ::i
      Integer ::j
      Integer ::Nr

      If ( mode .eq. 1 ) then
          Allocate( PPositions( 1:3 , 7 ) )
      ElseIf ( mode .eq. 2 ) then
          Write( * , * ) " ERRRRRRRRRRRRRROR  "
          Write( * , * ) "Tell Jonathan to program mode 2"
          Write( * , * ) " ERRRRRRRRRRRRRROR  "
      ElseIf ( mode .eq. 3 ) then
          Write( * , * ) " ERRRRRRRRRRRRRROR  "
          Write( * , * ) "Tell Jonathan to program mode 3"
          Write( * , * ) " ERRRRRRRRRRRRRROR  "
      End If

      PPositions( : , 1 ) = MatMul( lattice , Pos )
      If ( mode .eq. 1 ) then
        Nr = 2
        Do j = 1 , 3
          Do i = -1 , 1 , 2
             TempPos = Pos
             TempPos( j )  =  Pos( j ) + Dble( i )
             PPositions( : , Nr ) = MatMul( lattice , TempPos )
             Nr  =  Nr  +  1
          End Do
        End Do
      End If

   End Function PeriodicImageCartS





   !!
   !!
   !!  Supply a filename and a certain keyword
   !!  to the routine; the function returns
   !!  the number of occourences of the keyword in the
   !!  file
   !!
   Function count_words( fname , word ) Result( Nwords )

     Implicit None
     Character( len = * )  ::fname
     Character( len = * )  ::word
     Integer   ::Nwords

     Write( * , * ) "System call"
     Write( * , * ) 'grep '// '"'//Trim( word ) //'"'// ' '//Trim(fname) // ' | wc -l > wc.txt'
     Call execute_command_line( 'grep '//'"'// Trim( word ) //'"'// ' ' // Trim(fname) // ' | wc -l > wc.txt' )
     open( unit = 30 , action='read' , status='old' , file = "wc.txt" )
     Read( 30 , * ) Nwords
     close( 30 )
   End Function count_words

   Function get_periodic_images_and_link( data ) Result( images )

      Implicit None
      Real*8,dimension( : , : )  ::data
      Integer ::i
      Integer ::aa
      Integer ::bb
      Integer ::cc
      !Real*8,dimension( 1 : 27 * Size( data , 1 ) , 1 : Size( data , 2 ) + 1 )  ::images
      Real*8,dimension( 1 : Size( data , 1 ) + 1 ,  1 : 27 * Size( data , 2 ) )  ::images
      Integer :: row

        images = 0

        row = Size( data , 2 ) + 1
        Do aa = -1 , 1
          Do bb = -1 , 1
            Do cc = -1 , 1
              If ( aa .eq. 0 .and. bb .eq. 0 .and. cc .eq. 0 ) Cycle
              !!!$OMP Parallel Do Shared( aa , bb , cc , images , data ) Private( i , row )
              Do i = 1 , Size( data , 2 )
                images( 1 : 3 , i )   =   data( 1 : 3 , i )
                images( 4 , i )       =   Dble( i )
                images( 1 , row )     =   data( 1 , i ) + Dble( aa )
                images( 2 , row )     =   data( 2 , i ) + Dble( bb )
                images( 3 , row )     =   data( 3 , i ) + Dble( cc )
                images( 4 , row )     =   Dble( i )
                row = row + 1
              End Do
             !!!$OMP End Parallel Do
            End Do
          End Do
        End Do


   End Function get_periodic_images_and_link

   Function get_periodic_images_and_link_cart( data , lattice ) Result( images )

      Implicit None
      Real*8,dimension( : , : )  ::data
      Real*8,dimension( : , : )  ::lattice
      Integer ::i
      Integer ::aa
      Integer ::bb
      Integer ::cc
      !Real*8,dimension( 1 : 27 * Size( data , 1 ) , 1 : Size( data , 2 ) + 1 )  ::images
      Real*8,dimension( 1 : Size( data , 1 ) + 1 ,  1 : 27 * Size( data , 2 ) )  ::images
      Integer :: row

        images = 0

        row = Size( data , 2 ) + 1
        Do aa = -1 , 1
          Do bb = -1 , 1
            Do cc = -1 , 1
              If ( aa .eq. 0 .and. bb .eq. 0 .and. cc .eq. 0 ) Cycle
              !!!$OMP Parallel Do Shared( aa , bb , cc , images , data ) Private( i , row )
              Do i = 1 , Size( data , 2 )
                images( 1 : 3 , i )   =   data( 1 , i ) * lattice( : , 1 ) + &
                                          data( 2 , i ) * lattice( : , 2 ) + &
                                          data( 3 , i ) * lattice( : , 3 )
                images( 4 , i )       =   Dble( i )
                images( 1 , row )     =   data( 1 , i ) + Dble( aa )
                images( 2 , row )     =   data( 2 , i ) + Dble( bb )
                images( 3 , row )     =   data( 3 , i ) + Dble( cc )
                images( 4 , row )     =   Dble( i )
                images( 1:3 , row )     =  images( 1 , row ) * lattice( : , 1 ) + &
                                           images( 2 , row ) * lattice( : , 2 ) + &
                                           images( 3 , row ) * lattice( : , 3 )
                row = row + 1
              End Do
             !!!$OMP End Parallel Do
            End Do
          End Do
        End Do


   End Function get_periodic_images_and_link_cart


   Function get_periodic_images( data ) Result( images )

     Implicit None
     Real*8,dimension( : , : )  ::data
     Integer ::i
     Integer ::aa
     Integer ::bb
     Integer ::cc
     Real*8,dimension( 1:3 , 1 : 27 * Size( data , 2 ) )  ::images
     Integer :: row

     images = Huge( images( 1 , 1 ) )
     row = Size( data , 2 ) + 1
     Do aa = -1 , 1
       Do bb = -1 , 1
         Do cc = -1 , 1
           If ( aa .eq. 0 .and. bb .eq. 0 .and. cc .eq. 0 ) Cycle
              Do i = 1 , Size( data , 2 )
                images( 1:3 , i )   =  data( 1:3 , i )
                images( 1 , row )   =  data( 1 , i ) + Dble( aa )
                images( 2 , row )   =  data( 2 , i ) + Dble( bb )
                images( 3 , row )   =  data( 3 , i ) + Dble( cc )
                row = row + 1
              End Do
         End Do
       End Do
     End Do

   End Function get_periodic_images


   Function get_norm( in_vector )
      Implicit None
      Real*8,dimension( : )  ::in_vector
      Real*8                 ::get_norm

      get_norm  =   Dsqrt( Dot_Product( in_vector , in_vector ) )

   End Function get_norm


   Subroutine sort_I(tempdata,dimen,dd,num)
     implicit none
     Integer                           ::a,b
     Integer                           ::dimen,dd
     Real*8,dimension(1:dd)            ::temp1
     Real*8,dimension(1:dimen,1:dd)    ::tempdata
     Integer                           ::num

     Do a=1,dimen-1
       Do b=a+1,dimen
           If (tempdata(a,num).gt.tempdata(b,num)) then
             temp1(:)=tempdata(b,:)
             tempdata(b,:)=tempdata(a,:)
             tempdata(a,:)=temp1(:)
           End If
       End Do
     End Do

   End Subroutine


   Subroutine sort_II( tempdata , dimen , dd , num )
     implicit none
     Integer                               ::a , b
     Integer                               ::dimen , dd
     Real*8,dimension( 1 : dd )            ::temp1
     Real*8,dimension( 1 : dd , 1 : dimen )::tempdata
     Integer                               ::num

     Do a = 1 , dimen-1
       Do b = a + 1 , dimen
           If ( tempdata( num , a ) .gt. tempdata( num , b ) ) then
             temp1( : ) = tempdata( : , b )
             tempdata( : , b ) = tempdata( : , a )
             tempdata( : , a ) = temp1( : )
           End If
       End Do
     End Do

   End Subroutine

   Subroutine sort_III( tempdata , dd )
     implicit none
     Integer                    ::a , b
     Integer                    ::dd
     Real*8                     ::temp1
     Real*8,dimension( 1 : dd ) ::tempdata

     Do a = 1 , dd-1
       Do b = a + 1 , dd
         If ( tempdata( a ) .gt. tempdata( b ) ) then
           temp1  = tempdata( b )
           tempdata( b ) = tempdata( a )
           tempdata( a ) = temp1
         End If
       End Do
     End Do
   End Subroutine


   !Returns the nearest image of position B with respect to
   !position A
   !

   Function get_nearest_image( vectorA , vectorB )  Result(  vectorC )

     Implicit None
     Real*8,dimension( : )  ::vectorA
     Real*8,dimension( : )  ::vectorB
     Real*8,dimension( 1 : Size( vectorA , 1 ) ) ::vectorC
     Real*8    ::delta
     Integer   ::i

     Do i = 1 , Size( vectorA , 1 )
        delta         =   vectorA( i )  -  vectorB( i )
        vectorC( i )  =   vectorB( i )
        If ( Abs( delta ) .gt. 0.5d0 ) then
           vectorC( i )  =  vectorC( i )  +  Sign( 1d0 , delta )
        End If
     End Do

   End Function get_nearest_image



   !Function returns a vector which alters the supplied vector
   !in a way that all distances are agree with a minimum image
   !convention in agreement with periodic boundaries
   Function periodic_dist( vector_in ) Result( vector_out )

     Implicit None
     Real*8,dimension( : )  ::vector_in
     Real*8,dimension( 1 : Size( vector_in , 1 ) ) ::vector_out
     Integer ::i

     Do i = 1 , Size( vector_in , 1 )
       vector_out( i )  =  vector_in( i )
       If ( Abs( vector_out( i ) ) .gt. 0.5d0 ) then
         vector_out( i )   =   vector_out( i ) - Sign( 1d0 , vector_out( i ) )
       End If
     End Do

   End Function periodic_dist




   Function calc_spherical_theta( vector )

      use tool_interfaces, only :get_norm
      Implicit None
      Real*8,dimension( : )  ::vector
      Real*8                 ::r
      Real*8,parameter       ::pi  =  2d0 * DAcos( 0d0 )
      Real*8                 ::calc_spherical_theta
      Real*8                 ::theta

      If ( Size( vector , 1 ) .ne. 3 ) then
         Write( * , "(A)" , advance = 'No' ) "Error in calc theta"
         Write( * , "(A)"  ) "vector too large; only 3-vectors allowed "
      End If
      r  =  get_norm( vector )

      If ( vector( 3 ) .eq. 0d0 ) then
        theta  =  pi / 2d0
      Else
        theta  =  Datan( Dsqrt( ( vector( 1 ) / r )**2 + ( vector( 2 ) / r )**2 ) &
                              / ( vector( 3 ) / r ) )
      End If

      If ( vector( 3 ) .lt. 0d0 ) then
        theta  =  pi + theta
      End If

      If ( theta .lt. 0d0 ) then
        theta  =  Abs( theta )
      End If

      If ( theta .gt. pi ) then
        theta  =  theta - pi
      End If


      calc_spherical_theta  =  theta

   End Function calc_spherical_theta



   Function calc_spherical_phi( vector )

      use tool_interfaces, only :get_norm
      Implicit None
      Real*8,dimension( : )  ::vector
      Real*8                 ::r
      Real*8,parameter       ::pi  =  2d0 * DAcos( 0d0 )
      Real*8                 ::calc_spherical_phi
      Real*8                 ::phi

      If ( Size( vector , 1 ) .ne. 3 ) then
        Write( * , "(A)" , advance = 'No' ) "Error in calc phi"
        Write( * , "(A)" ) "vector too large; only 3-vectors allowed"
      End If
      r  =  get_norm( vector )
      phi   =  Datan( ( vector( 2 ) / r ) / ( vector( 1 ) / r ) )

      If ( vector( 1 ) .eq. 0d0 .and. vector( 2 ) .eq. 0d0 ) then
        phi  =  0d0
      Else If ( vector( 1 ) .gt. 0d0 .and. vector( 2 ) .lt. 0d0 ) then
        phi  =  2d0 * pi + phi
      Else If ( vector( 1 ) .eq. 0d0 .and. vector( 2 ) .gt. 0d0 ) then
        phi  =  pi/2d0
      Else If ( vector( 1 ) .eq. 0d0 .and. vector( 2 ) .lt. 0d0 ) then
        phi  =  3d0 * pi / 2d0
      Else If ( vector( 1 ) .le. 0d0 ) then
        phi  =  phi + pi
      End If

      If ( phi .lt. 0d0 ) then
        phi  =  2d0 * pi + phi
      End If
      If ( phi .gt. 2d0 * pi ) then
        phi  =  phi - 2d0 * pi
      End If

      calc_spherical_phi   =   phi

   End Function calc_spherical_phi


   Function calc_pearson_terms( xnew , ynew ) Result( data )

      Implicit None
      Real*8         ::xnew
      Real*8         ::ynew
      Real*8,dimension( 1 : 5 ) ::data

      data( 1 )   =  xnew  *  ynew
      data( 2 )   =  xnew
      data( 3 )   =  xnew**2
      data( 4 )   =  ynew
      data( 5 )   =  ynew**2

   End Function calc_pearson_terms



   Function mean_value_IIxII_array( in_array , N ) Result( out_array )

     Implicit None
     Integer   ::N
     Real*8,dimension( : , : )  ::in_array
     Real*8,dimension( 1 : Size( in_array , 1 ) , 1 : Size( in_array , 2 ) ) &
                         ::out_array

     out_array  =  in_array  /  Dble( N )

   End Function mean_value_IIxII_array




   Function calc_single_pearson_from_terms( in_data , N )

     Implicit None
     Integer       ::N
     Real*8,dimension( 1 : 5 )  ::in_data
     Real*8                     ::calc_single_pearson_from_terms
     Real*8 ::varx
     Real*8 ::vary
     Real*8 ::xmean
     Real*8 ::ymean
     Real*8 ::p

     xmean = in_data( 2 ) / Dble( N )
     ymean = in_data( 4 ) / Dble( N )


     varx = in_data( 3 ) / Dble( N ) - ( xmean )**2
     vary = in_data( 5 ) / Dble( N ) - ( ymean )**2

     varx = Dsqrt( varx )
     vary = Dsqrt( vary )


     p  =  ( in_data( 1 ) - xmean * in_data( 4 ) - ymean * in_data( 2 ) &
             + Dble( N ) * xmean * ymean ) / Dble( N )
     p  =  p / ( varx * vary )

     calc_single_pearson_from_terms = p
   End Function

   Function bin( theta , dtheta , Ntheta )

    implicit none
    Real*8            ::theta
    Real*8            ::dtheta
    Integer           ::Ntheta
    Integer           ::bin

    bin = Ceiling( theta / dtheta )
    If ( bin .lt. 1 )       bin  =  1
    If ( bin .gt. Ntheta )  bin  =  Ntheta
  End Function bin



  Function shift_atom_to_prim_cell( atom ) Result( atom_in_cell )

    Implicit None
    Real*8,dimension( : )  ::atom
    Real*8,dimension( 1 : Size( atom ) ) ::atom_in_cell
    Integer ::i

    atom_in_cell( : )  =  atom( : )
    Do i = 1 , Size( atom , 1 )
      If ( atom( i ) .lt. 0d0 ) atom_in_cell( i )  =  1d0 + atom_in_cell( i )
      If ( atom( i ) .gt. 1d0 ) atom_in_cell( i )  =  1d0 - atom_in_cell( i )
    End Do

  End Function shift_atom_to_prim_cell


  Subroutine write_poscar_format( cp , latc , species , Natoms , coord_type , coords )

    Implicit None
    Real*8    ::cp
    Real*8,dimension( : , : )  ::latc
    Character( len = * )   ::species
    Integer,dimension( : ) ::Natoms
    Character( len = * )   ::coord_type
    Real*8,dimension( : , : ) ::coords
    Integer ::i


    open( unit = 25 , action='write' , file ='POSCAR' , status='unknown' )
    Write( 25 , "(A)" ) "INSERT NAME HERE"
    Write( 25 , "(F15.8)" ) cp
    Do i = 1 , Size( latc , 1 )
      Write( 25 , "(3F10.4)" ) latc( i , : )
    End Do
    Write( 25 , "(A)" ) Trim( species )
    Write( 25 , "(10000I5)" ) Natoms( : )
    Write( 25 , "(A)" ) Trim( coord_type )
    Do i =  1 , Size( coords , 2 )
      Write( 25 , "(3F10.6)" ) coords( : , i )
    End Do
    close( 25 )

  End Subroutine write_poscar_format


  Function get_N_ats_unit_cell( Ntot , unit_size )

     Implicit None
     Integer                 ::Ntot
     Real*8,dimension( : )   ::unit_size
     Integer                 ::get_N_ats_unit_cell
     Integer                 ::i

     get_N_ats_unit_cell  =  Ntot

     Do i = 1 , Size( unit_size , 1 )
       get_N_ats_unit_cell  =  NINT( get_N_ats_unit_cell * unit_size( i ) )
     End Do

  End Function get_N_ats_unit_cell


  Function get_nearest_neighbour_no_periodicity( atoms , center , N , Nself ) Result( list )

     use tool_interfaces, only :get_norm
     Implicit None
     Real*8,dimension( : , : )  ::atoms
     Real*8,dimension( : )      ::center
     Integer                    ::N
     Integer,dimension( 1 : N ) ::list
     Real*8,allocatable,dimension( : , : ) ::temp_dat
     Integer  ::Nself
     Integer  ::i

     Allocate( temp_dat( 1 : Size( atoms , 2 ) , 1 : 2 ) )

     Do i = 1 , Size( atoms , 2 )
       If ( i .ne. Nself ) then
          temp_dat( i , 1 )  =  i
          temp_dat( i , 2 )  =  get_norm( atoms( : , i ) - center( : ) )
       Else
          temp_dat( i , : )  =  100000000d0
       End If
     End Do

     Call sort_I( temp_dat , Size( atoms , 2 ) , 2 , 2 )

     Do i = 1 , N
       list( i )  =  NINT( temp_dat( i , 1 ) )
     End Do


     Deallocate( temp_dat )

  End Function


  Function get_cart_from_spher( vector ) Result( cart_vec )

    Implicit None
    Real*8,dimension( : )  ::vector
    Real*8,dimension( 1 : 3 ) ::cart_vec

    If ( Size( vector , 1 ) .ne. 3 ) then
      Write( * , "(A)" , advance = 'No' ) "Vector passed to routine"
      Write( * , "(A)" ) "'get_cart_from_spher' has wrong size"
      Return
    End If

    cart_vec( 1 )  =  vector( 3 ) * Dcos( vector( 2 ) ) * Dsin( vector( 1 ) )
    cart_vec( 2 )  =  vector( 3 ) * Dsin( vector( 2 ) ) * Dsin( vector( 1 ) )
    cart_vec( 3 )  =  vector( 3 ) * Dcos( vector( 1 ) )

  End Function get_cart_from_spher



  Function spherical_sin_correction( in_data , dimen , u_limit) Result( out_data )

    Implicit None
    Integer   ::i
    Real*8    ::delta
    Real*8,dimension( : , : )  ::in_data
    Real*8,dimension( 1 : Size( in_data , 1 ) , &
                      1 : Size( in_data , 2 ) )  ::out_data
    Real*8  ::u_limit
    Integer ::dimen

    delta  =  u_limit  /  Dble( Size( in_data , dimen ) )


    If ( dimen .eq. 1 ) then
      Do i = 1 , Size( in_data , dimen )
        out_data( i , : )   =   in_data( i , : )  /  &
                             Dsin( Dble( i ) * delta - delta / 2d0 )
      End Do
    ElseIf ( dimen .eq. 2 ) then
      Do i = 1 , Size( in_data , dimen )
        out_data( : , i )   =   in_data( : , i )  /  &
                             Dsin( Dble( i ) * delta - delta / 2d0 )
      End Do
    Else
      Write( * , "(A)" ) "ERROR"
      Write( * , "(A)" ) "Sine correction failed wrong dimension given, original array returned"
      out_data( : , : )  =  in_data( : , : )
    End If

  End Function spherical_sin_correction



  Subroutine gen_file_number( Num , file_string )

    Implicit None
    Character( len = * ) ::file_string
    Integer  ::Num

    If ( Num .lt. 10 ) then
      Write( file_string , "(I1)" ) Num
    ElseIf ( Num .lt. 100 ) then
      Write( file_string , "(I2)" ) Num
    ElseIf ( Num .lt. 1000 ) then
      Write( file_string , "(I3)" ) Num
    ElseIf ( Num .lt. 10000 ) then
      Write( file_string , "(I4)" ) Num
    ElseIf ( Num .lt. 100000 ) then
      Write( file_string , "(I5)" ) Num
    End If

  End Subroutine gen_file_number


  Function histogram_2d( data , min1 , min2 , max1 , max2 , N1 , N2 ) Result( hist )

    Implicit None
    Real*8,dimension( : , : )  ::data
    Real*8     ::min1
    Real*8     ::min2
    Real*8     ::max1
    Real*8     ::max2
    Integer    ::N1
    Integer    ::N2
    Real*8     ::d1
    Real*8     ::d2
    Integer    ::i
    Integer,dimension( 1 : N1 , 1 : N2 )  ::hist
    Integer ::bin1
    Integer ::bin2



    hist = 0

    d1  =  ( max1 - min1 ) / Dble( N1 )
    d2  =  ( max2 - min2 ) / Dble( N2 )

    Do i = 1 , Size( data , 1 )
      bin1   =   Ceiling( data( 1 , i ) / d1 )
      bin2   =   Ceiling( data( 2 , i ) / d2 )
      If ( bin1 .le. 0 ) bin1  =  1
      If ( bin2 .le. 0 ) bin2  =  1
      If ( bin1 .gt. N1 ) bin1  =  N1
      If ( bin2 .gt. N2 ) bin2  =  N2
      hist( bin1 , bin2 ) = hist( bin1 , bin2 ) + 1
    End Do
  End Function histogram_2d



  Function Cross_Product( a , b ) Result( c )
    Implicit None
    Real*8,dimension( : ) ::a
    Real*8,dimension( : ) ::b
    Real*8,dimension( 1 : Size( a , 1 ) ) ::c

    c( 1 )  =  a( 2 ) * b( 3 ) - a( 3 ) * b( 2 )
    c( 2 )  =  a( 3 ) * b( 1 ) - a( 1 ) * b( 3 )
    c( 3 )  =  a( 1 ) * b( 2 ) - a( 2 ) * b( 1 )

  End Function Cross_Product


  Function get_2d_bin( theta , phi , dtheta , dphi , Ntheta , Nphi ) Result( bins )

     Implicit None
     Real*8   ::theta
     Real*8   ::phi
     Real*8   ::dtheta
     Real*8   ::dphi
     Integer  ::Ntheta
     Integer  ::Nphi
     Integer,dimension( 1 : 2 )  ::bins


    bins( 1 )  =  Ceiling( theta / dtheta )
    bins( 2 )  =  Ceiling( phi / dphi )


    If ( bins( 1 ) .le. 0 ) bins( 1 )  =  1
    If ( bins( 2 ) .le. 0 ) bins( 2 )  =  1
    If ( bins( 1 ) .gt. Ntheta ) bins( 1 )  =  Ntheta
    If ( bins( 2 ) .gt. Nphi ) bins( 2 )  =  Nphi

  End Function get_2d_bin

  !Function returns the a closest lower power of 2 of a given integer
  Function give_closest_lower_power_of_2( N1 ) Result( N2 )

     Implicit None
     Integer ::N1
     Integer ::N2

     N2 = N1
     If ( IAND ( N2 , N2 - 1 ) .ne. 0 ) then
           Do while ( N2 .ge. 0 )
              N2 = N2 - 1
              If ( IAND ( N2 , N2 - 1 ) .eq. 0 ) Exit
           End Do
     End If

     If ( N2 .eq. 0 ) then
       Write( * , * ) "ERROR in function give_closest_lower_power_of_2"
       Write( * , * ) "Number found return value is zero"
     End If


  End Function give_closest_lower_power_of_2


  !Pass vectors vectorA and vectorB
  !in direct coordinates , and the bravais matrix of the considered sysetm
  !Function returns vectorB in cartesian coordiantes closest to vectorA
  !Function get_nearest_image_cart( vectorA , vectorB , lattice )  Result(  vectorC )
  Function get_nearest_image_cartesian( vectorA , vectorB , lattice ) Result( vectorC )

    use tool_interfaces,only :get_norm
    Implicit None
    Real*8,dimension( : )  ::vectorA
    Real*8,dimension( : )  ::vectorB
    Real*8,dimension( : , : )  ::lattice
    Real*8,dimension( 1 : Size( vectorA , 1 ) ) ::vectorC
    Real*8,dimension( 1 : Size( vectorA , 1 ) ) ::tempA
    Real*8,dimension( 1 : Size( vectorA , 1 ) ) ::tempB
    Real*8   ::norm
    Real*8   ::normTemp
    Integer  ::i
    Integer  ::j
    Integer  ::k

    tempA =  Matmul( lattice , vectorA )
    norm  =  get_norm( Matmul( lattice , vectorA - vectorB ) )
    vectorC = Matmul( lattice , vectorB )
    !print*,tempA,norm,vectorC
    Do i = -1 , 1
       Do j = -1 , 1
          Do k = -1 , 1
             If ( i .eq. 0 .and. j .eq. 0 .and. k .eq. 0 ) Cycle
             tempB( 1 ) = vectorB( 1 ) + Dble( i )
             tempB( 2 ) = vectorB( 2 ) + Dble( j )
             tempB( 3 ) = vectorB( 3 ) + Dble( k )
             normTemp = get_norm( Matmul( lattice , vectorA - tempB ) )
             If ( normTemp .lt. norm ) then
               vectorC = Matmul( lattice , tempB )
               norm = normTemp
             End If
          End Do
       End Do
    End Do
    !print*,vectorC
  End function get_nearest_image_cartesian



  !Supply 2 arrays of atomic coordinates in direct coordinates
  !to the function. The two arrays have to be of equal size
  !At output the function will return an array which contains
  !to every atom from atomsB the nearest neighbour coordinates
  !of the array atomsA atoms

  !j.L 14.06.2018

  Function assign_closest_atoms_to_array( atomsA , atomsB ) Result( ClosestAtoms )

    use tool_interfaces,only :get_periodic_images,get_norm
    Real*8,dimension( : , : )  ::atomsA
    Real*8,dimension( : , : )  ::atomsB
    Real*8,dimension( 1:3 , 1:Size( atomsA , 2 ) )  ::ClosestAtoms
    Real*8,dimension( 1:3 , 1:27*Size( atomsA , 2 ) )  ::images
    Real*8   ::Sdist  !Shortest dist
    Real*8   ::Adist  !Actual dist
    Integer ::i
    Integer ::j

    If ( Size( atomsA , 2 ) .ne. Size( atomsB , 2 ) ) then
      Write( * , * ) "Error  in function assign_closest_atoms_to_array"
      Write( * , * ) "The dimensions of the supplied arrays do not match"
      Write( * , * ) "Return to calling routine with zeros on output"
      ClosestAtoms = 0d0
      Return
    End If

    images  =  get_periodic_images( atomsA )


    Do i = 1 , Size( atomsB , 2 )
         Sdist = Huge( Sdist )
         Do j = 1 , Size( images , 2 )
            Adist  =  get_norm( atomsB( : , i ) - images( : , j ) )
            If ( Adist .lt. Sdist ) then
               Sdist = Adist
               ClosestAtoms( : , i )  =  images( : , j )
            End If
         End Do
    End Do

  End Function assign_closest_atoms_to_array


  !Routine computes the gemetrical center of an atoms array and shifts
  !it to the point 0.5,0.5,0.5
  !On output the atomic positions are
  !defined with respect to 0.5,0.5,0.5
  !Only working for 3d vector positions
  !j.L. 15.06.2018
  Subroutine Drift_correction( atoms )

    Implicit None
    Real*8,dimension( : , : ) ::atoms
    Real*8,Save,dimension( 1:3 )  ::center = [0.5d0, 0.5d0, 0.5d0 ]
    Real*8,dimension( 1:3 ) :: COM
    Integer ::i

    COM = 0d0

    If ( Size( atoms , 1 ) .eq. 3 ) then

      Do i = 1 , Size( atoms , 2 )
         COM( : )  =  COM( : ) + atoms( : , i )
      End Do

      COM( : )  =  COM( : ) / Dble( Size( atoms , 2 ) )
      COM( : )  =  COM( : ) - center( : )
      Do i = 1 , Size( atoms , 2 )
        atoms( : , i )  =  atoms( : , i )  - COM( : )
      End Do
    Else if ( Size( atoms , 2 ) .eq. 3 ) then

      Do i = 1 , Size( atoms , 1 )
         COM( : )  =  COM( : ) + atoms( i , : )
      End Do

      COM( : )  =  COM( : ) / Dble( Size( atoms , 2 ) )
      COM( : )  =  COM( : ) - center( : )

      Do i = 1 , Size( atoms , 1 )
        atoms( i , : )  =  atoms( i , : )  - COM( : )
      End Do
    End If

  End Subroutine Drift_correction


  Subroutine transformation_matrix( latVa , latVb , tm )

!    use vector_interfaces , only       :vector_norm
    implicit none
    Real*8,dimension( : , : )  ::latVa
    Real*8,dimension( : , : )  ::latVb
    Real*8,dimension( : , : )  ::tm
    Integer                    ::i
    Integer                    ::j
    Integer                    ::N
    Integer                    ::M

    If ( Size( latVa( : , 1 ) ) .ne. Size( latVb( : , 1 ) ) ) then
       Write( * , "(A)" , advance = 'No') "Error in determination of"
       Write( *, "(A)" ) "transformation matrix unequal number of basis vectors"
       STOP
    End If
    If ( Size( latVa( 1 , : ) ) .ne. Size( latVb( 1 , : ) ) ) then
       Write( * , "(A)" , advance = 'No' ) "Error in determination of"
       Write( * , "(A)" ) "transformation matrix unequal number of basis vectors"
       STOP
    End If

    N  =  Size( latVa( : , 1 ) )
    M  =  Size( latVa( 1 , : ) )

    Do i = 1 , N
      Do j =  1 , M
         tm( i , j )  =  Dot_Product( latVa( : , i ) , latVb( : , j ) )  !/ &
                        !vector_norm( latVa( : , i ) ) / vector_norm( latVb( : , i ) )
      End Do
    End Do
  End Subroutine transformation_matrix
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !!!
  !!! inplace coordinate transformation
  !!! 
  Subroutine coordinate_transformationIP( A , x )


    Implicit None
    Real*8,intent( in ),dimension( : , : ) ::A
    Real*8,intent( inout ),dimension( : , : )  ::x

    Integer ::i
    Do i = 1 , Size( x , 2 )
      x( : , i )  =  Matmul( A , x( : , i ) )
    End Do
  End Subroutine coordinate_transformationIP



  Subroutine coordinate_transformation( A , x , xprime )

    !vector transformtaion
    !A supplied transformation matrix
    !x vectors to be transfromed
    !xprime transformed vectors

    implicit none
    Real*8,dimension( : , : )  ::A
    Real*8,dimension( : , : )  ::x
    Real*8,dimension( : , : )  ::xprime
    Integer                    ::i

    Do i =  1 , Size( x , 2 )
      !matrix has to be transposed because of fortran
      !like storing of the values
      xprime( : , i )  =  Matmul( A , x( : , i ) )
    End Do
  End Subroutine coordinate_transformation

  Subroutine lattice_consts( a , b , c , alpha , beta , gamma , bravais )

   !use vector_interfaces, only  :vector_norm
   use tool_interfaces,only :get_norm
   implicit none
   Real*8                    ::a,b,c
   Real*8                    ::alpha,beta,gamma
   Real*8,dimension( 3 , 3 ) ::bravais


   a  =  get_norm( bravais( : , 1 ) )
   b  =  get_norm( bravais( : , 2 ) )
   c  =  get_norm( bravais( : , 3 ) )


   gamma  =  Dot_Product( bravais( : , 1 ) , bravais( : , 2 ) ) / a / b
   gamma  =  Dacos( gamma )

   beta   =  Dot_Product( bravais( : , 1 ) , bravais( : , 3 ) ) / a / c
   beta   =  Dacos( beta )

   alpha  =  Dot_Product( bravais( : , 2 ) , bravais( : , 3 ) ) / b / c
   alpha  =  Dacos( alpha )

  End Subroutine lattice_consts
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Function make_lattice( a , b, c ) Result( lattice )
     !This is only a test version which will work for perfectly cubic systems
     Implicit None
     Real*8   ::a
     Real*8   ::b
     Real*8   ::c
     Real*8,dimension( 1 : 3 , 1 : 3 )  ::lattice


     lattice  =  0d0
     lattice( 1 , 1 )  =  a
     lattice( 2 , 2 )  =  b
     lattice( 3 , 3 )  =  c
  End Function make_lattice


  Function lattice_volume( Bravais ) Result( volume ) !j.L 20.07.2018
      !compute the volume of a given bravais lattice via
      !the scalar triple product

      use tool_interfaces,only :cross_product
      Implicit None
      Real*8,dimension( : , : )  ::Bravais
      Real*8  ::volume
      Real*8,dimension( 1:3 ) ::vec_prod

      If ( Size( Bravais , 1 ) .ne. 3 .or. Size( Bravais , 2 ) .ne. 3 ) then
         Write( * , * ) "         ERROR       "
         Write( * , * ) "Error in function lattice volume"
         Write( * , * ) "dimensions of bravais lattice not correct"
         Write( * , * ) "routine only working for dimension 3"
         Write( * , * ) "retrun value of 1"
         Write( * , * ) "         ERROR   "
         volume = 1d0
         Return
      End If

      vec_prod = cross_product( Bravais( : , 2 ) , Bravais( : , 3 ) )
      volume = Abs( Dot_Product( vec_prod , Bravais( : , 1 ) ) )

  End Function lattice_volume




  Function compute_angle( vecA , vecB ) Result( angle )
    use tool_interfaces, only :get_norm
    Implicit None
    Real*8,dimension( : )  ::vecA
    Real*8,dimension( : )  ::vecB
    Real*8  ::normA
    Real*8  ::normB
    Real*8  ::angle

    If ( Size( vecA , 1 ) .ne. Size( vecB , 1 ) ) then
      Write( * , * ) "Error in function compute_angles"
      Write( * , * ) "Vectors are of differemnt dimensions"
      Write( * , * ) "Returning zero"
    End If

    normA  =  get_norm( vecA )
    normB  =  get_norm( vecB )

    angle = Dacos( Dot_product( vecA , vecB ) / normA / normB )

  End Function compute_angle






  Module Downfolding

    Implicit None
    Interface direct_reduction
      Module Procedure Reductio1D
    End Interface direct_reduction

    Interface folded_reduction
      Module Procedure Downfold1D
    End Interface folded_reduction

    Contains

    !
    !  expanded array |1.2.3.4.5|6.7.8.9.10|11.12.13.14.15|
    !  reduced array  |1+6+11.2+7+12.3+8+13.4+9+14.5+10+15|
    !
    Function Reductio1D( data , Period ) Result( folded )

        Implicit None
        Real*8,dimension( : )  ::data
        Integer ::Period

        Real*8,allocatable,dimension( : )  ::folded
        Integer ::Nrest
        Integer ::Nbin
        Integer ::i
        Integer ::Bin_red

        Bin_red =  Size( data , 1 ) / Period
        Nrest = Modulo( Size( data , 1 ) , Period )

        Allocate( folded( 1:Bin_red ) )
        folded = 0d0

        If ( Nrest .ne. 0 ) then
          Write( * , * ) "Error in downfold1D the supplied"
          Write( * , * ) "period does not match the size"
          Write( * , * ) "of the supplied array"
          Return
        End If

        Do i = 1 , Size( data , 1 )
          Nbin = Modulo( i , Bin_red )
          If ( Nbin .eq. 0 ) Nbin = Bin_red
          folded( Nbin )  =  folded( Nbin )  +  data( i )
        End Do
    End Function Reductio1D



    Function Downfold1D( data , Period ) Result( folded )

        Implicit None
        Real*8,dimension( : ) ::data
        Integer   ::Period

        Real*8,allocatable,dimension( : )  ::folded
        Integer ::Nrest
        Integer ::Nbin
        Integer ::i
        Integer ::Bin_red

        Integer ::Folder

        Bin_red =  Size( data , 1 ) / Period
        Nrest = Modulo( Size( data , 1 ) , Period )


        Allocate( folded( 1:Bin_red ) )
        folded = 0d0

        If ( Nrest .ne. 0 ) then
          Write( * , * ) "Error in downfold1D the supplied"
          Write( * , * ) "period does not match the size"
          Write( * , * ) "of the supplied array"
          Return
        End If


        Do i = 1 , Size( data , 1 )
          Nbin = Modulo( i , Bin_red )
          If ( Nbin .eq. 0 ) Nbin = Bin_red
          Folder = i / Bin_red
          If ( i .eq. Size( data , 1 ) ) Folder  =  Period
          If ( Modulo( Folder , 2 ) .ne. 0 ) then
            If ( Nbin .eq. Bin_red ) then
              Nbin = Bin_red
            Else
              Nbin = Bin_red - Nbin
            End If
          Else
          End If
          folded( Nbin )  =  folded( Nbin )  +  data( i )
        End Do

    End Function Downfold1D

  End Module Downfolding

  Module IntegerGrid

    Interface ComputeIntegerGrid
      Module Procedure IntGrid,IntGridCartesian
    End Interface ComputeIntegerGrid

    Contains
      !Supply atoms in direct coordinates
      Function IntGrid( atoms , Nx , Ny , Nz ) Result( Grid )

        use tool_interfaces, only :get_periodic_images_and_link
        Implicit None
        Real*8,dimension( : , : )  ::atoms
        Integer ::Nx
        Integer ::Ny
        Integer ::Nz
        Integer,allocatable,dimension( : , : ) ::Grid

        Real*8,dimension( 1:2 ) ::bounds
        Integer ::k
        Real*8,allocatable,dimension( : , : )  ::Images
        Real*8,allocatable,dimension( : , : )  ::AtomsC
        Integer ::i
        Integer ::j
        Real*8,dimension( 1:3 )  ::dx
        Integer,dimension( 1:3 ) ::Ncells
        Logical ::done

        If ( .not. Allocated( Grid ) ) then
          Allocate( Grid( 1:Size( atoms , 1 ) , 1:Size( atoms , 2 ) ) )
        End If
        Allocate( Images( 1 : Size( atoms , 1  ) + 1 , 27 * Size( atoms , 2 ) ) )
        Allocate( AtomsC( 1 : Size( atoms , 1  ) + 1 , Size( atoms , 2 ) ) )

        Images  =  get_periodic_images_and_link( atoms )

        bounds( 1 ) = -1d0
        bounds( 2 ) =  0d0
        Ncells( 1 ) =  Nx
        Ncells( 2 ) =  Ny
        Ncells( 3 ) =  Nz
        Do j = 1 , 3
          dx( j ) = 1d0 / Dble( Ncells( j ) )
        End Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        AtomsC( 1:Size( atoms ,1 ) , : ) = atoms( : , : )
        k = 1
        Do
            If ( k .gt. 50 ) then
              Write( * , "(A)" ) "ERROR!!! Integer grid projection not found stopping program"
              Write( * , "(A)" ) "eval_tools routine"
              STOP
            End If

            Do i = 1 , Size( images , 2 )
               If ( images( 1 , i ) .gt. bounds( 1 ) .and. images( 1 , i ) .lt. bounds( 2 ) .and. &
                       images( 2 , i ) .gt. bounds( 1 ) .and. images( 2 , i ) .lt. bounds( 2 ) .and. &
                          images( 3 , i ) .gt. bounds( 1 ) .and. images( 3 , i ) .lt. bounds( 2 ) ) then
                        j  =  NINT( images( 4 , i ) )
                        AtomsC( 1:3 , j )   =   images( 1 : 3 , i )
               End If
            End Do

            Do i = 1 , Size( atomsC , 2 )
               atomsC( 1:3 , i ) = atomsC( 1:3 , i ) - bounds( 1 )
            End Do

            Do j = 1 , Size( atoms , 1 )   !xyz loop
               Do i = 1 , Size( atomsC , 2 ) !looping over atoms
                  Grid( j , i )  =  NINT( AtomsC( j , i ) / dx( j ) )
                  If ( Grid( j , i ) .gt. Ncells( j ) ) Grid( j , i ) = 1
                  If ( Grid( j , i ) .lt. 1 ) Grid( j , i ) = Ncells( j )
               End Do
            End Do
            Call check_projection( Grid , done )
            If( done ) then
                Write( * , "(A,I5,A)" ) "Integer grid projection converged after" , k , " steps"
                Deallocate( images , AtomsC )
                Return !relative orientation calculation
            End If
            bounds = bounds + 0.01d0
            k = k + 1
        End Do


        If ( done ) then
          Write( * , * )  "Integer grid projection successful"
        Else
          Write( * , * ) "Integer grid projection failed. Don't trust your results"
        End If

        Deallocate( images , AtomsC )
      End Function IntGrid


      Function IntGridCartesian( atoms , Nx , Ny , Nz , lattice ) Result( Grid )

        use tool_interfaces, only :get_periodic_images_and_link_cart,get_norm
        Implicit None
        Real*8,dimension( : , : )  ::atoms
        Integer ::Nx
        Integer ::Ny
        Integer ::Nz
        Real*8,dimension( : , : )  ::lattice
        Integer,allocatable,dimension( : , : ) ::Grid

        Real*8,dimension( 1:2 ) ::bounds
        Integer ::k
        Real*8,allocatable,dimension( : , : )  ::Images
        Real*8,allocatable,dimension( : , : )  ::AtomsC
        Integer ::i
        Integer ::j
        Real*8,dimension( 1:3 )  ::dx
        Integer,dimension( 1:3 ) ::Ncells
        Real*8 ::MaxNorm
        Real*8 ::MinNorm

        Logical ::done

        If ( .not. Allocated( Grid ) ) then
          Allocate( Grid( 1:Size( atoms , 1 ) , 1:Size( atoms , 2 ) ) )
        End If
        Allocate( Images( 1 : Size( atoms , 1  ) + 1 , 27 * Size( atoms , 2 ) ) )
        Allocate( AtomsC( 1 : Size( atoms , 1  ) + 1 , Size( atoms , 2 ) ) )

        Images  =  get_periodic_images_and_link_cart( atoms , lattice )

        Ncells( 1 ) =  Nx
        Ncells( 2 ) =  Ny
        Ncells( 3 ) =  Nz
        MaxNorm = 0d0
        MinNorm = Huge( MinNorm )
        Do j = 1 , 3
          dx( j ) = get_norm( lattice( : , j ) ) / Dble( Ncells( j ) )
          If ( get_norm( lattice( : , j ) ) .gt. MaxNorm ) then
            MaxNorm = get_norm( lattice( : , j ) )
          End If
          If ( get_norm( lattice( : , j ) ) .lt. MinNorm ) then
            MinNorm = get_norm( lattice( : , j ) )
          End If
        End Do
        MinNorm = MinNorm / 100d0   !Step size for projection
        bounds( 1 ) = -MaxNorm
        bounds( 2 ) =  0d0
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Do i = 1, Size( atoms , 2 )
           AtomsC( 1:Size( atoms ,1 ) , i )  =  MatMul( lattice , atoms( : , i ) )
           AtomsC( Size( atomsC ,1 ) , i ) = i
        End Do

        k = 1
        Do
            If ( k .gt. 110 ) then
              Write( * , "(A)" ) "ERROR!!! Integer grid projection not found stopping program"
              Write( * , "(A)" ) "eval_tools routine cartesian"
              STOP
            End If

            Do i = 1 , Size( images , 2 )
               If ( images( 1 , i ) .gt. bounds( 1 ) .and. images( 1 , i ) .lt. bounds( 2 ) .and. &
                       images( 2 , i ) .gt. bounds( 1 ) .and. images( 2 , i ) .lt. bounds( 2 ) .and. &
                          images( 3 , i ) .gt. bounds( 1 ) .and. images( 3 , i ) .lt. bounds( 2 ) ) then
                        j  =  NINT( images( 4 , i ) )
                        AtomsC( 1:3 , j )   =   images( 1 : 3 , i )
               End If
            End Do

            Do i = 1 , Size( atomsC , 2 )
               atomsC( 1:3 , i ) = atomsC( 1:3 , i ) - bounds( 1 )
            End Do

            Do j = 1 , Size( atoms , 1 )   !xyz loop
               Do i = 1 , Size( atomsC , 2 ) !looping over atoms
                  Grid( j , i )  =  NINT( AtomsC( j , i ) / dx( j ) )
                  If ( Grid( j , i ) .gt. Ncells( j ) ) Grid( j , i ) = 1
                  If ( Grid( j , i ) .lt. 1 ) Grid( j , i ) = Ncells( j )
               End Do
            End Do
            Call check_projection( Grid , done )
            If( done ) then
                Write( * , "(A,I5,A)" ) "Integer grid projection converged after" , k , " steps"
                Deallocate( images , AtomsC )
                Return !relative orientation calculation
            End If
            bounds = bounds + MinNorm
            k = k + 1
        End Do

      End Function IntGridCartesian

      Subroutine check_projection( Grid , done )

          Implicit None
          Integer,dimension( : , : ) ::Grid
          Logical ::done
          Integer ::i
          Integer ::j

          Do i = 1 , Size( Grid , 2 ) - 1
             Do j = i + 1 , Size( Grid , 2 )
                If ( Grid( 1 , i ) .eq. Grid( 1 , j ) .and. &
                       Grid( 2 , i ) .eq. Grid( 2 , j ) .and. &
                          Grid( 3 , i ) .eq. Grid( 3 , j ) ) then
                        done = .False.
                        Return
                End If
             End Do
          End Do
          done = .True.
      End Subroutine check_projection

      !Grid has to be a 3d integer grid
      Function MakeXYZlistPositive( Grid ) Result( list )

          Implicit None
          Integer,dimension( : , : )  ::Grid
          Integer,allocatable,dimension( : , : ) ::list
          Integer ::i
          Integer ::j



          If ( Allocated( list ) ) Deallocate( list )
          Allocate( list( 1:3 , 1 : Size( Grid , 2 ) ) )

          Do i = 1 , Size( Grid , 2 )
            Do j = 1 , Size( Grid , 2 )

            End Do
          End Do
      End Function MakeXYZlistPositive
  End Module IntegerGrid


  !!!
  !!! module is taking distributions on input
  !!! and will apply gaussian smearing
  !!! matrices to your data
  !!! 
  Module GaussianSmearing

    Implicit None
    Real*8,parameter ::pi = 2d0 * DAcos( 0d0 )
    Interface ComputeGaussian
      module procedure ComputeGaussian1d,ComputeGaussian2d,&
                       ComputeGaussian3d
    End Interface ComputeGaussian

    Interface ComputeGaussianSmearing
      module procedure ComputeGaussianSmearing3d
    End Interface ComputeGaussianSmearing

    Contains

      Function ComputeGaussian3d( Nx , Ny , Nz , sigmax , sigmay , sigmaz ) Result( Gauss )

        Implicit None
        Integer ::Nx
        Integer ::Ny
        Integer ::Nz
        Real*8  ::sigmaX
        Real*8  ::sigmaY
        Real*8  ::sigmaZ

        Real*8,allocatable,dimension( : , : , : )  ::Gauss  
      
        Integer ::i
        Integer ::j
        Integer ::k

        Real*8  ::normA
        Real*8  ::normB
        Real*8  ::normC

        Real*8  ::MeanX
        Real*8  ::MeanY
        Real*8  ::MeanZ

        normA = 1d0 / Dsqrt( 2d0 * pi * sigmaX**2 )
        normB = 1d0 / Dsqrt( 2d0 * pi * sigmaY**2 )
        normC = 1d0 / Dsqrt( 2d0 * pi * sigmaZ**2 )

        MeanX = Ceiling( Dble( Nx ) / 2d0 )
        MeanY = Ceiling( Dble( Ny ) / 2d0 )
        MeanZ = Ceiling( Dble( Nz ) / 2d0 )

        Allocate( Gauss( 1:Nx , 1:Ny , 1:Nz ) )

        Do k = 1 , Nz
          Do j = 1 , Ny
            Do i = 1 , Nx
              Gauss( i , j , k ) = normA * normB * normC * &
                        Dexp( -0.5d0 * ( Dble( i ) - MeanX )**2 / ( sigmaX**2 ) ) * &
                        Dexp( -0.5d0 * ( Dble( j ) - MeanY )**2 / ( sigmaY**2 ) ) * &
                        Dexp( -0.5d0 * ( Dble( k ) - MeanZ )**2 / ( sigmaZ**2 ) )
            End Do
          End Do
        End Do



      End Function ComputeGaussian3d


      Function ComputeGaussian2d( Nx , Ny , sigmaX , sigmaY ) Result( Gauss )

        Implicit None
        Integer ::Nx
        Integer ::Ny
        Real*8  ::sigmaX
        Real*8  ::sigmaY

        Real*8,allocatable,dimension( : , : ) ::Gauss

        Integer ::i
        Integer ::j

        Real*8  ::normA
        Real*8  ::normB

        Real*8  ::MeanX
        Real*8  ::MeanY

        Allocate( Gauss( 1:Nx , 1:Ny ) )


        normA = 1d0 / Dsqrt( 2d0 * pi * sigmaX**2 )
        normB = 1d0 / Dsqrt( 2d0 * pi * sigmaY**2 )

        MeanX = Ceiling( Dble( Nx ) / 2d0 )
        MeanY = Ceiling( Dble( Ny ) / 2d0 )


        Do j = 1 , Ny
          Do i = 1 , Nx
             Gauss( i , j ) = normA * normB * Dexp( -0.5d0 * ( Dble( i ) - MeanX )**2 / ( sigmaX**2 ) ) * &
                                              Dexp( -0.5d0 * ( Dble( j ) - MeanY )**2 / ( sigmaY**2 ) )
          End Do
        End Do
      End Function ComputeGaussian2d



      Function ComputeGaussian1d( Nx , sigmaX ) Result( Gauss )

        Implicit None
        Integer ::Nx
        Real*8  ::sigmaX

        Real*8,allocatable,dimension( : ) ::Gauss

        Integer ::i
        Real*8  ::normA
        Real*8  ::MeanX

        Allocate( Gauss( 1:Nx ) )

        MeanX = Dble( Nx ) / 2d0
        normA = Dsqrt( 1d0 / ( 2d0 * sigmaX**2 * pi ) )

        Do i = 1 , Nx
           Gauss( i ) = normA * Dexp( -0.5d0 * ( Dble( i ) - MeanX )**2 / ( sigmaX**2 ) )
        End Do
      End Function ComputeGaussian1d


      Subroutine ComputeGaussianSmearing3d( Dist3d , Nx , Ny , Nz , sigmaX , sigmaY , sigmaZ , work )

        Implicit None
        Real*8,allocatable,dimension( : , : , : ) ::Dist3d
        Real*8,allocatable,dimension( : , : , : ) ::work
        Integer ::Nx
        Integer ::Ny
        Integer ::Nz
        Real*8  ::sigmaX
        Real*8  ::sigmaY
        Real*8  ::sigmaZ
        Real*8,allocatable,dimension( : , : , : ) ::Gauss

        Integer ::i
        Integer ::j
        Integer ::k


        Integer ::Na
        Integer ::Nb
        Integer ::Nc

        Integer ::Nx2
        Integer ::Ny2
        Integer ::Nz2

        Real*8,allocatable,dimension( : , : , : ) ::tmp

        Allocate( Gauss( 1:Nx , 1:Ny , 1:Nz ) )
        Allocate( tmp( 1:Nx , 1:Ny , 1:Nz ) )
        Gauss = ComputeGaussian3d( Nx , Ny , Nz , sigmaX , sigmaY , sigmaZ )

        Na = Size( Dist3d , 1 )
        Nb = Size( Dist3d , 2 )
        Nc = Size( Dist3d , 3 )
       
        Nx2 = Floor( Nx / 2d0 )
        Ny2 = Floor( Ny / 2d0 )
        Nz2 = Floor( Nz / 2d0 )


        Allocate( work( 1:Na + Nx , 1:Nb + Ny , 1:Nc + Nz ) )


        work = 0d0

        Do k = 1 , Nc
          Do j = 1 , Nb
            Do i = 1 , Na
              work( i : i + Nx-1 , j : j + Ny-1 , k : k + Nz-1 )  =  &
                  work( i : i + Nx-1 , j : j + Ny-1 , k : k + Nz-1 ) + &
                        Gauss( : , : , : ) * Dist3d( i , j , k )
            End Do
          End Do
        End Do

        Deallocate( Gauss )


      End Subroutine ComputeGaussianSmearing3d



      !!!
      !!! routine adding array B to array A
      !!! not in use but was helpful for debugging
      Subroutine AddArray3d( A , B , Weight )

        Implicit None
        Real*8,dimension( : , : , : ) ::A
        Real*8,dimensioN( : , : , : ) ::B
        Real*8 ::Weight

        Integer ::i
        Integer ::j
        Integer ::k

        Do i = 1, Size( A , 3 )
          Do j = 1, Size( A , 2 )
            Do k = 1, Size( A , 1 )
              print*,A( k , j , i ),B( k , j , i ) ,Weight , B( k , j , i )*Weight
              A( k , j , i )  = A( k , j , i ) + B( k , j , i )*Weight
            End Do
          End Do
        End Do
      End Subroutine AddArray3d
  End Module GaussianSmearing



  !!
  !! this module should generate free file numbers
  !! it is not working yet
  !!
  Module FileIOTools

    Implicit None
    Contains
    Function GetFileNumber( start ) Result( fileUnit )

      Implicit None
      Integer,optional ::start
      Integer ::fileUnit
      Logical ::fileExist


      fileExist = .True.
      If ( present( start ) ) then
        fileUnit = start-1
      Else
        fileUnit = 20
      End If

      Do While ( .True. )
        fileUnit = fileUnit + 1
        INQUIRE( unit = fileUnit , exist = fileExist )
        If ( .not. fileExist ) Exit
      End Do
    End Function GetFileNumber


    Function GetNFileNumbers( N , start ) Result( fileUnit )

      Implicit None
      Integer ::N
      Integer,optional ::start
      Integer,dimension( 1:N ) ::fileUnit
      Logical,dimension( 1:N ) ::fileExist
      Integer ::fileUnitBegin

      Integer ::i

      fileExist = .False.
      If ( present( start ) ) then
        fileUnitBegin  =  start
      Else
        fileUnitBegin  =  20
      End If

      Do While ( .not. all( fileExist ) )
        fileExist = .False.
        Do i = 1 , N
          fileUnit( i ) = fileUnitBegin + ( i - 1 )
          INQUIRE( unit = fileUnit( i ) , exist = fileExist( i ) )
          print*,fileExist( i )
        End Do
        fileUnitBegin = FileUnitBegin + 1 
      End Do


    End Function GetNFileNumbers
  End Module FileIOTools




  !!! module can be used for signal processing
  !!! to obtain smoother signals 
  Module WindowingAnalysis

    Implicit None
    Interface GaussianFiltering
      module procedure GaussianFilter2d,GaussianFilter
    End Interface GaussianFiltering

    Real*8,parameter ::PI  =  2d0 * Dacos( 0d0 )
    Real*8,parameter ::SqrtPI  =  Dsqrt( 2d0 * Dacos( 0d0 ) )
    Real*8,parameter ::Half = 0.5d0

    Contains
    !!
    !!
    !! gaussian filter like carried out by
    !! routine GaussianFilter( Xaxis , Signal , omega0 )
    !! for a signal of the kind S(k,\omega)
    Subroutine GaussianFilter2d( dx , Signal , omega0 )

      Implicit None
      Real*8,intent( in )       ::dx
      Real*8,dimension( : , : ),intent( inout ) ::Signal
      Real*8,intent( in )  ::omega0
      Real*8,allocatable,dimension( : , : )  ::work

      Integer ::i
      Integer ::Nk

      Nk = Size( Signal , 1 )

      Allocate( work( 1:Size( Signal , 2 ) , 1:Nk ) )
      work = Transpose( Signal )
      Do i = 1 , Nk 
         Call GaussianFilter( dx , work( : , i ) , omega0 )
      End Do
      Signal = Transpose( work )
      Deallocate( work )
    End Subroutine GaussianFilter2d





    !!!
    !!! Apply gaussian filter
    !!! apply gaussian filter function 
    !!! of the form
    !!! \int_{\infty}^{\infty}R(\omega-\omega')f(\omega)d\omega'
    !!! R(\omega-\omega')=\frac{1}{\pi \omega_{0}}exp(-\frac{(\omega-\omega')^{2}}{\omega_{0}})
    !!! omega0 is a paramteer that has to be supplied
    Subroutine GaussianFilter( dx , Signal , omega0 )

      use trapez_int, only :trapezoid
      Implicit None
      Real*8,intent( in ) ::dx                             !! grid spacing of signal
      Real*8,dimension( : ),intent( inout )  ::Signal      !! Signal on 2d grid
      Real*8,intent( in )   ::omega0                       !! width of gaussian window
      Real*8,allocatable,dimension( : )  ::work            !! working array
      Real*8,allocatable,dimension( : )  ::GausOnGrid
      Real*8,allocatable,dimension( : )  ::Xaxis
      Integer ::N
      Real*8 ::DxHalf
      Real*8 ::x0
      Integer ::i 

      N =  Size( Signal ) 
      DxHalf = dx * Half

      Allocate( work( 1:N ) )
      Allocate( Xaxis( 1:N ) )
      Allocate( GausOnGrid( 1:N ) )

      Do i = 1 , N
        Xaxis( i )  =  dx * Dble( i ) - DxHalf
      End Do


      work = Signal
      Do i = 1 , N
         x0 = dx * Dble( i ) - DxHalf
         Call ComputeGauss( Xaxis , x0 , GausOnGrid , omega0 )
         Signal( i )  =  trapezoid( dx , work * GausOnGrid ) 
      End Do
      Contains

        Subroutine ComputeGauss( x , x0 , Func , omega0 )

          Implicit None
          Real*8,dimension( : )  ::x
          Real*8  ::x0
          Real*8,dimension( 1:Size( x ) )  ::Func
          Real*8  ::omega0
          Func  =  Dexp( -( ( x - x0 )**2 ) / omega0**2 ) / SqrtPI / omega0**2 
        End Subroutine ComputeGauss
    End Subroutine GaussianFilter
  End Module WindowingAnalysis



