

  Module SplitCelltoBox


          Type BoxList
               Integer,allocatable,dimension( : )  ::List
               Integer ::Natoms
          End Type

          Interface MakeBoxList
                  module procedure MakeLinkListAtomtoBox,&
                                   MakeLinkListAtomtoBoxCart
          End Interface




          Contains

          !! return link list assigining atoms to boxes
          !! do comutation in cartesian coords
          !!
          Function MakeLinkListAtomtoBoxCart( Boxes , Atoms , lattice ) Result( List )
             Implicit None
             Real*8,dimension( : , : ) ::Boxes
             Real*8,dimension( : , : ) ::Atoms
             Real*8,dimension( : , : ) ::lattice

             Integer ::i
             Integer ::j
             Integer ::counter

             Type( BoxList ),allocatable,dimension( : ) ::List

             Integer,allocatable,dimension( : ) ::TempList

             Real*8,dimensioN( 1:Size( Atoms , 1 ) ) ::TempAtom
             Real*8,dimensioN( 1:Size( Atoms , 1 ) ) ::TempBoxMin
             Real*8,dimensioN( 1:Size( Atoms , 1 ) ) ::TempBoxMax

             Integer ::BoxCounter
             
             
             
             Allocate( List( 1:Size( Boxes , 2 ) ) )
             
             Do i = 1 , Size( Boxes , 2 )
                Allocate( TempList( 1:Size( Atoms , 2 ) ) )
                BoxCounter = 1
                Do j = 1, Size( Boxes , 1 ) , 2
                   TempBoxMin( BoxCounter ) = Boxes( j , i )
                   TempBoxMax( BoxCounter ) = Boxes( j+1 , i )
                   BoxCounter = BoxCounter + 1
                End Do
                TempBoxMin = Matmul( lattice , TempBoxMin )
                TempBoxMax = Matmul( lattice , TempBoxMax )
                counter = 1
                Do j = 1 , Size( Atoms , 2 )
                   TempAtom = MatMul( lattice , Atoms( : , j ) )
                   If ( TempAtom( 1 ) .ge. TempBoxMin( 1 ) .and. &
                        TempAtom( 1 ) .lt. TempBoxMax( 1 ) .and. & 
                        TempAtom( 2 ) .ge. TempBoxMin( 2 ) .and. &
                        TempAtom( 2 ) .lt. TempBoxMax( 2 ) .and. &
                        TempAtom( 3 ) .ge. TempBoxMin( 3 ) .and. &
                        TempAtom( 3 ) .lt. TempBoxMax( 3 ) ) then
                   !----------------------------------------
                      TempList( counter )  = j
                      counter = counter + 1
                   End If
                End Do
                List( i )%Natoms  =  counter - 1
                Allocate( List( i )%List( 1 : List( i )%Natoms ) )
                Do j = 1 , List( i )%Natoms
                   List( i )%List( j ) = TempList( j )
                End Do
                Deallocate( TempList )
             End Do

          End Function MakeLinkListAtomtoBoxCart






          !! return link list assigining atoms to boxes
          !!
          !!
          Function MakeLinkListAtomtoBox( Boxes , Atoms ) Result( List )

             Implicit None
             Real*8,dimension( : , : ) ::Boxes
             Real*8,dimension( : , : ) ::Atoms

             Integer ::i
             Integer ::j
             Integer ::counter

             Type( BoxList ),allocatable,dimension( : ) ::List

             Integer,allocatable,dimension( : ) ::TempList

            
             Allocate( List( 1:Size( Boxes , 2 ) ) )


             Do i = 1 , Size( Boxes , 2 )
                Allocate( TempList( 1:Size( Atoms , 2 ) ) )
                counter = 1
                Do j = 1 , Size( Atoms , 2 )
                   If ( Atoms( 1 , j ) .ge. Boxes( 1 , i ) .and. &
                        Atoms( 1 , j ) .lt. Boxes( 2 , i ) .and. & 
                        Atoms( 2 , j ) .ge. Boxes( 3 , i ) .and. &
                        Atoms( 2 , j ) .lt. Boxes( 4 , i ) .and. &
                        Atoms( 3 , j ) .ge. Boxes( 5 , i ) .and. &
                        Atoms( 3 , j ) .lt. Boxes( 6 , i ) ) then
                   !----------------------------------------
                      TempList( counter )  = j
                      counter = counter + 1
                   End If
                End Do
                List( i )%Natoms  =  counter - 1
                Allocate( List( i )%List( 1 : List( i )%Natoms ) )
                Do j = 1 , List( i )%Natoms
                   List( i )%List( j ) = TempList( j )
                End Do
                Deallocate( TempList )
             End Do



         End Function MakeLinkListAtomtoBox



         !!!
         !!! returns array containing
         !!! minX , maxX , minY , maxY , minZ , maxZ
         !!! in the first dimension of the array
         !!! the second index numbers the boxes
         !!! Supply N1 , N2 , N3 as the number of boxes
         !!! in a certain direction
         !!! the box width is 1/Ni
         !!! ________________________________________
         !!! |         |         |         |         |
         !!! |         |         |         |         |
         !!! |         |         |         |         |
         !!! |---------|---------|---------|---------|
         !!! |         |         |         |         |
         !!! |         |         |         |         |
         !!! |_________|_________|_________|_________|
         !!!
         Function MakeBoxes( N1 , N2 , N3 ) Result( Boxes )

             Implicit None
             Integer ::N1
             Integer ::N2
             Integer ::N3

             Integer ::i
             Integer ::j
             Integer ::k
             
             Integer ::col

             Real*8,allocatable,dimension( : , : ) ::Boxes

             Real*8 ::dx
             Real*8 ::dy
             Real*8 ::dz


             dx = 1d0 / Dble( N1 )
             dy = 1d0 / Dble( N2 )
             dz = 1d0 / Dble( N3 )

             col = 1
             Allocate( Boxes( 1:6 , 1 : N1 * N2 * N3 ) )
             Do i = 1 , N1
               Do j = 1 , N2
                  Do k = 1 , N3
                     Boxes( 1 , col ) = Dble( i - 1 ) * dx
                     Boxes( 2 , col ) = Dble( i ) * dx
                     Boxes( 3 , col ) = Dble( j - 1 ) * dy
                     Boxes( 4 , col ) = Dble( j ) * dy
                     Boxes( 5 , col ) = Dble( k - 1 ) * dz
                     Boxes( 6 , col ) = Dble( k ) * dz
                     col = col +1
                  End Do
               End Do
             End Do
         End Function MakeBoxes

  End Module SplitCelltoBox



  Module link_lists

    Implicit None
    Interface link_list
      module procedure compute_link_list , &
                       compute_link_list_cart , &
                       compute_link_list_NN , &
                       compute_link_list_self_AVN,&
                       compute_link_list_self_PB_AVN,&
                       compute_link_list_to_point,&
                       compute_link_list_to_point_cart,&
                       compute_link_list_self_AVN_cartesian,&
                       ComputeLinkList_NNN_cart,&
                       LinkAtoBthentoC,&
                       LinkAtoBthentoCCart
    End Interface link_list
    Interface NeighborCounter
       module procedure CountNumberInteractionPartners
    End Interface NeighborCounter

    Interface UnrollList
       module procedure UnrollList2ColFrom
    End interface UnrollList

    Contains

    !Function returns array with the indices of the nearest negighbours
    !atomsA found in atomsB only the first nearest neighbour is used
    !Function compute_link_list( atomsA , atomsB ) Result( link_list )
    Function compute_link_list( atomsA , atomsB ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Integer,dimension( 1 : Size( atomsA , 2 ) ) ::link_list
      Integer  ::i
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NA
      Integer  ::NB

      NA = Size( atomsA , 2 )
      NB = Size( atomsB , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NB ) )
      Do i = 1 , NA
         Do j = 1 , NB
            dist_list( 1 , j )  =  get_norm( atomsA( : , i ) - get_nearest_image( atomsA( : , i ) , atomsB( : , j ) ) )
            dist_list( 2 , j )  =  j
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 )  !j.L.21.11.2018
         link_list( i )  =  NINT( dist_list( 2 , 1 ) )
      End Do
      Deallocate( dist_list )
    End Function compute_link_list

    !Function returns array with the indices of the nearest negighbours
    !atomsA found in atomsB only the first nearest neighbour is used
    !Function compute_link_list( atomsA , atomsB ) Result( link_list )
    Function compute_link_list_cart( atomsA , atomsB , lattice ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image_cartesian
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Real*8,dimension( : , : )  ::lattice
      Integer,dimension( 1 : Size( atomsA , 2 ) ) ::link_list

      Real*8,dimensioN( 1:Size( atomsA , 1 ) ) ::TempA
      Real*8,dimensioN( 1:Size( atomsA , 1 ) ) ::TempB
      Integer  ::i
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NA
      Integer  ::NB

      NA = Size( atomsA , 2 )
      NB = Size( atomsB , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NB ) )

      Do i = 1 , NA
         TempA = MatMul( lattice , atomsA( : , i ) )
         Do j = 1 , NB
            TempB =  get_nearest_image_cartesian( atomsA( : , i ) , &
                                             atomsB( : , j ) , lattice )
            dist_list( 1 , j )  =   get_norm( TempA - TempB )
            dist_list( 2 , j )  =  j
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 ) !j.L. 21.11.2018
         link_list( i )  =  NINT( dist_list( 2 , 1 ) )
      End Do
      Deallocate( dist_list )
    End Function compute_link_list_cart

    !Function returns array with the indices of the first NN nearest negighbours
    !atomsA found in atomsB only the first nearest neighbour is used
    !Function compute_link_list( atomsA , atomsB ) Result( link_list )
    Function compute_link_list_NN( atomsA , atomsB , NN ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image
      use MergeSorter , only :Sort
      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Integer  ::NN
      Integer,dimension( 1 : NN , 1 : Size( atomsA , 2 ) ) ::link_list
      Integer  ::i
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NA
      Integer  ::NB

      NA = Size( atomsA , 2 )
      NB = Size( atomsB , 2 )


      Allocate( dist_list( 1 : 2 , 1 : NB ) )
      Do i = 1 , NA
         Do j = 1 , NB
            dist_list( 1 , j )  =  get_norm( atomsA( : , i ) - get_nearest_image( atomsA( : , i ) , atomsB( : , j ) ) )
            dist_list( 2 , j )  =  j
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 )  !j.L. 21.11.2018
         Do j = 1 , NN
           link_list( j , i )  =  NINT( dist_list( 2 , j ) )
         End Do
      End Do
      Deallocate( dist_list )
    End Function compute_link_list_NN


    !Function returns array with the indices of the first NN nearest negighbours
    !atomsA found in atomsB with self avoidance 
    !Function compute_link_list( atomsA , atomsB ) Result( link_list )
    Function compute_link_list_self_AVN( atoms , NN ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image
      use MergeSorter, only:Sort
      Implicit None
      Real*8,dimension( : , : )  ::atoms
      Integer  ::NN
      Integer,dimension( 1 : NN , 1 : Size( atoms , 2 ) ) ::link_list
      Integer  ::i
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NA

      NA  =  Size( atoms , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NA * 27 ) )
      Do i = 1 , NA
         dist_list  =  Huge( dist_list( 1 , 1 ) )
         Do j = 1 , NA
            If ( i .eq. j ) then
              dist_list( 1 , j )  =  Huge( dist_list( 1 , j ) )
              Cycle
            End If
            dist_list( 1 , j )  =  get_norm( atoms( : , i ) - get_nearest_image( atoms( : , i ) , atoms( : , j ) ) )
            dist_list( 2 , j )  =  j
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 ) !j.L.21.11.2018
         Do j = 1 , NN
           link_list( j , i )  =  NINT( dist_list( 2 , j ) )
         End Do
      End Do
      Deallocate( dist_list )
    End Function compute_link_list_self_AVN


    !Function returns array with the indices of the first NN nearest negighbours
    !of atoms only the first NN nearest neighbours are used, self avoidance is included
    !Function compute_link_list( atomsA , atomsB ) Result( link_list )
    Function compute_link_list_self_AVN_cartesian( atoms , NN , lattice ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image_cartesian
      use MergeSorter, only:Sort
      Implicit None
      Real*8,dimension( : , : )  ::atoms
      Integer  ::NN
      Real*8,dimension( : , : ) ::lattice
      Integer,dimension( 1 : NN , 1 : Size( atoms , 2 ) ) ::link_list
      Integer  ::i
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NA
      Real*8,dimension( 1:3 ) ::tempA
      Real*8,dimension( 1:3 ) ::tempB

      NA  =  Size( atoms , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NA ) )
      Do i = 1 , NA
         tempA  =  Matmul( lattice , atoms( :, i ) )
         Do j = 1 , NA
            If ( i .eq. j ) then
              dist_list( 1 , j )  =  Huge( dist_list( 1 , j ) )
              Cycle
            End If
            tempB = get_nearest_image_cartesian( atoms( : , i ) , atoms( : , j ) , lattice )
            dist_list( 1 , j )  =  get_norm( tempA - tempB )
            dist_list( 2 , j )  =  j
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 )  !j.L. 21.11.2018
         Do j = 1 , NN
           link_list( j , i )  =  NINT( dist_list( 2 , j ) )
         End Do
      End Do

      Deallocate( dist_list )
    End Function compute_link_list_self_AVN_cartesian


    !Returns the first NN nearest neighbours out of list atomsB
    ! with respect to atomsA
    !                  atomsB3
    !                    |
    !    atomsB2 ----- atomA ---- atomsB1
    !                    |
    !                  atomsB4
    !The routine takes direct coordinates as an input
    Function ComputeLinkList_NNN_cart( atomsA , atomsB , lattice , NN ) Result( List )

      use tool_interfaces,only :get_nearest_image_cartesian,get_norm
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Real*8,dimension( : , : )  ::lattice
      Integer ::NN
      Real*8,allocatable,dimension( : , : ) ::distList
      Integer,allocatable,dimension( : , : )  ::List
      Real*8,dimension( 1:Size( atomsA , 1 ) ) ::TempPos
      Real*8,dimension( 1:Size( atomsA , 1 ) ) ::TempPosA
      Real*8 ::dist

      Integer ::i
      Integer ::j

      If ( .not. Allocated( List ) ) Allocate( List( 1:NN, 1:Size( atomsA , 2 ) ) )
      Allocate( distList( 1:2 , 1:Size( atomsB , 2 ) ) )

      Do i = 1 , Size( atomsA , 2 )
        TempPosA = MatMul( lattice , atomsA( : , i ) )
        Do j = 1, Size( atomsB , 2 )
           TempPos = get_nearest_image_cartesian( atomsA( : , i ) , atomsB( : , j ) , lattice )
           dist = get_norm( TempPos - TempPosA )
           distList( 1 , j ) = dist
           distList( 2 , j ) = Dble( j )
        End Do
        !Call Sort_II( distList , Size( distList , 2 ) , Size( distList , 1 ) , 1 )
        Call Sort( distList , 1 , Size( distList , 2 ) , 1 ) !j.L 21.11.2018
        Do j = 1, NN
           List( j , i ) = NINT( distList( 2 , j ) )
        End Do
      End Do

      Deallocate( distList )


    End Function





    !Function returns array with the indices of the first NN nearest negighbours
    !atoms, first nearest neigbhours are assigned to the first Nats atoms
    !
    Function compute_link_list_self_PB_AVN( atoms , Nats , NN ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image
      use MergeSorter, only:Sort
      Implicit None
      Real*8,dimension( : , : )  ::atoms
      Integer  ::NN
      Integer  ::Nats
      Integer,dimension( 1 : NN , 1 : Nats ) ::link_list
      Integer  ::i
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NA

      NA  =  Size( atoms , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NA ) )

      Do i = 1 , Nats
         Do j = 1 , NA
            If ( i .eq. j ) then
              dist_list( 1 , j )  =  Huge( dist_list( 1 , j ) )
              Cycle
            End If
            dist_list( 1 , j )  =  get_norm( atoms( 1 : 3 , i ) - atoms( 1 : 3 , j ) )
            dist_list( 2 , j )  =  atoms( 4 , j )
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 )  !j.L. 21.11.2018
         Do j = 1 , NN
           link_list( j , i )  =  NINT( dist_list( 2 , j ) )
         End Do
      End Do
      Deallocate( dist_list )
    End Function compute_link_list_self_PB_AVN


    !Function returns array with the indices of the first NN nearest negighbours
    !of atoms according to a given point called vector
    Function compute_link_list_to_point( center , atoms , NN ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image
      use MergeSorter, only:Sort
      Implicit None
      Real*8,dimension( : )  ::center
      Real*8,dimension( : , : )  ::atoms
      Integer  ::NN
      Integer,dimension( 1 : NN ) ::link_list
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NB

      NB = Size( atoms , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NB ) )
      Do j = 1 , NB
         dist_list( 1 , j )  =  get_norm( center( : ) - get_nearest_image( center( : ) , atoms( : , j ) ) )
         dist_list( 2 , j )  =  j
      End Do

      !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
      Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 ) !j.L. 21.11.2018
      Do j = 1 , NN
        link_list( j )  =  NINT( dist_list( 2 , j ) )
      End Do

      Deallocate( dist_list )
    End Function compute_link_list_to_point


    !Function returns array with the indices of the first NN nearest negighbours
    !of atoms according to a given point called vector
    Function compute_link_list_to_point_cart( center , atoms , NN , lattice ) Result( link_list )

      use tool_interfaces, only :get_norm,get_nearest_image_cartesian
      use MergeSorter,only :Sort
      Implicit None
      Real*8,dimension( : )  ::center
      Real*8,dimension( : , : )  ::atoms
      Integer  ::NN
      Integer,dimension( 1 : NN ) ::link_list
      Integer  ::j
      Real*8,allocatable,dimension( : , : ) ::dist_list
      Integer  ::NB
      Real*8,dimension( 1 : Size( center , 1 ) ) ::tempA
      Real*8,dimension( 1 : Size( center , 1 ) ) ::tempB
      Real*8,dimension( : , : ) ::lattice

      NB = Size( atoms , 2 )

      Allocate( dist_list( 1 : 2 , 1 : NB ) )
      tempA = Matmul( lattice , center )
      Do j = 1 , NB
         tempB  =  get_nearest_image_cartesian( center , atoms( : , j ) , lattice )
         dist_list( 1 , j )  =  get_norm( tempA - tempB )
         dist_list( 2 , j )  =  j
      End Do
      !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
      Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 )  !j.L. 21.11.2018
      Do j = 1 , NN
        link_list( j )  =  NINT( dist_list( 2 , j ) )
      End Do

      Deallocate( dist_list )
    End Function compute_link_list_to_point_cart



    ! Function returns a list containing the indices of the
    ! first NN nearest neighbours of atomsA with respect to atomsB
    ! no periodic boundary conditions are used
    ! atomsA 1:3 , 1:N array
    ! atomsB 1:3 , 1:N array
    ! lattice 3x3 bravais matrix
    ! NN number of nearest neighbours

    Function compute_link_list_noPB( atomsA , atomsB , NN ) Result( list )  !j.L 28.08.2018

      use tool_interfaces, only :get_norm
      use MergeSorter, only:Sort
      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Integer ::NN

      Real*8,allocatable,dimension( : , : )  ::dist_list
      Integer,allocatable,dimension( : , : )  ::list

      Integer ::i
      Integer ::j

      Allocate( dist_list( 1:2 , 1:size( atomsA , 2 ) ) )

      If ( .not. Allocated( list ) ) Allocate( list( 1:NN , 1:Size( atomsB , 2 ) ) )

      Do i = 1 , Size( atomsB , 2 )
         Do j = 1 , Size( atomsA , 2 )
            dist_list( 1 , j )  =  get_norm( atomsA( : , j ) - atomsB( : , i ) )
            dist_list( 2 , j )  =  j
         End Do
         !Call sort_II( dist_list , Size( dist_list , 2 ) , Size( dist_list , 1 ) , 1 )
         Call Sort( dist_list , 1 , Size( dist_list , 2 ) , 1 )  !j.L. 21.11.2018
         Do j = 1 , NN
            list( j , i )  =  NINT( dist_list( 2 , j ) )
         End Do
      End Do

      Deallocate( dist_list )
    End Function compute_link_list_noPB


    !
    ! Function counts the number of interaction partners
    ! within a supplied cutoff radius
    ! self interaction is avoided by a small number comparison
    ! routine uses direct atomic coordinates if the cutoff radius is supplied
    ! in direct coords supply unit matrix for lattice

    Function CountNumberInteractionPartners( atoms , Cutoff , lattice ) Result( NNN )

      use tool_interfaces, only :get_nearest_image,get_norm
      Implicit None
      Real*8,dimension( : , : )  ::atoms
      Real*8 ::cutoff
      Real*8,dimension( : , : ) ::lattice
      Integer ::i
      Real*8,dimension( 1:3 ) ::tempA
      Real*8,dimension( 1:3 ) ::tempB
      Real*8::dist
      Integer ::NNN

      !Measure is only taken for first atom
      NNN = 0
      tempA = Matmul( lattice , atoms( : , 1 ) )
      Do i = 2 , Size( atoms , 2 )
         tempB =  get_nearest_image( atoms( : , 1 ) , atoms( : , i ) )
         tempB =  Matmul( lattice , tempB )
         dist = get_norm( tempA - tempB )
         If ( dist .lt. Cutoff ) NNN = NNN + 1
      End Do
    End Function CountNumberInteractionPartners



    Function UnrollList2ColFrom( list_in ) Result( list_out )


       Implicit None
       Integer,dimension( : , : )  ::list_in
       Integer,Allocatable,dimension( : , : )  ::list_out
       Integer ::i
       Integer ::j
       Integer ::col

       If ( .not. Allocated( list_out ) ) then
          Allocate( list_out( 1:2 , 1 : Size( list_in , 1 ) * Size( list_in , 2 ) ) )
       End If

       col = 1
       Do i = 1 , Size( list_in , 2 )
         Do j = 1, Size( list_in , 1 )
            list_out( 1 , col ) = i
            list_out( 2 , col ) = list_in( j , i )
            col = col + 1
         End Do
       End Do
    End Function UnrollList2ColFrom


    Function LinkList_special( data , NN , MM ) Result( LinkList )

      use tool_interfaces, only :get_nearest_image,get_norm
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : )  ::data
      Integer ::NN
      Integer ::MM

      Integer ::i
      Integer ::j

      Integer,allocatable,dimension( : , : ) ::LinkList
      Real*8,allocatable,dimension( : , : ) ::TempData

      Allocate( LinkList( 1:NN , 1:MM ) )
      Allocate( TempData( 1:2 , 1:Size( data , 2 ) ) )

      Do i = 1, MM
        Do j = 1 , Size( data , 2 )
          If ( NINT( data( 4 , i ) ) .eq. NINT( data( 4 , j ) ) ) then
             TempData( 1 , j )  =  Huge( TempData( 1 , j ) )
             Cycle
          End If
          TempData( 1 , j ) = get_norm( get_nearest_image( data( 1:3 , i ) , data( 1:3 , j ) )&
                                  - data( 1:3 , i ) )
          TempData( 2 , j ) = j
        End Do
        Call Sort( TempData , 1 , Size( TempData , 2 ) , 1 )
        Do j = 1, NN
          LinkList( j , i ) = NINT( TempData( 2 , j ) )
        End Do
      End Do


    End Function



    !!
    !! function computes a link list.
    !! first it links the nearest NAB atoms of type A to B
    !! then NBC atoms of type B are linked to type C
    !! returned is a link list of NAB atoms linke to atom type C
    !! sounds silly ;-)
    !! an example where this is useful is the following molecule
    !!
    !!         H
    !!         |
    !!  H--N---C---N--H
    !!     |       |
    !!     H       H
    !!  and you want to have all nitrogen hydrogens that belong to the same
    !! FA molecule/ion
    !!
    Function LinkAtoBthentoC( atomsA , atomsB , atomsC , NAB , NBC ) Result( LinkList ) !j.L 07.11.2019

      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Real*8,dimension( : , : )  ::atomsC
      Integer ::NAB
      Integer ::NBC

      Integer,allocatable,dimension( : , : ) ::LinkList
      Integer,allocatable,dimension( : , : ) ::TempListAB
      Integer,allocatable,dimension( : , : ) ::TempListBC

      Integer ::i
      Integer ::j
      Integer ::k
      Integer ::col

      !! Linking atomA to atomB
      TempListAB = compute_link_list_NN( atomsB , atomsA , NAB )
      !! Linking atomsB to atomsC
      TempListBC = compute_link_list_NN( atomsC , atomsB , NBC )

      Allocate( LinkList( 1:NAB*NBC , Size( atomsC , 2 ) ) )
      Do i = 1, Size( atomsC , 2 )
        col = 1
        Do k = 1 , NBC
          Do j = 1, NAB
             LinkList( col , i ) = TempListAB( j , TempListBC( k , i ) )
             col = col + 1
          End Do
        End Do
      End Do
      Deallocate( TempListAB , TempListBC )
    End Function LinkAtoBthentoC

    !! routine dooes the same as LinkAtoBthentoC but with cartesian coords
    !! function computes a link list.
    !! first it links the nearest NAB atoms of type A to B
    !! then NBC atoms of type B are linked to type C
    !! returned is a link list of NAB atoms linke to atom type C
    !! sounds silly ;-)
    !! an example where this is useful is the following molecule
    !!
    !!         H
    !!         |
    !!  H--N---C---N--H
    !!     |       |
    !!     H       H
    !!  and you want to have all nitrogen hydrogens that belong to the same
    !! FA molecule/ion
    !!
    Function LinkAtoBthentoCCart( atomsA , atomsB , atomsC , NAB , NBC , lattice ) Result( LinkList ) !j.L 07.11.2019

      Implicit None
      Real*8,dimension( : , : )  ::atomsA
      Real*8,dimension( : , : )  ::atomsB
      Real*8,dimension( : , : )  ::atomsC
      Integer ::NAB
      Integer ::NBC
      Real*8,dimension( : , : ) ::lattice

      Integer,allocatable,dimension( : , : ) ::LinkList
      Integer,allocatable,dimension( : , : ) ::TempListAB
      Integer,allocatable,dimension( : , : ) ::TempListBC

      Integer ::i
      Integer ::j
      Integer ::k
      Integer ::col

      !! Linking atomA to atomB
      TempListAB = ComputeLinkList_NNN_cart( atomsB , atomsA , lattice , NAB )
      !! Linking atomsB to atomsC
      TempListBC = ComputeLinkList_NNN_cart( atomsC , atomsB , lattice , NBC )

      Allocate( LinkList( 1:NAB*NBC , Size( atomsC , 2 ) ) )
      Do i = 1, Size( atomsC , 2 )
        col = 1
        Do k = 1 , NBC
          Do j = 1, NAB
             LinkList( col , i ) = TempListAB( j , TempListBC( k , i ) )
             col = col + 1
          End Do
        End Do
      End Do
      Deallocate( TempListAB , TempListBC )
    End Function LinkAtoBthentoCCart




  End Module link_lists



  Module xyzList

      Interface ComputeXYZList
          module procedure Integer_grid_robust,XYZIntegerGrid,&
                           Integer_grid_robustAB
      End Interface ComputeXYZList
      Interface ComputeABCList
          module procedure Integer_grid_robustGD,&
                           Integer_grid_robustABGD
      End Interface ComputeABCList
      Contains

        Function XYZIntegerGrid( atoms , lattice , STD_routine , PrintList , Nx , Ny , Nz ) Result ( xyzList )

          use IntegerGrid, only :ComputeIntegerGrid
          Implicit None
          Real*8,dimension( : , : )  ::atoms
          Real*8,dimension( : , : )  ::lattice

          Logical ::STD_routine
          Logical ::PrintList

          Integer ::Nx
          Integer ::Ny
          Integer ::Nz

          Integer,allocatable,dimension( : , : , : )  ::xyzList
          Integer,allocatable,dimension( : , : )  ::TempList
          Integer,allocatable,dimension( : , : , : ) ::grid

          Integer ::i
          Integer ::j
          Integer ::k

          Integer ::x
          Integer ::y
          Integer ::z

          Integer ::Nidx1p
          Integer ::Nidx1m
          Integer ::Nidx2p
          Integer ::Nidx2m
          Integer ::Nidx3p
          Integer ::Nidx3m


          If ( .not. Allocated( xyzList ) ) &
                 Allocate( xyzList( 1:2 , 1:3 , 1:Size( atoms , 2 ) ) )

          !!!
          !!! computing integer grid by projection
          !!! divide direct coords by cells width in certain direction
          !!! and taking the NINT
          If ( STD_routine ) then
             Allocate( TempList( 1 : Size( atoms , 1 ) , 1 : Size( atoms , 2 ) ) )
             Write( * , * ) "Making XYZ List with Integer Grid Projection"
             Write( * , * ) "Routine should be used for Regular Grids"
             Write( * , * ) "Use this routine carefully"
             TempList = ComputeIntegerGrid( atoms , Nx , Ny , Nz , lattice )
             Allocate( grid( 1:Nx , 1:Ny , 1:Nz ) )
             Do i = 1, Size( atoms , 2 )
                x = TempList( 1 , i )
                y = TempList( 2 , i )
                z = TempList( 3 , i )
                grid( x , y , z )  = i
             End Do
             Deallocate( TempList )

             Do i = 1, Nx
                Do j = 1, Ny
                   Do k = 1 , Nz
                      If ( i .eq. 1 ) then
                        Nidx1p = grid( i + 1 , j , k )
                        Nidx1m = grid( Nx , j , k )
                      Else If ( i .eq. Nx ) then
                        Nidx1p = grid( 1 , j , k )
                        Nidx1m = grid( i-1 , j , k )
                      Else
                        Nidx1p = grid( i + 1 , j , k )
                        Nidx1m = grid( i-1 , j , k )
                      End If
                      xyzList( 1 , 1 , grid( i, j , k ) )  =  Nidx1p
                      xyzList( 2 , 1 , grid( i, j , k ) )  =  Nidx1m
                      If ( j .eq. 1 ) then
                        Nidx2p = grid( i , j + 1, k )
                        Nidx2m = grid( i , Ny , k )
                      Else If ( j .eq. Ny ) then
                        Nidx2p = grid( i , 1 , k )
                        Nidx2m = grid( i , j-1 , k )
                      Else
                        Nidx2p = grid( i , j+1 , k )
                        Nidx2m = grid( i , j-1 , k )
                      End If
                      xyzList( 1 , 2 , grid( i, j , k ) )  =  Nidx2p
                      xyzList( 2 , 2 , grid( i, j , k ) )  =  Nidx2m
                      If ( k .eq. 1 ) then
                        Nidx3p = grid( i , j , k+1 )
                        Nidx3m = grid( i , j , Nz )
                      Else If ( k .eq. Nz ) then
                        Nidx3p = grid( i , j , 1 )
                        Nidx3m = grid( i , j , k-1 )
                      Else
                        Nidx3p = grid( i , j , k+1 )
                        Nidx3m = grid( i , j , k-1 )
                      End If
                      xyzList( 1 , 3 , grid( i, j , k ) )  =  Nidx3p
                      xyzList( 2 , 3 , grid( i, j , k ) )  =  Nidx3m
                   End Do
                End Do
             End Do
             Deallocate( grid )
          Else
            !!! using a distance dependent approach
            !!! 1) Collect first nearest neighbors by sorting algorithm
            !!! 2) Assign those first 6 nearest neighbours to xyz
            !!!    directions by computing dot products with those directions and
            !!!    again exploiting a sort algorithm to assign bonds to directions
            Write( * , * ) "Making distance dependent XYZ List"
            Write( * , * ) "Sorting atoms according to distance"
            Write( * , * ) "Computing DotProducts with x,y,z direction and sorting"
            Write( * , * ) "to assign nearest neighbours to directions"
            xyzList  = Integer_grid_robust( atoms , lattice )
          End If

          If ( PrintList ) then
            Write( * , * ) "Integer List, can be used as a safety check"
            Write( * , "(8A7)" ) "Atom","x+" , "x-" , "y+" , "y-" , "z+" , "z-"
            Do i = 1 , Size( xyzList , 3 )
               Write( * , "(10I7)" )  i,xyzList( : , 1 , i ) , xyzList( : , 2 , i ) , xyzList( : , 3 , i )
            End Do
          End If
        End Function XYZIntegerGrid



        !!
        !!  Routine computes xyz nearest neighbour lists under the help of sorting alogirthms
        !!  first the routine assumes a coordination of your atoms by 6 nearest neighbours (nn)
        !!  those are computed according to a computing distances and using the 6 nearest atoms
        !!  next the dot products of the bonding vectors with the xyz directions are computed
        !!  next those dot products are again sorted according to their magnitude
        !!  to assign the bonds to certain directions
      Function Integer_grid_robust( atoms , lattice ) Result( TempListDir )

        use tool_interfaces, only:get_nearest_image_cartesian,get_norm
        use link_lists, only :link_list
        use MergeSorter, only :Sort
        Implicit None
        Real*8,dimension( : , : )  ::atoms
        Real*8,dimension( : , : )  ::lattice
        Real*8,Save,dimension( 1:3 , 1:3 ) ::Dirs = reshape( [ 1, 0, 0, 0,&
          1, 0, 0, 0, 1 ], shape( Dirs ) )

        Integer,allocatable,dimension( : , : , : )  ::TempListDir
        Integer,allocatable,dimension( : , : )  ::TempList
        Integer ::i
        Integer ::j
        Integer ::k
        Real*8,allocatable,dimension( : , : ) ::DotPs
        Real*8,dimension( 1:3 )  ::Image
        Real*8,dimension( 1:3 )  ::AtomC
        Integer ::Indx

        Allocate( TempList( 1:6 , 1: Size( atoms , 2 ) ) )
        TempList = link_list( atoms , 6 , lattice )
        Allocate( TempListDir( 1:2 , 1:3 , 1 : Size( atoms , 2 ) ) )

        Allocate( DotPs( 1:4 , 1:6 ) )

        Do i = 1, Size( atoms , 2 )
           AtomC = MatMul( lattice , atoms( : , i ) )
           Do j = 1 , 3
              Do k = 1 , 6
                 Indx = TempList( k , i )
                 Image = get_nearest_image_cartesian( atoms( : , i ) , atoms( : , Indx ) , lattice )
                 Image = Image - AtomC
                 Dotps( j , k )  =  Dot_product( Image , Dirs( : , j ) ) / get_norm( Image )
                 Dotps( 4 , k ) = Indx
              End Do
           End Do
           Do j = 1 , 3
              Call Sort( DotPs , 1 , 6 , j )
              TempListDir( 1 , j , i )  =  NINT( DotPs( 4 , 6 ) )
              TempListDir( 2 , j , i )  =  NINT( DotPs( 4 , 1 ) )
           End Do
        End Do
        Deallocate( TempList )

      End Function Integer_grid_robust


        !!
        !!  Routine computes xyz nearest neighbour lists under the help of sorting alogirthms
        !!  first the routine assumes a coordination of your atoms by 6 nearest neighbours (nn)
        !!  those are computed according to a computing distances and using the 6 nearest atoms
        !!  next the dot products of the bonding vectors with the xyz directions are computed
        !!  next those dot products are again sorted according to their magnitude
        !!  to assign the bonds to certain directions
        !!  the only difference to Integer_grid_robust is that one assumes a central atomA
        !!  and sourrounds them by atomsB
        !!               B
        !!               |
        !!               |
        !!       B ----- A ----- B
        !!               |
        !!               |
        !!               B
      Function Integer_grid_robustAB( atomsA , atomsB , lattice ) Result( TempListDir )

         use tool_interfaces, only:get_nearest_image_cartesian,get_norm
         use link_lists, only :link_list
         use MergeSorter, only :Sort
         Implicit None
         Real*8,dimension( : , : )  ::atomsA
         Real*8,dimension( : , : )  ::atomsB
         Real*8,dimension( : , : )  ::lattice
         Real*8,Save,dimension( 1:3 , 1:3 ) ::Dirs = reshape( [ 1, 0, 0, 0,&
           1, 0, 0, 0, 1 ], shape( Dirs ) )

         Integer,allocatable,dimension( : , : , : )  ::TempListDir
         Integer,allocatable,dimension( : , : )  ::TempList
         Integer ::i
         Integer ::j
         Integer ::k
         Real*8,allocatable,dimension( : , : ) ::DotPs
         Real*8,dimension( 1:3 )  ::Image
         Real*8,dimension( 1:3 )  ::AtomC
         Integer ::Indx

         Allocate( TempList( 1:6 , 1: Size( atomsA , 2 ) ) )
         TempList = link_list( atomsA , atomsB , lattice , 6 )
         Allocate( TempListDir( 1:2 , 1:3 , 1 : Size( atomsA , 2 ) ) )

         Allocate( DotPs( 1:4 , 1:6 ) )

         Do i = 1, Size( atomsA , 2 )
            AtomC = MatMul( lattice , atomsA( : , i ) )
            Do j = 1 , 3
               Do k = 1 , 6
                  Indx = TempList( k , i )
                  Image = get_nearest_image_cartesian( atomsA( : , i ) , atomsB( : , Indx ) , lattice )
                  Image = Image - AtomC
                  Dotps( j , k )  =  Dot_product( Image , Dirs( : , j ) ) / get_norm( Image )
                  Dotps( 4 , k ) = Indx
               End Do
            End Do
            Do j = 1 , 3
               Call Sort( DotPs , 1 , 6 , j )
               TempListDir( 1 , j , i )  =  NINT( DotPs( 4 , 6 ) )
               TempListDir( 2 , j , i )  =  NINT( DotPs( 4 , 1 ) )
            End Do
         End Do
         Deallocate( TempList )

      End Function Integer_grid_robustAB




        !!
        !!  Routine computes xyz nearest neighbour lists under the help of sorting alogirthms
        !!  first the routine assumes a coordination of your atoms by 6 nearest neighbours (nn)
        !!  those are computed according to a computing distances and using the 6 nearest atoms
        !!  next the dot products of the bonding vectors with the xyz directions are computed
        !!  next those dot products are again sorted according to their magnitude
        !!  to assign the bonds to certain directions
      Function Integer_grid_robustGD( atoms , lattice , Dirs ) Result( TempListDir )

        use tool_interfaces, only:get_nearest_image_cartesian,get_norm
        use link_lists, only :link_list
        use MergeSorter, only :Sort
        Implicit None
        Real*8,dimension( : , : )  ::atoms
        Real*8,dimension( : , : )  ::lattice
        Real*8,dimension( : , : )  ::Dirs

        Integer,allocatable,dimension( : , : , : )  ::TempListDir
        Integer,allocatable,dimension( : , : )  ::TempList
        Integer ::i
        Integer ::j
        Integer ::k
        Real*8,allocatable,dimension( : , : ) ::DotPs
        Real*8,dimension( 1:3 )  ::Image
        Real*8,dimension( 1:3 )  ::AtomC
        Integer ::Indx

        Allocate( TempList( 1:6 , 1: Size( atoms , 2 ) ) )

        TempList = link_list( atoms , 6 , lattice )
        
        If ( .not. Allocated( TempListDir ) ) &
           Allocate( TempListDir( 1:2 , 1:3 , 1 : Size( atoms , 2 ) ) )

        Allocate( DotPs( 1:4 , 1:6 ) )


        Do i = 1, Size( atoms , 2 )
           AtomC = MatMul( lattice , atoms( : , i ) )
           Do j = 1 , Size( Dirs , 2 )
              Do k = 1 , 6
                 Indx   =  TempList( k , i )
                 Image  =  get_nearest_image_cartesian( atoms( : , i ) , atoms( : , Indx ) , lattice )
                 Image  =  Image - AtomC
                 Dotps( j , k )  =  Dot_product( Image , Dirs( : , j ) ) / get_norm( Image )
                 Dotps( 4 , k ) = Indx
              End Do
           End Do
           Do j = 1 , 3
              Call Sort( DotPs , 1 , 6 , j )
              TempListDir( 1 , j , i )  =  NINT( DotPs( 4 , 6 ) )
              TempListDir( 2 , j , i )  =  NINT( DotPs( 4 , 1 ) )
           End Do
        End Do
        Deallocate( TempList )

      End Function Integer_grid_robustGD


        !!
        !!  Routine computes xyz nearest neighbour lists under the help of sorting alogirthms
        !!  first the routine assumes a coordination of your atoms by 6 nearest neighbours (nn)
        !!  those are computed according to a computing distances and using the 6 nearest atoms
        !!  next the dot products of the bonding vectors with the xyz directions are computed
        !!  next those dot products are again sorted according to their magnitude
        !!  to assign the bonds to certain directions
        !!  the only difference to Integer_grid_robust is that one assumes a central atomA
        !!  and sourrounds them by atomsB
        !!               B
        !!               |
        !!               |
        !!       B ----- A ----- B
        !!               |
        !!               |
        !!               B
      Function Integer_grid_robustABGD( atomsA , atomsB , lattice , Dirs ) Result( TempListDir )

         use tool_interfaces, only:get_nearest_image_cartesian,get_norm
         use link_lists, only :link_list
         use MergeSorter, only :Sort
         Implicit None
         Real*8,dimension( : , : )  ::atomsA
         Real*8,dimension( : , : )  ::atomsB
         Real*8,dimension( : , : )  ::lattice
         Real*8,dimension( : , : )  ::Dirs

         Integer,allocatable,dimension( : , : , : )  ::TempListDir
         Integer,allocatable,dimension( : , : )  ::TempList
         Integer ::i
         Integer ::j
         Integer ::k
         Real*8,allocatable,dimension( : , : ) ::DotPs
         Real*8,dimension( 1:3 )  ::Image
         Real*8,dimension( 1:3 )  ::AtomC
         Integer ::Indx

         Allocate( TempList( 1:6 , 1: Size( atomsA , 2 ) ) )
         TempList = link_list( atomsA , atomsB , lattice , 6 )
         Allocate( TempListDir( 1:2 , 1:3 , 1 : Size( atomsA , 2 ) ) )

         Allocate( DotPs( 1:4 , 1:6 ) )


         Do i = 1, Size( atomsA , 2 )
            AtomC = MatMul( lattice , atomsA( : , i ) )
            Do j = 1 , Size( Dirs , 2 )
               Do k = 1 , 6
                  Indx = TempList( k , i )
                  Image = get_nearest_image_cartesian( atomsA( : , i ) , atomsB( : , Indx ) , lattice )
                  Image = Image - AtomC
                  Dotps( j , k )  =  Dot_product( Image , Dirs( : , j ) ) / get_norm( Image )
                  Dotps( 4 , k ) = Indx
               End Do
            End Do
            Do j = 1 , 3
               Call Sort( DotPs , 1 , 6 , j )
               TempListDir( 1 , j , i )  =  NINT( DotPs( 4 , 6 ) )
               TempListDir( 2 , j , i )  =  NINT( DotPs( 4 , 1 ) )
            End Do
         End Do
         Deallocate( TempList )

      End Function Integer_grid_robustABGD

  End Module xyzList


  Module GenBondTypes
    !!
    !! module contains general routines
    !! to determine atoms
    !! belonging to certain bond types
    !!

    Implicit None
    Interface GenBonds
      module procedure ListOfABBonds,ListOfAABonds,&
                       ListOfABABonds
    End Interface GenBonds

    Interface BondPosition
      module procedure ABBondPositions,BABBondPositions
    End Interface BondPosition

    Interface ComputeBond
      module procedure ComputeABBond,ComputeBABBond
    End Interface ComputeBond

    Contains

    !!
    !! routine checks which atoms of atomsA
    !! are connected to exactly 1 atom of atomsB
    !! within the cutoff criterion
    !! function returns a list
    !! first dimension contains index of atomsA
    !! and second index gives you the index pointing to atomB
    !! ato which A is connected
    !!
    !!
    !!
    !! initializes A-A bonds
    !! link list containing List( 1 , j ) -> index A
    !! and List( 2 , j ) -> index B
    !!
    Function ListOfABBonds( atomsA , atomsB , lattice ,cutoff ) Result( List )
      use tool_interfaces, only :get_nearest_image_cartesian,get_norm
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : ) ::atomsA
      Real*8,dimension( : , : ) ::atomsB
      Real*8,dimension( : , : ) ::lattice
      Real*8 ::cutoff

      Integer ::NN
      Integer ::NNN

      Integer ::i
      Integer ::j

      Real*8,allocatable,dimension( : , : ) ::DistList
      Real*8,dimension( 1:3 ) ::PosA
      Real*8,dimension( 1:3 ) ::PosB

      Integer,allocatable,dimension( : , : ) ::List
      Integer,allocatable,dimension( : , : ) ::ListTmp

      Allocate( DistList( 1:2 , 1:Size( atomsB , 2 ) ) )
      Allocate( ListTmp( 1:2 , 1:Size( atomsB , 2 ) ) )


      NNN = 0   !counter number of bonds satisifying crit


      Do i = 1, Size( atomsA , 2 )
        PosA  =  MatMul( lattice , atomsA( : , i ) )
        Do j = 1, Size( atomsB , 2 )
          PosB = get_nearest_image_cartesian( atomsA( : , i ) , atomsB( : , j ) , lattice )
          DistList( 1 , j ) = get_norm( PosA - PosB )
          DistList( 2 , j ) = j
        End Do
        Call Sort( DistList , 1 , Size( DistList, 2 ) , 1 )
        NN = 0
        Do j = 1 , Size( DistList , 2 )
          If ( DistList( 1 , j ) .gt. cutoff ) Exit
          NN = NN + 1
        End Do
        If ( NN .eq. 1 ) then
          NNN = NNN + 1
          ListTmp( 1 , NNN ) = i
          ListTmp( 2 , NNN ) = NINT( DistList( 2 , 1 )  )
        End If
      End Do
      Allocate( List( 1:2 , 1:NNN ) )
      List = ListTmp( : , 1:NNN )

      Deallocate( ListTmp )




    End Function ListOfABBonds

    !!
    !! routine checks which atoms of atoms
    !! are connected to exactly 1 atom of atoms excluding itself
    !! and avoiding double counting
    !! within the cutoff criterion
    !! function returns a list
 !!
 !! initializes A-A bonds
 !! link list containing List( 1 , j ) -> index A
 !! and List( 2 , j ) -> index B
 !!
    Function ListOfAABonds( atoms , lattice ,cutoff ) Result( List )
      use tool_interfaces, only :get_nearest_image_cartesian,get_norm
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : ) ::atoms
      Real*8,dimension( : , : ) ::lattice
      Real*8 ::cutoff

      Integer ::NN
      Integer ::NNN

      Integer ::i
      Integer ::j

      Real*8,allocatable,dimension( : , : ) ::DistList
      Real*8,dimension( 1:3 ) ::PosA
      Real*8,dimension( 1:3 ) ::PosB

      Integer,allocatable,dimension( : , : ) ::List
      Integer,allocatable,dimension( : , : ) ::ListTmp

      Allocate( DistList( 1:2 , 1:Size( atoms , 2 ) ) )
      Allocate( ListTmp( 1:2 , 1:Size( atoms , 2 ) ) )


      NNN = 0   !counter number of bonds satisifying crit


      MainLoop: Do i = 1, Size( atoms , 2 )
        PosA  =  MatMul( lattice , atoms( : , i ) )
        Do j = 1, NNN
          If ( ListTmp( 2 , j ) .eq. i ) Cycle MainLoop
        End Do
        Do j = 1, Size( atoms , 2 )
          If ( i .eq. j ) then
            DistList( 1 , j ) = Huge( DistList( 1 , j ) )
            Cycle
          End If
          PosB = get_nearest_image_cartesian( atoms( : , i ) , atoms( : , j ) , lattice )
          DistList( 1 , j ) = get_norm( PosA - PosB )
          DistList( 2 , j ) = j
        End Do
        Call Sort( DistList , 1 , Size( DistList, 2 ) , 1 )
        NN = 0
        Do j = 1 , Size( DistList , 2 )
          If ( DistList( 1 , j ) .gt. cutoff ) Exit
          NN = NN + 1
        End Do
        If ( NN .eq. 1 ) then
          NNN = NNN + 1
          ListTmp( 1 , NNN ) = i
          ListTmp( 2 , NNN ) = NINT( DistList( 2 , 1 )  )
        End If
      End Do MainLoop
      Allocate( List( 1:2 , 1:NNN ) )
      List = ListTmp( : , 1:NNN )
      Deallocate( ListTmp )
    End Function ListOfAABonds


    !!
    !! routine checks which atoms of atomsA
    !! are connected to exactly 1 atom of atomsB
    !! within the cutoff criterion
    !! function returns a list
    !! first dimension contains index of atomsA
    !! and second index gives you the index pointing to atomB
    !! ato which A is connected
    !!
    !!
    !!
    !! initializes A-A bonds
    !! link list containing List( 1 , j ) -> index A
    !! and List( 2 , j ) -> index B
    !!
    Function ListOfABABonds( atomsA , atomsB , lattice , cutoff , dummy ) Result( List )
      use tool_interfaces, only :get_nearest_image_cartesian,get_norm
      use MergeSorter, only :Sort
      Implicit None
      Real*8,dimension( : , : ) ::atomsA
      Real*8,dimension( : , : ) ::atomsB
      Real*8,dimension( : , : ) ::lattice
      Real*8 ::cutoff
      Integer ::dummy

      Integer ::NN
      Integer ::NNN

      Integer ::i
      Integer ::j

      Real*8,allocatable,dimension( : , : ) ::DistList
      Real*8,dimension( 1:3 ) ::PosA
      Real*8,dimension( 1:3 ) ::PosB

      Integer,allocatable,dimension( : , : ) ::List
      Integer,allocatable,dimension( : , : ) ::ListTmp

      Allocate( DistList( 1:2 , 1:Size( atomsB , 2 ) ) )
      Allocate( ListTmp( 1:3 , 1:Size( atomsB , 2 ) ) )

      dummy = dummy + 2

      NNN = 0   !counter number of bonds satisifying crit


      Do i = 1, Size( atomsA , 2 )
        PosA  =  MatMul( lattice , atomsA( : , i ) )
        Do j = 1, Size( atomsB , 2 )
          PosB = get_nearest_image_cartesian( atomsA( : , i ) , atomsB( : , j ) , lattice )
          DistList( 1 , j ) = get_norm( PosA - PosB )
          DistList( 2 , j ) = j
        End Do
        Call Sort( DistList , 1 , Size( DistList, 2 ) , 1 )
        NN = 0
        Do j = 1 , Size( DistList , 2 )
          If ( DistList( 1 , j ) .gt. cutoff ) Exit
          NN = NN + 1
        End Do
        If ( NN .eq. 2 ) then
          NNN = NNN + 1
          ListTmp( 2 , NNN ) = i
          ListTmp( 1 , NNN ) = NINT( DistList( 2 , 1 )  )
          ListTmp( 3 , NNN ) = NINT( DistList( 2 , 2 )  )
        End If
      End Do
      Allocate( List( 1:3 , 1:NNN ) )
      List = ListTmp( : , 1:NNN )
      Deallocate( ListTmp )
    End Function ListOfABABonds



    !!! Functions should be changed to 
    !!! subroutines otherwise
    !!! allocation will be done always again

    Function ABBondPositions( atomsA , atomsB , LinkList ) Result( Pos )

       use tool_interfaces,only :get_nearest_image
       Implicit None
       Real*8,dimension( : , : ) ::atomsA
       Real*8,dimension( : , : ) ::atomsB
       Integer,dimension( : , : ) ::LinkList
       Real*8,allocatable,dimension( : , : ) ::Pos

       Real*8,dimension( 1:3 ) ::Pos1
       Real*8,dimension( 1:3 ) ::Pos2

       Integer ::i

       If ( .not. Allocated( Pos ) ) Allocate( Pos( 1:3 , 1:Size( LinkList , 2 ) ) )

       Do i = 1, Size( LinkList , 2 )
         Pos1 = atomsA( : , LinkList( 1 , i ) )
         Pos2 = atomsB( : , LinkList( 2 , i ) )
         Pos2 = get_nearest_image( Pos1 , Pos2 )
         Pos( : , i ) = ( Pos1 + Pos2 ) * 0.5d0
       End Do
    End Function ABBondPositions



    Function BABBondPositions( atomsA , atomsB , atomsC , LinkList ) Result( Pos )

       use tool_interfaces,only :get_nearest_image
       Implicit None
       Real*8,dimension( : , : ) ::atomsA
       Real*8,dimension( : , : ) ::atomsB
       Real*8,dimension( : , : ) ::atomsC
       Integer,dimension( : , : ) ::LinkList
       Real*8,allocatable,dimension( : , : ) ::Pos

       Real*8,dimension( 1:3 ) ::Pos1
       Real*8,dimension( 1:3 ) ::Pos2
       Real*8,dimension( 1:3 ) ::Pos3

       Integer ::i

       If ( .not. Allocated( Pos ) ) Allocate( Pos( 1:3 , 1:Size( LinkList , 2 ) ) )

       Do i = 1, Size( LinkList , 2 )
         Pos1 = atomsA( : , LinkList( 2 , i ) )
         Pos2 = atomsB( : , LinkList( 1 , i ) )
         Pos3 = atomsC( : , LinkList( 3 , i ) )
         Pos2 = get_nearest_image( Pos1 , Pos2 )
         Pos3 = get_nearest_image( Pos1 , Pos3 )
         Pos( : , i ) = ( Pos1 + Pos2 + Pos3 ) / 3d0
       End Do
    End Function BABBondPositions


    Function ComputeABBond( atomsA , atomsB , LinkList ) Result( Vecs )

       use tool_interfaces,only :get_nearest_image,get_norm
       Implicit None
       Real*8,dimension( : , : ) ::atomsA
       Real*8,dimension( : , : ) ::atomsB
       Integer,dimension( : , : ) ::LinkList
       Real*8,allocatable,dimension( : , : ) ::Vecs

       Real*8,dimension( 1:3 ) ::Pos1
       Real*8,dimension( 1:3 ) ::Pos2

       Integer ::i

       If( .not. Allocated( Vecs ) ) Allocate( Vecs( 1:3 , 1:Size( LinkList , 2 ) ) )

       Do i = 1, Size( LinkList , 2 )
         Pos1 = atomsA( : , LinkList( 1 , i ) )
         Pos2 = atomsB( : , LinkList( 2 , i ) )
         Pos2 = get_nearest_image( Pos1 , Pos2 )
         Vecs( : , i ) = Pos1 - Pos2
         Vecs( : , i ) = Vecs( : , i ) / get_norm( Vecs( : , i ) )
       End Do
    End Function ComputeABBond


    Function ComputeBABBond( atomsA , atomsB , atomsC , LinkList ) Result( Vecs )

       use tool_interfaces,only :get_nearest_image,get_norm
       Implicit None
       Real*8,dimension( : , : ) ::atomsA
       Real*8,dimension( : , : ) ::atomsB
       Real*8,dimension( : , : ) ::atomsC
       Integer,dimension( : , : ) ::LinkList
       Real*8,allocatable,dimension( : , : ) ::Vecs

       Real*8,dimension( 1:3 ) ::Pos1
       Real*8,dimension( 1:3 ) ::Pos2
       Real*8,dimension( 1:3 ) ::Pos3

       Integer ::i

       If ( .not. Allocated( Vecs ) ) Allocate( Vecs( 1:3 , 1:Size( LinkList , 2 ) ) )

       Do i = 1, Size( LinkList , 2 )
         Pos1 = atomsA( : , LinkList( 2 , i ) )
         Pos2 = atomsB( : , LinkList( 1 , i ) )
         Pos3 = atomsC( : , LinkList( 3 , i ) )
         Pos2 = get_nearest_image( Pos1 , Pos2 )
         Pos3 = get_nearest_image( Pos1 , Pos3 )
         Vecs( : , i ) = Pos1 - ( Pos2 + Pos3 ) * 0.5d0
         Vecs( : , i ) = Vecs( : , i ) / get_norm( Vecs( : , i ) )
       End Do
    End Function ComputeBABBond


  End Module GenBondTypes
