
 Module MergeSorter

   !
   ! by calling Sort one will obtain an ordered array
   ! for cvalling Sort pass array as a vector and the
   ! minimum and maximum index for the considered array
   ! in standard Fortran 1 and N

   Implicit None
   Interface Sort
     Module Procedure MergeSort1D,MergeSort2D
   End Interface Sort

   Contains
     !!!
     !!!
     !!! This part of the routine is taken for vectors or two dimensional arrays
     !!!
     Subroutine Merge2D( x , left , middle , right , Col ) !j.L. 21.11.2018

       !! 
       !! Merge part of the sort algo; x is total data array and integers
       !! left, middle and right are giving the bounds for the array x
       !! to be splitted into its part; Ldata and Rdata contain the 
       !! two parts of the original array x
       !! 

       Implicit None
       Real*8,dimension( : , : )  ::x
       Integer ::left
       Integer ::right
       Integer ::middle
       Integer ::Col               !Coloumn of array according to which is sorted
       Real*8,dimension( 1:Size( x , 1 ) , 1 : middle - left + 1 )  ::Ldata
       Real*8,dimension( 1:Size( x , 1 ) , 1 : right - middle )  ::Rdata
       Integer ::N1
       Integer ::N2

       Integer ::i
       Integer ::j
       Integer ::k

       N1  =  middle - left + 1  ! Size of temp array Ldata
       N2  =  right - middle     ! Size of temp array Rdata

       !print*,"Left"
       !print*,left,middle
       !print*,Size( Ldata , 1 ),middle-left + 1,N1
       !print*,x
       !print*,"Right"
       !print*,right,middle
       !print*,Size( Rdata , 1 ),right - middle,N2


       Ldata = x( : , left : left + N1 )
       Rdata = x( : , middle + 1 : right + N2 )

       i = 1
       j = 1
       k = left
       Do While ( i .le. N1 .and. j .le. N2 )
          If ( Ldata( Col , i ) .le. Rdata( Col , j ) ) then
            x( : , k )  =  Ldata( : , i )
            i = i + 1
          Else
            x( : , k )  =  Rdata( : , j )
            j = j + 1
          End If
          k = k + 1
       End Do
       
       !Copy remaining elements of left hand side
       Do While( i .le. N1 )
           x( : , k )  =  Ldata( : , i )
           i = i + 1
           k = k + 1
       End Do

       !Copy remaining elemets of right hand side
       Do While( j  .le. N2 ) 
         x( : , k )  =  Rdata( : , j )
         j = j + 1
         k = k + 1
       End Do

       !print*,x


     End Subroutine Merge2D

     Recursive Subroutine MergeSort2D( x , left , right , Col )

       !! MergeSort1d computes the index where the array has to be split ->
       !! -> midddle; next calls itself again until base case is reached 
       !! determined by the if clause;
       !!

       Implicit None
       Real*8,dimension( : , : )  ::x
       Integer ::left
       Integer ::right
       Integer ::middle
       Integer ::Col

       !print*,"MAIN",left,right
       If ( left .lt. right ) then
          middle = ( left + ( right - 1 ) ) / 2
          Call MergeSort2D( x , left , middle , Col )
          Call MergeSort2D( x , middle + 1 , right , Col )
          Call Merge2D( x , left , middle , right , Col )
       End If

     End Subroutine MergeSort2D

     !!!
     !!!
     !!! This part of the routine is taken for vectors or one dimensional arrays
     !!!
     Subroutine Merge1D( x , left , middle , right ) !j.L. 21.11.2018

       !! 
       !! Merge part of the sort algo; x is total data array and integers
       !! left, middle and right are giving the bounds for the array x
       !! to be splitted into its part; Ldata and Rdata contain the 
       !! two parts of the original array x
       !! 

       Implicit None
       Real*8,dimension( : )  ::x
       Integer ::left
       Integer ::right
       Integer ::middle
       Real*8,dimension( 1 : middle - left + 1 )  ::Ldata
       Real*8,dimension( 1 : right - middle )  ::Rdata
       Integer ::N1
       Integer ::N2

       Integer ::i
       Integer ::j
       Integer ::k

       N1  =  middle - left + 1  ! Size of temp array Ldata
       N2  =  right - middle     ! Size of temp array Rdata

       !print*,"Left"
       !print*,left,middle
       !print*,Size( Ldata , 1 ),middle-left + 1,N1
       !print*,x
       !print*,"Right"
       !print*,right,middle
       !print*,Size( Rdata , 1 ),right - middle,N2


       Ldata = x( left : left + N1 )
       Rdata = x( middle + 1 : right + N2 )

       i = 1
       j = 1
       k = left
       Do While ( i .le. N1 .and. j .le. N2 )
          If ( Ldata( i ) .le. Rdata( j ) ) then
            x( k )  =  Ldata( i )
            i = i + 1
          Else
            x( k )  =  Rdata( j )
            j = j + 1
          End If
          k = k + 1
       End Do
       
       !Copy remaining elements of left hand side
       Do While( i .le. N1 )
           x( k )  =  Ldata( i )
           i = i + 1
           k = k + 1
       End Do

       !Copy remaining elemets of right hand side
       Do While( j  .le. N2 ) 
         x( k )  =  Rdata( j )
         j = j + 1
         k = k + 1
       End Do

       !print*,x


     End Subroutine Merge1D

     Recursive Subroutine MergeSort1D( x , left , right )

       !! MergeSort1d computes the index where the array has to be split ->
       !! -> midddle; next calls itself again until base case is reached 
       !! determined by the if clause;
       !!

       Implicit None
       Real*8,dimension( : )  ::x
       Integer ::left
       Integer ::right
       Integer ::middle

       !print*,"MAIN",left,right
       If ( left .lt. right ) then
          middle = ( left + ( right - 1 ) ) / 2
          Call MergeSort1D( x , left , middle )
          Call MergeSort1D( x , middle + 1 , right )
          Call Merge1D( x , left , middle , right )
       End If

     End Subroutine MergeSort1D

 End Module MergeSorter



! Program test
!
!   use MergeSorter, only: Sort
!   Implicit None
!   Real*8,allocatable,dimension( : ) ::x
!   Real*8,allocatable,dimension( : ) ::y
!   Real*8,allocatable,dimension( : , : ) ::xx
!   Real*8,allocatable,dimension( : , : ) ::yy
!   Integer ::i
!   Integer ::j
!   Integer ::N
!   Integer ::N1
!   Integer ::N2
!
!   N = 5
!
!   N1 = 10
!   N2 = 5
!
!   Allocate( x( 1:N) )
!   Allocate( y( 1:N) )
!   Allocate( xx( 1:N2 , 1:N1 ) )
!   Allocate( yy( 1:N2 , 1:N1 ) )
!
!   !Call Random_Seed( Size = i )
!   !Do i = 1, N
!   !  Call Random_Number( x( i ) )
!   !End Do
!
!   !Do i = 1, N
!   !  print*,x( i )
!   !  y( i )  =  x( i )
!   !End Do
!
!   !! Sort 1d array
!   !Call Sort( x , 1 , N )
!
!   !print*,
!   !Do i = 1, N
!   !  print*,x( i )
!   !End Do
!
!   !Call check( x , y , N )
!
!   Do i = 1 , N1
!     Do j = 1 , N2
!        Call Random_number( xx( j , i ) )
!     End Do
!     print*,xx( : , i )
!   End Do
!
!   !! Sort 2d array according to 1st column
!   Call Sort( xx , 1 , N1 , 1 )
!   print*,
!   Do i = 1 , N1
!     print*,xx( : , i )
!   End Do
! End Program
!
! Subroutine check( x , y , N )
!
!   Implicit None
!   Integer ::N
!   Real*8,dimension( 1:N ) ::x
!   Real*8,dimension( 1:N ) ::y
!
!   Integer ::i
!   Integer ::j
!   Logical ::Found
!
!   
!   Do i = 1, N
!      Found = .False.
!      Do j = 1 , N
!        If ( y( i ) .eq. x( j ) ) then
!          print*,"Found" , y( i )
!          Found = .True.
!        End If
!      End Do
!      If ( .not. Found ) then
!        print*,"Not found ",y( i ),i
!      End If
!   End Do
! End Subroutine
! 
! Subroutine check( x , y , N )
!
!   Implicit None
!   Integer ::N
!   Real*8,dimension( 1:N ) ::x
!   Real*8,dimension( 1:N ) ::y
!
!   Integer ::i
!   Integer ::j
!   Logical ::Found
!
!   
!   Do i = 1, N
!      Found = .False.
!      Do j = 1 , N
!        If ( y( i ) .eq. x( j ) ) then
!          print*,"Found" , y( i )
!          Found = .True.
!        End If
!      End Do
!      If ( .not. Found ) then
!        print*,"Not found ",y( i ),i
!      End If
!   End Do
! End Subroutine
