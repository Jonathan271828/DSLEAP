


  Subroutine lu_inversion( A , N )                            !Jonathan Lahnsteiner   05.August 2016
   !This routine calls an LU decomposition, so
   !the matrix will be first written as LU=A
   !next the routine calls a backsubstitution algorithm
   !Which gives the inverse of your matrix
   !Your initial input matrix will be overwritten with the inverse of your matrix
   
   implicit none
   Integer                               ::N
   Real*8,dimension( 1 : N , 1 : N )     ::A
   Real*8,allocatable,dimension( : , : ) ::II
   Integer                               ::i
   Integer                               ::j
   Integer,allocatable,dimension( : )    ::list
   Integer                               ::D
   
   Allocate( II( 1 : N , 1 : N ) )
   Allocate( list( 1 : N ) )
  
   Do i = 1 , N
     Do j = 1 , N
       II( i , j )  =  0d0
     End Do
     II( i , i )  =  1d0
   End Do 
  
   call lu_decomp( A , N , list , D )
  
   Do i = 1 , N
     call lu_bksb( A , N , list , II( : , i ) )
   End Do
  
   A = II
  
   Deallocate( II , list )
  
  End Subroutine
  
  Subroutine lu_decomp( A , N , list , D )
   !This routine perfroms an LU decomposition
   !of an NxN matrix A supplied to the routine
   !This routine was programmed after a routine shown in Numercial Recipes by 
   !William H. Press
   !Brian P. Flannery
   !Saul A. Teukolsky
   !William T. Vetterling
   !Cambridge University Press
   !1986
   !page 24
   implicit none
   Integer                           ::N
   Real*8,dimension( 1 : N , 1 : N ) ::A
   Integer,dimension( 1 : N )        ::list
   Real*8,allocatable,dimension( : ) ::scaler       
   Integer                           ::i
   Integer                           ::j
   Integer                           ::k
   Integer                           ::D
   Real*8                            ::maxpiv              !Maximum pivot element of every row
   Real*8                            ::sumer
   Real*8                            ::test
   Integer                           ::maxi
   Real*8,parameter                  ::mini = 1d-20
  
   Allocate( scaler( 1 : N ) )  !Stroes the scaling factors for the individual rows
   D = 1 !stores if number of permutations is even or odd; start with 1 no permutations performed
  
   !Nested loop for finding the largest elements in every row
   Do i = 1 , N
      maxpiv = 0d0
      Do j = 1 , N
         If ( Abs( A( i , j ) ) .gt. maxpiv  ) &
              maxpiv = Abs( A( i , j ) ) 
      End Do
      If ( maxpiv .eq. 0d0 ) then
        Write( * , "(A)" )
        Write( * , * ) "!!!Error, LU factorization ERROR!!!"
        Write( * , * ) "!!!Error, your matrix is singular ERROR!!!"
        Write( * , * ) "!!!Error, LU factorization ERROR!!!"
        Write( * , "(A)" )
        STOP
      End If
      scaler( i )  =  1d0 / maxpiv
   End Do
  
   Do j = 1 , N
     If ( j .gt. 1 ) then
       !Calculating betas except for diagonals
       Do i = 1 , j - 1
          sumer  =  A( i , j )
          If ( i .gt. 1 ) then
             Do k = 1 , i - 1
                sumer  =   sumer - A( i , k ) * A( k , j )
             End Do
             A( i , j )   =  sumer
          End If
       End Do
     End If
     maxpiv = 0d0
     
     !Calculating alphas and also diagonal betas 
     !simultanously the max pivot element is searched
     !The division is not written to the matrix
     Do i = j , N
       sumer  =  A( i , j )
       If ( j .gt. 1 ) then
         Do k = 1 , j - 1
           sumer  =  sumer  - A( i , k ) * A( k , j )
         End Do
         A( i , j )  =  sumer
       End If
       test  =  scaler( i ) * ABS( sumer )
       If ( test .ge. maxpiv ) then
         maxpiv  =  test
         maxi    =  i 
       End If
     End Do
  
     !Checking if pivoting is necessary
     If ( j .ne. maxi ) then
       !Interchanging two rows to have max element on diagonal
       Do k = 1 , N 
         test  =  A( maxi , k )
         A( maxi , k )  =   A( j , k )
         A( j , k )     =   test
       End Do
       D  =  -D
       scaler( maxi )   =   scaler( j )
     End If
  
     list( j )  =  maxi
  
     !Here the division is carried out
     !If the pivot element would be zero now
     !it would be set to a minium value
     !the minimum value can be neglected if desired
     If ( j .ne. N ) then
        If ( A( j , j ) .eq. 0 ) A( j , j )  =  mini
        test = 1d0 / A( j , j )
        Do i = j + 1 , N
          A( i , j )  =  A( i , j ) * test
        End Do
     End If
     If ( A( N , N ) .eq. 0 ) A( N , N ) = mini
   End Do
  
   Deallocate( scaler )
  End Subroutine
  
  
  Subroutine lu_bksb( A , N , list , b )
    !This routine performs the backsubstituion
    !after the LU decomposition
    !this routine can be used either for
    !the solution of equations where b is the right hand side of
    !the considered equation or one is able to invert matrices
    !if b are the coloumns of the unitary matrix
    !This routine was programmed according to
    !William H. Press
    !Brian P. Flannery
    !Saul A. Teukolsky
    !William T. Vetterling
    !Cambridge University Press
    !1986
    !page 24
    implicit none
    Integer                           ::N            !dimension of matrix
    Real*8,dimension( 1 : N , 1 : N ) ::A            !NxN matrix -->> LU decomposed
    Real*8,dimension( 1 : N )         ::b            !Right hand side of equation
    Integer,dimension( 1 : N )        ::list         !pivoting list returned by the lu_decomp routine
    Integer                           ::i
    Integer                           ::j          
    Integer                           ::fnz          !fist non zero element in b
    Real*8                            ::sumer        !variable used for suming
  
    fnz = 0
  
    !The first backsubstituion Ly=b
    Do i = 1 , N
      sumer           =  b( list( i ) )
      b( list( i ) )  =  b( i )
      If ( fnz .ne. 0 ) then
        Do j = fnz , i - 1
           sumer = sumer - A( i , j ) * b( j )
        End Do
      ElseIf ( sumer .ne. 0 ) then
        fnz  =  i
      End If
      b( i )  =  sumer
    End Do
  
    !Second backsubtitution Ux=y
    Do i = N , 1 , -1
      sumer  =  b( i )
      If ( i .lt. N ) then
        Do j = i + 1 , N      
          sumer =  sumer - A( i , j ) * b( j )
        End Do
      End If    
      b( i )  =  sumer  /  A( i , i )
    End Do
  
  End Subroutine
