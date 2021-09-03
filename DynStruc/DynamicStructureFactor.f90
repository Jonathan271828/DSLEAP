 
!!!
!!!
!!! Module contains routines that generate Qgrids for you
!!!
!!!

Module QgridCollection
  

  Implicit None

  Contains
     !!!
     !!!
     !!!  comensurate qgrid with respect to the used supercell
     !!!
     Subroutine ComputeQgridDirect( Qgrid , Nx , Ny , Nz , WhoCalls )

       Implicit None
       Real*8,allocatable,dimension( : , : )  ::Qgrid
       Integer ::Nx
       Integer ::Ny
       Integer ::Nz
       Character( len = * ),optional ::WhoCalls
       Integer ::i
       Integer ::j
       Integer ::k
       Integer ::col

       If ( Allocated( Qgrid ) ) Deallocate( Qgrid )

       Allocate( Qgrid( 1:3 , 1:Nx*Ny*Nz ) )
       Write( * , * ) "Sampling full Brillouin zone for"
       If ( present( WhoCalls ) ) then
           Write( * , * ) Trim( WhoCalls )
       End If

       col = 1
       Do i = 0 , Nx-1
         Do j = 0 , Ny-1
           Do k = 0 , Nz-1
              Qgrid( 1 , col )  = Modulo( i + Nx / 2 , Nx ) -Nx / 2 + 1
              Qgrid( 2 , col )  = Modulo( j + Ny / 2 , Ny ) -Ny / 2 + 1
              Qgrid( 3 , col )  = Modulo( k + Nz / 2 , Nz ) -Nz / 2 + 1
              Write( * , "(I4,3F15.8)" ) col,Qgrid( : , col )
              col = col + 1
           End Do
         End Do
       End Do
     End Subroutine ComputeQgridDirect


     !!!
     !!! sample first brillouin zone
     !!! to consider harmonic crystal test case
     !!!
     Subroutine SampleFirstBrillouinZone( Qgrid , N , WhoCalls )

       Implicit None
       Real*8,allocatable,dimension( : , : ) ::Qgrid( : , : )
       Character( len = * ),optional ::WhoCalls
       Integer ::N
       Integer ::i
       Integer ::col

       If ( Allocated( Qgrid ) ) Deallocate( Qgrid )

       Allocate( Qgrid( 1:3 , 1:N ) )
       Write( * , * ) "Sampling Brillouin zone along x"
       If ( present( WhoCalls ) ) then
           Write( * , * ) Trim( WhoCalls )
       End If
       Qgrid = 0d0
       col = 1
       Do i = 0 , N-1
          Qgrid( 1 , col )  = MODULO( i+ N/2-1, N)-N/2+1
          Write( * , "(I4,3F15.8)" ) col,Qgrid( : , col )
          col = col  +  1
       End Do
     End Subroutine SampleFirstBrillouinZone


     Subroutine SampleFirstBrillouinZoneCubic( Qgrid , Nx , Ny , Nz , WhoCalls )
      !!!! sample the qpoints in the first irreducible brillouin zone
      !!!! 
      
       Implicit None
       Real*8,allocatable,dimension( : , : ) ::Qgrid
       Integer  ::Nx
       Integer  ::Ny
       Integer  ::Nz
       Character( len = * ),optional ::WhoCalls

       Integer  ::Nx2
       Integer  ::Ny2
       Integer  ::Nz2
       Integer ::i
       Integer ::j
       Integer ::k
       Integer ::col

       Nx2 =  Nx / 2 
       Ny2 =  Ny / 2 
       Nz2 =  Nz / 2 
       
       
       If ( Allocated( Qgrid ) ) Deallocate( Qgrid )

       Allocate( Qgrid( 1:3 , 1:Nx2*Ny2*Nz2 + Nx2*Ny2 + Nx2*Nz2 + Ny2*Nz2 ) )
       Write( * , * ) "Sampling irreducible Brillouin zone for"
       If ( present( WhoCalls ) ) then
           Write( * , * ) Trim( WhoCalls )
       End If
       
       col = 1
       Do i = 0 , Nx2-1
         Do j = 0 , Ny2-1
           Do k = 0 , Nz2-1
              Qgrid( 1 , col )  = Modulo( i + Nx / 2 , Nx ) -Nx / 2 + 1
              Qgrid( 2 , col )  = Modulo( j + Ny / 2 , Ny ) -Ny / 2 + 1
              Qgrid( 3 , col )  = Modulo( k + Nz / 2 , Nz ) -Nz / 2 + 1
              Write( * , "(I4,3F15.8)" ) col,Qgrid( : , col )
              col = col + 1
           End Do
         End Do
       End Do

       !!Ny2*Nz2 term
       i = Nx - 1
       Do j = 0 , Ny2-1
         Do k = 0 , Nz2-1
            Qgrid( 1 , col )  = Modulo( i + Nx / 2 , Nx ) -Nx / 2 + 1
            Qgrid( 2 , col )  = Modulo( j + Ny / 2 , Ny ) -Ny / 2 + 1
            Qgrid( 3 , col )  = Modulo( k + Nz / 2 , Nz ) -Nz / 2 + 1
            Write( * , "(I4,3F15.8)" ) col,Qgrid( : , col )
            col = col + 1
         End Do
       End Do

       !! Nx2*Nz2 term
       j = Ny - 1
       Do i = 0 , Nx2-1
          Do k = 0 , Nz2-1
             Qgrid( 1 , col )  = Modulo( i + Nx / 2 , Nx ) -Nx / 2 + 1
             Qgrid( 2 , col )  = Modulo( j + Ny / 2 , Ny ) -Ny / 2 + 1
             Qgrid( 3 , col )  = Modulo( k + Nz / 2 , Nz ) -Nz / 2 + 1
             Write( * , "(I4,3F15.8)" ) col,Qgrid( : , col )
             col = col + 1
          End Do
       End Do

       !! Ny2*Nz2 term
       k = Nz - 1
       Do i = 0 , Nx2-1
         Do j = 0 , Ny2-1
             Qgrid( 1 , col )  = Modulo( i + Nx / 2 , Nx ) -Nx / 2 + 1
             Qgrid( 2 , col )  = Modulo( j + Ny / 2 , Ny ) -Ny / 2 + 1
             Qgrid( 3 , col )  = Modulo( k + Nz / 2 , Nz ) -Nz / 2 + 1
             Write( * , "(I4,3F15.8)" ) col,Qgrid( : , col )
             col = col + 1
         End Do
       End Do

     End Subroutine SampleFirstBrillouinZoneCubic

     !!
     !!
     !! cubic brillouin zone sampling
     !!
     Subroutine BrillouinZoneCubicGXMGRX( Qgrid , Nx , Ny , Nz , WhoCalls )
       Implicit None
       Real*8,allocatable,dimension( : , : ) ::Qgrid
       Integer  ::Nx
       Integer  ::Ny
       Integer  ::Nz
       Character( len = * ),optional ::WhoCalls

       Integer ::i
       Integer ::col

       Integer ::N
       Real*8,allocatable,dimension( : , : ) ::QGridTemp


       N = Nx * Ny * Nz

       Allocate( QGridTemp( 1:3 , N ) )

       Write( * , * ) "Sampling Brillouin zone along" 
       Write( * , * ) "Gamma-X-M-Gamma-R"
       If ( present( WhoCalls ) ) then
               Write( * , "(A)" ) Trim( WhoCalls )
       End If
       Write( * , * ) "in cubic box"

       col = 1
       QGridTemp = 0d0
       !Gamma -> X
       Do i = 0 , Ny / 2
          QgridTemp( 2 , col ) = i
          col = col + 1
       End Do

       !X -> M
       Do i = 1 , Nx / 2
          QgridTemp( 2 , col ) = Ny / 2
          QgridTemp( 1 , col ) = i
          col = col + 1
       End Do

       !M -> Gamma
       Do i = Nx / 2 - 1 , 0 , -1
          QgridTemp( 2 , col ) = i
          QgridTemp( 1 , col ) = i
          col = col + 1
       End Do
       !Gamma -> R
       Do i = 1, Nx/2
          QgridTemp( 1 , col ) = i
          QgridTemp( 2 , col ) = i
          QgridTemp( 3 , col ) = i
          col = col + 1
       End Do

       !R -> X
       Do i = Nx / 2 - 1 , 0 , -1 
          QgridTemp( 1 , col ) = i
          QgridTemp( 2 , col ) = Ny / 2
          QgridTemp( 3 , col ) = i
          col = col + 1
       End Do


       If ( Allocated( Qgrid ) ) Deallocate( Qgrid ); Allocate( Qgrid( 3 , col-1 ) )
       Write( * , * ) "Qpoints"
       Do i = 1, col-1
         Qgrid( : , i ) = QgridTemp( : , i ) !+ shift
         Write( * , "(I4,3F15.8)" ) i,QgridTemp( : , i )
       End Do
       Deallocate( QgridTemp )
     End Subroutine BrillouinZoneCubicGXMGRX


     !!
     !! orthorhombic structure
     !!
     Subroutine BrillouinZoneGXUZGS( Qgrid , Nx , Ny , Nz , WhoCalls )

       Implicit None
       Real*8,allocatable,dimension( : , : ) ::Qgrid
       Integer  ::Nx
       Integer  ::Ny
       Integer  ::Nz
       Character( len = * ),optional ::WhoCalls

       Integer ::i
       Integer ::col

       Integer ::N
       Real*8,allocatable,dimension( : , : ) ::QGridTemp


       N = Nx * Ny * Nz

       Allocate( QGridTemp( 1:3 , N ) )



       Write( * , * ) "Sampling Brillouin zone along" 
       Write( * , * ) "Gamma-X-S-Y-Gamma-Z-U-R-T-Y"
       If ( present( WhoCalls ) ) then
           Write( * , "(A)" ) Trim( WhoCalls )
       End If

       col = 1
       QGridTemp = 0d0
       !! gamma -> X
       Do i = 0 , Nx / 2
          QgridTemp( 1 , col )  =  i
          col = col + 1
       End Do

       !! X -> S
       Do i = 1 , Ny / 2
          QgridTemp( 1 , col )  =  Nx / 2
          QgridTemp( 2 , col )  =  i 
          col = col + 1
       End Do

       !! S -> Y
       Do i = Nx / 2 - 1 , 0 , -1
          QgridTemp( 1 , col )  =  i
          QgridTemp( 2 , col )  =  Ny / 2 
          col = col + 1
       End Do

       !! Y -> Gamma
       Do i = Ny / 2 -1, 0 , -1
          QgridTemp( 2 , col )  =  i
          col = col + 1
       End Do


       !! Gamma -> Z
       Do i = 1 , Nz / 2
          QgridTemp( 3 , col )  =  i
          col = col + 1
       End Do
       
       !! Z -> U
       Do i = 1 , Nx / 2
          QgridTemp( 3 , col )  =  Nz / 2
          QgridTemp( 1 , col )  =  i
          col = col + 1
       End Do

       !! U -> R
       Do i = 1 , Ny / 2
          QgridTemp( 1 , col )  =  Nx / 2
          QgridTemp( 2 , col )  =  i
          QgridTemp( 3 , col )  =  Nz / 2
          col = col + 1
       End Do

       !! R -> T
       Do i = Nx / 2 - 1 , 0 , -1
          QgridTemp( 1 , col )  =  i
          QgridTemp( 2 , col )  =  Ny / 2 
          QgridTemp( 3 , col )  =  Nz / 2 
          col = col + 1
       End Do

       !! T -> Z
       Do i = Ny / 2 -1, 0 , -1
          QgridTemp( 2 , col )  =  i
          QgridTemp( 3 , col )  =  Nz / 2
          col = col + 1
       End Do


       If ( Allocated( Qgrid ) ) Deallocate( Qgrid ); Allocate( Qgrid( 3 , col-1 ) )
       Write( * , * ) "Qpoints"
       Do i = 1, col-1
         Qgrid( : , i ) = QgridTemp( : , i ) !+ shift
         Write( * , "(I4,3F15.8)" ) i,QgridTemp( : , i )
       End Do
       Deallocate( QgridTemp )

     End Subroutine BrillouinZoneGXUZGS

     
     !!Compute qgrid in cartesian coordinates
     !! including the factor 2pi
     Subroutine ComputeCartesianQgrid( Dir , Cart , lattice )

       Implicit None
       Real*8,dimension( : , : ) ::Dir
       Real*8,dimensioN( : , : ) ::Cart
       Real*8,dimension( : , : ) ::lattice
       Real*8,parameter ::PI2 = 4d0 * Dacos( 0d0 )                          !!! factor 2pi

       Integer ::i

       Do i = 1 , Size( Dir , 2 )
         Cart( : , i )  = PI2 * MatMul( lattice , Dir( : , i ) )
       End Do
     End Subroutine ComputeCartesianQgrid


End Module QgridCollection


 Module PhononHelperRoutines
         Implicit None
         Contains
         !!!
         Subroutine GetMasses( Masses , fname )

            Implicit None
            Real*8,allocatable,dimension( : ) ::Masses
            Character( len = * ) ::fname
            Integer ::N
            Integer ::i
            
            open( unit = 33 , status='old' , file = Trim( fname ) , action='read' )
            Read( 33 , * ) N
              Allocate( Masses( 1:N ) ) 
              Do i = 1, N
                 Read( 33 , * ) Masses( i )
                 Masses( i ) = Dsqrt( Masses( i ) )
              End Do
            close( 33 )


         End Subroutine GetMasses

      

     


         !
         ! read vectors from MDRun code
         ! phonon finite differences approach
         ! makes a projection
         ! on phonon eigenvectors possible
         !
         !

         Subroutine ReadEigenvectorsMDRun( BasisVectors ,  fname )

           Implicit None
           Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
           Character( len=* )  ::fname
           Integer ::Nq
           Integer ::N1
           Integer ::N2
           Integer ::i
           Integer ::j
           Integer ::k
           Real*8  ::dummy1 , dummy2 , dummy3


           open( unit = 31 , action='read' , status='old' , file = Trim( fname ) )

           Read( 31 , * ) Nq , N1
           Read( 31 , * )
           N2 = N1 / 3

           Allocate( BasisVectors( 1:3 , 1:N2 , 1:N1 , 1:Nq ) )

           Write( * , * ) "Using phonon polarization vectors"
           !! Using 
           Do i = 1, Nq 
              Do j = 1 , N1
                 Read( 31 , * )
                 Do k = 1 , N2
                    Read( 31 , * ) BasisVectors( 1 , k , j , i ) , dummy1 , &
                                   BasisVectors( 2 , k , j , i ) , dummy2 , &
                                   BasisVectors( 3 , k , j , i ) , dummy3
                    Write( * , "(3F15.8)" ) BasisVectors( : , k , j , i )
                 End Do
              End Do
              If ( i .ne. Nq ) then
                Read( 31 , * )
                Read( 31 , * )
                Read( 31 , * )
              End If
              Write( * , * )
           End Do
           close( 31 )

         End Subroutine ReadEigenvectorsMDRun
     


         Subroutine ReadEigenvectorsImproved( BasisVectors , NQgrid , fname )

           Implicit None
           Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
           Integer ::NQgrid
           Character( len = * ) ::fname
           Integer ::i
           Integer ::j
           Integer ::N1
           Integer ::N2
           Real*8,allocatable,dimension( : , : , : ) ::TempVectors

           open( file = Trim( fname ) , status = 'old' , action = 'read' ,&
                 unit = 31 )
           Read( 31 , * ) N1 , N2

           Allocate( TempVectors( 1:3 , 1:N1 , 1:N2 ) )
           Do i = 1 , N2
              Read( 31 , * )
              Do j = 1 , N1
                 Read( 31 , * ) TempVectors( : , j , i )
              End Do
           End Do

           Allocate( BasisVectors( 1:3 , 1:N1 , 1:N2 , 1:NQgrid ) )

           Do i = 1 , NQgrid
              BasisVectors( : , : , : , i )  =  TempVectors
           End Do
           

           close( 31 )


         End Subroutine ReadEigenvectorsImproved

         Subroutine ReadEigenvectorsModeB( BasisVectors , fname )

           Implicit None
           Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
           Character( len = * )  ::fname
           
           Integer ::N1
           Integer ::N2
           Integer ::N3

           Integer ::i
           Integer ::j
           Integer ::k


          
           open( file = Trim( fname ) , status = 'old' , action = 'read' ,&
                unit = 31 )

            Read( 31 , * ) N1 , N2 , N3
            Allocate( BasisVectors( 1:3 , 1:N3 , 1:N2 , 1:N1 ) )
            Do i = 1 , N1
               Do j = 1, N2
                  Do k = 1, N3
                     Read( 31 , * ) BasisVectors( : , k , j , i )
                  End Do
                  Read( 31 , * )
               End Do
               Read( 31 , * )
            End Do
            close( 31 )


         End Subroutine ReadEigenVectorsModeB
     

         Subroutine ReadEigenVectorsModeC( BasisVectors , QGrid  , fname )

                 Implicit None
                 Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
                 Real*8,allocatable,dimension( : , :  )  ::QGrid
                 Character( len = * )  ::fname

                 Integer ::qq
                 Integer ::branch
                 Integer ::atom
                 Integer ::xyz

                 Write( * , * ) "Using eigenmode C"
                 Write( * , * ) "Qvectors are overwritten"
                 open( file = Trim( fname ) , status = 'old' , action = 'read' ,&
                     unit = 31 )

                 Read( 31 , * ) qq , branch , atom , xyz

                 If ( allocated( BasisVectors ) ) Deallocate( BasisVectors )
                 If ( allocated( Qgrid ) ) Deallocate( QGrid )

                 Allocate( BasisVectors( 1:xyz , 1:atom , 1:branch , 1:qq ) )
                 Allocate( QGrid( 1:xyz , 1:qq ) )

                 Write( * , * ) "Using Q-vectors instead"
                 Do qq = 1 , Size( BasisVectors , 4 )
                    Read( 31 , * ) QGrid( : , qq )
                    Write( * , "(3F15.8)" ) QGrid( : , qq )
                    Do branch = 1, Size( BasisVectors , 3 )
                       Do atom = 1 , Size( BasisVectors , 2 )
                           Read( 31 , * )  BasisVectors( : , atom , branch , qq )
                       End Do
                       Read( 31 , * )
                    End Do
                    Read( 31 , * )
                 End Do
                 close( 31 )
         End Subroutine ReadEigenVectorsModeC

         Function ComputeNorm( velos ) Result( Average )

                 Implicit None
                 Real*8,dimension( : , : )  ::Velos
                 Real*8  ::Average

                 Integer ::i

                 Average  =  0d0

                 Do i = 1, Size( Velos , 2 )
                    Average = Average + Dot_Product( Velos( : , i ),Velos( : , i ) )
                 End Do

                 Average = Average / Size( Velos , 2 )


         End Function




         !!!
         !!! compute dot product of 3xN matrix
         !!!
         Function Compute2x2Dot( A , B ) Result( Dotp )

           Implicit None
           Real*8,dimension( : , : )  ::A
           Real*8,dimension( : , : )  ::B

           Integer ::i

           Real*8  ::Dotp
           

           Dotp = 0d0
           Do i = 1 , Size( A , 2 )
              Dotp = Dotp + A( 1 , i ) * B( 1 , i )
              Dotp = Dotp + A( 2 , i ) * B( 2 , i )
              Dotp = Dotp + A( 3 , i ) * B( 3 , i )
           End Do


         End Function Compute2x2Dot
     


         Subroutine MakeDiagonalBasis( SpaceDim , Natoms , NQ , BasisVectors )

             Implicit None
             Integer,intent( in ) ::SpaceDim
             Integer,intent( in ) ::Natoms
             Integer,intent( in ) ::NQ

             Real*8,intent(inout),allocatable,dimension( : , : , : , : ) ::BasisVectors


             Integer ::i
             Real*8,allocatable,dimension( : , : , : ) ::Temp

             Integer ::col
             Integer ::row


             Allocate( BasisVectors( 1:SpaceDim , 1:Natoms , 1:SpaceDim*Natoms , NQ ) )
             Allocate( Temp( 1:SpaceDim , 1:Natoms , 1:SpaceDim*Natoms ) )

             Temp = 0d0
             col = 1
             row = 1
             Do i = 1 , Natoms*SpaceDim
                Temp( row , col , i ) = 1d0
                row = row + 1
                If ( row .gt. SpaceDim ) then
                        row = 1
                        col = col + 1
                End If
             End Do


             Do i = 1 , SpaceDim*Natoms
                BasisVectors( : , : , : , i )  =  Temp
             End Do


         End Subroutine MakeDiagonalBasis



         Subroutine ReadQpointsFromFile( QGrid , fname )

             Implicit None
             Real*8,allocatable,dimension( : , : ) ::QGrid
             Character( len = * ) ::fname
             Integer ::N
             Integer ::i
             open( unit = 31 , status = 'unknown' , file = Trim( fname ) , action='read' )
             Read( 31 , * ) N
             If ( Allocated( QGrid ) ) Deallocate( QGrid )
             Allocate( QGrid( 1:3 , 1:N ) )
             Write( * , * ) "Reading QVectors from file " , Trim( fname )
             Do i = 1 , N
                Read( 31 , * ) QGrid( : , i ) 
                Write( * , "(10000F15.8)" ) QGrid( : , i )
             End Do
             close( 31 )
         End Subroutine ReadQpointsFromFile



         Subroutine ReadScatteringLengths( Array , fname )

             Implicit None
             Real*8,allocatable,dimension( : )  ::Array
             Character( len = * ) ::fname
             Integer ::N
             Integer ::i

             open( unit = 31 , status = 'unknown' , file = Trim( fname ) , action='read' )
             Read( 31 , * ) N
             If ( Allocated( Array ) ) Deallocate( Array )
             Allocate( Array( 1:N ) )
             Do i = 1 , N
                Read( 31 , * ) Array( i ) 
                Write( * , "(10000F15.8)" ) Array( i )
             End Do
             close( 31 )

         End Subroutine ReadScatteringLengths


 End Module PhononHelperRoutines



!!!
!!!                   ~~
!!!                 ~~~~~~
!!!              ~~~~~~~~~~~~
!!!            ~~~~~~~~~~~~~~~~
!!!  Projected velocity autocorrelation function
!!!            ~~~~~~~~~~~~~~~~
!!!              ~~~~~~~~~~~~
!!!                 ~~~~~~
!!!                   ~~

 Module ComputePVACF

   use CrossCorrComplex, only :CrossCorrVar
   use SplitCelltoBox , only :BoxList
   Implicit None
   Type PVACFVar 
     Real*8,allocatable,dimension( : , : )  ::QgridDir                  !!! Qgrid in direct coordinates
     Integer ::Nqgrid                                                   !!! size qgrid
     Integer ::ActStep                                                  !!! actual MD step
     Integer ::MaxStep                                                  !!! total number of MD steps
     Real*8  ::Tstep                                                    !!! time step
     Real*8,dimensioN( 1:3 , 1:3 )  ::InvLattice                        !!! Inverse lattice
     Real*8,allocatable,dimension( : ) ::norm
     Integer ::QPathNr
     Integer,allocatable,dimension( : )  ::Qpath                        !!! indices containing qpath
     Type( CrossCorrVar ) ::TimeCorrelation
     Logical ::TimeCorrOnOff
     Integer ::Natoms

     !!! projected stuff
     Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
     Type( CrossCorrVar ),allocatable,dimension( : ) ::TimeCorrelationProj
     Real*8,allocatable,dimension( : , : )  ::PosOld
     Real*8,allocatable,dimensioN( : , : , : )  ::Velocity 
     Complex*16,allocatable,dimension( : , : , : )  ::PVACF   !!! dynamic structure factor in time domain
     Type( BoxList ),allocatable,dimension( : )  ::LinkUC     !!! Link atoms to unit cells
     Real*8,allocatable,dimension( : , : )  ::Boxes           !!! centers for fourier transform
     Real*8,allocatable,dimension( : ) ::Masses
     Logical :: FilterFunction                                !!! checks if exponential filter function is applied
   End Type PVACFVar

   Real*8,parameter ::PI2 = 4d0 * Dacos( 0d0 )                          !!! factor 2pi
   Complex*16,parameter ::ImagI = Cmplx( 0d0 , 1d0 )                    !!! imaginary number i


   Contains



     Subroutine ComputePVACFInit( data , Nsteps , dt , Nx , Ny , Nz , atoms , lattice , fname )

       use FlagReader, only :ReadFlags
       use CrossCorrComplex, only :CrossCorrelationInit
       use PhononHelperRoutines , only :GetMasses,&
                                        ReadEigenVectorsModeC,&
                                        MakeDiagonalBasis,&
                                        ReadQpointsFromFile
       use QgridCollection, only :SampleFirstBrillouinZoneCubic

       Implicit None
       Type( PVACFVar ) ::data
       Integer ::NSteps
       Real*8  ::dt
       Integer ::Nx
       Integer ::Ny
       Integer ::Nz
       Real*8,allocatable,dimension( : , : )  ::atoms
       Real*8,dimension( : , : ) ::lattice
       Character( len = * )  ::fname


       Logical ::QFound
       Real*8  ::temp
       Integer ::i


       data%MaxStep     =  NSteps
       data%Tstep       =  dt
       data%InvLattice  =  lattice
       Call lu_inversion( data%InvLattice , 3 )


       Qfound = .False.
       Call ReadFlags( "Qpath" , temp , QFound , fname )
       If ( .not. Qfound ) then
         data%QPathNr = 1
       Else
         data%QPathNr = NINT( temp )
       End If

       data%Natoms = Size( atoms , 2 )


       Call ReadBoxList( data , Nx , Ny , Nz , Size( atoms , 2 ) , "BoxList.in" )
       Call GetMasses( data%Masses , "masses.in" )





       Inquire( File = "BasisVector.in" , exist = QFound )
       If ( QFound ) then
          Write( * , * ) "ComputePVACFInit is reading BasisVectors and Qpoints"
          Write( * , * ) "from file BasisVector.in"
          Call ReadEigenvectorsModeC( data%BasisVectors , data%QGridDir , "BasisVector.in" )
          data%NQgrid = Size( data%QGridDir , 2 )
     
       Else
          Write( * , * ) "ComputePVACFInit No BasisVector.in file found"
          Write( * , * ) "Checking for a QVectorsFile"
          Inquire( File="QVectors.in" , exist = Qfound )
          If ( Qfound ) then
             Write( * , * ) "In routine ComputePVACFInit Q-vectors are read"
             Write( * , * ) "QVectors.in"
             Call ReadQpointsFromFile( data%QGridDir , "QVectors.in" )
             data%NQgrid  =  Size( data%QgridDir , 2 )
          Else
             Write( * , * ) "In routine ComputePVACFInit no QVectors.in"
             Write( * , * ) "supplied. Using default Q-mesh of first irreducible "
             Write( * , * ) "Brillouin zone"
             Call SampleFirstBrillouinZoneCubic( data%QgridDir , Nx , Ny , Nz )
             data%NQgrid  =  Size( data%QgridDir , 2 )
          End If
          Write( * , * ) "Using a diagonal basis since no file was supplied"
          Write( * , * ) "Diagonal basis gives some decomposition but does"
          Write( * , * ) "make sense physically"
          Call MakeDiagonalBasis( 3 , data%LinkUC( 1 )%Natoms , data%Nqgrid , data%BasisVectors )
       End If
       !!
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       data%QgridDir = PI2 * data%QgridDir
       data%PosOld = atoms
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       

       data%TimeCorrOnOff = .True.
       !! supply window size
       Call ReadFlags( "DSTC" , temp , data%TimeCorrOnOff , fname )
       If ( .not. data%TimeCorrOnOff ) then
          temp = Real( NSteps - 1 )
          data%TimeCorrOnOff = .True.
       End If

       If ( data%TimeCorrOnOff ) then
          Write( * , * ) "Time Averaging of the projected velocity autocorrelation"
          Write( * , * ) "is switched on. The window size was set to"
          Write( * , * ) NINT( temp )
          Allocate( data%TimeCorrelationProj( 1:Size( data%BasisVectors , 3 ) ) )

          !! Projected stuff
          Do i = 1, Size( data%BasisVectors , 3 )
              !! -1 because velocities are needed
              Call CrossCorrelationInit( data%TimeCorrelationProj( i ) ,&
                   NINT( temp ) , NSteps-1 , Size( data%QgridDir , 2 ) )
          End Do
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       End If



     End Subroutine ComputePVACFInit

     Subroutine ReadBoxList( data , N1 , N2 , N3 , Natoms , fname )

        use SplitCelltoBox, only :MakeBoxes,MakeBoxList
        Implicit None
        Type( PVACFVar ) ::data
        Integer  ::N1
        Integer  ::N2
        Integer  ::N3
        Integer  ::Natoms
        Character( len = * ) ::fname
        Integer  ::Npc
        Real*8,allocatable,dimension( : , : )  ::Boxes

        Integer ::i
        

        Npc = Natoms / ( N1 * N2 * N3 )
       
        Boxes = MakeBoxes( N1 , N2 , N3 )


        Allocate( data%LinkUC( 1:N1*N2*N3 ) )
        open( unit = 33 , action='read' , file = Trim( fname ) , status='old' )
        Do i = 1 , Size( data%LinkUC , 1 )
           Allocate( data%LinkUC( i )%List( 1:Npc ) )
           data%LinkUC( i )%Natoms = Npc
           Read( 33 , * ) data%LinkUC( i )%List( : )
        End Do
        close( 33 )


        
        Allocate( data%Velocity( 1:3 , 1:NPc , 1:N1*N2*N3 ) )
        Allocate( data%Boxes( 1:3 , 1 : N1 * N2 * N3 ) )
        Do i = 1 , Size( data%Boxes , 2 )
           data%Boxes( 1 , i )  =  0.5d0 * ( Boxes( 1 , i ) + Boxes( 2 , i ) )
           data%Boxes( 2 , i )  =  0.5d0 * ( Boxes( 3 , i ) + Boxes( 4 , i ) )
           data%Boxes( 3 , i )  =  0.5d0 * ( Boxes( 5 , i ) + Boxes( 6 , i ) )
        End Do

     End Subroutine ReadBoxList


     Subroutine ComputeVelocity( data , PosNew , lattice )

       use tool_interfaces, only :get_nearest_image,get_norm,&
                                  get_norm
       Implicit None
       Type( PVACFVar ) ::data
       Real*8,dimension( : , : )     ::PosNew
       Real*8,dimension( : , : ) ::lattice

       Integer ::i
       Integer ::j
       Integer ::Indx

       Do i = 1 , Size( data%LinkUC , 1 )
         Do j = 1 , data%LinkUC( i )%Natoms
            Indx = data%LinkUC( i )%List( j )
            data%Velocity( : , j , i )  =  MatMul( lattice , get_nearest_image( data%PosOld( : , Indx ) ,&
                                           PosNew( : , Indx ) ) - data%PosOld( : , Indx ) )
            data%Velocity( : , j , i )  =  data%Velocity( : , j , i ) * data%Masses( j )
            data%PosOld( : , Indx )     =  PosNew( : , Indx )
         End Do
       End Do

     End Subroutine ComputeVelocity



     !!
     !! comuting the spatial fourier transform for obtaining the
     !! dynamic structure factor, nearest-neighbor convention used
     !! S( q , t ) Sum_{R,R'} exp(iq(R(t)-R'(t)))
     !! the modes are projected on eigenvectors that have to be supplied as 
     !! input
     !!
     Subroutine ComputePVACFMain( data , atoms , lattice )

       use tool_interfaces, only:get_nearest_image,get_norm
       use CrossCorrComplex, only :CrossCorrelationMain
       Implicit None
       Type( PVACFVar ) ::data
       Real*8,dimension( : , : )  ::atoms
       Real*8,dimension( : , : )  ::lattice
       Integer ::qq
       Integer ::i
       Integer ::j
       Integer ::k
       Integer ::Indx
       Real*8 ::dpx
       Real*8 ::dpy
       Real*8 ::dpz

       Real*8 ::Dotp

       Complex*16,allocatable,dimension( : , : )  ::Proj



       If ( data%ActStep .gt. 1 ) then
         Call ComputeVelocity( data , atoms , lattice )
       End If



       Allocate( Proj( 1:data%NQgrid , 1:Size( data%BasisVectors , 3 ) ) )
       Proj = 0d0
       !$OMP Parallel Do Private(qq,j,i,k,indx,dpx,dpy,dpz,Dotp) Shared(atoms,data) &
       !$OMP& Reduction (+:Proj)
       Do qq = 1 , data%Nqgrid
         Do j = 1 , Size( data%BasisVectors , 3 )
            !!looping over unit cells
            Do i = 1 , Size( data%LinkUC , 1 )
               !!! looping over atoms in unit cell
               Do k = 1, Size( data%BasisVectors , 2 )
                  Indx =  data%LinkUC( i )%List( k )
                  dpx  =  data%QgridDir( 1 , qq ) * atoms( 1 , Indx ) 
                  dpy  =  data%QgridDir( 2 , qq ) * atoms( 2 , Indx ) 
                  dpz  =  data%QgridDir( 3 , qq ) * atoms( 3 , Indx ) 
                  !dpx  =  data%QgridDir( 1 , qq ) * data%Boxes( 1 , i )
                  !dpy  =  data%QgridDir( 2 , qq ) * data%Boxes( 2 , i )
                  !dpz  =  data%QgridDir( 3 , qq ) * data%Boxes( 3 , i )
                  Dotp = dpx + dpy + dpz
                  Proj( qq , j )  =  Proj( qq , j ) + &
                              Dot_product( data%BasisVectors( : , k , j , qq ) , &
                                             data%Velocity( : , k , i )  ) * CDexp( ImagI * Dotp )
               End Do
            End Do
         End Do
       End Do
       !$OMP End Parallel Do


       If ( data%ActStep .gt. 1 ) then
         Do i = 1 , Size( data%BasisVectors , 3 )
             Call CrossCorrelationMain( data%TimeCorrelationProj( i ) , Proj( : , i ) , & 
                                        Proj( : , i ) )
         End Do
       End If

       Deallocate( Proj )
       data%ActStep =  data%ActStep + 1

     End Subroutine ComputePVACFMain


     Subroutine ComputePVACFFinalize( data , fname )

       use tool_interfaces,only :give_closest_lower_power_of_2, get_norm,&
                                 gen_file_number
       use fourier_interface, only :fourier_dp
       use CrossCorrComplex,only:FinalizeCrossCorr
       use WindowingAnalysis, only:GaussianFiltering
       Implicit None
       Type( PVACFVar ) ::data
       Character( len = * ) ::fname
       Character( len = 3 ) ::fnum
       Complex*16,allocatable,dimension( : , : ) ::FourierArray

       Integer ::N2
       Integer ::qq
       Integer ::i
       Integer ::j
       Real*8 ::deltaFrequ
       Real*8 ::Qpath
       Real*8,allocatable,dimension( : ) ::Qdist


       !! compute average over time correltion functions
       If ( data%TimeCorrOnOff ) then
         Do i = 1 , Size( data%BasisVectors , 3 )
              Call FinalizeCrossCorr( data%TimeCorrelationProj( i ) )
              If ( i .eq. 1 ) then
                  Allocate( data%PVACF( 1:Size( data%TimeCorrelationProj( i )%Average , 1 ),&
                                                  1:Size( data%TimeCorrelationProj( i )%Average , 2 ),&
                                                  1:Size( data%BasisVectors , 3 ) ) )
              End If
              data%PVACF( : , : , i ) = data%PVACF( : , : , i ) +&
                                               data%TimeCorrelationProj( i )%Average( : , : )
         End Do
       End If

       !! apply filtering
       N2 =  give_closest_lower_power_of_2( Size( data%PVACF , 2 ) )
       If ( data%FilterFunction ) then
               Do i = 1 , N2
                  Do j = 1, Size( data%BasisVectors , 3 )
                     data%PVACF( : , i , j ) = &
                          data%PVACF( : , i , j ) * Dexp( -Dble( i ) * 6d0 / Dble( N2 ) )
                  End Do
               End Do
       End If


       !!! write time-dependent projected velovity autocorrelation
       Do j = 1, Size( data%BasisVectors, 3 )
         Call gen_file_number( j , fnum )
         open( unit = 302 , status='unknown' , file = "PVACF_vs_t.out"//Trim(fnum) , action='write' )
          Do i = 1 , N2
            Write( 302 , * ) Dble( i )*data%Tstep - data%Tstep / 2d0,&
                                          ( Abs( data%PVACF( qq , i , j ) )&
                                          , qq = 1, Size( data%PVACF( : , i , j ) ) )
          End Do
         close( 302 )
       End Do

       !!! fourier transform dynamic structure factor
       Allocate( FourierArray( 1:N2 , data%Nqgrid ) )
       Allocate( Qdist( 1 : Size( data%QgridDir , 2 ) ) )
       Do i = 1, Size( data%QgridDir , 2 )
          Qdist( i )  =  get_norm( MatMul( data%InvLattice , data%QgridDir( : , i ) ) )
       End Do

       Do j = 1 , Size( data%BasisVectors , 3 )
          FourierArray = Transpose( data%PVACF( 1:data%Nqgrid , 1:N2 , j ) )
          Do qq = 1 , data%Nqgrid
            Call fourier_dp( FourierArray( : , qq ) , -1 )
          End Do
          deltafrequ = 1d0 / ( 2d0 * data%Tstep ) / ( Dble( N2 ) / 2d0 )

          !! Write fourier transform without gaussian smoothing
          Call gen_file_number( j , fnum )
          open( unit = 302 , status='unknown' , file = Trim( fname )//Trim( fnum ) , action='write' )
          Do i = 1 , N2 / 2 
             Write( 302 , * ) ( deltafrequ * Dble( i ) - deltafrequ / 2d0 ),&
                                Abs( FourierArray( i , : ) ) / data%Natoms / PI2 !/ Qdist
          End Do
          close( 302 )
       End Do


       Deallocate( FourierArray )


       open( unit = 302 , action='write' , status='unknown' , file = "QpathPVACF.out" )
         Qpath = 0d0
         Write( 302  , "(F15.8)" ) Qpath
         Do i = 2 , Size( data%QgridDir , 2 )
            Qpath = Qpath + get_norm( MatMul( data%InvLattice ,&
                                              data%QgridDir( : , i )/PI2-data%QgridDir( : , i - 1 )/PI2 ) )
            Write( 302 , "(F15.8)" ) Qpath 
         End Do
       close( 302 )
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       Deallocate( data%PVACF )

     End Subroutine ComputePVACFFinalize

 End Module ComputePVACF



!!!
!!!                   ~~
!!!                 ~~~~~~
!!!              ~~~~~~~~~~~~
!!!            ~~~~~~~~~~~~~~~~
!!!  Projected Kspace velocity autocorrelation function
!!!  no basis vectors are used in this routine
!!!            ~~~~~~~~~~~~~~~~
!!!              ~~~~~~~~~~~~
!!!                   ~~

 Module ComputePVACFKProjOnly

   use CrossCorrComplexND, only :CrossCorrVarND
   use SplitCelltoBox , only :BoxList
   Implicit None
   Type PVACFVarKProj 
     Real*8,allocatable,dimension( : , : )  ::QgridDir                  !!! Qgrid in direct coordinates
     Integer ::Nqgrid                                                   !!! size qgrid
     Integer ::ActStep                                                  !!! actual MD step
     Integer ::MaxStep                                                  !!! total number of MD steps
     Real*8  ::Tstep                                                    !!! time step
     Real*8,dimensioN( 1:3 , 1:3 )  ::InvLattice                        !!! Inverse lattice
     Integer ::QPathNr
     Integer,allocatable,dimension( : )  ::Qpath                        !!! indices containing qpath
     Type( CrossCorrVarND ) ::TimeCorrelation
     Logical ::TimeCorrOnOff
     Integer ::Natoms
     Logical ::RotateCell   !! enables sqrt2 representation for a cubic cell

     !!! projected stuff
     Type( CrossCorrVarND ),allocatable,dimension( : ) ::TimeCorrelationProj
     Real*8,allocatable,dimension( : , : )  ::PosOld
     Real*8,allocatable,dimensioN( : , : , : )  ::Velocity 
     Complex*16,allocatable,dimension( : , : , : )  ::VACFKProj !!! VACF K space projection
     Type( BoxList ),allocatable,dimension( : )  ::LinkUC     !!! Link atoms to unit cells
     Real*8,allocatable,dimension( : , : )  ::Boxes           !!! centers for fourier transform
     Real*8,allocatable,dimension( : ) ::Masses
     Logical :: FilterFunction                                !!! checks if exponential filter function is applied
   End Type PVACFVarKProj

   Real*8,parameter ::PI2 = 4d0 * Dacos( 0d0 )                          !!! factor 2pi
   Complex*16,parameter ::ImagI = Cmplx( 0d0 , 1d0 )                    !!! imaginary number i
   Contains

     Subroutine ComputePVACFInitKProjOnly( data , Nsteps , dt , Nx ,&
                                           Ny , Nz , atoms , lattice , fname )

       use FlagReader, only :ReadFlags
       use tool_interfaces, only :get_norm
       use CrossCorrComplexND, only :CrossCorrelationInitND
       use QgridCollection, only :SampleFirstBrillouinZoneCubic
       use PhononHelperRoutines, only :GetMasses,&
                                       ReadQpointsFromFile

       Implicit None
       Type( PVACFVarKProj ) ::data
       Integer ::NSteps
       Real*8  ::dt
       Integer ::Nx
       Integer ::Ny
       Integer ::Nz
       Real*8,allocatable,dimension( : , : )  ::atoms
       Real*8,dimension( : , : ) ::lattice
       Character( len = * )      ::fname


       Logical ::QFound
       Real*8  ::temp

       Integer ::i


       data%MaxStep =  NSteps
       data%Tstep   =  dt
       data%InvLattice  =  lattice
       Call lu_inversion( data%InvLattice , 3 )

       Qfound = .False.
       Call ReadFlags( "Qpath" , temp , QFound , fname )
       If ( .not. Qfound ) then
         data%QPathNr = 1
       Else
         data%QPathNr = NINT( temp )
       End If

       data%Natoms = Size( atoms , 2 )


       Call ReadBoxList( data , Nx , Ny , Nz , Size( atoms , 2 ) , "BoxList.in" )
       Call GetMasses( data%Masses , "masses.in" )



       Inquire( File="QVectors.in" , exist = Qfound )
       If ( Qfound ) then
          Write( * , * ) "In routine ComputePVACFInitKProjOnly Q-vectors are read"
          Call ReadQpointsFromFile( data%QGridDir , "QVectors.in" )
          data%NQgrid  =  Size( data%QgridDir , 2 )
       Else
          Write( * , * ) "In routine ComputePVACFInitKProjOnly no QVectors.in"
          Write( * , * ) "supplied. Using default Q-mesh of first irreducible "
          Write( * , * ) "Brillouin zone"
          Call SampleFirstBrillouinZoneCubic( data%QgridDir , Nx , Ny , Nz )
          data%NQgrid  =  Size( data%QgridDir , 2 )
       End If


       data%QgridDir = PI2 * data%QgridDir 
       data%PosOld  =  atoms
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       

       data%TimeCorrOnOff = .False.
       !! supply window size
       Call ReadFlags( "DSTC" , temp , data%TimeCorrOnOff , fname )
       If ( .not. data%TimeCorrOnOff ) then
          temp  =  Real( NSteps - 1 )
          data%TimeCorrOnOff = .True.
       End If

       If ( data%TimeCorrOnOff ) then
          Write( * , * ) "Time Averaging of the K-Space velocity autocorrelation"
          Write( * , * ) "is switched on. The window size was set to"
          Write( * , * ) NINT( temp )
          !! Projected stuff
          !! -1 because velocities are needed
          Allocate( data%TimeCorrelationProj( 1:Size( data%LinkUC( 1 )%List , 1  ) ) )
          Do i = 1 , Size( data%LinkUC( 1 )%List , 1 )
             Call CrossCorrelationInitND( data%TimeCorrelationProj( i ) ,&
                       NINT( temp ) , NSteps-1 , Size( data%QgridDir , 2 ) )
          End Do
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       End If


     End Subroutine ComputePVACFInitKProjOnly
         


     Subroutine ReadBoxList( data , N1 , N2 , N3 , Natoms , fname )

       use SplitCelltoBox, only :MakeBoxes,MakeBoxList
       Implicit None
       Type( PVACFVarKProj ) ::data
       Integer  ::N1
       Integer  ::N2
       Integer  ::N3
       Integer  ::Natoms
       Character( len = * )  ::fname
       Integer  ::Npc
       Real*8,allocatable,dimension( : , : )  ::Boxes

       Integer ::i


       Npc = Natoms / ( N1 * N2 * N3 )
       
       
       Boxes = MakeBoxes( N1 , N2 , N3 )


       Allocate( data%LinkUC( 1:N1*N2*N3 ) )


       open( unit = 33 , action='read' , file = Trim( fname ) , status='old' )
       Do i = 1 , Size( data%LinkUC , 1 )
          Allocate( data%LinkUC( i )%List( 1:Npc ) )
          data%LinkUC( i )%Natoms = Npc
          Read( 33 , * ) data%LinkUC( i )%List( : )
       End Do
       close( 33 )


       
       Allocate( data%Velocity( 1:3 , 1:NPc , 1:N1*N2*N3 ) )
       Allocate( data%Boxes( 1:3 , 1 : N1 * N2 * N3 ) )
       Do i = 1 , Size( data%Boxes , 2 )
          data%Boxes( 1 , i )  =  0.5d0 * ( Boxes( 1 , i ) + Boxes( 2 , i ) )
          data%Boxes( 2 , i )  =  0.5d0 * ( Boxes( 3 , i ) + Boxes( 4 , i ) )
          data%Boxes( 3 , i )  =  0.5d0 * ( Boxes( 5 , i ) + Boxes( 6 , i ) )
       End Do


     End Subroutine ReadBoxList
         


     Subroutine ComputeVelocity( data , PosNew , lattice )

       use tool_interfaces, only :get_nearest_image,get_norm,&
                                  get_norm
       Implicit None
       Type( PVACFVarKProj ) ::data
       Real*8,dimension( : , : )     ::PosNew
       Real*8,dimension( : , : ) ::lattice

       Integer ::i
       Integer ::j
       Integer ::Indx

       Do i = 1 , Size( data%LinkUC , 1 )
         Do j = 1 , data%LinkUC( i )%Natoms
            Indx = data%LinkUC( i )%List( j )
            data%Velocity( : , j , i )  =  MatMul( lattice , get_nearest_image( data%PosOld( : , Indx ) ,&
                                           PosNew( : , Indx ) ) - data%PosOld( : , Indx ) )
            data%Velocity( : , j , i )  =  data%Velocity( : , j , i ) * data%Masses( j )
            data%PosOld( : , Indx )     =  PosNew( : , Indx )
         End Do
       End Do

     End Subroutine ComputeVelocity


     !!
     !! comuting the spatial fourier transform for obtaining the
     !! dynamic structure factor, nearest-neighbor convention used
     !! S( q , t ) Sum_{R,R'} exp(iq(R(t)-R'(t)))
     !! the modes are projected on eigenvectors that have to be supplied as 
     !! input
     !!
     Subroutine ComputePVACFMainKProjOnly( data , atoms , lattice )

       use tool_interfaces, only:get_nearest_image,get_norm
       use CrossCorrComplexND, only :CrossCorrelationMainND
       Implicit None
       Type( PVACFVarKProj ) ::data
       Real*8,dimension( : , : )  ::atoms
       Real*8,dimension( : , : )  ::lattice
       Integer ::qq
       Integer ::i
       Integer ::k
       Integer ::Indx
       Real*8 ::dpx
       Real*8 ::dpy
       Real*8 ::dpz

       Real*8 ::Dotp

       Complex*16,allocatable,dimension( : , : , : )::Proj





       If ( data%ActStep .gt. 1 ) then
         Call ComputeVelocity( data , atoms , lattice )
       End If


       Allocate( Proj( 1:3 , 1:data%NQgrid , Size( data%LinkUC( 1 )%List , 1 ) ) )
       Proj  = 0d0
       !$OMP Parallel Do Private(qq,i,k,indx,dpx,dpy,dpz,Dotp) Shared(atoms,data) &
       !$OMP& Reduction (+:Proj)
       Do qq = 1 , data%Nqgrid
          !!looping over unit cells
          Do i = 1 , Size( data%LinkUC , 1 )
             !!! looping over atoms in unit cell
             Do k = 1, Size( data%LinkUC( i )%List( : ) , 1 )
                Indx =  data%LinkUC( i )%List( k )
                dpx  =  data%QgridDir( 1 , qq ) * atoms( 1 , Indx )
                dpy  =  data%QgridDir( 2 , qq ) * atoms( 2 , Indx )
                dpz  =  data%QgridDir( 3 , qq ) * atoms( 3 , Indx )
                !dpx   =   data%QgridDir( 1 , qq ) * data%Boxes( 1 , i )
                !dpy   =   data%QgridDir( 2 , qq ) * data%Boxes( 2 , i )
                !dpz   =   data%QgridDir( 3 , qq ) * data%Boxes( 3 , i )
                Dotp =  dpx + dpy + dpz
                Proj( : , qq , k )  =  Proj( : , qq , k ) + data%Velocity( : , k , i ) * CDexp( ImagI * Dotp )
             End Do
          End Do
       End Do
       !$OMP End Parallel Do



       If ( data%ActStep .gt. 1 ) then
         Do i = 1 , Size( data%LinkUC( 1 )%List , 1 )  !! loop over atoms in UC -> non projected branches
            Call CrossCorrelationMainND( data%TimeCorrelationProj( i ) , Proj( : , : , i ) , & 
                                         Proj( : , : , i ) )
         End Do
       End If
       data%ActStep =  data%ActStep + 1

     End Subroutine ComputePVACFMainKProjOnly

     Subroutine ComputePVACFFinalizeKProjOnly( data , fname )

       use tool_interfaces,only :give_closest_lower_power_of_2, get_norm,&
                                 gen_file_number
       use fourier_interface, only :fourier_dp
       use CrossCorrComplexND,only:FinalizeCrossCorrND
       use WindowingAnalysis, only:GaussianFiltering
       Implicit None
       Type( PVACFVarKProj ) ::data
       Character( len = * ) ::fname
       Character( len = 3 ) ::fnum
       Complex*16,allocatable,dimension( : , : ) ::FourierArray

       Integer ::N2
       Integer ::qq
       Integer ::i
       Integer ::j
       Real*8 ::deltaFrequ
       Real*8 ::Qpath
       Real*8,allocatable,dimension( : ) ::Qdist


       !! compute average over time correltion functions

       If ( data%TimeCorrOnOff ) then
         Do i = 1 , Size( data%LinkUC( 1 )%List , 1 )
            Call FinalizeCrossCorrND( data%TimeCorrelationProj(i) )
            If ( i .eq. 1 ) then
               Allocate( data%VACFKProj( 1:Size( data%TimeCorrelationProj( i )%Average , 1 ),&
                                               1:Size( data%TimeCorrelationProj( i )%Average , 2 ),&
                                               1:Size( data%LinkUC( 1 )%List , 1 ) ) )
            End If
            data%VACFKProj( : , : , i ) = data%VACFKProj( : , : , i ) +&
                                                   data%TimeCorrelationProj(i)%Average( : , : )
         End Do
       End If

       N2 =  give_closest_lower_power_of_2( Size( data%VACFKProj , 2 ) )
       If ( data%FilterFunction ) then
               Do j = 1 , Size( data%LinkUC( 1 )%List , 1 )
                  Do i = 1 , N2
                        data%VACFKProj( : , i , j ) = &
                             data%VACFKProj( : , i , j ) * Dexp( -Dble( i ) * 6d0 / Dble( N2 ) )
                  End Do
               End Do
        End If



       !!! write time-dependent projected velocity autocorrelation S(q,t)
       Do j = 1, Size( data%LinkUC( 1 )%List , 1 )
         Call gen_file_number( j , fnum )
         open( unit = 302 , status='unknown' , file = "PVACF_KO_vs_t.out"//Trim(fnum) , action='write' )
          Do i = 1 , N2
            Write( 302 , * ) Dble( i ) * data%Tstep - data%Tstep/2d0,&
                                           ( Abs( data%VACFKProj( qq , i , j ) ) ,&
                                             qq = 1, Size( data%VACFKProj( : , i , j ) ) )
          End Do
         close( 302 )
       End Do

       !!! fourier transform dynamic structure factor
       Allocate( FourierArray( 1:N2 , data%Nqgrid ) )
       Allocate( Qdist( 1 : Size( data%QgridDir , 2 ) ) )
       Do j = 1, Size( data%LinkUC( 1 )%List , 1 )
           Do i = 1, Size( data%QgridDir , 2 )
              Qdist( i )  =  get_norm( MatMul( data%InvLattice , data%QgridDir( : , i ) ) )
           End Do

           FourierArray = Transpose( data%VACFKProj( 1:data%Nqgrid , 1:N2 , j ) )
           Do qq = 1 , data%Nqgrid
             Call fourier_dp( FourierArray( : , qq ) , -1 )
           End Do
           deltafrequ = 1d0 / ( 2d0 * data%Tstep ) / ( Dble( N2 ) / 2d0 )

           !! Write fourier transform without gaussian smoothing
           Call gen_file_number( j , fnum )
           open( unit = 302 , status='unknown' , file = Trim( fname ) // Trim( fnum ) , action='write' )
           Do i = 1 , N2 / 2 
              Write( 302 , * ) ( deltafrequ * Dble( i ) - deltafrequ / 2d0 ),&
                                 Abs( FourierArray( i , : ) ) / data%Natoms / PI2
           End Do
           close( 302 )
       End Do
       Deallocate( FourierArray )


       open( unit = 302 , action='write' , status='unknown' , file = "QpathPVACF_KO.out" )
         Qpath = 0d0
         Write( 302  , "(F15.8)" ) Qpath
         Do i = 2 , Size( data%QgridDir , 2 )
            Qpath = Qpath + get_norm( MatMul( data%InvLattice ,&
                                              data%QgridDir( : , i )/PI2-data%QgridDir( : , i - 1 )/PI2 ) )
            Write( 302 , "(F15.8)" ) Qpath 
         End Do
       close( 302 )
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       Deallocate( data%VACFKProj )

     End Subroutine ComputePVACFFinalizeKProjOnly

 End Module ComputePVACFKProjOnly



!!!
!!!                   ~~
!!!                 ~~~~~~
!!!              ~~~~~~~~~~~~
!!!            ~~~~~~~~~~~~~~~~
!!!  Dynamic structure factor computation
!!!            ~~~~~~~~~~~~~~~~
!!!              ~~~~~~~~~~~~
!!!                 ~~~~~~
!!!                   ~~

 Module ComputeDynamicStructureFactor


   use CrossCorrComplex, only :CrossCorrVar
   use SplitCelltoBox , only :BoxList
   Implicit None
   Type DynStrucVar 
     Real*8,allocatable,dimension( : , : )  ::QgridDir                  !!! Qgrid in direct coordinates
     Complex*16,allocatable,dimension( : , : )  ::DynStrucFac           !!! dynamic structure factor in time domain
     Integer ::Nqgrid                                                   !!! size qgrid
     Integer ::ActStep                                                  !!! actual MD step
     Integer ::MaxStep                                                  !!! total number of MD steps
     Real*8  ::Tstep                                                    !!! time step
     Real*8,dimensioN( 1:3 , 1:3 )  ::InvLattice                        !!! Inverse lattice
     Integer ::QPathNr
     Integer,allocatable,dimension( : )  ::Qpath                        !!! indices containing qpath
     Type( CrossCorrVar ) ::TimeCorrelation
     Logical ::TimeCorrOnOff
     Integer ::Natoms
     Logical :: FilterFunction                                !!! checks if exponential filter function is applied
     Real*8,allocatable,dimension( : ) ::ScatteringLength     !!! scattering length
     Type( BoxList ),allocatable,dimension( : )  ::LinkUC     !!! Link atoms to unit cells
     Real*8,allocatable,dimension( : , : )  ::Boxes           !!! centers for fourier transform
   End Type DynStrucVar


   Real*8,parameter ::PI2 = 4d0 * Dacos( 0d0 )                          !!! factor 2pi
   Complex*16,parameter ::ImagI = Cmplx( 0d0 , 1d0 )                    !!! imaginary number i
   Contains

     Subroutine ComputeDynamicStructureFactorInit( data , Nsteps , dt , Nx , Ny ,&
                                                       Nz , atoms , lattice , fname )

       use FlagReader, only :ReadFlags
       use CrossCorrComplex, only :CrossCorrelationInit
       use QgridCollection , only :SampleFirstBrillouinZoneCubic
       use PhononHelperRoutines, only :ReadQpointsFromFile,&
                                       ReadScatteringLengths
       Implicit None
       Type( DynStrucVar ) ::data
       Integer ::NSteps
       Real*8  ::dt
       Integer ::Nx
       Integer ::Ny
       Integer ::Nz
       Real*8,allocatable,dimension( : , : ) ::atoms
       Real*8,dimension( : , : ) ::lattice
       Character( len = * ) ::fname
       Logical ::CheckQpointsFile


       Logical ::ScatterFound
       Logical ::QFound
       Real*8  ::temp


       data%MaxStep    =  NSteps
       data%Tstep      =  dt
       data%InvLattice =  lattice
       Call lu_inversion( data%InvLattice , 3 )

       ScatterFound  =  .False.
       Qfound = .False.
       Call ReadFlags( "Qpath" , temp , QFound , fname )
       If ( .not. Qfound ) then
         data%QPathNr = 1
       Else
         data%QPathNr = NINT( temp )
       End If

       data%Natoms = Size( atoms , 2 )
       
       !Call ReadBoxList( data , Nx , Ny , Nz , Size( atoms , 2 ) , "BoxList.in" )

       INQUIRE( File = "QVectors.in" , Exist = CheckQpointsFile )
       If ( CheckQpointsFile ) then
         Write( * , * ) "In routine ComputeDynamicStructureFactorInit Q-vectors are read"
         Call ReadQpointsFromFile( data%QGridDir , "QVectors.in" )
         data%NQgrid  =  Size( data%QgridDir , 2 )
       Else
         Write( * , * ) "In routine ComputeDynamicStructureFactorInit no QVectors.in"
         Write( * , * ) "supplied. Using default Q-mesh of first irreducible "
         Write( * , * ) "Brillouin zone"
         Call SampleFirstBrillouinZoneCubic( data%QgridDir , Nx , Ny , Nz )
         data%NQgrid  =  Size( data%QgridDir , 2 )
       End If
       

       Inquire( File = "ScatteringL.in" , exist = ScatterFound )
       If ( ScatterFound ) then
               Call ReadScatteringLengths( data%ScatteringLength , "ScatteringL.in" )
       Else
               Write( * , * ) "No ScateringL.in file supplied"
               Write( * , * ) "Setting all scattering lengths to 1"
               Allocate( data%ScatteringLength( &
                         1:Size( Atoms , 2 ) / ( Nx * Ny * Nz ) ) )
       End If

       data%QgridDir = PI2 * data%QgridDir


       data%TimeCorrOnOff = .True.
       !! supply window size
       Call ReadFlags( "DSTC" , temp , data%TimeCorrOnOff , fname )
       If ( .not. data%TimeCorrOnOff ) then
          temp = Real( NSteps )
          data%TimeCorrOnOff = .True.
       End If

       If ( data%TimeCorrOnOff ) then
          Write( * , * ) "Time Averaging of the dynamic structure factor"
          Write( * , * ) "is switched on. The window size was set to"
          Write( * , * ) NINT( temp )
          Call CrossCorrelationInit( data%TimeCorrelation , NINT( temp ) ,&
                                     NSteps , Size( data%QgridDir , 2 ) )
       End If


     End Subroutine ComputeDynamicStructureFactorInit
     

     Subroutine ReadBoxList( data , N1 , N2 , N3 , Natoms , fname )

        use SplitCelltoBox, only :MakeBoxes,MakeBoxList
        Implicit None
        Type( DynStrucVar ) ::data
        Integer  ::N1
        Integer  ::N2
        Integer  ::N3
        Integer  ::Natoms
        Character( len = * ) ::fname
        Integer  ::Npc
        Real*8,allocatable,dimension( : , : )  ::Boxes

        Integer ::i
        

        Npc = Natoms / ( N1 * N2 * N3 )
       
        Boxes = MakeBoxes( N1 , N2 , N3 )


        Allocate( data%LinkUC( 1:N1*N2*N3 ) )
        open( unit = 33 , action='read' , file = Trim( fname ) , status='old' )
        Do i = 1 , Size( data%LinkUC , 1 )
           Allocate( data%LinkUC( i )%List( 1:Npc ) )
           data%LinkUC( i )%Natoms = Npc
           Read( 33 , * ) data%LinkUC( i )%List( : )
        End Do
        close( 33 )


        
        Allocate( data%Boxes( 1:3 , 1 : N1 * N2 * N3 ) )
        Do i = 1 , Size( data%Boxes , 2 )
           data%Boxes( 1 , i )  =  0.5d0 * ( Boxes( 1 , i ) + Boxes( 2 , i ) )
           data%Boxes( 2 , i )  =  0.5d0 * ( Boxes( 3 , i ) + Boxes( 4 , i ) )
           data%Boxes( 3 , i )  =  0.5d0 * ( Boxes( 5 , i ) + Boxes( 6 , i ) )
        End Do

     End Subroutine ReadBoxList


     !!
     !! computing the spatial fourier transform for obtaining the
     !! dynamic structure factor, nearest-neighbor convention used
     !! S( q , t ) Sum_{R,R'} exp(iq(R(t)-R'(t)))
     !! this does the same like routine above; but one loop is contracted
     !! here what makes it much faster
     Subroutine ComputeDynamicStructureFactorMain( data , atoms , lattice )

       use tool_interfaces, only:get_nearest_image,get_norm
       use CrossCorrComplex, only :CrossCorrelationMain
       Implicit None
       Type( DynStrucVar ) ::data
       Real*8,allocatable,dimension( : , : )  ::atoms
       Real*8,dimension( : , : )  ::lattice
       Integer ::qq
       Integer ::i
       Integer ::k
       Real*8 ::dpx
       Real*8 ::dpy
       Real*8 ::dpz
       Real*8 ::Dotp
       Integer ::AIndx

       Complex*16,allocatable,dimension( : )  ::P


       Allocate( P( 1:data%Nqgrid ) )
       P = 0d0
       !$OMP Parallel Do Private(qq,i,dpx,dpy,dpz,Dotp) Shared(atoms,data) &
       !$OMP& Reduction (+:P)
       Do qq = 1 , data%Nqgrid
          Do i = 1 , Size( atoms , 2 )
             dpx    =  data%QgridDir( 1 , qq ) * atoms( 1 , i ) 
             dpy    =  data%QgridDir( 2 , qq ) * atoms( 2 , i ) 
             dpz    =  data%QgridDir( 3 , qq ) * atoms( 3 , i ) 
             Dotp   =  dpx + dpy + dpz
             P( qq )  =  P( qq ) + &
                         !data%ScatteringLength( k ) * CDexp( ImagI * Dotp )
                         CDexp( ImagI * Dotp )
          End Do
       End Do
       !$OMP End Parallel Do


       Call CrossCorrelationMain( data%TimeCorrelation , P )
       data%ActStep =  data%ActStep + 1
       Deallocate( P )

     End Subroutine ComputeDynamicStructureFactorMain


     Subroutine ComputeDynamicStructureFactorFinalize( data , fname )

       use tool_interfaces,only :give_closest_lower_power_of_2, get_norm
       use fourier_interface, only :fourier_dp
       use CrossCorrComplex,only:FinalizeCrossCorr
       use BaseStatistics, only:Mean
       use WindowingAnalysis, only:GaussianFiltering
       Implicit None
       Type( DynStrucVar ) ::data
       Character( len = * ) ::fname
       Complex*16,allocatable,dimension( : , : ) ::FourierArray

       Integer ::N2
       Integer ::qq
       Integer ::i
       Real*8 ::deltaFrequ
       Real*8,allocatable,dimension( : )  ::Average
       Real*8 ::Qpath
       Real*8,allocatable,dimension( : ) ::Qdist


       !! compute average over time correltion functions
       If ( data%TimeCorrOnOff ) then
         Call FinalizeCrossCorr( data%TimeCorrelation )
         Allocate( data%DynStrucFac( 1:Size( data%TimeCorrelation%Average , 1 ),&
                                     1:Size( data%TimeCorrelation%Average , 2 ) ) )
         data%DynStrucFac = data%TimeCorrelation%Average
       End If

       N2 =  give_closest_lower_power_of_2( Size( data%DynStrucFac , 2 ) )
       !! apply filter
       If ( data%FilterFunction ) then
               Do i = 1 , N2
                  data%DynStrucFac( : , i )   =  data%DynStrucFac( : , i ) * Dexp( -Dble( i ) * 6d0 / Dble( N2 ) )
               End Do
       End If


       !! store time average of Dynamic Structure factor
       !! -> Static structure factor
       Allocate( Average( 1:data%NQgrid ) )
       Average  =  Mean( Abs( data%DynStrucFac ) )

       !!! write time-dependent structure factor S(q,t)
       open( unit = 302 , status='unknown' , file = "StructureFactor_vs_t.out" , action='write' )
         Do i = 1 , N2
            Write( 302 , * ) Dble( i )*data%Tstep-data%Tstep/2d0 , &
                                           ( Abs( data%DynStrucFac( qq , i ) ) &
                                           , qq = 1, Size( data%DynStrucFac( : , i ) ) )
         End Do
       close( 302 )
       

       open( unit = 302 , status='unknown' , file = "StructureFactorComplex_vs_t.out" , action='write' )
         Do i = 1 , N2
            Write( 302 , * ) Dble( i )*data%Tstep-data%Tstep/2d0 , &
                                           ( Dble( data%DynStrucFac( qq , i ) ) , &
                                             DImag( data%DynStrucFac( qq , i ) ) &
                                           , qq = 1, Size( data%DynStrucFac( : , i ) ) )
         End Do
       close( 302 )

       !!! fourier transform dynamic structure factor
       Allocate( FourierArray( 1:N2 , data%Nqgrid ) )
       FourierArray = Transpose( data%DynStrucFac( 1:data%Nqgrid , 1:N2 ) )


       !!! Write the fourier array to a file
       !Do i = 1 , N2
       !     Write( 302 , * ) Dble( i )*data%Tstep-data%Tstep/2d0 , &
       !                                    ( Abs( ) &
       !                                    , qq = 1, Size( data%DynStrucFac( : , i ) ) )
       !End Do





       Do qq = 1 , data%Nqgrid
         Call fourier_dp( FourierArray( : , qq ) , -1 )
       End Do
       deltafrequ = 1d0 / ( 2d0 * data%Tstep ) / ( Dble( N2 ) / 2d0 )


       !!! compute q grid
       Allocate( Qdist( 1 : Size( data%QgridDir , 2 ) ) )
       Do i = 1, Size( data%QgridDir , 2 )
          Qdist( i )  =  get_norm( MatMul( data%InvLattice , data%QgridDir( : , i ) ) )
       End Do

       !! Write fourier transform without gaussian smoothing
       open( unit = 302 , status='unknown' , file = Trim( fname ) , action='write' )
       Do i = 1 , N2 / 2 
          Write( 302 , * ) ( deltafrequ * Dble( i ) - deltafrequ / 2d0 ),&
                             Abs( FourierArray( i , : ) ) / data%Natoms / PI2 
       End Do
       close( 302 )
       Deallocate( FourierArray )


       open( unit = 302 , action='write' , status='unknown' , file = "DSFQpath.out" )
       open( unit = 303 , action='write' , status ='unknown' , file = "StructureFactor.out" )
         Qpath = 0d0
         Write( 302 , * ) Qpath
         Write( 303 , * ) Qpath , Average( 1 )
         Do i = 2 , Size( data%QgridDir , 2 )
            Qpath = Qpath + get_norm( MatMul( data%InvLattice ,&
                                              data%QgridDir( : , i )/PI2-data%QgridDir( : , i - 1 )/PI2 ) )
            Write( 302 , * ) Qpath 
            Write( 303 , * ) Qpath , Average( i )
         End Do
       close( 302 )
       close( 303 )
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


       Deallocate( data%DynStrucFac )

     End Subroutine ComputeDynamicStructureFactorFinalize

 End Module ComputeDynamicStructureFactor


!!!
!!!                   ~~
!!!                 ~~~~~~
!!!              ~~~~~~~~~~~~
!!!            ~~~~~~~~~~~~~~~~
!!!     Projected Dynamic Structure factor
!!!            ~~~~~~~~~~~~~~~~
!!!              ~~~~~~~~~~~~
!!!                 ~~~~~~
!!!                   ~~
 

 Module ComputeProjectedDynamicStructureFactor

   use CrossCorrComplex, only :CrossCorrVar
   use SplitCelltoBox , only :BoxList
   Implicit None
   Type ProjDynFacVar 
     Real*8,allocatable,dimension( : , : )  ::QgridDir                  !!! Qgrid in direct coordinates
     Real*8,allocatable,dimension( : , : )  ::QgridCart                 !!! Qgrid in cartesian coordinates
     Integer ::Nqgrid                                                   !!! size qgrid
     Integer ::ActStep                                                  !!! actual MD step
     Integer ::MaxStep                                                  !!! total number of MD steps
     Real*8  ::Tstep                                                    !!! time step
     Real*8,dimension( 1:3 , 1:3 )  ::lattice                           !!! lattice
     Real*8,dimensioN( 1:3 , 1:3 )  ::InvLattice                        !!! Inverse lattice
     Complex*16,allocatable,dimension( : ) ::Pinit                      !!! spec
     Complex*16,allocatable,dimension( : ) ::P
     Real*8,allocatable,dimension( : ) ::norm
     Integer ::QPathNr
     Integer,allocatable,dimension( : )  ::Qpath                        !!! indices containing qpath
     Type( CrossCorrVar ) ::TimeCorrelation
     Logical ::TimeCorrOnOff
     Integer ::Natoms

     !!! projected stuff
     Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
     Type( CrossCorrVar ),allocatable,dimension( : ) ::TimeCorrelationProj
     Complex*16,allocatable,dimension( : , : , : )  ::DynStrucFacProj
     Real*8,allocatable,dimension( : , : )  ::PosOld
     Real*8,allocatable,dimensioN( : , : , : )  ::Velocity 
     Type( BoxList ),allocatable,dimension( : )  ::LinkUC     !!! Link atoms to unit cells
     Real*8,allocatable,dimension( : , : )  ::Boxes           !!! centers for fourier transform
     Real*8,allocatable,dimension( : ) ::Masses
     Logical :: FilterFunction                                !!! checks if exponential filter function is applied
   End Type ProjDynFacVar

   Real*8,parameter ::PI2 = 4d0 * Dacos( 0d0 )                          !!! factor 2pi
   Complex*16,parameter ::ImagI = Cmplx( 0d0 , 1d0 )                    !!! imaginary number i
   Contains

     Subroutine ComputeProjectedDynamicStructureFactorInit( data , Nsteps , dt , Nx ,&
                                 Ny , Nz , atoms , lattice , fname )

       use FlagReader, only :ReadFlags
       use CrossCorrComplex, only :CrossCorrelationInit
       use QgridCollection, only :SampleFirstBrillouinZoneCubic
       use PhononHelperRoutines, only :GetMasses,&
                                       ReadEigenvectorsModeC,&
                                       MakeDiagonalBasis,&
                                       ReadQpointsFromFile



       Implicit None
       Type( ProjDynFacVar ) ::data
       Integer ::NSteps
       Real*8  ::dt
       Integer ::Nx
       Integer ::Ny
       Integer ::Nz
       Real*8,allocatable,dimension( : , : )  ::atoms
       Real*8,dimension( : , : ) ::lattice
       Character( len = * ) ::fname


       Logical ::QFound
       Real*8  ::temp
       Integer ::i


       data%MaxStep =  NSteps
       data%Tstep   =  dt
       data%lattice = lattice
       data%InvLattice =  data%lattice
       Call lu_inversion( data%InvLattice , 3 )

       Qfound = .False.
       Call ReadFlags( "Qpath" , temp , QFound , fname )
       If ( .not. Qfound ) then
         data%QPathNr = 1
       Else
         data%QPathNr = NINT( temp )
       End If

       data%Natoms = Size( atoms , 2 )

       Call ReadBoxList( data , Nx , Ny , Nz , Size( atoms , 2 ) , "BoxList.in" )
       Call GetMasses( data%Masses , "masses.in" )

       !Inquire( File = "BasisVector.in" , exist = QFound )
       !If ( QFound ) then
       !   Write( * , * ) "ComputeProjectedDynamicStructureFactorInit"
       !   Write( * , * ) "is reading BasisVectors and Qpoints"
       !   Write( * , * ) "from file BasisVector.in"
       !   Call ReadEigenvectorsModeC( data%BasisVectors , data%QGridDir , "BasisVector.in" )
       !   data%NQgrid = Size( data%QGridDir , 2 )
       !Else
       !   Write( * , * ) "ComputeProjectedDynamicStructureFactorInit"
       !   Write( * , * ) "No BasisVector.in file found"
       !   Write( * , * ) "Checking for a QVectorsFile"
       Inquire( File="QVectors.in" , exist = Qfound )
       If ( Qfound ) then
          Write( * , * ) "In routine ComputeProjectedDynamicStructureFactorInit"
          Write( * , * ) "Q-vectors are read"
          Write( * , * ) "QVectors.in"
          Call ReadQpointsFromFile( data%QGridDir , "QVectors.in" )
          data%NQgrid  =  Size( data%QgridDir , 2 )
       Else
          Write( * , * ) "In routine ComputeProjectedDynamicStructureFactorInit"
          Write( * , * ) "no QVectors.in"
          Write( * , * ) "supplied. Using default Q-mesh of first irreducible "
          Write( * , * ) "Brillouin zone"
          Call SampleFirstBrillouinZoneCubic( data%QgridDir , Nx , Ny , Nz )
          data%NQgrid  =  Size( data%QgridDir , 2 )
       End If
       !   Write( * , * ) "Using a diagonal basis since no file was supplied"
       !   Write( * , * ) "Diagonal basis gives some decomposition but does"
       !   Write( * , * ) "make sense physically"
       !   !Call MakeDiagonalBasis( 3 , data%LinkUC( 1 )%Natoms , data%Nqgrid , data%BasisVectors )
       !End If
       !!
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       data%QgridDir = PI2 * data%QgridDir
       !! projected dynamic structure factor
       data%PosOld = atoms
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       data%TimeCorrOnOff = .False.
       !! supply window size
       Call ReadFlags( "DSTC" , temp , data%TimeCorrOnOff , fname )
       If ( .not. data%TimeCorrOnOff ) then
          temp  =  Real( NSteps - 1 )
          data%TimeCorrOnOff = .True.
       End If
       If ( data%TimeCorrOnOff ) then
          Write( * , * ) "Time Averaging of the projected dynamic structure factor"
          Write( * , * ) "is switched on. The window size was set to"
          Write( * , * ) NINT( temp )
          Allocate( data%TimeCorrelationProj( 1:6 ) )  !! allocate for longitudinal and transverse part
          !! Projected stuff
          Do i = 1, 6
             !! -1 because velocities are needed
             Call CrossCorrelationInit( data%TimeCorrelationProj( i ) ,&
                 NINT( temp ) , NSteps-1 , Size( data%QgridDir , 2 ) )
          End Do
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       End If


     End Subroutine ComputeProjectedDynamicStructureFactorInit

     Subroutine ReadBoxList( data , N1 , N2 , N3 , Natoms , fname )

       use SplitCelltoBox, only :MakeBoxes,MakeBoxList
       Implicit None
       Type( ProjDynFacVar ) ::data
       Integer  ::N1
       Integer  ::N2
       Integer  ::N3
       Integer  ::Natoms
       Character( len = * )  ::fname
       Integer  ::Npc
       Real*8,allocatable,dimension( : , : )  ::Boxes

       Integer ::i


       Npc = Natoms / ( N1 * N2 * N3 )
       
       Boxes = MakeBoxes( N1 , N2 , N3 )


       Allocate( data%LinkUC( 1:N1*N2*N3 ) )
       open( unit = 33 , action='read' , file = Trim( fname ) , status='old' )
       Do i = 1 , Size( data%LinkUC , 1 )
          Allocate( data%LinkUC( i )%List( 1:Npc ) )
          data%LinkUC( i )%Natoms = Npc
          Read( 33 , * ) data%LinkUC( i )%List( : )
       End Do
       close( 33 )


       
       Allocate( data%Velocity( 1:3 , 1:NPc , 1:N1*N2*N3 ) )
       Allocate( data%Boxes( 1:3 , 1 : N1 * N2 * N3 ) )
       Do i = 1 , Size( data%Boxes , 2 )
          data%Boxes( 1 , i )  =  0.5d0 * ( Boxes( 1 , i ) + Boxes( 2 , i ) )
          data%Boxes( 2 , i )  =  0.5d0 * ( Boxes( 3 , i ) + Boxes( 4 , i ) )
          data%Boxes( 3 , i )  =  0.5d0 * ( Boxes( 5 , i ) + Boxes( 6 , i ) )
       End Do


     End Subroutine ReadBoxList
       

     Subroutine ComputeVelocity( data , PosNew , lattice )
             use tool_interfaces, only :get_nearest_image,get_norm,&
                             get_norm
             Implicit None
             Type( ProjDynFacVar ) ::data
             Real*8,dimension( : , : )     ::PosNew
             Real*8,dimension( : , : ) ::lattice

             Integer ::i
             Integer ::j
             Integer ::Indx

             Do i = 1 , Size( data%LinkUC , 1 )
               Do j = 1 , data%LinkUC( i )%Natoms
                  Indx = data%LinkUC( i )%List( j )
                  data%Velocity( : , j , i )  =  MatMul( lattice , get_nearest_image( data%PosOld( : , Indx ) ,&
                                                 PosNew( : , Indx ) ) - data%PosOld( : , Indx ) )
                  data%Velocity( : , j , i )  =  data%Velocity( : , j , i ) * data%Masses( j )
                  data%PosOld( : , Indx )     =  PosNew( : , Indx )
               End Do
             End Do

     End Subroutine ComputeVelocity



     !!
     !! comuting the spatial fourier transform for obtaining the
     !! dynamic structure factor, nearest-neighbor convention used
     !! S( q , t ) Sum_{R,R'} exp(iq(R(t)-R'(t)))
     !! the modes are projected on eigenvectors that have to be supplied as 
     !! input
     !!
     Subroutine ComputeProjectedDynamicStructureFactorMain( data , atoms , lattice )

       use tool_interfaces, only:get_nearest_image,get_norm
       use CrossCorrComplex, only :CrossCorrelationMain
       Implicit None
       Type( ProjDynFacVar ) ::data
       Real*8,dimension( : , : )  ::atoms
       Real*8,dimension( : , : )  ::lattice
       Integer ::qq
       Integer ::i
       Integer ::j
       Integer ::k
       Integer ::Indx
       Real*8  ::dpx
       Real*8  ::dpy
       Real*8  ::dpz

       Real*8  ::Dotp
       Real*8,dimension( 1:3 )  ::UnitQ


       Complex*16,allocatable,dimension( : , : )  ::DynProj


       If ( data%ActStep .gt. 1 ) then
         Call ComputeVelocity( data , atoms , lattice )
       End If

       Allocate( DynProj( 1:data%NQGrid , 1:6 ) )
       DynProj  = 0d0
       !$OMP Parallel Do Private(qq,i,j,k,indx,dpx,dpy,dpz,Dotp) Shared(atoms,data) &
       !$OMP& Reduction (+:DynProj)
       Do qq = 1 , data%Nqgrid
            !!looping over unit cells
            UnitQ  =  data%QgridDir( : , qq ) / get_norm( data%QgridDir( : , qq ) )
            Do i = 1 , Size( data%LinkUC , 1 )
               !!! looping over atoms in unit cell
               Do k = 1, Size( data%Velocity , 2 )
                  Indx  =  data%LinkUC( i )%List( k )
                  dpx   =  data%QgridDir( 1 , qq ) * atoms( 1 , Indx ) 
                  dpy   =  data%QgridDir( 2 , qq ) * atoms( 2 , Indx ) 
                  dpz   =  data%QgridDir( 3 , qq ) * atoms( 3 , Indx ) 
                  Dotp  =  dpx + dpy + dpz

                  !! Longitudinal part
                  DynProj( qq , 1:3 )  =  DynProj( qq , 1:3 ) + &
                            Dot_Product( data%Velocity( : , k , i ) , UnitQ ) * & 
                            UnitQ * CDexp( ImagI * Dotp ) / data%Masses( k )

                  !! transverse part
                  DynProj( qq , 4:6 )  =  DynProj( qq , 4:6 ) + &
                                       ( data%Velocity( : , k , i ) - &
                                  Dot_Product( data%Velocity( : , k , i ) , UnitQ ) * UnitQ ) * &
                                           UnitQ * CDexp( ImagI * Dotp ) / data%Masses( k )
               End Do
            End Do
       End Do
       !$OMP End Parallel Do


       If ( data%ActStep .gt. 1 ) then
         Do j = 1 , Size( DynProj , 2 )
            Call CrossCorrelationMain( data%TimeCorrelationProj( j ) , DynProj( : , j ) )
         End Do
       End If
       data%ActStep =  data%ActStep + 1

     End Subroutine ComputeProjectedDynamicStructureFactorMain


     Subroutine ComputeProjectedDynamicStructureFactorFinalize( data , fname )

       use tool_interfaces,only :give_closest_lower_power_of_2, get_norm,&
                                 gen_file_number
       use fourier_interface, only :fourier_dp
       use CrossCorrComplex,only:FinalizeCrossCorr
       use WindowingAnalysis, only:GaussianFiltering
       Implicit None
       Type( ProjDynFacVar ) ::data
       Character( len = * ) ::fname
       Character( len = 3 ) ::fnum
       Complex*16,allocatable,dimension( : , : ) ::FourierArray

       Integer ::N2
       Integer ::qq
       Integer ::i
       Integer ::j
       Real*8 ::deltaFrequ
       Real*8 ::Qpath


       !! compute average over time correltion functions
       If ( data%TimeCorrOnOff ) then
         Do j = 1, Size( data%TimeCorrelationProj , 1 )
            Call FinalizeCrossCorr( data%TimeCorrelationProj( j ) )
            If ( j .eq. 1 ) then
               Allocate( data%DynStrucFacProj( 1:Size( data%TimeCorrelationProj( j )%Average , 1 ),&
                                       1:Size( data%TimeCorrelationProj( j )%Average , 2 ),&
                                       1:6 ) )
            End If
            data%DynStrucFacProj( : , : , j )  =  data%TimeCorrelationProj( j )%Average( : , : )
         End Do
       End If

       N2 =  give_closest_lower_power_of_2( Size( data%DynStrucFacProj , 2 ) )
       !!! apply filter function 
       If ( data%FilterFunction ) then
               Do i = 1 , N2
                  Do j = 1, Size( data%DynStrucFacProj , 3 )
                     data%DynStrucFacProj( : , i , j ) = &
                          data%DynStrucFacProj( : , i , j ) * Dexp( -Dble( i ) * 6d0 / Dble( N2 ) )
                  End Do
               End Do
       End If



       !!! write time-dependent structure factor S(q,t)
       Do j = 1, Size( data%DynStrucFacProj , 3 )
         Call gen_file_number( j , fnum )
         open( unit = 302 , status='unknown' , file = &
                       "StructureFactorProj_vs_t.out"//Trim(fnum) , action='write' )
          Do i = 1 , N2
             Write( 302 , * ) Dble( i )*data%Tstep-data%Tstep/2d0,&
                              ( Abs( data%DynStrucFacProj( qq , i , j ) ) ,&
                              qq = 1, Size( data%DynStrucFacProj( : , i , j ) ) )
          End Do
         close( 302 )
       End Do

       !!! fourier transform dynamic structure factor
       Allocate( FourierArray( 1:N2 , data%Nqgrid ) )
       Do j = 1 , Size( data%DynStrucFacProj , 3 )
          FourierArray = Transpose( data%DynStrucFacProj( 1:data%Nqgrid , 1:N2 , j ) )
          Do qq = 1 , data%Nqgrid
            Call fourier_dp( FourierArray( : , qq ) , -1 )
          End Do
          deltafrequ = 1d0 / ( 2d0 * data%Tstep ) / ( Dble( N2 ) / 2d0 )

          !! Write fourier transform without gaussian smoothing
          Call gen_file_number( j , fnum )
          open( unit = 302 , status='unknown' , file = Trim( fname )//Trim( fnum ) , action='write' )
          Do i = 1 , N2 / 2 
             Write( 302 , * ) ( deltafrequ * Dble( i ) - deltafrequ / 2d0 ) ,&
                                Abs( FourierArray( i , : ) ) / data%Natoms / PI2
          End Do
          close( 302 )
       End Do


       Deallocate( FourierArray )


       open( unit = 302 , action='write' , status='unknown' , file = "QpathDSFProj.out" )
         Qpath = 0d0
         Write( 302  , * ) Qpath
         Do i = 2 , Size( data%QgridDir , 2 )
            Qpath = Qpath + get_norm( MatMul( data%InvLattice ,&
                                              data%QgridDir( : , i )/PI2-data%QgridDir( : , i - 1 )/PI2 ) )
            Write( 302 , * ) Qpath 
         End Do
       close( 302 )
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       Deallocate( data%DynStrucFacProj )

     End Subroutine ComputeProjectedDynamicStructureFactorFinalize

 End Module ComputeProjectedDynamicStructureFactor




!!!
!!!                   ~~
!!!                 ~~~~~~
!!!              ~~~~~~~~~~~~
!!!            ~~~~~~~~~~~~~~~~
!!!     Compute Single Phonon scattering
!!!            ~~~~~~~~~~~~~~~~
!!!              ~~~~~~~~~~~~
!!!                 ~~~~~~
!!!                   ~~
 

 Module QResolvedSinglePhononScattering

   use CrossCorrComplex, only :CrossCorrVar
   use SplitCelltoBox , only :BoxList
   Implicit None
   Type SinglePhononVar 
     Real*8,allocatable,dimension( : , : )  ::QgridDir                  !!! Qgrid in direct coordinates
     Real*8,allocatable,dimension( : , : )  ::QgridCart                 !!! Qgrid in cartesian coordinates
     Integer ::Nqgrid                                                   !!! size qgrid
     Integer ::ActStep                                                  !!! actual MD step
     Integer ::MaxStep                                                  !!! total number of MD steps
     Real*8  ::Tstep                                                    !!! time step
     Real*8,dimension( 1:3 , 1:3 )  ::lattice                           !!! lattice
     Real*8,dimensioN( 1:3 , 1:3 )  ::InvLattice                        !!! Inverse lattice
     Complex*16,allocatable,dimension( : ) ::Pinit                      !!! spec
     Complex*16,allocatable,dimension( : ) ::P
     Real*8,allocatable,dimension( : ) ::norm
     Integer ::QPathNr
     Integer,allocatable,dimension( : )  ::Qpath                        !!! indices containing qpath
     Type( CrossCorrVar ) ::TimeCorrelation
     Logical ::TimeCorrOnOff
     Integer ::Natoms

     !!! projected stuff
     Real*8,allocatable,dimension( : , : , : , : )  ::BasisVectors
     Complex*16,allocatable,dimension( : , : )  ::DynStrucFac
     Real*8,allocatable,dimension( : , : )  ::Equilibrium     !!! equilibrium positions
     Real*8,allocatable,dimensioN( : , : , : )  ::Displacement 
     Type( BoxList ),allocatable,dimension( : )  ::LinkUC     !!! Link atoms to unit cells
     Real*8,allocatable,dimension( : , : )  ::Boxes           !!! centers for fourier transform
     Real*8,allocatable,dimension( : ) ::Masses
     Logical :: FilterFunction                                !!! checks if exponential filter function is applied
   End Type SinglePhononVar

   Real*8,parameter ::PI2 = 4d0 * Dacos( 0d0 )                          !!! factor 2pi
   Complex*16,parameter ::ImagI = Cmplx( 0d0 , 1d0 )                    !!! imaginary number i
   Contains

     Subroutine QResolvedSinglePhononScatteringInit( data , Nsteps , dt , Nx ,&
                                 Ny , Nz , atoms , lattice , fname )

       use FlagReader, only :ReadFlags
       use CrossCorrComplex, only :CrossCorrelationInit
       use QgridCollection, only :SampleFirstBrillouinZoneCubic
       use PhononHelperRoutines, only :GetMasses,&
                                       ReadEigenvectorsModeC,&
                                       MakeDiagonalBasis,&
                                       ReadQpointsFromFile
       use ReadingInput, only :INPUT,OpenPOSCARFile,ReadingCoordinatesATO



       Implicit None
       Type( SinglePhononVar ) ::data
       Integer ::NSteps
       Real*8  ::dt
       Integer ::Nx
       Integer ::Ny
       Integer ::Nz
       Real*8,allocatable,dimension( : , : )  ::atoms
       Real*8,dimension( : , : ) ::lattice
       Character( len = * ) ::fname


       Logical ::QFound
       Real*8  ::temp
       Integer ::i

       Type( Input ) ::Equi
       Logical       ::EquiFound


       data%MaxStep =  NSteps
       data%Tstep   =  dt
       data%lattice = lattice
       data%InvLattice =  data%lattice
       Call lu_inversion( data%InvLattice , 3 )

       Qfound = .False.
       Call ReadFlags( "Qpath" , temp , QFound , fname )
       If ( .not. Qfound ) then
         data%QPathNr = 1
       Else
         data%QPathNr = NINT( temp )
       End If

       data%Natoms = Size( atoms , 2 )

       Call ReadBoxList( data , Nx , Ny , Nz , Size( atoms , 2 ) , "BoxList.in" )
       !Call GetMasses( data%Masses , "masses.in" )

       !Inquire( File = "BasisVector.in" , exist = QFound )
       !If ( QFound ) then
       !   Write( * , * ) "ComputeProjectedDynamicStructureFactorInit"
       !   Write( * , * ) "is reading BasisVectors and Qpoints"
       !   Write( * , * ) "from file BasisVector.in"
       !   Call ReadEigenvectorsModeC( data%BasisVectors , data%QGridDir , "BasisVector.in" )
       !   data%NQgrid = Size( data%QGridDir , 2 )
       !Else
       !   Write( * , * ) "ComputeProjectedDynamicStructureFactorInit"
       !   Write( * , * ) "No BasisVector.in file found"
       !   Write( * , * ) "Checking for a QVectorsFile"
       Inquire( File="QVectors.in" , exist = Qfound )
       If ( Qfound ) then
          Write( * , * ) "In routine QResolvedSinglePhononScatteringInit"
          Write( * , * ) "Q-vectors are read"
          Write( * , * ) "QVectors.in"
          Call ReadQpointsFromFile( data%QGridDir , "QVectors.in" )
          data%NQgrid  =  Size( data%QgridDir , 2 )
       Else
          Write( * , * ) "In routine QResolvedSinglePhononScatteringInit"
          Write( * , * ) "no QVectors.in"
          Write( * , * ) "supplied. Using default Q-mesh of first irreducible "
          Write( * , * ) "Brillouin zone"
          Call SampleFirstBrillouinZoneCubic( data%QgridDir , Nx , Ny , Nz )
          data%NQgrid  =  Size( data%QgridDir , 2 )
       End If
       !   Write( * , * ) "Using a diagonal basis since no file was supplied"
       !   Write( * , * ) "Diagonal basis gives some decomposition but does"
       !   Write( * , * ) "make sense physically"
       !   !Call MakeDiagonalBasis( 3 , data%LinkUC( 1 )%Natoms , data%Nqgrid , data%BasisVectors )
       !End If
       !!
       !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       data%QgridDir = PI2 * data%QgridDir
       !! projected dynamic structure factor

       Equi%fname  =  "EquiPos.in"
       Inquire( File= Trim( Equi%fname ) , exist = EquiFound )
       If ( EquiFound ) then
            !! get equlibrium structure
            Equi%FileNr = 111
            Call OpenPOSCARFile( Equi )
            Allocate( data%Equilibrium( 1:3 , 1:Equi%Natoms ) )
            Call ReadingCoordinatesATO( Equi , data%Equilibrium )
            close( Equi%FileNr )
       Else
            Write( * , * ) "#########  WARNING FILE EquiPos.in not found"
            Write( * , * ) "Taking first configuration of XDATCAR instead as "
            Write( * , * ) "equilibrium position"
            data%Equilibrium = atoms
       End If
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       data%TimeCorrOnOff = .False.
       !! supply window size
       Call ReadFlags( "DSTC" , temp , data%TimeCorrOnOff , fname )
       If ( .not. data%TimeCorrOnOff ) then
          temp  =  Real( NSteps - 1 )
          data%TimeCorrOnOff = .True.
       End If
       If ( data%TimeCorrOnOff ) then
          Write( * , * ) "Time Averaging of the projected dynamic structure factor"
          Write( * , * ) "is switched on. The window size was set to"
          Write( * , * ) NINT( temp )
          !! Projected stuff
          Call CrossCorrelationInit( data%TimeCorrelation , &
                 NINT( temp ) , NSteps-1 , Size( data%QgridDir , 2 ) )
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       End If


     End Subroutine QResolvedSinglePhononScatteringInit



     Subroutine ReadBoxList( data , N1 , N2 , N3 , Natoms , fname )

       use SplitCelltoBox, only :MakeBoxes,MakeBoxList
       Implicit None
       Type( SinglePhononVar ) ::data
       Integer  ::N1
       Integer  ::N2
       Integer  ::N3
       Integer  ::Natoms
       Character( len = * )  ::fname
       Integer  ::Npc
       Real*8,allocatable,dimension( : , : )  ::Boxes

       Integer ::i


       Npc = Natoms / ( N1 * N2 * N3 )
       
       Boxes = MakeBoxes( N1 , N2 , N3 )


       Allocate( data%LinkUC( 1:N1*N2*N3 ) )
       open( unit = 33 , action='read' , file = Trim( fname ) , status='old' )
       Do i = 1 , Size( data%LinkUC , 1 )
          Allocate( data%LinkUC( i )%List( 1:Npc ) )
          data%LinkUC( i )%Natoms = Npc
          Read( 33 , * ) data%LinkUC( i )%List( : )
       End Do
       close( 33 )


       
       Allocate( data%Displacement( 1:3 , 1:NPc , 1:N1*N2*N3 ) )
       Allocate( data%Boxes( 1:3 , 1 : N1 * N2 * N3 ) )
       Do i = 1 , Size( data%Boxes , 2 )
          data%Boxes( 1 , i )  =  0.5d0 * ( Boxes( 1 , i ) + Boxes( 2 , i ) )
          data%Boxes( 2 , i )  =  0.5d0 * ( Boxes( 3 , i ) + Boxes( 4 , i ) )
          data%Boxes( 3 , i )  =  0.5d0 * ( Boxes( 5 , i ) + Boxes( 6 , i ) )
       End Do


     End Subroutine ReadBoxList
       


     Subroutine ComputeDisplacement( data , PosNew , lattice )
             use tool_interfaces, only :get_nearest_image,get_norm,&
                             get_norm
             Implicit None
             Type( SinglePhononVar ) ::data
             Real*8,dimension( : , : )     ::PosNew
             Real*8,dimension( : , : ) ::lattice

             Integer ::i
             Integer ::j
             Integer ::Indx

             Do i = 1 , Size( data%LinkUC , 1 )
               Do j = 1 , data%LinkUC( i )%Natoms
                  Indx = data%LinkUC( i )%List( j )
                  data%Displacement( : , j , i )  =  MatMul( lattice , &
                                                 get_nearest_image( data%Equilibrium( : , Indx ) ,&
                                                 PosNew( : , Indx ) ) - data%Equilibrium( : , Indx ) )
               End Do
             End Do

     End Subroutine ComputeDisplacement



     !!
     !! comuting the spatial fourier transform for obtaining the
     !! dynamic structure factor, nearest-neighbor convention used
     !! S( q , t ) Sum_{R,R'} exp(iq(R(t)-R'(t)))
     !! the modes are projected on eigenvectors that have to be supplied as 
     !! input
     !!
     Subroutine QResolvedSinglePhononScatteringMain( data , atoms , lattice )

       use tool_interfaces, only:get_nearest_image,get_norm
       use CrossCorrComplex, only :CrossCorrelationMain
       Implicit None
       Type( SinglePhononVar ) ::data
       Real*8,dimension( : , : )  ::atoms
       Real*8,dimension( : , : )  ::lattice
       Integer ::qq
       Integer ::i
       Integer ::j
       Integer ::k
       Integer ::Indx
       Real*8 ::dpx
       Real*8 ::dpy
       Real*8 ::dpz

       Real*8 ::Dotp
       Real*8,dimension( 1:3 )  ::UnitQ
       Complex*16,dimension( 1:data%NQGrid )  ::DynProj


       If ( data%ActStep .gt. 1 ) then
         Call ComputeDisplacement( data , atoms , lattice )
       End If

       DynProj  = 0d0
       !$OMP Parallel Do Private(qq,i,j,k,indx,dpx,dpy,dpz,Dotp) Shared(atoms,data) &
       !$OMP& Reduction (+:DynProj)
       Do qq = 1 , data%Nqgrid
            !!looping over unit cells
            UnitQ  =  data%QgridDir( : , qq ) / get_norm( data%QgridDir( : , qq ) )
            UnitQ  =  data%QgridDir( : , qq )
            Do i = 1 , Size( data%LinkUC , 1 )
               !!! looping over atoms in unit cell
               Do k = 1, Size( data%Displacement , 2 )
                  Indx  =  data%LinkUC( i )%List( k )
                  !dpx   =  data%QgridDir( 1 , qq ) * atoms( 1 , Indx ) 
                  !dpy   =  data%QgridDir( 2 , qq ) * atoms( 2 , Indx ) 
                  !dpz   =  data%QgridDir( 3 , qq ) * atoms( 3 , Indx ) 
                  dpx   =  data%QgridDir( 1 , qq ) * data%Boxes( 1 , i ) 
                  dpy   =  data%QgridDir( 2 , qq ) * data%Boxes( 2 , i )
                  dpz   =  data%QgridDir( 3 , qq ) * data%Boxes( 3 , i )
                  Dotp  =  dpx + dpy + dpz

                  DynProj( qq )  =  DynProj( qq ) + &
                            Dot_Product( data%Displacement( : , k , i ) , UnitQ ) * & 
                                         CDexp( ImagI * Dotp )

               End Do
            End Do
       End Do
       !$OMP End Parallel Do


       If ( data%ActStep .gt. 1 ) then
           Call CrossCorrelationMain( data%TimeCorrelation , DynProj( : ) )
       End If
       data%ActStep =  data%ActStep + 1

     End Subroutine QResolvedSinglePhononScatteringMain


     Subroutine QResolvedSinglePhononScatteringFinalize( data , fname )

       use tool_interfaces,only :give_closest_lower_power_of_2, get_norm,&
                                 gen_file_number
       use fourier_interface, only :fourier_dp
       use CrossCorrComplex,only:FinalizeCrossCorr
       use WindowingAnalysis, only:GaussianFiltering
       Implicit None
       Type( SinglePhononVar ) ::data
       Character( len = * ) ::fname
       Character( len = 3 ) ::fnum
       Complex*16,allocatable,dimension( : , : ) ::FourierArray

       Integer ::N2
       Integer ::qq
       Integer ::i
       Integer ::j
       Real*8 ::deltaFrequ
       Real*8 ::Qpath


       !! compute average over time correltion functions
       If ( data%TimeCorrOnOff ) then
            Call FinalizeCrossCorr( data%TimeCorrelation )
            Allocate( data%DynStrucFac( 1:Size( data%TimeCorrelation%Average , 1 ),&
                                        1:Size( data%TimeCorrelation%Average , 2 ) ) )
            data%DynStrucFac( : , : )  =  data%TimeCorrelation%Average( : , : )
       End If

       N2 =  give_closest_lower_power_of_2( Size( data%DynStrucFac , 2 ) )
       !!! apply filter function 
       If ( data%FilterFunction ) then
               Do i = 1 , N2
                  data%DynStrucFac( : , i ) = &
                       data%DynStrucFac( : , i ) * Dexp( -Dble( i ) * 6d0 / Dble( N2 ) )
               End Do
       End If



       !!! write time-dependent structure factor S(q,t)
       open( unit = 302 , status='unknown' , file = &
                     "SinglePhononScattering_vs_t.out" , action='write' )
        Do i = 1 , N2
           Write( 302 , * ) Dble( i )*data%Tstep-data%Tstep/2d0,&
                            ( Abs( data%DynStrucFac( qq , i ) ) ,&
                            qq = 1, Size( data%DynStrucFac( : , i ) ) )
        End Do
       close( 302 )

       !!! fourier transform dynamic structure factor
       Allocate( FourierArray( 1:N2 , data%Nqgrid ) )
       FourierArray = Transpose( data%DynStrucFac( 1:data%Nqgrid , 1:N2 ) )
       Do qq = 1 , data%Nqgrid
         Call fourier_dp( FourierArray( : , qq ) , -1 )
       End Do
       deltafrequ = 1d0 / ( 2d0 * data%Tstep ) / ( Dble( N2 ) / 2d0 )

       !! Write fourier transform without gaussian smoothing
       open( unit = 302 , status='unknown' , file = Trim( fname ) , action='write' )
       Do i = 1 , N2 / 2 
          Write( 302 , * ) ( deltafrequ * Dble( i ) - deltafrequ / 2d0 ) ,&
                             Abs( FourierArray( i , : ) ) / data%Natoms / PI2
       End Do
       close( 302 )


       Deallocate( FourierArray )


       open( unit = 302 , action='write' , status='unknown' , file = "QpathSPS.out" )
         Qpath = 0d0
         Write( 302  , * ) Qpath
         Do i = 2 , Size( data%QgridDir , 2 )
            Qpath = Qpath + get_norm( MatMul( data%InvLattice ,&
                                              data%QgridDir( : , i )/PI2-data%QgridDir( : , i - 1 )/PI2 ) )
            Write( 302 , * ) Qpath 
         End Do
       close( 302 )
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       Deallocate( data%DynStrucFac )

     End Subroutine QResolvedSinglePhononScatteringFinalize

 End Module QResolvedSinglePhononScattering


 Module PhononInterface

   use ComputePVACF, only:PVACFVar,&
                     ComputePVACFInit,&
                     ComputePVACFMain,&
                     ComputePVACFFinalize

   use ComputePVACFKProjOnly,only:PVACFVarKProj,&
                     ComputePVACFInitKProjOnly,&
                     ComputePVACFMainKProjOnly,&
                     ComputePVACFFinalizeKProjOnly

   use ComputeDynamicStructureFactor, only:DynStrucVar,&
                     ComputeDynamicStructureFactorInit,&
                     ComputeDynamicStructureFactorMain,&
                     ComputeDynamicStructureFactorFinalize
   
   use ComputeProjectedDynamicStructureFactor, only:ProjDynFacVar,&
                         ComputeProjectedDynamicStructureFactorInit,&
                         ComputeProjectedDynamicStructureFactorMain,&
                         ComputeProjectedDynamicStructureFactorFinalize
   
   use QResolvedSinglePhononScattering, only:SinglePhononVar,&
                         QResolvedSinglePhononScatteringInit,&
                         QResolvedSinglePhononScatteringMain,&
                         QResolvedSinglePhononScatteringFinalize
   
   Implicit None
   Logical :: PVACFOnOff
   Logical :: PVACFKProjOnOff
   Logical :: DynStrucFacOnOff
   Logical :: ProjectedDynStrucFacOnOff
   Logical :: SinglePhononOnOff
   Type ( PVACFVar ) ::PVACFData
   Type ( PVACFVarKProj ) ::PVACFDataKProj
   Type ( DynStrucVar ) ::DynStrucData
   Type ( ProjDynFacVar ) ::ProjectedDynStrucData
   Type ( SinglePhononVar ) ::SinglePhononData
   Real*8,allocatable,dimension( : , : )  ::Positions   !! atomic positions in direct coordinates
   Integer ::Natoms                  !! total number of atoms

   Contains

   !!! setting up the lattice oscillations
   !!! or phonons computation routines
   !!! possibilities are
   !!! projected velocity autocorrelation function 
   !!! eigenvectors and Kspace -> PVACF
   !!! projected velocity autocorrelation function -> PVACF
   !!! eigenvectors and k-space
   !!! projected velocity autocorrelation function
   !!! K space -> PVACK
   !!! projected velocity autocorrelation function -> PVACF
   Subroutine PhononInterfaceInit( MDParams , Atoms , lattice )

     use FlagReader, only :ReadFlags
     use tool_interfaces, only :GetSelfCorrSteps
     use ReadingInput, only:atom_type
     use ReadInputFile, only:SimParams
     Implicit None
     Type( SimParams )     ::MDParams
     Type( atom_type ),dimension( : ) ::atoms
     Real*8,dimension( 1:3 , 1:3 )  ::lattice

     Call ReadFlags( "PVACF" , PVACFOnOff , MDParams%InputFileName )
     Call ReadFlags( "PVACK" , PVACFKProjOnOff , MDParams%InputFileName )
     Call ReadFlags( "DYNSTRUC" , DynStrucFacOnOff , MDParams%InputFileName )
     Call ReadFlags( "PROJFAC" , ProjectedDynStrucFacOnOff , MDParams%InputFileName )
     Call ReadFlags( "SPHON" , SinglePhononOnOff , MDParams%InputFileName )

     If ( PVACFOnOff .or. DynStrucFacOnOff .or. ProjectedDynStrucFacOnOff &
                     .or. PVACFKProjOnOff .or. SinglePhononOnOff ) then
       Natoms = Sum( atoms( : )%Natoms )
       Allocate( Positions( 1:3 , 1:Natoms ) )
       Call PrepareAtoms( atoms )
     End If

     If ( PVACFOnOff ) then
       PVACFData%FilterFunction  =  MDParams%FilterFunction
       Call ComputePVACFInit( PVACFData ,&
                        GetSelfCorrSteps( MDParams%Nstart , MDParams%Nend,&
                                          MDParams%BlockSample ) , MDParams%timestep , &
                               MDParams%Nx , MDParams%Ny , MDParams%Nz ,&
                               Positions , lattice , MDParams%InputFileName )
     End If

     If ( PVACFKProjOnOff ) then
       PVACFDataKProj%FilterFunction  =  MDParams%FilterFunction
       Call ComputePVACFInitKProjOnly( PVACFDataKProj ,&
                        GetSelfCorrSteps( MDParams%Nstart , MDParams%Nend, &
                                          MDParams%BlockSample ) , MDParams%timestep , &
                               MDParams%Nx , MDParams%Ny , MDParams%Nz ,&
                               Positions , lattice , MDParams%InputFileName )
     End If
     If ( DynStrucFacOnOff ) then
       DynStrucData%FilterFunction  =  MDParams%FilterFunction
       Call ComputeDynamicStructureFactorInit( DynStrucData ,&
                        GetSelfCorrSteps( MDParams%Nstart , MDParams%Nend,&
                                          MDParams%BlockSample ) , MDParams%timestep , &
                               MDParams%Nx , MDParams%Ny , MDParams%Nz ,&
                               Positions , lattice , MDParams%InputFileName )
     End If
     If ( ProjectedDynStrucFacOnOff ) then
       ProjectedDynStrucData%FilterFunction  =  MDParams%FilterFunction
       Call ComputeProjectedDynamicStructureFactorInit( ProjectedDynStrucData ,&
                        GetSelfCorrSteps( MDParams%Nstart , MDParams%Nend, &
                                          MDParams%BlockSample ) , MDParams%timestep , &
                               MDParams%Nx , MDParams%Ny , MDParams%Nz ,&
                               Positions , lattice , MDParams%InputFileName )
     End If
     

     If ( SinglePhononOnOff ) then
       SinglePhononData%FilterFunction  =  MDParams%FilterFunction
       Call QResolvedSinglePhononScatteringInit( SinglePhononData ,&
                        GetSelfCorrSteps( MDParams%Nstart , MDParams%Nend, &
                                          MDParams%BlockSample ) , MDParams%timestep , &
                               MDParams%Nx , MDParams%Ny , MDParams%Nz ,&
                               Positions , lattice , MDParams%InputFileName )
     End If

   End Subroutine PhononInterfaceInit








   Subroutine PhononInterfaceMain( MDParams , Atoms , lattice )

     use ReadingInput, only:atom_type
     use ReadInputFile, only:SimParams
     Implicit None
     Type( SimParams )     ::MDParams
     Type( atom_type ),dimension( : ) ::atoms
     Real*8,dimension( 1:3 , 1:3 )  ::lattice

     If ( PVACFOnOff .or. DynStrucFacOnOff .or. ProjectedDynStrucFacOnOff &
            .or. PVACFKProjOnOff .or. SinglePhononOnOff ) then
       Call PrepareAtoms( atoms )
     End If

     If ( PVACFOnOff ) then
       Call ComputePVACFMain( PVACFData , Positions , lattice )
     End If

     If ( PVACFKProjOnOff ) then
       Call ComputePVACFMainKProjOnly( PVACFDataKProj , Positions , lattice )
     End If

     If ( DynStrucFacOnOff ) then
       Call ComputeDynamicStructureFactorMain( DynStrucData , Positions , lattice )
     End If

     If ( ProjectedDynStrucFacOnOff ) then
       Call ComputeProjectedDynamicStructureFactorMain( ProjectedDynStrucData , Positions , lattice )
     End If
     
     If ( SinglePhononOnOff ) then
       Call QResolvedSinglePhononScatteringMain( SinglePhononData , Positions , lattice )
     End If

   End Subroutine PhononInterfaceMain


   Subroutine PrepareAtoms( atoms )

     use ReadingInput, only:atom_type
     Implicit None
     Integer ::i
     Integer ::j
     Type( atom_type ),dimension( : ) ::atoms
     Integer ::col

     col = 1
     Do i = 1, atoms( 1 )%Nspecies
       Do j = 1 , atoms( i )%Natoms
         Positions( : , col ) = atoms( i )%direct_coords( : , j )
         col = col + 1
       End Do
     End Do
   End Subroutine PrepareAtoms

   Subroutine PhononInterfaceFinal

     Implicit None
  
     If ( PVACFOnOff ) then
       Call ComputePVACFFinalize( PVACFData , "PVACF.out" )
     End If
     
     If ( PVACFKProjOnOff ) then
       Call ComputePVACFFinalizeKProjOnly( PVACFDataKProj , "PVACF_KO.out" )
     End If
     
     If ( DynStrucFacOnOff ) then
       Call ComputeDynamicStructureFactorFinalize( DynStrucData , "DynamicStrucFac.out" )
     End If

     If ( ProjectedDynStrucFacOnOff ) then
       Call ComputeProjectedDynamicStructureFactorFinalize( ProjectedDynStrucData , "ProjectedDSF.out" )
     End If
     
     If ( SinglePhononOnOff ) then
       Call QResolvedSinglePhononScatteringFinalize( SinglePhononData , "SinglePhononScattering.out" )
     End If

   End Subroutine PhononInterfaceFinal

 End Module PhononInterface
