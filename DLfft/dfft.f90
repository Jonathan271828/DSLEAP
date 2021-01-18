
  Module fourier_interface

      Implicit None
      Contains

        Subroutine fourier_sp( data , isign )
      
          !data is a input array and will be replaced
          !by its fourier or N times fourier back transform
          !depending on the isign flag which is -1 or 1
          !respectively.
          !N the length of the array has to be integer power of 2
      
          Implicit None
          Complex*8,dimension( : )  ::data
          Integer                   ::isign
      
          Integer   ::N       !Size of input array
          Integer   ::i       !Loop counter
          Integer   ::istep   !counts steps during Davidson Lanczos
          Integer   ::m       !Loop counter
          Integer   ::j       !index
          Integer   ::mmax
          Integer   ::N2      !Half size of the data array
          Real*8    ::theta   !factor in exponential function
          Complex*8 ::temp    !variable to store temporary stuf
          Complex*16 ::w      !trigonometric reccurence
          Complex*16 ::wp
          Complex*8  ::ws
          Real*8,parameter ::pi = 2d0 * Dacos( 0d0 )
      
      
          N = Size( data , 1 )
      
          If ( Iand( N , N - 1 ) .ne. 0 ) then
            Write( * , "(A)" ) "N is not a power of 2; the program has to be stopped"
            STOP
          End If
      
          N2  =  N / 2
          j = N2
          !This section is usually caled the bit reversal part
          !This section sorts the first half of the array with odd
          !array indices
          !the second part of the array conatains only even indexed
          !quantities
          Do i = 1 , N - 2
            If ( j .gt. i ) then
               temp  =  data( i + 1 )
               data( i + 1 )  =  data( j + 1 )
               data( j + 1 )  =  temp
            End If
            m = N2
            Do
              If ( m .lt. 2 .or. j .lt. m ) Exit
              j  =  j - m 
              m  =  m / 2
            End Do
            j = j + m
          End Do
      
          mmax = 1
      
          !Danielson-Lanczos section
          !-------------------------
          Do    !Loop will be executed log_2(N) times
            If ( N .le. mmax ) Exit
            istep =  2 * mmax
            theta =  pi / ( isign * mmax )    !Initialization for trigonometric reccurence
            !wp   =          cos( theta ) - 1 , Sin( theta )
            wp    =  Cmplx( -2.0 * Sin( 0.5 * theta )**2 , Sin( theta ) )
            w     =  Cmplx( 1.0 , 0.0 )
            Do m = 1 , mmax
               ws = w
               Do i = m , N , istep
                 j = i + mmax
                 temp       =  ws * data( j )
                 data( j )  =  data( i )  -  temp
                 data( i )  =  data( i )  +  temp
               End Do
               w = w * wp + w
            End Do
            mmax = istep
          End Do
      
      
        End Subroutine fourier_sp
      
      
        
        Subroutine fourier_dp( data , isign )
      
          !data is a input array and will be replaced
          !by its fourier or N times fourier back transform
          !depending on the isign flag which is -1 or 1
          !respectively.
          !N the length of the array has to be integer power of 2
      
          Implicit None
          Complex*16,dimension( : ),intent( inout ) ::data
          Integer,intent( in )                      ::isign
      
          Integer   ::N       !Size of input array
          Integer   ::i       !Loop counter
          Integer   ::istep   !counts steps during Davidson Lanczos
          Integer   ::m       !Loop counter
          Integer   ::j       !index
          Integer   ::mmax
          Integer   ::N2      !Half size of the data array
          Real*8    ::theta   !factor in exponential function
          Complex*16 ::temp    !variable to store temporary stuf
          Complex*16 ::w      !trigonometric reccurence
          Complex*16 ::wp
          Complex*16  ::ws
          Real*8,parameter ::pi = 2d0 * Dacos( 0d0 )
      
      
          N = Size( data , 1 )

          If ( Iand( N , N - 1 ) .ne. 0 ) then
            Write( * , "(A)" ) "N is not a power of 2; the program has to be stopped"
            STOP
          End If
      
          N2  =  N / 2
          j = N2
          !This section is usually caled the bit reversal part
          !This section sorts the first have of the array with odd
          !array indices
          !the second part of the array conatains only even indexed
          !quantities
          Do i = 1 , N - 2
            If ( j .gt. i ) then
               temp  =  data( i + 1 )
               data( i + 1 )  =  data( j + 1 )
               data( j + 1 )  =  temp
            End If
            m = N2
            Do
              If ( m .lt. 2 .or. j .lt. m ) Exit
              j  =  j - m 
              m  =  m / 2
            End Do
            j = j + m
          End Do
      
          mmax = 1
      
          !Danielson-Lanczos section
          !-------------------------
          Do    !Loop will be executed log_2(N) times
            If ( N .le. mmax ) Exit
            istep =  2 * mmax
            theta =  pi / ( isign * mmax )    !Initialization for trigonometric reccurence
            !wp   =          cos( theta ) - 1 , Sin( theta )
            wp    =  Cmplx( -2d0 * DSin( 0.5d0 * theta )**2 , DSin( theta ) , kind = 16 )
            w     =  Cmplx( 1d0 , 0d0 , kind = 16 )
            Do m = 1 , mmax
               ws = w
               Do i = m , N , istep
                 j = i + mmax
                 temp       =  ws * data( j )
                 data( j )  =  data( i )  -  temp
                 data( i )  =  data( i )  +  temp
               End Do
               w = w * wp + w
            End Do
            mmax = istep
          End Do
        
        End Subroutine fourier_dp




        !!!
        !!!
        !!! routine does a discrete fourier transform
        !!! in the brute force N^2 scaling
        !!! way; on input supply the data array containing your signal
        !!! and the sign of the exponential
        !!! should not be used
        Subroutine DiscreteFourierTransform( data , isign )

          Implicit None
          Complex*16,dimension( : ),intent( inout )  ::data
          Integer,intent( in )  ::isign
          Real*8,parameter  ::pi2 = 4d0 *Dacos( 0d0 )
          Integer ::n
          Integer ::k
          Complex*16 ::W
          Complex*16 ::Wnk
          Integer ::Nsig
          Complex*16,dimension( 1:Size( data ) ) ::temp

          Nsig = Size( data )


          W  =  CDexp( DCmplx( 0d0 , pi2 / Nsig ) )
          temp = 0d0
          !! could we written as a matrix multiplication
          Do n = 0 , Nsig - 1
            Do k = 0 , Nsig - 1
              Wnk = W**( n*k*isign )
              temp( n+1 )  =  temp( n+1 ) + data( k+1 )*Wnk
            End Do
          End Do
          data = temp
        End Subroutine DiscreteFourierTransform



        Function determine_frequency_grid( Nsteps , dstep ) Result( frequency )


          Implicit None
          Integer,intent( in ) ::Nsteps
          Real*8,intent( in )  ::dstep

          Real*8  ::deltaFrequ
          Real*8,allocatable,dimension( : ) ::frequency
          Integer ::i
          Integer ::ArraySize

          ArraySize = Nsteps / 2


          deltaFrequ = 1d0 / ( 2d0 * dstep ) / ( Dble( Nsteps ) / 2d0 )
          Allocate( frequency( 1:ArraySize ) )

          Do i = 1 , ArraySize
             frequency( i )  =  ( deltafrequ * Dble( i ) - deltafrequ / 2d0 )
          End Do
        End Function determine_frequency_grid

  End Module fourier_interface
