

      Module ReadInputFile


              !!! this data structure contains your MD specifications
              Type SimParams
                   Integer ::Nstart         !! starting geometry
                   Integer ::Nend           !! geometry to stop sampling
                   Real*8  ::timeStep       !! time step of your MD
                   Integer ::BlockSample    !! Sample every BlockSample th structure
                   Integer ::Nx             !! number of unit cells in x
                   Integer ::Ny             !! number of unit cells in y
                   Integer ::Nz             !! number of unit cells in z
                   Logical ::FilterFunction !! apply filter function or not
                   Character( len=: ),allocatable ::InputFileName
              End Type
              Contains

                      !!! I am the routine reading your input file
                      !!! if you don't like me feel free to adapt me
                      !!! to your desires ( only by adppating the flag
                      !!! names ;-) )
                      !!! don't forget to recompile afterwards
                Subroutine ReadInput( data , fname )

                  use FlagReader, only: ReadFlags
                  
                  Implicit None
                  Type( SimParams )  ::data      !!! your MD
                                                 !!! specification 
                                                 !!! structure it knows
                                                 !!! everyting in the end ;-)
                  Character( len = * )  ::fname  !!! input file name where flags are specified

                  Real*8 :: temp    !!! temp variable to read values to

                  Logical ::found   !!! initial value is set in flag routines 

                  data%InputFileName  =  fname
                  !!! check for supplied start structure
                  Call ReadFlags( "NSTART" , temp , found , fname )
                  If ( found ) then
                     data%Nstart = NINT( temp )
                     Write( * , "(A)" , advance = 'No' ) 'You decided '
                     Write( * , "(A)" ) 'to start sampling at structure'
                     Write( * ,  "(I5)" ) data%Nstart
                  Else
                     data%Nstart = 1
                     Write( * , "(A)" ) 'Start sampling at'
                     Write( * ,  "(I5)" ) data%Nstart
                  End If

                  !!! check for supplied end structure
                  Call ReadFlags( "NEND" , temp , found , fname )
                  If ( found ) then
                     data%Nend = NINT( temp )
                     Write( * , "(A)" , advance = 'No' ) 'You decided '
                     Write( * , "(A)" ) 'to end sampling at structure'
                     Write( * ,  "(I5)" ) data%Nend
                  Else
                     Write( * , "(A)" , advance = "No" ) 'I am sorry'
                     Write( * , "(A)" , advance = "No" ) ' but I have'
                     Write( * , "(A)" , advance = "No" ) ' to know '
                     Write( * , "(A)" , advance = "No" ) ' how many '
                     Write( * , "(A)" ) 'strcutures there are?'
                     Write( * , "(A)" ) 'Stopping code'
                     STOP
                  End If

                  !!! choose down smapling step
                  Call ReadFlags( "TSTEP" , temp , found , fname )
                  If ( found ) then
                     data%timestep = temp
                     Write( * , "(A)" , advance = 'No' ) 'You decided '
                     Write( * , "(A)" ) 'to have a time step of'
                     Write( * ,  "(F5.2)" ) data%timestep
                  Else
                     data%timestep = 1
                     Write( * , "(A)" ) 'Time step is'
                     Write( * ,  "(F6.3)" ) data%timestep
                  End If

                  
                  Call ReadFlags( "DSTEP" , temp , found , fname )
                  If ( found ) then
                     data%BlockSample = NINT( temp )
                     Write( * , "(A)" , advance = 'No' ) 'You decided '
                     Write( * , "(A)" ) 'down sample with every'
                     Write( * ,  "(I5)" ) data%BlockSample
                  Else
                     data%BlockSample = 1
                     Write( * , "(A)" ) 'No down sampling used'
                     Write( * ,  "(I5)" ) data%BlockSample
                  End If

                  !!! correct for the time step when downsampling is
                  !!! used
                  data%timestep  =  data%timestep * data%BlockSample
                  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                  !!! checking cell dimensions in x direction
                  Call ReadFlags( "NX" , temp , found , fname )
                  If ( found ) then
                     data%Nx = NINT( temp )
                     Write( * , "(A)" , advance = 'No' ) 'Number of'
                     Write( * , "(A)" , advance = 'No' ) ' cells in '
                     Write( * , "(A)" ) 'x'
                     Write( * ,  "(I5)" ) data%Nx
                  Else
                     Write( * , "(A)" ) 'I really need cell dimensions'
                     Write( * , "(A)" ) 'I have to stop sorry'
                     Write( * , "(A)" ) 'Supply them with '
                     Write( * , "(A)" ) "NX=XXX"
                     STOP
                  End If

                  !!! checking cell dimensions in y direction
                  Call ReadFlags( "NY" , temp , found , fname )
                  If ( found ) then
                     data%Ny = NINT( temp )
                     Write( * , "(A)" , advance = 'No' ) 'Number of'
                     Write( * , "(A)" , advance = 'No' ) ' cells in '
                     Write( * , "(A)" ) 'y'
                     Write( * ,  "(I5)" ) data%Ny
                  Else
                     Write( * , "(A)" ) 'I really need cell dimensions'
                     Write( * , "(A)" ) 'I have to stop sorry'
                     Write( * , "(A)" ) 'Supply them with '
                     Write( * , "(A)" ) "NY=YYY"
                     STOP
                  End If
                  

                  !!! checking cell dimensions in z direction
                  Call ReadFlags( "NZ" , temp , found , fname )
                  If ( found ) then
                     data%Nz = NINT( temp )
                     Write( * , "(A)" , advance = 'No' ) 'Number of'
                     Write( * , "(A)" , advance = 'No' ) ' cells in '
                     Write( * , "(A)" ) 'z'
                     Write( * ,  "(I5)" ) data%Nz
                  Else
                     Write( * , "(A)" ) 'I really need cell dimensions'
                     Write( * , "(A)" ) 'I have to stop sorry'
                     Write( * , "(A)" ) 'Supply them with '
                     Write( * , "(A)" ) "NZ=ZZZ"
                     STOP
                  End If 
                  


                  !!! checking filter function
                  !!! filter function is an exp(-K*x)
                  data%FilterFunction  =  .False.
                  Call ReadFlags( "FILTER" , data%FilterFunction&
                                           , fname )
                  If ( found ) then
                     Write( * , "(A)" , advance = 'No' ) 'Apply Exp'
                     Write( * , "(A)" , advance = 'No' ) ' filter'
                     Write( * , "(A)" ) ' function'
                  Else
                     Write( * , "(A)" ) 'Filter Function switched off'
                  End If 

                  Write( * , "(A)" ) "Reading input parameters done" 
                  Write( * , "(A)" ) "-----------------------------" 

                End Subroutine ReadInput
      End Module ReadInputFile

