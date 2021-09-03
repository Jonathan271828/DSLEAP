




      !!!! genereal input reader for POSCAR format
      Module ReadingInput

              Implicit None
              Type INPUT
                   Character ( len = : ),allocatable ::fname   !file name POSCAR
                   Character ( len = : ),allocatable ::Comment !comment line
                   Character( len=2 ),allocatable    ::AtomTypes(:) !!atom types array containing names
                   Integer,allocatable    ::AtomNumbers(:) !! containing atom numbers
                   Integer ::NTypes   !! number of atom types
                   Integer ::Natoms   !! total number of atoms
                   Real*8    :: latticeparameter  !!lattice parameter second line
                   Real*8,dimension( 1:3 , 1:3 )     ::lattice !!lattice constants
                   Integer ::FileNr   !! input file number
              End Type

              !!! use this to store atoms for ever type indicidually
              Type atom_type
                Real*8,allocatable,dimension( : , : )  ::direct_coords
                Real*8,allocatable,dimension( : , : )  ::cartesian_coords
                Integer     ::Natoms
                Integer  ::Nspecies
                Character( len = 2 ) ::Atype
              End Type atom_type



              Interface ReadCoordinates
                      module procedure ReadingCoordinatesAtomType,&
                                       ReadingCoordinatesATO
              End Interface





              Contains
                 !!!
                 !!!
                 !!! routine will open your poscar file
                 !!! then reads the comment line
                 !!! the lattice parameter
                 !!! the lattice
                 !!! AtomTypes
                 !!! AtomNumbers
                 !!!
                 Subroutine OpenPOSCARFile( data )


                         Type( INPUT ) ::data
                         Character( len = 100 ) ::buffer
                         Integer ::i


                       open( data%FileNr , status='old' ,file =&
                               data%fname , action = 'read')

                       Read( data%FileNr , * ) buffer
                       data%Comment = Trim( buffer )

                       !!! Read lattice constant
                       Read( data%FileNr , * ) data%latticeparameter

                       !!!Reading Lattice from poscar file
                       Write( * , * ) 'Lattice times lattice parameter'
                       Do i = 1 , 3
                          Read( data%FileNr , * ) data%lattice( : , i )
                          data%lattice( : , i ) = data%lattice( : , i )&
                                  * data%latticeparameter
                          Write( * , "(3F15.8)" ) data%lattice( : , i )
                       End Do

                       !!! determine atom types and their number
                       Call ReadLineUL( buffer , data%fileNr )
                       Call PrepareAtomTypes( buffer , data%AtomTypes )
                       data%NTypes = Size( data%AtomTypes , 1 )
                       Allocate( data%AtomNumbers( 1:data%NTypes ) )
                       Read( data%FileNr , * ) data%AtomNumbers

                       Write( * , "(A6,I2,A15,A)" ) 'Found ' , &
                               data%NTypes ,&
                                  ' Atom types in ', data%fname
                       Do i = 1, data%NTypes
                          Write( * , * ) data%AtomTypes( i ) , &
                                         data%AtomNumbers( i )
                       End Do

                       data%Natoms  =  Sum( data%AtomNumbers )

                 End Subroutine OpenPOSCARFile

                 !!!
                 !!!
                 !!! reading character line of unknown length
                 !!!

                 Subroutine ReadLineUL( line , fileNr )

                       Implicit None
                       Character( len = * )  ::line
                       Integer ::fileNr

                       Integer ::ioerr


                       ioerr = 0

                       Read( fileNr , "(A)", ADVANCE='NO', iostat =&
                              ioerr ) line
                 End Subroutine ReadLineUL


                 !!!! routine takes line with atom types as input and
                 !!!! counts the number of atom types
                 !!!! returns array with length number of atom types and 
                 !!!! the atom types
                 Subroutine PrepareAtomTypes( line , AtomTypes )
                       
                       Implicit None
                       Character( len = * ) ::line
                       Character( len = 2 ),allocatable  ::AtomTypes( : )

                       Integer ::i
                       Integer ::Natoms

                       Logical ::sampling

                       Character( len = 2 )  ::Atom

                       Natoms = 0

                       sampling = .False.

                       !!!
                       !!!
                       !!!  count the number of atoms
                       !!!
                       Do i = 1 , 100
                          !! check for space
                          If ( Sampling ) then
                             If ( line(i:i) .ne. ' ' ) then
                                  Atom = Atom // line( i:i )
                             Else
                                  Sampling =  .False.
                                  Natoms = Natoms + 1
                             End If
                          Else
                             If ( line(i:i) .ne. ' ' ) then
                                     Sampling = .True.
                                     Atom = line( i:i )
                             End If
                          End If
                       End Do

                       Allocate( AtomTypes( 1:Natoms ) )


                       Natoms = 0
                       Do i = 1 , 100
                          !! check for space
                          If ( Sampling ) then
                             If ( line(i:i) .ne. ' ' ) then
                                  Atom = Trim(Atom) // line( i:i )
                             Else
                                  Sampling =  .False.
                                  Natoms = Natoms + 1
                                  AtomTypes( Natoms ) = Atom
                             End If
                          Else
                             If ( line(i:i) .ne. ' ' ) then
                                     Sampling = .True.
                                     Atom = line( i:i )
                             End If
                          End If
                       End Do


                End Subroutine PrepareAtomTypes



                !!!
                !!!
                !!! reading all atoms in a single array called array
                !!!
                Subroutine ReadingCoordinatesATO( data , array )

                        Implicit None
                        Type( INPUT ) ::data
                        Real*8,dimension( : , : )  ::array

                        Integer ::i
                        Integer ::ioerr

                        ioerr = 0

                        Read( data%FileNr , * ) !! read direct configuration line
                                      !! in XDATCAR
                        Do i = 1 , data%Natoms
                           Read( data%fileNr , * , iostat = ioerr ) &
                                   array( : , i )
                           If ( ioerr .ne. 0 ) then
                              Write( * , "(A,A)" ) 'Error in input'&
                                    ,'reader '
                              Write( * , * ) 'Terminateing read'
                              Exit
                           End If
                        End Do


                End Subroutine ReadingCoordinatesATO



                !!!
                !!! reading atom into a atom type structure
                !!! containing for every atom type an index.
                !!! within this index all coords in direct and
                !!! cartesian are stored for a single type
                !!!
                Subroutine ReadingCoordinatesAtomType( data , atoms ,&
                                       OK )

                     Implicit None
                     Type( INPUT ) :: data
                     Type( atom_type ),allocatable,dimension( : )::atoms

                     Integer ::i
                     Integer ::j
                     Integer ::ioerr
                     Logical,optional  ::OK

                     ioerr  =  0
                     If ( present( OK ) ) then
                       OK = .True.
                     End If

                     Read( data%FileNr , * , iostat = ioerr )
                     If ( ioerr .ne. 0 ) then
                         Write( * , "(A,A)" ) 'Error in input'&
                               ,'reader '
                         Write( * , * ) 'Terminateing read'
                         If ( present( OK ) ) then
                            OK = .False.
                         End If
                         Return
                     End If


                     Do i = 1 , data%NTypes
                       Do j = 1 , data%AtomNumbers( i )

                           Read( data%fileNr , * , iostat = ioerr ) &
                                   atoms( i )%direct_coords( : , j )
                           atoms( i )%cartesian_coords( : , j ) = &
                                   MatMul( data%lattice , &
                                      atoms( i )%direct_coords( : , j ))
                           If ( ioerr .ne. 0 ) then
                              Write( * , "(A,A)" ) 'Error in input'&
                                    ,'reader '
                              Write( * , * ) 'Terminateing read'
                              If ( present( OK ) ) then
                                 OK = .False.
                              End If
                              Return
                           End If
                       End Do
                     End Do
                     Return

                End Subroutine ReadingCoordinatesAtomType


                !!!
                !!! initializing data atom type data staructure
                !!! has to be called after atom types 
                !!! were determined by OpenPOSCARFile
                !!!
                Subroutine AtomTypeInit( InputParams , atoms )
                    
                    Implicit None
                    Type( INPUT ) ::InputParams
                    Type( atom_type ),allocatable,dimension( : )::atoms

                    Integer ::i

                    Write( * , * ) "Initializing atom structure"
                    Allocate( atoms( 1:InputParams%NTypes ) )

                    Do i = 1 , InputParams%NTypes
                       Allocate( atoms( i )%direct_coords( 1:3 , &
                                 1:InputParams%AtomNumbers( i ) ) )
                       Allocate( atoms( i )%cartesian_coords( 1:3 , &
                                 1:InputParams%AtomNumbers( i ) ) )
                       atoms( i )%Natoms  =  InputParams%AtomNumbers( i )
                       atoms( i )%Nspecies =  InputParams%NTypes
                       atoms( i )%Atype  =  InputParams%AtomTypes( i )
                       atoms( i )%direct_coords = 66d0
                    End Do

                End Subroutine AtomTypeInit


      End Module ReadingInput


