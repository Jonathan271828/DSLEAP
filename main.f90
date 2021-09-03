



      Program MainPhononCollection
          !!!
          !!!
          !!!
          !!! this is the main code of the Phonon collection
          !!!
          !!!
          
          !use PhononInterface, only :PhononInterfaceMain
          use ReadInputFile, only :ReadInput,SimParams
          use ReadingInput, only  :INPUT,OpenPOSCARFile,atom_type,&
                                  AtomTypeInit,ReadCoordinates
          use PhononInterface

          Implicit None
          Type( SimParams ) ::MDParameters
          Type( INPUT ) ::StrucReader
          Type( atom_type ),allocatable,dimension( : ) ::Atoms
          Character( len = : ),allocatable  ::InputFile
          Logical ::ReadOK

          !!!
          !local variables only needed in main routine
          Integer ::Nstruc   !! actual structure to be analyzed


          !!! file containing input flags
          InputFile           =  'Phonon.in'
          StrucReader%FileNr  =  666

          !! open input parameter file and read data from it
          !! 
          Call ReadInput( MDParameters , InputFile )
          !!
          !! here we read the data of the XDATCAR file
          !! but could also be a POSCAR
          !! since you are doing phonon comuptations
          !! I hope you did them in the microcannoical
          !! ensemble otherise it does not make sense
          !! because you compute dynamic properties from
          !! autocorrelation functions
          !! variable cell shapes are not supported
          StrucReader%fname = './XDATCAR'
          Call OpenPOSCARFile( StrucReader )
          !!! initializing atom array
          Call AtomTypeInit( StrucReader , Atoms )

          !!! read first structure
          Call ReadCoordinates( StrucReader , Atoms )
          !!! do initialization
          Call PhononInterfaceInit( MDParameters , Atoms , &
                                    StrucReader%Lattice )

          Do Nstruc = 1 , MDParameters%Nend
             If ( Nstruc .ge. MDParameters%Nstart .and. &
                   Modulo( Nstruc , MDParameters%BlockSample ) .eq. 0 ) then
                   !! Call relevant phonon routines
                   Call PhononInterfaceMain( MDParameters , Atoms , &
                                             StrucReader%Lattice )
                   If ( Modulo( Nstruc , 10 ) .eq. 0 ) then
                           Write( * , * ) "Analysed " , Nstruc , &
                                          "structures"
                   End If
             End If
             Call ReadCoordinates( StrucReader , Atoms , ReadOK )
             If ( .not. ReadOK ) Exit
          End Do

          Call PhononInterfaceFinal



      End Program MainPhononCollection

