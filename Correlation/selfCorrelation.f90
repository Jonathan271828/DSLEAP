  


  !!
  !!
  !! Module is used to compute a time averaged
  !! autocorrelation for scalar quantities
  !! 
  !!
  !!
  Module SelfCorrScalar  !j.L 26.06.2019

    Implicit None
    Type selfCorrVar
      Integer ::TotLen                    !! Total number of MD steps
      Integer ::WinShift                  !! number after which a new sample is taken
      Real*8,allocatable,dimension( : , : , : )  ::startConf  !! storing starting configurations
      Real*8,allocatable,dimension( : , : ) ::CorrFuncs       !! correlation functions
      Integer  ::Nstart                   !! counting number of samples to start
      Integer,allocatable,dimension( : )  ::StepCounter!! Counts steps of correlation functions independently 
      Real*8,allocatable,dimension( : )  ::Average
      Integer ::MDstep
      Integer ::NNStartConfs
    End Type selfCorrVar

    Contains

      Subroutine SelfCorrelationInit( data , Window , MDSteps , Nbonds , dimen )
        !!
        !! initializes stuff, how many starting configurations are needed
        !! 
        !!

        Implicit None
        Type( selfCorrVar )  ::data
        Integer ::Window
        Integer ::MDSteps
        Integer ::Nbonds     !! number of contributing entities
        Integer,Save ::NCorrs = 1
        Integer,optional ::dimen

        data%WinShift =  Window
        data%TotLen   =  MDSteps
        data%Nstart   =  Ceiling( 0.5 * MDSteps / Window )
        
        Write( * , * ) "Self correlation is averaging over ",data%Nstart
        Write( * , * ) "starting configurations in configuration" , NCorrs
        
        If ( present( dimen ) ) then
           Allocate( data%startConf( 1:dimen , 1:Nbonds , 1:data%Nstart ) )
        Else
           Allocate( data%startConf( 1:3 , 1:Nbonds , 1:data%Nstart ) )
        End If

        data%startConf = 0d0
        Allocate( data%CorrFuncs( 1:MDSteps , 1:data%Nstart ) )
        data%CorrFuncs = 0d0
        Allocate( data%StepCounter( 1:data%Nstart ) )
        data%StepCounter = 0

        data%NNStartConfs = 0
        data%MDstep = 1

        NCorrs  =  NCorrs  + 1
      End Subroutine SelfCorrelationInit



      Subroutine SelfCorrelationMain( data , Vectors )

        !!
        !!
        !!  Main routine that has to be called at every MD step
        !!  that should be analyzed; MD starts new trajectories as
        !!  the step is lower than the total length of the simulation
        !!  routine initializes new self correlation after WinShift
        !!  steps; the function ComputeCorrelationValue can be exchanged for any other
        !!  desired function that maps the initial and the current porerty that
        !!  has to be correlated to a scalar 
        !!

        Implicit None
        Type( selfCorrVar ) ::data
        Real*8,dimension( : , : ) ::Vectors
        Integer ::i
        Integer ::j



        If ( Modulo( data%MDstep , data%WinShift ) .eq. 0 .or. data%MDstep .eq. 1 ) then
          If ( data%MDstep .lt. data%TotLen * 0.5 ) then
             data%NNStartConfs  =  data%NNStartConfs  +  1
             data%startConf( : , : , data%NNStartConfs )  =  Vectors
          End If
        End If

        Do i = 1 , data%NNStartConfs
           data%StepCounter( i ) = data%StepCounter( i ) + 1
           Do  j = 1, Size( Vectors , 2 )
              data%CorrFuncs( data%StepCounter( i ) , i ) =  data%CorrFuncs( data%StepCounter( i ) , i ) + &
                 ComputeCorrelationValue( Vectors( : , j ) , data%startConf( : , j , i ) )
           End Do
        End Do
        data%MDstep = data%MDstep + 1

      End Subroutine SelfCorrelationMain


      Function ComputeCorrelationValue( VectorA , VectorB ) Result( dotP )

        use tool_interfaces , only :get_norm
        Implicit None
        Real*8,dimension( : )  ::VectorA
        Real*8,dimension( : )  ::VectorB
        Real*8  ::dotP

        dotP  =  Dot_product( vectorA , vectorB ) 


      End Function ComputeCorrelationValue



      Subroutine FinalizeSelfCorr( data )

        !!!
        !!
        !! Function computes the time average of the different starting points
        !! that were taken during the analysis
        !!

        Implicit None
        Type( SelfCorrVar ) ::data
        Integer ::i
        Integer ::Nbonds

        Nbonds = Size( data%startConf , 2 )


        Allocate( data%Average( 1:data%StepCounter( data%Nstart ) ) )
        data%Average  =  0d0

        Do i = 1 , data%Nstart
           data%Average( : ) = data%Average( : ) + data%CorrFuncs( 1:data%StepCounter( data%Nstart ) , i ) &
                            / Dble( data%Nstart ) / Dble( Nbonds )
        End Do

        Deallocate( data%CorrFuncs )
        Deallocate( data%StepCounter )
        Deallocate( data%startConf )


      End Subroutine FinalizeSelfCorr



      Subroutine WriteSelfCorrelation( data , Fname )

        !!!
        !! Subroutine writes time average to desired file
        !!

        Implicit None
        Type( SelfCorrVar ) ::data
        Character( len = * )  ::Fname
        Integer ::i

        open( unit = 3000000 ,action = 'write' , file = Trim( Fname ) , status='unknown' )
          Do i = 1, Size( data%Average , 1 )
             Write( 3000000 , "(I10,F15.8)" ) i , data%Average( i )
          End Do
        close( 3000000 )

        Deallocate( data%Average )
      End Subroutine WriteSelfCorrelation

  End Module SelfCorrScalar


  !!
  !!
  !! Module is used to compute a time averaged
  !! autocorrelation for 3 dimensional quantities
  !! 
  !!
  !!
  Module SelfCorr  !j.L 12.12.2018

    Implicit None
    Type selfCorrVar
      Integer ::TotLen                    !! Total number of MD steps
      Integer ::WinShift                  !! number after which a new sample is taken
      Real*8,allocatable,dimension( : , : , : )  ::startConf  !! storing starting configurations
      Real*8,allocatable,dimension( : , : ) ::CorrFuncs       !! correlation functions
      Integer  ::Nstart                   !! counting number of samples to start
      Integer,allocatable,dimension( : )  ::StepCounter!! Counts steps of correlation functions independently 
      Real*8,allocatable,dimension( : )  ::Average
      Integer ::MDstep
      Integer ::NNStartConfs
    End Type selfCorrVar

    Contains

      Subroutine SelfCorrelationInit( data , Window , MDSteps , Nbonds , dimen )
        !!
        !! initializes stuff, how many starting configurations are needed
        !! 
        !!

        Implicit None
        Type( selfCorrVar )  ::data
        Integer ::Window
        Integer ::MDSteps
        Integer ::Nbonds     !! number of contributing entities
        Integer,Save ::NCorrs = 1
        Integer,optional ::dimen

        data%WinShift =  Window
        data%TotLen   =  MDSteps
        data%Nstart   =  Ceiling( 0.5 * MDSteps / Window )
        
        Write( * , * ) "Self correlation is averaging over ",data%Nstart
        Write( * , * ) "starting configurations in configuration" , NCorrs
        
        If ( present( dimen ) ) then
           Allocate( data%startConf( 1:dimen , 1:Nbonds , 1:data%Nstart ) )
        Else
           Allocate( data%startConf( 1:3 , 1:Nbonds , 1:data%Nstart ) )
        End If

        data%startConf = 0d0
        Allocate( data%CorrFuncs( 1:MDSteps , 1:data%Nstart ) )
        data%CorrFuncs = 0d0
        Allocate( data%StepCounter( 1:data%Nstart ) )
        data%StepCounter = 1

        data%NNStartConfs = 0
        data%MDstep = 1

        NCorrs  =  NCorrs  + 1
      End Subroutine SelfCorrelationInit



      Subroutine SelfCorrelationMain( data , Vectors )

        !!
        !!
        !!  Main routine that has to be called at every MD step
        !!  that should be analyzed; MD starts new trajectories as
        !!  the step is lower than the total length of the simulation
        !!  routine initializes new self correlation after WinShift
        !!  steps; the function ComputeCorrelationValue can be exchanged for any other
        !!  desired function that maps the initial and the current porerty that
        !!  has to be correlated to a scalar 
        !!

        Implicit None
        Type( selfCorrVar ) ::data
        Real*8,dimension( : , : ) ::Vectors
        Integer ::i
        Integer ::j



        If ( Modulo( data%MDstep , data%WinShift ) .eq. 0 .or. data%MDstep .eq. 1 ) then
          If ( data%MDstep .lt. data%TotLen * 0.5 ) then
             data%NNStartConfs  =  data%NNStartConfs  +  1
             data%startConf( : , : , data%NNStartConfs )  =  Vectors
          End If
        End If

        Do i = 1 , data%NNStartConfs
           Do  j = 1, Size( Vectors , 2 )
              data%CorrFuncs( data%StepCounter( i ) , i ) =  data%CorrFuncs( data%StepCounter( i ) , i ) + &
                 ComputeCorrelationValue( Vectors( : , j ) , data%startConf( : , j , i ) )
           End Do
           data%StepCounter( i ) = data%StepCounter( i ) + 1
        End Do
        data%MDstep = data%MDstep + 1

      End Subroutine SelfCorrelationMain


      Function ComputeCorrelationValue( VectorA , VectorB ) Result( dotP )

        use tool_interfaces , only :get_norm
        Implicit None
        Real*8,dimension( : )  ::VectorA
        Real*8,dimension( : )  ::VectorB
        Real*8  ::dotP

        dotP  =  Dot_product( vectorA , vectorB ) &
                / get_norm( VectorA ) / get_norm( VectorB )

      End Function ComputeCorrelationValue



      Subroutine FinalizeSelfCorr( data )

        !!!
        !!
        !! Function computes the time average of the different starting points
        !! that were taken during the analysis
        !!

        Implicit None
        Type( SelfCorrVar ) ::data
        Integer ::i
        Integer ::Nbonds

        Nbonds = Size( data%startConf , 2 )


        data%StepCounter  =  data%StepCounter - 1
        Allocate( data%Average( 1:data%StepCounter( data%Nstart ) ) )
        data%Average  =  0d0

        Do i = 1 , data%Nstart
           data%Average( : ) = data%Average( : ) + data%CorrFuncs( 1:data%StepCounter( data%Nstart ) , i ) &
                            / Dble( data%Nstart ) / Dble( Nbonds )
        End Do

        Deallocate( data%CorrFuncs )
        Deallocate( data%StepCounter )
        Deallocate( data%startConf )


      End Subroutine FinalizeSelfCorr



      Subroutine WriteSelfCorrelation( data , Fname )

        !!!
        !! Subroutine writes time average to desired file
        !!

        Implicit None
        Type( SelfCorrVar ) ::data
        Character( len = * )  ::Fname
        Integer ::i

        open( unit = 3000000 ,action = 'write' , file = Trim( Fname ) , status='unknown' )
          Do i = 1, Size( data%Average , 1 )
             Write( 3000000 , "(I10,F15.8)" ) i , data%Average( i )
          End Do
        close( 3000000 )

        Deallocate( data%Average )
      End Subroutine WriteSelfCorrelation

  End Module SelfCorr










  
  !!
  !!
  !! Module is used to compute a time averaged
  !! autocorrelation for 3 dimensional quantities
  !! 
  !!
  !!
  Module CrossCorr  !j.L 12.12.2018

    Implicit None
    Type CrossCorrVar
      Integer ::TotLen                    !! Total number of MD steps
      Integer ::WinShift                  !! number after which a new sample is taken
      Real*8,allocatable,dimension( : , : )  ::startConf  !! storing starting configurations
      Real*8,allocatable,dimension( : , : , : ) ::CorrFuncs       !! correlation functions
      Integer  ::Nstart                   !! counting number of samples to start
      Integer,allocatable,dimension( : )  ::StepCounter!! Counts steps of correlation functions independently 
      Real*8,allocatable,dimension( : , : )  ::Average
      Integer ::MDstep
      Integer ::NNStartConfs
      Integer ::Ntypes
    End Type CrossCorrVar

    Contains

      Subroutine CrossCorrelationInit( data , Window , MDSteps , Ntypes )
        !!
        !! initializes stuff, how many starting configurations are needed
        !! 
        !!

        Implicit None
        Type( CrossCorrVar )  ::data
        Integer ::Window
        Integer ::MDSteps
        Integer,Save ::NCorrs = 1
        Integer :: Ntypes 

        data%WinShift =  Window
        data%TotLen   =  MDSteps
        data%Nstart   =  Ceiling( 0.5 * MDSteps / Window )
        data%Ntypes   =  Ntypes
        
        Write( * , * ) "Cross correlation is averaging over ",data%Nstart
        Write( * , * ) "starting configuration correlation nr" , NCorrs
        
        Allocate( data%startConf( 1:Ntypes , 1:data%Nstart ) )
        data%startConf = 0d0
        Allocate( data%CorrFuncs( 1:Ntypes , 1:MDSteps , 1:data%Nstart ) )
        data%CorrFuncs = 0d0
        Allocate( data%StepCounter( 1:data%Nstart ) )
        data%StepCounter = 1

        data%NNStartConfs = 0
        data%MDstep = 1


        NCorrs  =  NCorrs  + 1
      End Subroutine CrossCorrelationInit



      Subroutine CrossCorrelationMain( data , Vectors )

        !!
        !!
        !!  Main routine that has to be called at every MD step
        !!  that should be analyzed; MD starts new trajectories as
        !!  the step is lower than the total length of the simulation
        !!  routine initializes new self correlation after WinShift
        !!  steps; the function ComputeCorrelationValue can be exchanged for any other
        !!  desired function that maps the initial and the current porerty that
        !!  has to be correlated to a scalar 
        !!

        Implicit None
        Type( CrossCorrVar ) ::data
        Real*8,dimension( : ) ::Vectors
        Integer ::i



        If ( Modulo( data%MDstep , data%WinShift ) .eq. 0 .or. data%MDstep .eq. 1 ) then
          If ( data%MDstep .lt. data%TotLen * 0.5 ) then
             data%NNStartConfs  =  data%NNStartConfs  +  1
             data%startConf( : , data%NNStartConfs )  =  Vectors
          End If
        End If

        Do i = 1 , data%NNStartConfs
           data%CorrFuncs( : , data%StepCounter( i ) , i ) =  data%CorrFuncs( : , data%StepCounter( i ) , i ) + &
                                                               Vectors( : ) * data%startConf( : , i )
           data%StepCounter( i ) = data%StepCounter( i ) + 1
        End Do
        data%MDstep = data%MDstep + 1

      End Subroutine CrossCorrelationMain


      Function ComputeCrossCorrelationValue( DatA , DatB ) Result( dotP )

        use tool_interfaces , only :get_norm
        Implicit None
        Real*8  ::DatA
        Real*8  ::DatB
        Real*8  ::dotP

        dotp = DatA * DatB 

      End Function ComputeCrossCorrelationValue



      Subroutine FinalizeCrossCorr( data )

        !!!
        !!
        !! Function computes the time average of the different starting points
        !! that were taken during the analysis
        !!

        Implicit None
        Type( CrossCorrVar ) ::data
        Integer ::i

        data%StepCounter  =  data%StepCounter - 1
        Allocate( data%Average( 1:data%Ntypes , 1:data%StepCounter( data%Nstart ) ) )
        data%Average  =  0d0

        Do i = 1 , data%Nstart
           data%Average( : , : ) = data%Average( : , : ) + data%CorrFuncs( : , 1:data%StepCounter( data%Nstart ) , i ) &
                            / Dble( data%Nstart )
        End Do

        Deallocate( data%CorrFuncs )
        Deallocate( data%StepCounter )
        Deallocate( data%startConf )


      End Subroutine FinalizeCrossCorr



      Subroutine WriteCrossCorrelation( data , Fname )

        !!!
        !! Subroutine writes time average to desired file
        !!

        Implicit None
        Type( CrossCorrVar ) ::data
        Character( len = * )  ::Fname
        Integer ::i

        open( unit = 3000000 ,action = 'write' , file = Trim( Fname ) , status='unknown' )
          Do i = 1, Size( data%Average , 1 )
             Write( 3000000 , "(I10,30000F15.8)" ) i , data%Average( : , i )
          End Do
        close( 3000000 )

        Deallocate( data%Average )
      End Subroutine WriteCrossCorrelation

  End Module CrossCorr


  !!
  !!
  !! Module is used to compute a time averaged
  !! autocorrelation for 1 dimensional quantities
  !! 
  !!
  !!
  Module CrossCorrComplex  !j.L 02.11.2019

    Implicit None
    Type CrossCorrVar
      Integer ::TotLen                    !! Total number of MD steps
      Integer ::WinShift                  !! number after which a new sample is taken
      Complex*16,allocatable,dimension( : , : )  ::startConf  !! storing starting configurations
      Complex*16,allocatable,dimension( : , : , : ) ::CorrFuncs       !! correlation functions
      Integer  ::Nstart                   !! counting number of samples to start
      Integer,allocatable,dimension( : )  ::StepCounter!! Counts steps of correlation functions independently 
      Complex*16,allocatable,dimension( : , : )  ::Average
      Integer ::MDstep
      Integer ::NNStartConfs
      Integer ::Ntypes
    End Type CrossCorrVar

    Contains

      Subroutine CrossCorrelationInit( data , Window , MDSteps , Ntypes )
        !!
        !! initializes stuff, how many starting configurations are needed
        !! 
        !!

        Implicit None
        Type( CrossCorrVar )  ::data
        Integer ::Window
        Integer ::MDSteps
        Integer,Save ::NCorrs = 1
        Integer :: Ntypes 

        data%WinShift =  Window
        data%TotLen   =  MDSteps
        data%Nstart   =  Ceiling( 0.5 * MDSteps / Window )
        data%Ntypes   =  Ntypes
        
        Write( * , * ) "Cross correlation is averaging over ",data%Nstart
        Write( * , * ) "starting configuration correlation nr" , NCorrs
        
        Allocate( data%startConf( 1:Ntypes , 1:data%Nstart ) )
        data%startConf = 0d0
        Allocate( data%CorrFuncs( 1:Ntypes , 1:MDSteps , 1:data%Nstart ) )
        data%CorrFuncs = 0d0
        Allocate( data%StepCounter( 1:data%Nstart ) )
        data%StepCounter = 1

        data%NNStartConfs = 0
        data%MDstep = 1


        NCorrs  =  NCorrs  + 1
      End Subroutine CrossCorrelationInit



      Subroutine CrossCorrelationMain( data , Vectors , VectorsB , Norm )

        !!
        !!
        !!  Main routine that has to be called at every MD step
        !!  that should be analyzed; MD starts new trajectories as
        !!  the step is lower than the total length of the simulation
        !!  routine initializes new self correlation after WinShift
        !!  steps; the function ComputeCorrelationValue can be exchanged for any other
        !!  desired function that maps the initial and the current porerty that
        !!  has to be correlated to a scalar 
        !!

        Implicit None
        Type( CrossCorrVar ) ::data
        Complex*16,dimension( : ) ::Vectors
        Complex*16,optional,dimension( : ) ::VectorsB
        Real*8,optional :: Norm
        Integer ::i



        If ( Modulo( data%MDstep , data%WinShift ) .eq. 0 .or. data%MDstep .eq. 1 ) then
          If ( data%MDstep .lt. data%TotLen * 0.5 ) then
             data%NNStartConfs  =  data%NNStartConfs  +  1
             if ( present( VectorsB  ) ) then
                     data%startConf( : , data%NNStartConfs )  =  DConjg( VectorsB )
             Else
                     data%startConf( : , data%NNStartConfs )  =  DConjg( Vectors )
             End If
          End If
        End If

        Do i = 1 , data%NNStartConfs
           data%CorrFuncs( : , data%StepCounter( i ) , i ) =  data%CorrFuncs( : , data%StepCounter( i ) , i ) + &
                                                               Vectors( : ) * data%startConf( : , i )
           data%StepCounter( i ) = data%StepCounter( i ) + 1
        End Do
        data%MDstep = data%MDstep + 1

      End Subroutine CrossCorrelationMain


      Function ComputeCrossCorrelationValue( DatA , DatB ) Result( dotP )

        use tool_interfaces , only :get_norm
        Implicit None
        Real*8  ::DatA
        Real*8  ::DatB
        Real*8  ::dotP

        dotp = DatA * DatB 

      End Function ComputeCrossCorrelationValue



      Subroutine FinalizeCrossCorr( data )

        !!!
        !!
        !! Function computes the time average of the different starting points
        !! that were taken during the analysis
        !!

        Implicit None
        Type( CrossCorrVar ) ::data
        Integer ::i


        data%StepCounter  =  data%StepCounter - 1
        Allocate( data%Average( 1:data%Ntypes , 1:data%StepCounter( data%Nstart ) ) )
        data%Average  =  0d0
           

        Do i = 1 , data%Nstart
           data%Average( : , : ) = data%Average( : , : ) + data%CorrFuncs( : , 1:data%StepCounter( data%Nstart ) , i ) &
                            / Dble( data%Nstart )
        End Do

        Deallocate( data%CorrFuncs )
        Deallocate( data%StepCounter )
        Deallocate( data%startConf )


      End Subroutine FinalizeCrossCorr



      Subroutine WriteCrossCorrelation( data , Fname )

        !!!
        !! Subroutine writes time average to desired file
        !!

        Implicit None
        Type( CrossCorrVar ) ::data
        Character( len = * )  ::Fname
        Integer ::i

        open( unit = 3000000 ,action = 'write' , file = Trim( Fname ) , status='unknown' )
          Do i = 1, Size( data%Average , 1 )
             Write( 3000000 , "(I10,30000F15.8)" ) i , Abs( data%Average( : , i ) )
          End Do
        close( 3000000 )

        Deallocate( data%Average )
      End Subroutine WriteCrossCorrelation

  End Module CrossCorrComplex

  !!
  !!
  !! Module is used to compute a time averaged
  !! autocorrelation for 1 dimensional quantities
  !! 
  !!
  !!
  Module CrossCorrComplexND  !j.L 07.09.2020

    Implicit None
    Type CrossCorrVarND
      Integer ::TotLen                    !! Total number of MD steps
      Integer ::WinShift                  !! number after which a new sample is taken
      Complex*16,allocatable,dimension( : , : , : ) ::startConf  !! storing starting configurations
      Complex*16,allocatable,dimension( : , : , : ) ::CorrFuncs       !! correlation functions
      Integer  ::Nstart                   !! counting number of samples to start
      Integer,allocatable,dimension( : )  ::StepCounter!! Counts steps of correlation functions independently 
      Complex*16,allocatable,dimension( : , : )  ::Average
      Integer ::MDstep
      Integer ::NNStartConfs
      Integer ::Ntypes
    End Type CrossCorrVarND

    Contains

      Subroutine CrossCorrelationInitND( data , Window , MDSteps , Ntypes , SpD )
        !!
        !! initializes stuff, how many starting configurations are needed
        !! 
        !!

        Implicit None
        Type( CrossCorrVarND )  ::data
        Integer ::Window
        Integer ::MDSteps
        Integer,Save ::NCorrs = 1
        Integer :: Ntypes
        Integer,optional ::SpD
        Integer ::DimVec

        data%WinShift =  Window
        data%TotLen   =  MDSteps
        data%Nstart   =  Ceiling( 0.5 * MDSteps / Window )
        data%Ntypes   =  Ntypes






        If ( present( SpD ) ) then
                DimVec = SpD
        Else
                DimVec = 3
        End If

        
        Write( * , * ) "Cross correlation is averaging over ",data%Nstart
        Write( * , * ) "starting configuration correlation nr" , NCorrs
        
        Allocate( data%startConf( 1:DimVec , 1:Ntypes , 1:data%Nstart ) )
        data%startConf = 0d0
        Allocate( data%CorrFuncs( 1:Ntypes , 1:MDSteps , 1:data%Nstart ) )
        data%CorrFuncs = 0d0
        Allocate( data%StepCounter( 1:data%Nstart ) )
        data%StepCounter = 1

        data%NNStartConfs = 0
        data%MDstep = 1


        NCorrs  =  NCorrs  + 1
      End Subroutine CrossCorrelationInitND



      Subroutine CrossCorrelationMainND( data , Vectors , VectorsB , Norm )

        !!
        !!
        !!  Main routine that has to be called at every MD step
        !!  that should be analyzed; MD starts new trajectories as
        !!  the step is lower than the total length of the simulation
        !!  routine initializes new self correlation after WinShift
        !!  steps; the function ComputeCorrelationValue can be exchanged for any other
        !!  desired function that maps the initial and the current porerty that
        !!  has to be correlated to a scalar 
        !!

        Implicit None
        Type( CrossCorrVarND ) ::data
        Complex*16,dimension( : , : ) ::Vectors
        Complex*16,optional,dimension( : , : ) ::VectorsB
        Real*8,optional :: Norm
        Integer ::i
        Integer ::j



        If ( Modulo( data%MDstep , data%WinShift ) .eq. 0 .or. data%MDstep .eq. 1 ) then
          If ( data%MDstep .lt. data%TotLen * 0.5 ) then
             data%NNStartConfs  =  data%NNStartConfs  +  1
             if ( present( VectorsB  ) ) then
                     data%startConf( : , : , data%NNStartConfs )  =  DConjg( VectorsB )
             Else
                     data%startConf( : , : , data%NNStartConfs )  =  DConjg( Vectors )
             End If
          End If
        End If

        Do i = 1 , data%NNStartConfs
           Do j = 1 , Size( Vectors , 2 )
               data%CorrFuncs( j , data%StepCounter( i ) , i ) =  data%CorrFuncs( j , data%StepCounter( i ) , i ) + &
                                                                     Dot_Product( Vectors( : , j ) , data%startConf( : , j , i ) )
           End Do
           data%StepCounter( i ) = data%StepCounter( i ) + 1
        End Do
        data%MDstep = data%MDstep + 1

      End Subroutine CrossCorrelationMainND


      Function ComputeCrossCorrelationValue( DatA , DatB ) Result( dotP )

        use tool_interfaces , only :get_norm
        Implicit None
        Real*8  ::DatA
        Real*8  ::DatB
        Real*8  ::dotP

        dotp = DatA * DatB 

      End Function ComputeCrossCorrelationValue



      Subroutine FinalizeCrossCorrND( data )

        !!!
        !!
        !! Function computes the time average of the different starting points
        !! that were taken during the analysis
        !!

        Implicit None
        Type( CrossCorrVarND ) ::data
        Integer ::i


        data%StepCounter  =  data%StepCounter - 1
        Allocate( data%Average( 1:data%Ntypes , 1:data%StepCounter( data%Nstart ) ) )
        data%Average  =  0d0
           

        Do i = 1 , data%Nstart
           data%Average( : , : ) = data%Average( : , : ) + data%CorrFuncs( : , 1:data%StepCounter( data%Nstart ) , i ) &
                            / Dble( data%Nstart )
        End Do

        Deallocate( data%CorrFuncs )
        Deallocate( data%StepCounter )
        Deallocate( data%startConf )


      End Subroutine FinalizeCrossCorrND



      Subroutine WriteCrossCorrelationND( data , Fname )

        !!!
        !! Subroutine writes time average to desired file
        !!

        Implicit None
        Type( CrossCorrVarND ) ::data
        Character( len = * )  ::Fname
        Integer ::i

        open( unit = 3000000 ,action = 'write' , file = Trim( Fname ) , status='unknown' )
          Do i = 1, Size( data%Average , 1 )
             Write( 3000000 , "(I10,30000F15.8)" ) i , Abs( data%Average( : , i ) )
          End Do
        close( 3000000 )

        Deallocate( data%Average )
      End Subroutine WriteCrossCorrelationND

  End Module CrossCorrComplexND


  !!
  !!
  !! Module is used to compute a time averaged
  !! autocorrelation for 3 dimensional quantities
  !! 
  !!
  !!
  Module SelfCorrProjectedVAC  !j.L 12.12.2018

    Implicit None
    Type selfCorrVarCMPLX
      Integer ::TotLen                    !! Total number of MD steps
      Integer ::WinShift                  !! number after which a new sample is taken
      Complex*16,allocatable,dimension( : , : , : )  ::startConf  !! storing starting configurations
      Complex*16,allocatable,dimension( : , : , : , : ) ::CorrFuncs       !! correlation functions
      Integer  ::Nstart                   !! counting number of samples to start
      Integer,allocatable,dimension( : )  ::StepCounter!! Counts steps of correlation functions independently 
      Complex*16,allocatable,dimension( : , : , : )  ::Average
      Integer ::MDstep
      Integer ::NNStartConfs
    End Type selfCorrVarCMPLX

    Contains

      Subroutine SelfCorrelationInitVAC( data , Window , MDSteps , Nbonds , dimen )
        !!
        !! initializes stuff, how many starting configurations are needed
        !! 
        !!

        Implicit None
        Type( selfCorrVarCMPLX )  ::data
        Integer ::Window
        Integer ::MDSteps
        Integer ::Nbonds     !! number of contributing entities
        Integer,Save ::NCorrs = 1
        Integer,optional ::dimen

        data%WinShift =  Window
        data%TotLen   =  MDSteps
        data%Nstart   =  Ceiling( 0.5 * MDSteps / Window )
        
        Write( * , * ) "Self correlation complex VAC is averaging over ",data%Nstart
        Write( * , * ) "starting configurations in configuration" , NCorrs
        
        If ( present( dimen ) ) then
           Allocate( data%startConf( 1:dimen , 1:Nbonds , 1:data%Nstart ) )
           Allocate( data%CorrFuncs( 1:dimen , 1:Nbonds , 1:MDSteps , 1:data%Nstart ) )
        Else
           Allocate( data%startConf( 1:3 , 1:Nbonds , 1:data%Nstart ) )
           Allocate( data%CorrFuncs( 1:3 , 1:Nbonds , 1:MDSteps , 1:data%Nstart ) )
        End If

        data%startConf = 0d0
        data%CorrFuncs = 0d0
        Allocate( data%StepCounter( 1:data%Nstart ) )
        data%StepCounter = 1

        data%NNStartConfs = 0
        data%MDstep = 1

        NCorrs  =  NCorrs  + 1
      End Subroutine SelfCorrelationInitVAC



      Subroutine SelfCorrelationMainVAC( data , Vectors )

        !!
        !!
        !!  Main routine that has to be called at every MD step
        !!  that should be analyzed; MD starts new trajectories as
        !!  the step is lower than the total length of the simulation
        !!  routine initializes new self correlation after WinShift
        !!  steps; the function ComputeCorrelationValue can be exchanged for any other
        !!  desired function that maps the initial and the current porerty that
        !!  has to be correlated to a scalar 
        !!

        Implicit None
        Type( selfCorrVarCMPLX ) ::data
        Complex*16,dimension( : , : ) ::Vectors
        Integer ::i

        If ( Modulo( data%MDstep , data%WinShift ) .eq. 0 .or. data%MDstep .eq. 1 ) then
          If ( data%MDstep .lt. data%TotLen * 0.5 ) then
             data%NNStartConfs  =  data%NNStartConfs  +  1
             data%startConf( : , : , data%NNStartConfs )  =  Vectors
          End If
        End If

        Do i = 1 , data%NNStartConfs
           !Do  j = 1, Size( Vectors , 2 )
              data%CorrFuncs( : , : , data%StepCounter( i ) , i ) = Vectors( : , : ) * DConjG( data%startConf( : , : , i ) )
           !End Do
           data%StepCounter( i ) = data%StepCounter( i ) + 1
        End Do
        data%MDstep = data%MDstep + 1

      End Subroutine SelfCorrelationMainVAC


      Function ComputeCorrelationValue( VectorA , VectorB ) Result( dotP )

        use tool_interfaces , only :get_norm
        Implicit None
        Real*8,dimension( : )  ::VectorA
        Real*8,dimension( : )  ::VectorB
        Real*8  ::dotP

        dotP  =  Dot_product( vectorA , vectorB ) &
                / get_norm( VectorA ) / get_norm( VectorB )

      End Function ComputeCorrelationValue



      Subroutine FinalizeSelfCorrVAC( data )

        !!!
        !!
        !! Function computes the time average of the different starting points
        !! that were taken during the analysis
        !!

        Implicit None
        Type( SelfCorrVarCMPLX ) ::data
        Integer ::i
        Integer ::Nbonds

        Nbonds = Size( data%startConf , 2 )


        data%StepCounter  =  data%StepCounter - 1
        Allocate( data%Average( 1:Size( data%CorrFuncs , 1 ),&
                                1:Size( data%CorrFuncs , 2 ),&
                                1:data%StepCounter( data%Nstart ) ) )
        data%Average  =  0d0

        Do i = 1 , data%Nstart
           data%Average( : , : , : ) = data%Average( : ,  : , : ) + &
                       data%CorrFuncs( : , : , 1:data%StepCounter( data%Nstart ) , i ) &
                            / Dble( data%Nstart )
        End Do

        Deallocate( data%CorrFuncs )
        Deallocate( data%StepCounter )
        Deallocate( data%startConf )


      End Subroutine FinalizeSelfCorrVAC



      Subroutine WriteSelfCorrelationVAC( data , Fname )

        !!!
        !! Subroutine writes time average to desired file
        !!

        Implicit None
        Type( SelfCorrVarCMPLX ) ::data
        Character( len = * )  ::Fname
        Integer ::i
        Integer ::j
        Integer ::k

        open( unit = 3000000 ,action = 'write' , file = Trim( Fname ) , status='unknown' )
          Do i = 1, Size( data%Average , 3 )
             !!! I only write out the real part because I checked for several simulations
             !!! that the imaginary part is numerical noise for most of the computations
             !!! also if the imagianry part is not zero the crystal would not be stable
             !!! this would not make sense in an MD
             Write( 3000000 , "(I10,300000F15.8)" ) i , ( ( Real( data%Average( j , k , i ) ) , &
                                                    j = 1, Size( data%Average , 1 ) ) ,&
                                                    k = 1, Size( data%Average , 2 ) ) 
          End Do
        close( 3000000 )

        Deallocate( data%Average )
      End Subroutine WriteSelfCorrelationVAC

  End Module SelfCorrProjectedVAC
