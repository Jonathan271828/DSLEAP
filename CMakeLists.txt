




project(PC Fortran)

#enable_language(Fortran)
SET(CMAKE_FC_COMPILE_NAME gfortran)

#binary director
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/build/modules)


if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set( CMAKE_Fortran_FLAGS         "-O3 -fopenmp -Wunused" )
    set( CMAKE_Fortran_FLAGS_DEBUG   "-O3 -fopenmp -g3 -fbounds-check -Wmaybe-uninitialized -Wuninitialized" )
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp" )
endif()





message( ${CMAKE_Fortran_FLAGS} )

#building numerical integration
add_subdirectory(integrate)
## building eval_tools
add_subdirectory(eval_tools)
## build Input Reader
add_subdirectory(ReadInputFiles)
#Matrix Inversion
add_subdirectory(MatInv)
## add MergeSort
add_subdirectory(Sort)
## correlation library
add_subdirectory(Correlation)
## fourier transform
add_subdirectory(DLfft)
## small statistics library
add_subdirectory(stat)
## add linkLists
add_subdirectory(linkList)
## Phonon Stuff
add_subdirectory(DynStruc)



## building Input File Reader
add_dependencies(eval_tools integrate)
## linking Input File Reader
add_dependencies(ReadInputFiles eval_tools)
## linking Correlation
add_dependencies(Correlation eval_tools)
## linking linkLists
add_dependencies(linkList 
	eval_tools 
	Sort)
## linking DynamicStructure factor
add_dependencies(DynStruc
	Correlation
	stat
	linkList
	ReadInputFiles
	StructureReader
	stat
	eval_tools
	Sort
	MatInv
	DLfft
	)







add_executable(PC main.f90)

target_link_libraries(PC
	DynStruc
	ReadInputFiles
	StructureReader
	Correlation
	linkList
	stat
	DLfft
	eval_tools
	integrate
	Sort
	MatInv
	)



