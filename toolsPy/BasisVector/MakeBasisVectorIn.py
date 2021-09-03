#!/usr/bin/env python3





import sys
import argparse
import numpy as np
import yaml



def getOptions( args = sys.argv[1:] ):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "-yaml" , "--phonoPy", type = str , help = "The band.yaml file of PhononPy" )
    parser.add_argument( "-Nx" , "--NcellsX" , type=int , help = "Number of replicated unit cells in x direction" )
    parser.add_argument( "-Ny" , "--NcellsY" , type=int , help = "Number of replicated unit cells in y direction" )
    parser.add_argument( "-Nz" , "--NcellsZ" , type=int , help = "Number of replicated unit cells in z direction" )
    options = parser.parse_args( args )
    return options














def ExtractEigenVectors( fname ):
    """
    extract phonon eigenvectors from phonpy
    band.yaml file
    indexes of return array
    0th index -> qpoint
    1st index -> branch
    2nd index -> atom
    3rd index -> cartesian components 
    """
    eigenvector = []
    with open( fname , 'r' ) as stream:
        data = yaml.safe_load( stream )
        for qpoint in  data["phonon"][:]:
            temp = [ bands["eigenvector"] for bands in qpoint["band"] ]
            eigenvector.append( temp )
    eigenvector = np.asarray( eigenvector )
    return eigenvector[ : , : , :  , : , 0 ]



def ExtractQpositions( fname ):
    """
    extracting all q-vectors of phonopy
    """
    qPos = []
    with open( fname , 'r' ) as stream:
        data = yaml.safe_load( stream )
        for qpoint in  data["phonon"][:]:
            qPos.append( qpoint[ "q-position"] )

    qPos = np.asarray( qPos )
    return qPos







def MakeVecToString( x ):
    line = "  "
    for a in x:
        aout = "{:15.8f}".format( a )
        line += aout
    return line





class LoadEigenVectors:
    """
    class uses a PhonoPy yaml file to extract
    QVectors and the harmonic eigenvectors
    These are then written to a BasisVector.in file
    in the correct format for the PC code
    fnamePhonon --> contains the name of the phonopy yaml file
    Nx          --> contains the number of replicated uni cells in x
    Ny          --> contains the number of replicated uni cells in y
    Nz          --> contains the number of replicated uni cells in z
    """
    def __init__( self , fnamePhono , Nx , Ny , Nz ):
        if ( fnamePhono == None ):
            print( "No phonopy band.yaml file supplied" )
            print( "please supply the name with -yaml FileName" )
            sys.exit()
        if ( Nx == None or Ny == None or Nz == None ):
            print( "No cell dimensions are given" )
            print( "Please supply with -Nx [number of cells in x] " )
            print( "Please supply with -Ny [number of cells in y] " )
            print( "Please supply with -Nz [number of cells in z] " )
            sys.exit()
        self.Phonopy = ExtractEigenVectors( fnamePhono )
        self.Qpoints = ExtractQpositions( fnamePhono )
        self.NN = np.asarray( [ Nx , Ny , Nz ] )
        self.MakeIntegerIndices()


    def MakeIntegerIndices( self ):
        for i in range( self.Qpoints.shape[0] ):
            self.Qpoints[ i , : ] = np.asarray( [ self.NN[j]*self.Qpoints[ i , j ] 
                                        for j in range( self.NN.shape[0] ) ] )


    def WriteOutput( self ):
        """
        Writing the output in the DSLEAP format
        can be directly used as input
        """
        with open( "BasisVector.in" , "w" ) as outfile:
            line = "   " + str( self.Phonopy.shape[0] ) + "   " \
                         + str( self.Phonopy.shape[1] ) + "   " + str( self.Phonopy.shape[2] ) + \
                   "   " + str( self.Phonopy.shape[3] )
            outfile.write( line + "\n" )
            for qq in range( self.Phonopy.shape[0] ):
                line = MakeVecToString( self.Qpoints[ qq , : ] ) 
                outfile.write( line + "\n" )
                for branch in range( self.Phonopy.shape[1] ):
                    for atom in range( self.Phonopy.shape[2] ):
                        line = MakeVecToString( self.Phonopy[ qq , branch, atom , : ] )
                        outfile.write( line + "\n" )
                    outfile.write( '\n' )
                outfile.write( '\n' )












if __name__=="__main__":
    options = getOptions( sys.argv[1:] )
    x = LoadEigenVectors( options.phonoPy , options.NcellsX , 
                                            options.NcellsY , 
                                            options.NcellsZ )
    x.WriteOutput()
