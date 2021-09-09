#!/usr/bin/env python3



import numpy as np
import yaml
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import orthogonal_procrustes
import argparse



__author__ = "Jonathan Lahnsteiner"           
__copyright__ = ""                            
__license__ = "GPL"                           
__version__ = "1.0.0"                         
__maintainer__ = "Jonathan Lahnsteiner"       
__email__ = "jonathan.lahnsteiner@gmx.at"     
__status__ = "Production"                     








def getOptions( args = sys.argv[1:] ):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "-f" , "--basisVectorfile", 
             help = "Your Phonopy input file from which you want to extract data qpoints.yaml or band.yaml" ,
             default = "qpoints.yaml" )
    options = parser.parse_args( args )
    return options





def ExtractStructure( fname ):
    """
    extract structure from a band.yaml phonopy
    or qpoints.yaml
    file
    return lattice contains
    bravais matrix
    structure is an array from 0 to Natoms-1
    containing dictionaries
    coordinates -> direct coordinates
    mass -> atomic mass
    symbol -> atom type
    """
    with open( fname ) as infile:
        data = yaml.safe_load( infile )
        lattice =  np.asarray( data["lattice"] )
        struc   =  data["points"]
    return lattice, struc



def ExtractFrequencies( fname ):
    """
    extract band structure from phonopy yaml
    file.
    return bands first dimension
    is qpath ,second is branch index
    Qdist contains the qpath
    as a 1D array
    """
    bands = []
    Qdist = []
    with open( fname , 'r' ) as stream:
        data = yaml.safe_load( stream )
        for qpoint in  data["phonon"][:]:
            Qdist.append( qpoint["distance"] )
            temp = [ bands["frequency"] for bands in qpoint["band"] ]
            bands.append( temp )
    bands = np.asarray( bands )
    Qdist = np.asarray( Qdist )
    return bands , Qdist



def ExtractEigenVectors( fname ):
    """
    extract phonon eigenvectors from phonpy
    band.yaml or qpoint.yaml file
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



def ExtractQpath( fname ):
    qpath = []
    with open( fname , 'r' ) as stream:
        data = yaml.safe_load( stream )
        for qpoint in  data["phonon"][:]:
            qpath.append( qpoint["distance"] )
    return np.asarray( qpath )



def ExtractQpositions( fname ):
    qPos = []
    with open( fname , 'r' ) as stream:
        data = yaml.safe_load( stream )
        for qpoint in  data["phonon"][:]:
            qPos.append( qpoint[ "q-position"] )

    qPos = np.asarray( qPos )
    return qPos



def ExtractAtomTypes( structure ):
    """
    extract different types from
    atom structure used in phonopy
    atoms will be counted also
    return types , numbertypes
    """
    types = []
    for i in range( len( structure ) ):
        act = structure[i][ "symbol" ]
        include = True
        for comp in types:
            if ( act == comp ):
                include = False
        if ( include ):
            types.append( act )
    Ntypes = np.zeros( [ len( types ) ] )
    for i in range( len( types ) ):
        for j in range( len( structure ) ):
            if ( types[i] == structure[j]["symbol"] ):
                Ntypes[i] += 1
    return types , Ntypes



class MakeEigenVectorFile:
    def __init__( self , fname ):
        """
        fname -> name of band.yaml or qpoints.yaml file
        0th index -> qpoint
        1st index -> branch
        2nd index -> atom
        3rd index -> cartesian components 
        """
        print( "Reading input" )
        self.EigenVectors  =  ExtractEigenVectors( fname ) 
        self.QVectors      =  ExtractQpositions( fname )

        print( self.QVectors )




    def WriteEigenVectorsDSLEAPFormat( self ):
        """
        write all eigenvetors to file 
        BasisVector.in for the DSLEAP routine
        to compute projected VACF correlation function
        """
        with open( "BasisVector.in" , 'w' ) as outfile:
            outfile.write( "  " + str( self.QVectors.shape[0] ) + 
                           "  " + str( self.EigenVectors.shape[1] ) + 
                           "  " + str( self.EigenVectors.shape[2] ) + "   3 \n" )

            qq  =  0
            for i in range( self.QVectors.shape[0] ):
                txt = "  ".join( [ "{:15.8f}".format( x ) for x in self.QVectors[ qq , : ] ] ) + "\n"
                outfile.write( txt )
                for j in range( self.EigenVectors.shape[1] ):
                    data = np.reshape( self.EigenVectors[ i , j , : , : ] , 
                                    self.EigenVectors.shape[2]*self.EigenVectors.shape[3] )
                    for k in range( self.EigenVectors.shape[2] ):
                        data = "  ".join( [ "{:15.8f}".format( 
                               self.EigenVectors[ i , j , k , l ] )
                                  for l in range( self.EigenVectors.shape[3] ) ] )
                        outfile.write( data + "\n" )
                    outfile.write( "\n" )
                outfile.write( "\n" )
                qq  +=  1



    def WriteQVectorsFile( self ):
        """
        writing the Qvectors in DSLEAP formtat
        to QVectors.in
        """
        with open( "QVectors.in" , "w" ) as outfile:
            outfile.write( "  " + str( self.QVectors.shape[0] ) + "\n" )
            for qq in range( self.QVectors.shape[0] ):
                outfile.write( ' '.join( [ "{:15.8f}".format( x ) 
                                       for x in self.QVectors[ qq , : ] ] ) + "\n" )



        




                     
if __name__=="__main__":
    options = getOptions( sys.argv[1:] )
    x = MakeEigenVectorFile( options.basisVectorfile )
    x.WriteEigenVectorsDSLEAPFormat()
    x.WriteQVectorsFile()
