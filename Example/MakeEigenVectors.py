#!/usr/bin/env python3


import sys
import numpy as np




class CreateEigenVectorsExample:
    def __init__( self ):
        """
        class creates a BasisVector.in file
        for the Morse potential exmample
        self.NQ        -> number of Qpoints
        self.Nbranches -> number of branches 
                          of phonon dispersion
        self.NatomsUC  -> number of atoms in the
                          unit cell
        self.dim       -> space dimensionality
        """
        self.NQ = 16
        self.Nbranches = 3
        self.NatomsUC  = 1
        self.dim       = 3
        self.QVectors = np.zeros( ( self.NQ , 3 ) )
        self.MakeQVectors()
        self.EigenVectors = np.zeros( [ self.NQ , 
                                        self.Nbranches ,
                                        self.NatomsUC , 
                                        self.dim ] )
        self.MakeEigenVectors()



    def MakeQVectors( self ):
        """
        Making Q vectors along Gamma - X
        for example Morse system
        """
        N = self.QVectors.shape[0]
        for i in range( N//2 + 1 ):
            self.QVectors[ i , 0 ] = i

        for i in range( 1 , N//2 ):
            self.QVectors[ N - i, 0 ] = i


    def MakeEigenVectors( self ):
        """
        making eigenvectors for the defined qpath
        """ 
        sqrt2 = np.sqrt(2)
        Isqrt2 =  1.0 / sqrt2
        EVectors = np.asarray( [ [ Isqrt2  , Isqrt2 , 0 ] ,
                                 [ Isqrt2  ,-Isqrt2 , 0 ] , 
                                 [ 0 , 0 , 1 ] ] )
        for i in range( self.NQ ):
            for j in range( self.Nbranches ):
                self.EigenVectors[ i , j , 0 , : ] =  EVectors[ j , : ]



    def WriteFile( self ):
        """
        this routine writes the phonon eigenvectors
        in the BasisVector.in file format.
        Not after every branch there is an empty line
        and after every qpoint is an empty line
        """
        with open( "BasisVector.in" , "w" ) as outfile:
            firstLine = "  " + str( self.NQ ) + \
                        "  " + str( self.Nbranches ) + \
                        "  " + str( self.NatomsUC ) + \
                        "  " + str( self.dim ) + "\n"
            outfile.write( firstLine )
            for qq in range( self.NQ ):  ## loop over Q vectors
                lineQ = [ "{:15.8f}".format( x ) for x in 
                          self.QVectors[ qq , : ] ]
                lineQ = "".join( lineQ )
                outfile.write( lineQ + "\n" )
                for branch in range( self.Nbranches ):  ## loop over branches
                    for atom in range( self.NatomsUC ): ## loop over atoms in unit cell
                        line = [ "{:15.8f}".format( x ) for x in 
                                     self.EigenVectors[ qq , branch , atom , : ] ]
                        line = "".join( line )
                        outfile.write( line + "\n" )
                    outfile.write( "\n" )
                outfile.write( "\n" )

























if __name__=="__main__":
    EigenVectors  =  CreateEigenVectorsExample()
    EigenVectors.WriteFile()

