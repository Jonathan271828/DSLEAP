#!/usr/bin/env python3



import sys



class CreateBoxListMonoAtom:
    """
    object creates the BoxList.in file
    for the DSLEAP code example
    """
    def __init__( self , N ):
        """
        N is the total number of atoms
        in the system. 
        there is a single atom
        """
        self.N = N
        self.MakeList()


    def MakeList( self ):
        """
        generates the box list for the
        supplied example
        """
        self.BoxList = []
        for i in range( 1 , self.N + 1 ):
            self.BoxList.append( i )



    def WriteFile( self ):
        """
        In every unit cell there is only a single atom
        therefore there is a single number per line
        """
         with open( "BoxList.in" , "w" ) as outfile:
             for i in range( self.N ):
                 outfile.write( "  " + str( i ) + "\n" )
















if __name__=="__main__":
    MakeBoxList = CreateBoxListMonoAtom( int( sys.argv[1] ) )
    MakeBoxList.WriteFile()
