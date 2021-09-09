


import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import yaml
from matplotlib.widgets import Slider




__author__ = "Jonathan Lahnsteiner"           
__copyright__ = ""                            
__license__ = "GPL"                           
__version__ = "1.0.0"                         
__maintainer__ = "Jonathan Lahnsteiner"       
__email__ = "jonathan.lahnsteiner@gmx.at"     
__status__ = "Production"                     






def getOptions( args = sys.argv[1:] ):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "-df" , "--dynfile", help = "Your input file of q resolved power spectrum" )
    parser.add_argument( "-T" , "--Temperature" , type = float , help = 
                         "Insert temperature of MD if reweighing by Planck distribution is desired " , default = None )
    options = parser.parse_args( args )
    return options




class GenBandPlot:
    def __init__( self , fname1 , T = None ):
        """
        Input:
              fname1 -> file containing the power spectrum
              the file structure is assumed to be 
              frequency axis | qpoint-1 | qpoint-2 | .....

              T -> temperature at which the power spectrum was sampled
                   if the temperature is supplied the ppower spectrum will
                   be placnk reweighed

        The module can interactively visulaize data computed with the DSLEAP
        code
        """
        self.PowerSpecFile  = fname1
        self.PowerSpectrum = np.loadtxt( self.PowerSpecFile )
        self.T = T
        self.kb =  0.02083671182786661
        
        self.ReweighPowerSpectrum()
        


    def ReweighPowerSpectrum( self ):
        """
        if a temperature self.T was supplied the power spectrum will be
        reweighed by Planck distribution in the classical limit
        """
        if ( self.T != None ):
            self.PowerSpectrum[ : , 1: ] = np.transpose( 
                    np.asarray( [ self.PowerSpectrum[ : , i ] * 
                    self.PowerSpectrum[ : , 0 ] / ( self.T * self.kb ) 
                    for i in range( 1 , self.PowerSpectrum.shape[1] ) ] ) )



    def PlotDSLEAPdata( self , Nbands = None ):
        """
        Plots a colormap of the power spectrum
        computed with the DSLEAP 
        """
        xmin = 0
        xmax = self.PowerSpectrum.shape[ 1 ] 
        ymin = self.PowerSpectrum[ 0 , 0 ]
        ymax = self.PowerSpectrum[ -1 , 0 ]    #here i have to take min and max of power spec otherwise it is wrongly scaled
        delta = ( xmax - xmin ) / ( self.PowerSpectrum.shape[1]-1 )
        fig , ax = plt.subplots()
        plt.subplots_adjust( left=0.25, bottom=0.25 )
        #define size of sliders in the plot
        minax =  plt.axes( [0.2,0.15,0.65,0.03] )
        maxax =  plt.axes( [0.2,0.1,0.65,0.03] )
        # define sliders with min max dstep values
        maxi  =  np.nanmax( np.log( self.PowerSpectrum[ : , 1: ] ) )
        mini  =  np.nanmin( np.log( self.PowerSpectrum[ : , 1: ] ) )

        mean  =  np.mean( [ maxi , mini ] )
        var   =  np.var( [ maxi , mini ] )
        initA =  mean - var
        initB =  mean + var 

        smin  =  Slider( minax , "Min" , mini , maxi , valinit = initA , valstep = 0.02 )
        smax  =  Slider( maxax , "Max" , mini , maxi , valinit = initB , valstep = 0.02 )
        #plot graph
        im = ax.imshow( np.log( self.PowerSpectrum[ ::-1 , 1: ] ) , interpolation = "nearest" , 
                        aspect="auto" , extent=(xmin,xmax,ymin,ymax) , cmap="gnuplot2" ,vmin=-5 , vmax=0 ) 
        #update values adjusted by sliders
        def update( val ):
            a = smin.val
            b = smax.val
            im.set_clim( vmin = a , vmax = b )  ## reset clim according to slider movement

        ## update stuff here
        fig.colorbar( im , ax=ax )
        smin.on_changed( update )
        smax.on_changed( update )
        ax.set_ylabel( "Frequency [inverse time unit]" )
        ax.set_xlabel( "qpath" )
        plt.show()


if __name__=="__main__":
    options = getOptions( sys.argv[1:] )
    x = GenBandPlot( options.dynfile , options.Temperature )
    x.PlotDSLEAPdata( )
