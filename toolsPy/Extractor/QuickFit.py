#!/usr/bin/env python3


import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy.stats import skewnorm
from scipy.special import erf
from scipy.signal import savgol_filter
from scipy.integrate import simps

from scipy.signal import find_peaks
from scipy.signal import peak_widths




__author__ = "Jonathan Lahnsteiner"
__copyright__ = ""
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jonathan Lahnsteiner"
__email__ = "jonathan.lahnsteiner@gmx.at"
__status__ = "Production"







def getOptions( args = sys.argv[1:] ):
    """
    parse user input when the file is executed as main file
    """
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument( "-f" , "--SignalFile", 
             help = "Supply filename of the file you want to fit" )
    parser.add_argument( "-c" , "--column", 
             help = "Supply the coulmn of the from the file you want to fit" , type = int )
    options = parser.parse_args( args )
    return options







#### part where functions that can be used for fitting are defined


def Lorentz( x , gamma , I  , x0 ):
    """
    Lorentz distribution 
    gamma controls width
    I controls height
    x0 is the position of the Lorentzian
    """
    gamma2  =  gamma * gamma
    return I * ( gamma2 / ( ( x - x0 )**2 + gamma2 ) )



def MaxwellBoltz( x , K1 , a ):
    """
    compute maxwell boltzmann distribution
    K1 determines the height of the distribution
    and a is the variance parameter
    """
    x2  =  x*x
    return K1 * np.sqrt( 2 / np.pi ) * np.exp( -x2 / ( 2.0 * a * a ) ) / ( a*a*a )



def Gaussian( x , K1 , K2 , x0 ):
    """
    use Gaussian for fitting where the Gaussian
    is non normalized. K1 determines height
    K2 the width of the Gaussian and x0 it's position
    """
    x2 = ( x - x0 )**2
    return K1 * np.exp( -K2 * x2 )



def SkewNormal( x , K , sigma , alpha , x0 ):
    """
    K  -> height of distribution
    sigma is the width of the distribution
    alpha is the skewness of the distribution
    x0 is the position of the distribution
    """
    xp = ( x - x0 ) / sigma   # convert to normal distribution
    N = np.pi * sigma
    N = 1.0 / N
    return K * N * np.exp( -0.5 * xp**2 ) * ( 1.0 + erf( alpha*xp / np.sqrt(2) ) )



def RayleighDist( x , sigma ): 
    """
    return Raleigh distribution
    with sigma as scale parameter
    """
    sigma2 = sigma * sigma             
    xp = x / sigma2                    
    return xp * np.exp(- 0.5 * xp**2 )



def LinearFunc( x , k , d ):
    return k * x + d



def ComputeMSD( x , y ):
    """
    compute mean square displacement between
    two data sets
    """
    return np.sum( ( x - y )**2 ) / x.shape[0]







class MakeQuickFit:
    def __init__( self , fname , col , tresh = 0.85 ):
        """
        supply a filename a or a data array
        as first argument
        the first column is assumed to be 
        the x-axis, the remaining columns are
        the signal
        supply col to determine the column to be fittedin your array
        tresh is the treshold value when self.MakeBestFit is used 
        """
        if ( type( fname ) == str ):
            data   =   np.loadtxt( fname )
        else:
            # assuming data was supplied
            data  = fname
        self.N =  1000 #data.shape[0]
        self.data = data[ 0:self.N , : ]
        self.col = col
        self.tresh = tresh

        self.Fits        =  np.zeros( [ self.data.shape[0] , 5 ] )
        self.Correlation =  np.zeros( 5 )
        self.FitParams   =  np.zeros( [ 5 , 2 ] )
        self.labels = [ 'Lorentz' , 'Gauss' , 'MaxwellBoltzmann' , 'SkewNorm' , 'Rayleigh' ]
        #self.CorrectForBaseLine()




    def CorrectForBaseLine( self ):
        """
        function can be used when signal contains a linear drift
        linear function will be fitted to the signal. This 
        fitted linear function is substracted from the signal
        as a base line correction
        """
        popt , pcov = curve_fit( LinearFunc , self.data[ : , 0 ] ,
                                             self.data[ : , self.col ] )
        y = LinearFunc( self.data[ : , 0 ] , popt[0] , popt[1] )
        self.data[ : , self.col ]   -=   y





    def MakeFitLorentz( self ):
        """
        Fit a Lorentzian function to your data sets
        the function will determine the fit parameters
        and the correlation according to pearsonr

        the fitting parameters are stored in self.FitParams[ : , 1 ]
        and the pearsonr is returned by the function

        function does not return the Lorentz Gamma parameter but
        the peaks FHWM

        the data is stored in the order
        frequency , FWHM
        """
        idx  =  np.argmax( self.data[ : , self.col ] )
        x0   =  self.data[ idx , 0 ]
        I0   =  self.data[ idx , self.col ]
        dx   =  self.data[ 1 , 0 ] - self.data[ 0 , 0 ]
        try:
            popt , pcov = curve_fit( Lorentz , self.data[ : , 0 ] ,
                                           self.data[ : , self.col ] ,
                                           p0 = [ 1.0 , I0 , x0 ] )
            self.Fits[ : , 0 ]    =  Lorentz( self.data[ : , 0 ] , popt[0] , popt[1] , popt[2] )
            idxes  =  find_peaks( self.Fits[ : , 0 ] )[0]
            width  =  peak_widths( self.Fits[ : , 0 ] , idxes , rel_height = 0.5 )[0][0]
            MSD = ComputeMSD( self.Fits[ : , 0 ] , self.data[ : , self.col ] )
        except:
            print( "Fit Parameters in Lorentz could not be found" )
            popt = [ 0 , 0 , 0 ]
            width = 0
            self.Fits[ : , 0 ] = 0

        self.Correlation[ 0 ] =  pearsonr( self.data[ : , self.col ] , self.Fits[ : , 0 ] )[0]**2
        self.FitParams[ 0 , 0 ]  =  popt[-1]
        self.FitParams[ 0 , 1 ]  =  width * dx
        print( "Lorentz correlation" , self.Correlation[0] )
        return self.Correlation[0]




    def MakeFitGauss( self ):
        """
        Fitting a Gaussian function to the supplied signal

        the data is stored in the self.FitParams[ 1 , : ]
        order of the data is again
        frequency , FWHM 
        """
        idx  =  np.argmax( self.data[ : , self.col ] )
        x0   =  self.data[ idx , 0 ]
        I0   =  self.data[ idx , self.col ]
        dx   =  self.data[ 1 , 0 ] - self.data[ 0 , 0 ]
        try:
            popt , pcov = curve_fit( Gaussian , self.data[ : , 0 ] ,
                                            self.data[ : , self.col ] ,
                                            p0 = [ 1.0 , I0 , x0 ] )
            self.Fits[ : , 1 ]    =  Gaussian( self.data[ : , 0 ] , popt[0] , popt[1] , popt[2] )
            idxes = find_peaks( self.Fits[ : , 0 ] )[0]
            width = peak_widths( self.Fits[ : , 0 ] , idxes , rel_height = 0.5 )[0][0]
        except:
            print( "Fit Parameters in Gaussian could not be found" )
            popt = [ 0 , 0 , 0 ]
            width = 0
            self.Fits[ : , 1 ]    =  0
        self.Correlation[ 1 ] =  pearsonr( self.data[ : , self.col ] , self.Fits[ : , 1 ] )[0]**2
        self.FitParams[ 1 , : ] = popt[-1]
        #determine width
        self.FitParams[ 1 , 1 ]  =  width * dx
        print( "Gaussian correlation" , self.Correlation[1] )




    def MakeFitMB( self ):
        """
        fit Maxwell-Boltzman distribution
        to the selected signal
        
        the data is stored in the self.FitParams[ 2 , : ]
        order of the data is again
        frequency , FWHM 
        """
        idx  =  np.argmax( self.data[ : , self.col ] )
        x0   =  self.data[ idx , 0 ]
        I0   =  self.data[ idx , self.col ]
        dx   =  self.data[ 1 , 0 ] - self.data[ 0 , 0 ]
        try:
            popt , pcov = curve_fit( MaxwellBoltz , self.data[ : , 0 ] ,
                                           self.data[ : , self.col ] ,
                                           p0 = [ I0 , 0.1 ] )
            self.Fits[ : , 2 ]    =  MaxwellBoltz( self.data[ : , 0 ] , popt[0] , popt[1] )
            #determine width
            idxes = find_peaks( self.Fits[ : , 2 ] )[0]
            width = peak_widths( self.Fits[ : , 2 ] , idxes , rel_height = 0.5 )[0][0]
        except:
            print( "Fit Parameters in Maxwell Boltzmann could not be found" )
            popt = [ 100000000 , 100000000 ]
            self.Fits[ : , 2 ]    = 0
            width = 0
        self.Correlation[ 2 ] =  pearsonr( self.data[ : , self.col ] , self.Fits[ : , 2 ] )[0]**2
        self.FitParams[ 2 , : ] = 2.0 * popt[-1] * np.sqrt( 2.0 / np.pi )
        self.FitParams[ 2 , 1 ]  =  width * dx
        print( "Maxwell Boltzmann correlation" , self.Correlation[2] )




    def MakeFitSkew( self ):
        """
        fit Skewed -Gaussian distribution
        to the selected signal
        
        the data is stored in the self.FitParams[ 3 , : ]
        order of the data is again
        frequency , FWHM 
        """
        idx  =  np.argmax( self.data[ : , self.col ] )
        x0   =  self.data[ idx , 0 ]
        I0   =  self.data[ idx , self.col ]
        dx   =  self.data[ 1 , 0 ] - self.data[ 0 , 0 ]
        try:
            popt , pcov = curve_fit( SkewNormal , self.data[ : , 0 ] ,
                                           self.data[ : , self.col ] ,
                                           p0 = [ I0 , 1.0 , 1.0 , x0 ] )
            self.Fits[ : , 3 ]    =  SkewNormal( self.data[ : , 0 ] , popt[0] , popt[1] , popt[2] , popt[3] )
            #determine width
            idxes = find_peaks( self.Fits[ : , 3 ] )[0]
            width = peak_widths( self.Fits[ : , 3 ] , idxes , rel_height = 0.5 )[0][0]
        except:
            print( "Fit Parameters in SkewGauss could not be found" )
            popt = [ 0 , 1.0 , 0 , 0 ]
            width = 0
        self.Correlation[ 3 ] =  pearsonr( self.data[ : , self.col ] , self.Fits[ : , 3 ] )[0]**2
        delta  =  popt[-2] / np.sqrt( 1+popt[-2]**2 )
        self.FitParams[ 3 , 0 ] = popt[-1] + popt[1] * delta * np.sqrt( 2.0 / np.pi )
        self.FitParams[ 3 , 1 ]  =  width * dx
        print( "Skew Normal correlation" , self.Correlation[3] )



    def MakeFitRayleigh( self ):
        """
        fit Raleigh distribution
        to the selected signal
        
        the data is stored in the self.FitParams[ 4 , : ]
        order of the data is again
        frequency , FWHM 
        """
        idx  =  np.argmax( self.data[ : , self.col ] )
        x0   =  self.data[ idx , 0 ]
        I0   =  self.data[ idx , self.col ]
        dx   =  self.data[ 1 , 0 ] - self.data[ 0 , 0 ]
        try:
            popt , pcov = curve_fit( RayleighDist , self.data[ : , 0 ] ,
                                           self.data[ : , self.col ] ,
                                           p0 = [ x0 ] )
            self.Fits[ : , 4 ]      =   RayleighDist( self.data[ : , 0 ] , popt[0] )
            #determine width
            idxes = find_peaks( self.Fits[ : , 4 ] )[0]
            width = peak_widths( self.Fits[ : , 4 ] , idxes , rel_height = 0.5 )[0][0]
        except:
            print( "Fit Parameters in Rayleigh could not be found" )
            popt   =  [ 0 ]
            width  =  0.0


        self.Correlation[ 4 ]   =   pearsonr( self.data[ : , self.col ] , self.Fits[ : , 4 ] )[0]**2
        self.FitParams[ 4 , : ] =   popt[0] * np.sqrt( np.pi / 2.0 )
        self.FitParams[ 4 , 1 ]  =  width* dx

        print( "Rayleigh correlation" , self.Correlation[4] )



    def MakeFit( self ):
        """
        fitting all supported distributions
        the fit parameters are found in self.FitParams
        the order of the parameters is
        as the functions are called below
        self.FitParams[ 0 , : ]  ->  Lorentz
        self.FitParams[ 1 , : ]  ->  Gaussian
        self.FitParams[ 2 , : ]  ->  MaxwellBoltzman
        .
        .
        .
        .

        """
        self.MakeFitLorentz()
        self.MakeFitGauss()
        self.MakeFitMB()
        self.MakeFitSkew()
        self.MakeFitRayleigh()
        self.PlotData()



    def LorentzOnly( self ):
        """
        Fit Lorentzian and return FitParams and 
        fit quality
        """
        pcov  =  self.MakeFitLorentz()
        return self.FitParams[ 0 , : ] , pcov

    def GaussOnly( self ):
        """
        Fit Gauss and return FitParams and 
        fit quality
        """
        self.MakeFitGauss()
        return self.FitParams[ 1 , : ]

    def SkewNormOnly( self ):
        """
        Fit Skewed Gaussian and return FitParams and 
        fit quality
        """
        self.MakeFitSkew()
        return self.FitParams[ 3 , : ]
    

    def MBOnly( self ):
        """
        Fit MaxwellBoltzmann and return FitParams and 
        fit quality
        """
        self.MakeFitMB()
        return self.FitParams[ 2 , : ]
    
    def RayleighOnly( self ):
        """
        Fit Raleigh distribution and return FitParams and 
        fit quality
        """
        self.MakeFitSkew()
        return self.FitParams[ 3 , : ]



    def PlotData( self ):
        """
        plot original data 
        the data smoothed with a Savitzky-Golay filter
        and show the Fitted signals obtained from the
        different distributions
        """
        plt.title( 'q = ' + str( self.col ) )
        plt.plot( self.data[ : , 0 ] , self.data[ : , self.col ] , label = 'orig' )
        func = savgol_filter( self.data[ : , self.col ] , 51 , 3 )
        plt.plot( self.data[ : , 0 ] , func , label = 'smooth' )
        for i in range( self.Fits.shape[1] ):
            plt.plot( self.data[ : , 0 ] , self.Fits[ : , i ] , label = self.labels[ i ] )


        plt.scatter( self.FitParams[ self.indx , 0 ] , np.max( self.Fits ) )

    
        plt.legend()
        plt.show()



    def MakeBestFit( self , plot = False ):
        """
        Fit all supported distributions
        and check wich one describes the data best
        the routine prints the best distribution
        to the secreen 
        and shows a plot when plot = True
        """
        self.MakeFitLorentz()
        self.MakeFitGauss()
        self.MakeFitMB()
        self.MakeFitSkew()
        #self.MakeFitRayleigh()
        for i in range( self.Correlation.shape[0] ):
            if ( np.isnan( self.Correlation[i] ) ):
                self.Correlation[i] = 0.0
        self.indx = np.argmax( self.Correlation )
        indx = self.indx
        print( 'Best description found by ' , self.labels[ indx ] )
        print( 'Fit Parameters ' , self.FitParams[ indx , : ] )
        print( 'Correlation cofficient r^2' , self.Correlation[ indx ] )
        if ( plot ):
            self.PlotData()
        if ( self.Correlation[ indx ] < self.tresh ):
            self.FitParams[ indx , : ] = [ -10 ] * len( self.FitParams[ indx , : ] )
        return self.FitParams[ indx , : ] , self.Correlation[ indx ]





if __name__=="__main__":
    options  =  getOptions( sys.argv[1:] )
    x = MakeQuickFit( options.SignalFile , options.column )

    ## if this file is directly executed the BestFit will be called
    x.MakeBestFit( True )
