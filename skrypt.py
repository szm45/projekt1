# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 17:13:57 2024
"""

import numpy as np
from argparse import ArgumentParser
import os as os
import tkinter as tk


class Transformacje:
    
    
    def __init__(self, model: str="WGS84"):
        """
        Parametry elipsoidy definiujące jej kształt obejmują:
            a- długość głównej półosi elipsoidy, czyli jej równikowy promień,
            e2-kwadrat mimośrodu, obliczany jako ((kwadrat równikowego promienia + kwadrat biegunowego promienia) / kwadrat równikowego promienia).
        """
        
        if model == "WGS84":
            self.a = 6378137.000
            self.e2 = 0.00669437999013
        elif model == "GRS80":
            self.a = 6378137.000
            self.e2 = 0.00669438002290
        elif model == "Krasowski":
            self.a = 6378245.000
            self.e2 = 0.00669342162296
        else:
            raise NotImplementedError(f"{model} zła elipsoida - możliwe elipsoidy WGS84, GRS80, Krasowski.")

    def dms(self, x):
        '''
        Funkcja dms zamienienia radiany na stopnie minuty i sekundy.   

        Parametry
        ----------
        x : FLOAT
            [radiany].

        Returns
        x : STR
            [dms] - stopnie, minuty, sekundy
        '''
        sig = ' '
        if x < 0:
            sig = '-'
            x = abs(x)
        x = x * 180/np.pi
        d = int(x)
        m = int(60 * (x - d))
        s = (x - d - m/60)*3600
        if s > 59.999995:
            s = 0
            m = m + 1
        if m == 60:
            m = 0
            d = d + 1
        
        d = str(d)
        if len(d) == 1:
            d = "  " + d
        elif len(d) == 2:
            d = " " + d
        elif len(d) == 3:
            d = d
            
        if m < 10:
            m = str(m)
            m = "0" + m
            
        if s < 10:
            s = "%.5f"%s
            s = str(s)
            s= "0" + s
        else:
            s = ("%.5f"%s)
            
        x1=(f'{d}°{m}′{s}″')  
        return(x1)
        
        
    def NP(self, f):
        '''
       Funkcja oblicza promień przekroju w pierwszej pionowej płaszczyźnie, który jest potrzebny m.in. do zastosowania algorytmu Hirvonena
        
        Parameters
        ----------
        f : FLOAT
            [radiany] - szerokość geodezyjna

        Returns
        -------
        N : float
            [metry] - promień przekroju w pierwszym wertykale

        '''
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        return(N)
    
    
    def algorytm_hirvonena(self, X, Y, Z, output="dec_degree"):
        '''
        Algorytm Hirvonena służy do konwersji współrzędnych prostokątnych (x, y, z) na współrzędne geodezyjne (fi, Lambda, h). Jest to proces iteracyjny, który umożliwia osiągnięcie dokładności rzędu 1 milimetra po kilku powtórzeniach procedury.

         Parametry
         ----------
         X, Y, Z : FLOAT
              współrzędne w układzie orto-kartezjańskim, 

         Returns
         -------
         fi : FLOAT
         [st. dz.] - szerokość geodezyjna.
         lam : FLOAT
         [st. dz.] - długośc geodezyjna.
         h : FLOAT
         [m] - wysokość elipsoidalna
         output [STR] - opcjonalne, domylne 
         -dec_degree - st. dz
         -dms - st, min, sek
         -radiany - rad 
         """
        '''

        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z/(p * (1 - self.e2)))
        while True:
            N = Transformacje.NP(self, f)
            h = (p / np.cos(f)) - N
            fp = f
            f = np.arctan(Z / (p * (1 - self.e2 * (N / (N+h)))))
            if np.abs(fp - f) < (0.000001/206265):
                break
        l = np.arctan2(Y, X)
        if output == "dec_degree":
            fi=(f*180/np.pi)
            lam=(l*180/np.pi)
            return (fi, lam, h)
        elif output == "dms":
            fi = Transformacje.dms(self, f)
            lam = Transformacje.dms(self, l)
            return (fi,lam,h) 
        elif output == 'radiany':
            fi=f
            lam=l
            return(f,l,h)
        else:
            raise NotImplementedError(f"{output} - output format not defined")
