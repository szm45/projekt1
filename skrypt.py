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

    def odwrotny_hirvonen(self, f, l, h):
       '''
      Algorytm odwrotny do algorytmu Hirvonena służy do przekształcania współrzędnych geodezyjnych (B, L, H) 
      na współrzędne ortokartezjańskie (x, y, z). To proces umożliwiający przejście z opisu geodezyjnego punktu 
      na powierzchni ziemi do odpowiadającej mu lokalizacji w trójwymiarowym układzie współrzędnych kartezjańskich.
       Parametry
       ----------
       f : FLOAT
           [st. dz] - szerokość
       l : FLOAT
           [st. dz] - długośc 
       h : FLOAT
           [m] - wysokość 
       Returns
       -------
        X, Y, Z : FLOAT
             [m] - współrzędne orto-kartezjańskie

       '''
       f=np.radians(f)
       l=np.radians(l)
       
       N = Transformacje.NP(self, f)
       X = (N + h) * np.cos(f) * np.cos(l)
       Y = (N + h) * np.cos(f) * np.sin(l)
       Z = (N *(1-self.e2) + h) * np.sin(f) 
       
       return(X,Y,Z)
   
   
    def flh2PL1992(self, f, l):
       '''
       Układ współrzędnych 1992 (PUWG-92) to system płaskich współrzędnych prostokątnych,
       który używa odwzorowania Gaussa-Krügera dla elipsoidy GRS80 w ramach pojedynczej dziesięciostopniowej strefy.

       Parametry
       ----------
       f : FLOAT
           [st. dz.] - szerokość
       l : FLOAT
           [st. dz] - długośc

       Returns
       -------
        X1992, Y1992 : FLOAT
             [m] - współrzędne (1992)

       '''
       
       if l > 25.5 or l < 13.5:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(l))} ten południk nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
           
       if f > 55 or f < 48.9:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(f))} ten równoleżnik nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
           
       f = np.radians(f)
       l = np.radians(l)
       a2 = self.a**2
       b2 = a2 * (1 - self.e2)
       e_2 = (a2 - b2)/b2
       l0 = np.radians(19)
       dl = l - l0
       dl2 = dl**2
       dl4 = dl**4
       t = np.tan(f)
       t2 = t**2
       t4 = t**4
       n2 = e_2 * (np.cos(f)**2)
       n4 = n2 ** 2
       N = Transformacje.NP(self, f)
       e4 = self.e2**2
       e6 = self.e2**3
       A0 = 1 - (self.e2/4) - ((3*e4)/64) - ((5*e6)/256)
       A2 = (3/8) * (self.e2 + e4/4 + (15*e6)/128)
       A4 = (15/256) * (e4 + (3*e6)/4)
       A6 = (35*e6)/3072
       sigma = self.a * ((A0 * f) - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
       xgk = sigma + ((dl**2)/2) * N * np.sin(f) * np.cos(f) * (1 + ((dl**2)/12)*(np.cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
       ygk = dl * N * np.cos(f) * (1 + (dl2/6) * (np.cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
       x92 = xgk * 0.9993 - 5300000
       y92 = ygk * 0.9993 + 500000
       return(x92,y92)
   
   
    def flh2PL2000(self, f, l):
       '''
      Układ współrzędnych 2000 to system prostych współrzędnych płaskich.
      Wykorzystuje on odwzorowanie Gaussa-Krügera dla elipsoidy GRS 80 w czterech określonych strefach, na południkach
      15°E, 18°E, 21°E i 24°E.

       Parametry
       ----------
       f : FLOAT
           [st. dz.] - szerokość 
       l : FLOAT
           [st. dz.] - długośc 

       Returns
       -------
        X2000, Y2000 : FLOAT
             [m] - współrzędne 

       '''
         
       if l >= 13.5 and l < 16.5:
           l0 = np.radians(15)
       elif l >= 16.5 and l < 19.5:
           l0 = np.radians(18)
       elif l >= 19.5 and l < 22.5:
           l0 = np.radians(21)
       elif l >= 22.5 and l <= 25.5:
           l0 = np.radians(24)
       else:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(l))} ten południk nie mieści się w zakresie")
       
       if f > 55 or f < 48.9:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(f))} ten równoleżnik nie mieści się w zakresie")
           
       f = np.radians(f)
       l = np.radians(l)
       a2 = self.a**2
       b2 = a2 * (1 - self.e2)
       e_2 = (a2 - b2)/b2
       dl = l - l0
       dl2 = dl**2
       dl4 = dl**4
       t = np.tan(f)
       t2 = t**2
       t4 = t**4
       n2 = e_2 * (np.cos(f)**2)
       n4 = n2 ** 2
       N = Transformacje.NP(self, f)
       e4 = self.e2**2
       e6 = self.e2**3
       A0 = 1 - (self.e2/4) - ((3*e4)/64) - ((5*e6)/256)
       A2 = (3/8) * (self.e2 + e4/4 + (15*e6)/128)
       A4 = (15/256) * (e4 + (3*e6)/4)
       A6 = (35*e6)/3072
       sigma = self.a * ((A0 * f) - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
       xgk = sigma + ((dl**2)/2) * N * np.sin(f) * np.cos(f) * (1 + ((dl**2)/12)*(np.cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
       ygk = dl * N * np.cos(f) * (1 + (dl2/6) * (np.cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
       strefa = round(l0 * 180/np.pi)/3
       x00 = xgk * 0.999923
       y00 = ygk * 0.999923 + strefa * 1000000 + 500000
       return(x00,y00)
   
    def dXYZ(self, xa, ya, za, xb, yb, zb):
        '''
       Funkcja oblicza różnicę współrzędnych między punktami A i B, należy jej użyć do obliczeń macierzy NEU.

        Parametry
        ----------
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne orto-kartezjańskie, 

        Returns
        -------
        dXYZ : ARRAY
            macierz różnicy współrzędnych

        '''
        dXYZ = np.array([xb-xa, yb-ya, zb-za])
        return(dXYZ)
    
    
    def rneu(self, f, l):
        '''
        Funkcja generuje macierz obrotu R, niezbędna jest ona do stworzenia macierzy NEU.

        Parametry
        ----------
        f : FLOAT
            [st. dz.] - szerokość 
        l : FLOAT
            [st. dz.] - długośc 

        Returns
        -------
        R ARRAY
            macierz obrotu R
             
        '''
        f=np.radians(f)
        l=np.radians(l)
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                      [-np.sin(f)*np.sin(l),  np.cos(l), np.cos(f)*np.sin(l)],
                      [np.cos(f),             0,         np.sin(f)          ]])
        return(R)
    
    
    def xyz2neu(self, f, l, xa, ya, za, xb, yb, zb):
        '''
        Układ współrzędnych horyzontalnych opisuje położenie obiektów astronomicznych względem
        lokalnego horyzontu. Zenit i nadir, są ważnymi punktami w tym układzie. Współrzędne horyzontalne
        zmieniają się wraz z ruchem obserwatora i czasem, co pozwala określić chwilową pozycję obiektów na niebie.

        Parametry
        ----------
        f : FLOAT
            [st. dz.] - szerokość 
        l : FLOAT
            [st. dz.] - długośc
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne orto-kartezjańskie

        Returns
        -------
        n , e, u : STR
            wsp. horyz.
            

        '''
        dX = Transformacje.dXYZ(self, xa, ya, za, xb, yb, zb)
        R = Transformacje.rneu(self, f,l)
        neu = R.T @ dX
        n = neu[0];   e = neu[1];   u = neu[2]
        n = "%.16f"%n; e = "%.16f"%e; u="%.16f"%u
        dlugosc = []
        xx = len(n); dlugosc.append(xx)
        yy = len(e); dlugosc.append(yy)
        zz = len(u); dlugosc.append(zz)
        P = 50
        
        while xx < P:
            n = str(" ") + n
            xx += 1
        
        while yy < P:
            e = str(" ") + e
            yy += 1
            
        while zz < P:
            u = str(" ") + u
            zz +=1
            
        return(n, e, u)
    
    eoeoeoeoeoeoeo