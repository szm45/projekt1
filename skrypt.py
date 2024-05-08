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
        
        
    def NP(self, fi):
        '''
       Funkcja oblicza promień przekroju w pierwszej pionowej płaszczyźnie, który jest potrzebny m.in. do zastosowania algorytmu Hirvonena
        
        Parameters
        ----------
        fi : FLOAT
            [radiany] - szerokość geodezyjna

        Returns
        -------
        N : float
            [metry] - promień przekroju w pierwszym wertykale

        '''
        N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
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
        fi = np.arctan(Z/(p * (1 - self.e2)))
        while True:
            N = Transformacje.NP(self, fi)
            h = (p / np.cos(fi)) - N
            fip = fi
            fi = np.arctan(Z / (p * (1 - self.e2 * (N / (N+h)))))
            if np.abs(fip - fi) < (0.000001/206265):
                break
        lam = np.arctan2(Y, X)
        if output == "dec_degree":
            fi=(fi*180/np.pi)
            lam=(lam*180/np.pi)
            return (fi, lam, h)
        elif output == "dms":
            fi = Transformacje.dms(self, fi)
            lam = Transformacje.dms(self, lam)
            return (fi,lam,h) 
        elif output == 'radiany':
            
            
            return(fi,lam,h)
        else:
            raise NotImplementedError(f"{output} - output format not defined")

    def odwrotny_hirvonen(self, fi, lam, h):
       '''
      Algorytm odwrotny do algorytmu Hirvonena służy do przekształcania współrzędnych geodezyjnych (B, L, H) 
      na współrzędne ortokartezjańskie (x, y, z). To proces umożliwiający przejście z opisu geodezyjnego punktu 
      na powierzchni ziemi do odpowiadającej mu lokalizacji w trójwymiarowym układzie współrzędnych kartezjańskich.
       Parametry
       ----------
       fi : FLOAT
           [st. dz] - szerokość
       lam : FLOAT
           [st. dz] - długośc 
       h : FLOAT
           [m] - wysokość 
       Returns
       -------
        X, Y, Z : FLOAT
             [m] - współrzędne orto-kartezjańskie

       '''
       fi=np.radians(fi)
       lam=np.radians(lam)
       
       N = Transformacje.NP(self, fi)
       X = (N + h) * np.cos(fi) * np.cos(lam)
       Y = (N + h) * np.cos(fi) * np.sin(lam)
       Z = (N *(1-self.e2) + h) * np.sin(fi) 
       
       return(X,Y,Z)
   
   
    def flh2PL1992(self, fi, lam):
       '''
       Układ współrzędnych 1992 (PUWG-92) to system płaskich współrzędnych prostokątnych,
       który używa odwzorowania Gaussa-Krügera dla elipsoidy GRS80 w ramach pojedynczej dziesięciostopniowej strefy.

       Parametry
       ----------
       fi : FLOAT
           [st. dz.] - szerokość
       lam : FLOAT
           [st. dz] - długośc

       Returns
       -------
        X1992, Y1992 : FLOAT
             [m] - współrzędne (1992)

       '''
       
       if lam > 25.5 or lam < 13.5:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(lam))} ten południk nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
           
       if fi > 55 or fi < 48.9:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(fi))} ten równoleżnik nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
           
       fi = np.radians(fi)
       lam = np.radians(lam)
       a2 = self.a**2
       b2 = a2 * (1 - self.e2)
       e_2 = (a2 - b2)/b2
       l0 = np.radians(19)
       dl = lam - l0
       dl2 = dl**2
       dl4 = dl**4
       t = np.tan(fi)
       t2 = t**2
       t4 = t**4
       n2 = e_2 * (np.cos(fi)**2)
       n4 = n2 ** 2
       N = Transformacje.NP(self, fi)
       e4 = self.e2**2
       e6 = self.e2**3
       A0 = 1 - (self.e2/4) - ((3*e4)/64) - ((5*e6)/256)
       A2 = (3/8) * (self.e2 + e4/4 + (15*e6)/128)
       A4 = (15/256) * (e4 + (3*e6)/4)
       A6 = (35*e6)/3072
       sigma = self.a * ((A0 * fi) - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi))
       xgk = sigma + ((dl**2)/2) * N * np.sin(fi) * np.cos(fi) * (1 + ((dl**2)/12)*(np.cos(fi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(fi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
       ygk = dl * N * np.cos(fi) * (1 + (dl2/6) * (np.cos(fi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(fi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
       x92 = xgk * 0.9993 - 5300000
       y92 = ygk * 0.9993 + 500000
       return(x92,y92)
   
   
    def flh2PL2000(self, fi, lam):
       '''
      Układ współrzędnych 2000 to system prostych współrzędnych płaskich.
      Wykorzystuje on odwzorowanie Gaussa-Krügera dla elipsoidy GRS 80 w czterech określonych strefach, na południkach
      15°E, 18°E, 21°E i 24°E.

       Parametry
       ----------
       fi : FLOAT
           [st. dz.] - szerokość 
       lam : FLOAT
           [st. dz.] - długośc 

       Returns
       -------
        X2000, Y2000 : FLOAT
             [m] - współrzędne 

       '''
         
       if lam >= 13.5 and lam < 16.5:
           l0 = np.radians(15)
       elif lam >= 16.5 and lam < 19.5:
           l0 = np.radians(18)
       elif lam >= 19.5 and lam < 22.5:
           l0 = np.radians(21)
       elif lam >= 22.5 and lam <= 25.5:
           l0 = np.radians(24)
       else:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(lam))} ten południk nie mieści się w zakresie")
       
       if fi > 55 or fi < 48.9:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(fi))} ten równoleżnik nie mieści się w zakresie")
           
       fi = np.radians(fi)
       lam = np.radians(lam)
       a2 = self.a**2
       b2 = a2 * (1 - self.e2)
       e_2 = (a2 - b2)/b2
       dl = lam - l0
       dl2 = dl**2
       dl4 = dl**4
       t = np.tan(fi)
       t2 = t**2
       t4 = t**4
       n2 = e_2 * (np.cos(fi)**2)
       n4 = n2 ** 2
       N = Transformacje.NP(self, fi)
       e4 = self.e2**2
       e6 = self.e2**3
       A0 = 1 - (self.e2/4) - ((3*e4)/64) - ((5*e6)/256)
       A2 = (3/8) * (self.e2 + e4/4 + (15*e6)/128)
       A4 = (15/256) * (e4 + (3*e6)/4)
       A6 = (35*e6)/3072
       sigma = self.a * ((A0 * fi) - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi))
       xgk = sigma + ((dl**2)/2) * N * np.sin(fi) * np.cos(fi) * (1 + ((dl**2)/12)*(np.cos(fi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(fi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
       ygk = dl * N * np.cos(fi) * (1 + (dl2/6) * (np.cos(fi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(fi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
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
    
    
    def rneu(self, fi, lam):
        '''
        Funkcja generuje macierz obrotu R, niezbędna jest ona do stworzenia macierzy NEU.

        Parametry
        ----------
        fi : FLOAT
            [st. dz.] - szerokość 
        lam : FLOAT
            [st. dz.] - długośc 

        Returns
        -------
        R ARRAY
            macierz obrotu R
             
        '''
        fi=np.radians(fi)
        lam=np.radians(lam)
        R = np.array([[-np.sin(fi)*np.cos(lam), -np.sin(lam), np.cos(fi)*np.cos(lam)],
                      [-np.sin(fi)*np.sin(lam),  np.cos(lam), np.cos(fi)*np.sin(lam)],
                      [np.cos(fi),             0,         np.sin(fi)          ]])
        return(R)
    
    
    def xyz2neu(self, fi, lam, xa, ya, za, xb, yb, zb):
        '''
        Układ współrzędnych horyzontalnych opisuje położenie obiektów astronomicznych względem
        lokalnego horyzontu. Zenit i nadir, są ważnymi punktami w tym układzie. Współrzędne horyzontalne
        zmieniają się wraz z ruchem obserwatora i czasem, co pozwala określić chwilową pozycję obiektów na niebie.

        Parametry
        ----------
        fi : FLOAT
            [st. dz.] - szerokość 
        lam : FLOAT
            [st. dz.] - długośc
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne orto-kartezjańskie

        Returns
        -------
        n , e, u : STR
            wsp. horyz.
            

        '''
        dX = Transformacje.dXYZ(self, xa, ya, za, xb, yb, zb)
        R = Transformacje.rneu(self, fi,lam)
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
    
    def wczytanie(self, Dane):
        with open (Dane,"r") as plik:
            tip = np.genformtxt(plik, delimiter=".",dtype = '<U20', skip_header = 4)
            X=[]
            Y=[]
            Z=[]
            for i in tip:
                x=i[0]
                X.append(float(x))
                y=i[1]
                Y.append(float(y))
                z=i[2]
                Z.append(float(z))
            wielkosc = len(X)
        return(X, Y, Z, wielkosc)
    
def zapisanie(self, X, Y, Z, f, l, h, x92, y92, x00, y00, N, E, U, xyz_txt, neu_txt ): 
    '''
    funkcja zapisuje wyniki obliczeń (x, y, z, f, l, h, x92, y92, x1992, y1992, x2000, y2000 ,neu).
    Tworzy z nich tabele.

    Parametry
    ----------
    X, Y, Z : LIST
         [metry] - współrzędne w układzie orto-kartezjańskim, 
     f : LIST
         [dms] - szerokość geodezyjna..
     l : LIST
         [dms] - długośc geodezyjna.
     h : LIST
         [metry] - wysokość elipsoidalna
    X1992, Y1992 : LIST
         [metry] - współrzędne w układzie 1992
     X2000, Y2000 : LIST
         [metry] - współrzędne w układzie 2000
    n, e, u : str
        współrzędne horyzontalne

    Returns
    -------
    PLIK TXT

    '''
    for i in range(len(X)):
        X[i] = Transformacje.zamiana_float2string(self, X[i])
        Y[i] = Transformacje.zamiana_float2string(self, Y[i])
        Z[i] = Transformacje.zamiana_float2string(self, Z[i])
    
    with open(xyz_txt , "w",  encoding="utf-8") as plik:
        plik.write(f"Wyniki_obliczen_Geodezyjnych; X, Y, Z, fi, lambda, h, x1992, y1992, x2000, y2000.\n")
        plik.write(f"Znak '-' w koordynatach; x1992, y1992, x2000, y2000 oznacza, że dla podanych współrzędnych ortokartezjańskich (X, Y, Z) po obliczeniu współrzędnych geodezyjnych fi i lambda. fi i lambda nie należą do dozwolonych współrzędnych \ngeodezyjnych układów PL1992, PL2000.\n")
        plik.write("-"*221)
        plik.write(f"\n")
        plik.write(f"|          X          |          Y          |          Z          |          fi         |        lambda       |          h          |        x1992        |        y1992        |        x2000        |        y2000        |")
        plik.write(f"\n")
        plik.write("-"*221)
        plik.write(f"\n")
        for x, y, z, f, l, h, x92, y92, x00, y00 in zip(X, Y, Z, f, l, h, x92, y92, x00, y00):
            plik.write(f"|{x}|{y}|{z}|     {f}|     {l}|{h}|{x92}|{y92}|{x00}|{y00}|")
            plik.write(f"\n")
        plik.write("-"*221)
    
    with open(neu_txt , "w", encoding="utf-8") as plik1:
        plik1.write(f"Wyniki_obliczen_Geodezyjnych; n, e, u.\n")
        plik1.write("-"*154)
        plik1.write(f"\n")
        plik1.write(f"|                        n                         |                        e                         |                        u                         |")
        plik1.write(f"\n")
        plik1.write("-"*154)
        plik1.write(f"\n")
        for n, e, u in zip(N, E, U):
            plik1.write(f"|{n}|{e}|{u}|")
            plik1.write(f"\n")
        plik1.write("-"*154)
        
if __name__ == "__main__":
    geo = Transformacje("GRS80")
    geo.wczytanie_zapisanie_pliku("wsp_inp.txt")