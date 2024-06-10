# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:30:01 2024
"""

from argparse import ArgumentParser
import skrypt
from skrypt import Transformacje

parser = ArgumentParser()
parser.add_argument('-m', '--m', type=str, help="należy podać elipsoide (Krasowski, GRS80, WGS84")
   
parser.add_argument('-x', '--x', type=float)
parser.add_argument('-y', '--y', type=float)
parser.add_argument('-z', '--z', type=float)
args = parser.parse_args()
   
    
geo = Transformacje(model = args.m)
    
    
f, l, h = geo.algorytm_hirvonena(args.x, args.y, args.z)
fi, lam, ha = geo.algorytm_hirvonena(args.x, args.y, args.z, output="dms")

print("")
print("")
print("Elipsoida:", args.m)
print(f"Wyniki_hirvonen'; fi = {fi}, lam = {lam}, ha = {ha:^.3f}[m]")
  
fi, lam, ha = geo.algorytm_hirvonena(args.x, args.y, args.z)
if lam >= 13.5 and lam <= 25.5 and fi <= 55.0 and fi >= 48.9:
    x92, y92 = geo.flh2PL1992(fi,lam)
    x00, y00 = geo.flh2PL2000(fi,lam)
    print(f"Wyniki_z_transformacji_1992_oraz_2000; X1992 = {x92:^.3f}[m], Y1992 = {y92:^.3f}[m], X2000 = {x00:^.3f}[m], Y2000 = {y00:^.3f}[m]")
else:
    x92 = " '-' " 
    y92 = " '-' " 
    x00 = " '-' " 
    y00 = " '-' " 
    print(f"Wyniki_z_transformacji_1992_i_2000; X1992 = {x92}[m], Y1992 = {y92}[m], X2000 = {x00}[m], Y2000 = {y00}[m]")
    print("niewłaściwe położenie")
  
print("")
print("")