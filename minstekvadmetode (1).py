# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 11:24:12 2021

@author: Simon
"""

import numpy as np
import math
import matplotlib.pyplot as plt

c = 299792458
e = 1.602176565*10**-19
h_planck = 6.62607015*10**(-34)

bølgelengde = [365*10**-9, 405*10**-9, 436*10**-9, 546*10**-9, 577*10**-9]
frekvens = []

for i in bølgelengde:
    frekvens.append(c/i)

Array = np.array([[frekvens[0], -1.884],[frekvens[1], -1.473],
                  [frekvens[2], -1.332],[frekvens[3], -0.753],
                  [frekvens[4], -0.712]])

def Areg(Array):

    sumx_i = Array.sum(axis=0)[0]
    sumy_i = Array.sum(axis=0)[1]
    sumx_isquared = np.sum(Array**2, axis=0)[0]
    sumx_i_y_i = 0
    N = Array.shape[0]
    
    for i in Array:
        sumx_i_y_i += i[0]*i[1]

    A = (sumx_isquared * sumy_i - sumx_i * sumx_i_y_i)/(N*sumx_isquared - sumx_i**2)
    
    B = (N * sumx_i_y_i - sumx_i * sumy_i)/(N * sumx_isquared - sumx_i**2)
    
    AB = [A,B]
    
    return AB


def detteersigmay(Array):
    lmao = 0
    A = Areg(Array)[0]
    B = Areg(Array)[1]
    N = Array.shape[0]
    for i in range(N):
        lmao += (Array[i,1] - A - B * Array[i,0]  )**2    
        
    sigmay = math.sqrt(lmao/(N - 2))
    return sigmay


def sigmaA(Array):
    
    N = Array.shape[0]
    sumx_isquared = np.sum(Array**2, axis=0)[0]
    sumx_i = Array.sum(axis=0)[0]
    
    sigmaA =detteersigmay(Array) * math.sqrt(sumx_isquared/(N*sumx_isquared - sumx_i**2))
    
    return sigmaA
    

def sigmaB(Array):
    
    N = Array.shape[0]
    sumx_isquared = np.sum(Array**2, axis=0)[0]
    sumx_i = Array.sum(axis=0)[0]
    
    sigmaB = detteersigmay(Array) * math.sqrt(N/(N*sumx_isquared - sumx_i**2))
    
    return sigmaB


def ytilplot(Array):
    
    AB = Areg(Array)
    y= list()
    
    a = Array[:,0]
    for i in range(Array.shape[0]):
        y.append(AB[0] + AB[1] * a[i])

    return y

    
def plotreg(Array):
    
    plt.style.use('seaborn-whitegrid')
    plt.errorbar(Array[:,0], Array[:,1], yerr=detteersigmay(Array), fmt='o', color='black',
             ecolor='b', elinewidth=3, capsize=0)
    plt.plot(Array[:,0],ytilplot(Array), 'r')
    plt.xlabel('Frekvens [Hz]')
    plt.ylabel('Stoppotensiale [V]')
    plt.show()

    
    
#plotter minste kvadratsmetode for Array 2
plotreg(Array)
    
sigmay = detteersigmay(Array)

B = Areg(Array)[1]

h_estimat = -B*e

delta_h = e*sigmaB(Array)


h_estimat_øvre = h_estimat + delta_h

# #finner y0, A og B
# y0 = Areg(Array)[0]
# A = Areg(Array)[0]
# B = Areg(Array)[1]

# #finner g i likning: B = g/k
# k = 9.81/B
    
# #finner usikkerhet i y, y0 og B
# usikkerhet = lmao(Array2)
    
    
# # print(sigmaB(Array2))
    




