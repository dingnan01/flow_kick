# -*- coding: utf-8 -*-
"""
functions taking the derivatives of the five variables
"""

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


#constants
a_n = 1
a_e = 2
N_n = 1
N_e = 2
b_nn = 1
b_en = 2
b_ne = .5
b_ee = .5
m_n = .1
m_e = .1
k_n = .1
k_e = .05
sigma = 1
gamma_n = 1
gamma_e = 1
c = 2.25

#kicks
tau = 10
H = 0.5
delta = .1

#haying
def kick(X):
    const = -H/(delta*(X[0]+X[1])+X[2]+X[3])
    
    kick = (const*delta*X[0], 
            const*delta*X[1],
            const*X[2],
            const*X[3],
            0)
    
    return kick


def run_model():
    
    x0 = (1, 1, 1, 1, 1)
    
    X=[]
    X = np.append(X, odeint(f, x0, range(1,tau)))
    X= np.append(X, X[-5:]+kick(X[-5:]))
    X=np.reshape(X, (tau, 5))
    
    for i in range(10):
        X = np.append(X, odeint(f, X[-1,:], range(1,tau)))
        X= np.append(X, X[-5:]+kick(X[-5:]))
        X=np.reshape(X, (tau*(i+2), 5))

    return X


def f(X,t):
    print(X)
    f_1=(a_n*(X[4]-N_n)-b_nn*X[2]-b_en*X[3])*X[0]
    f_2=(a_e*(X[4]-N_e)-b_ee*X[3]-b_ne*X[2])*X[1]
    
    return (
            f_1 -m_n*X[0],
            f_2-m_e*X[1],
            m_n*X[0] - k_n*X[2],
            m_e*X[1] - k_e*X[3],
            c - sigma*X[4] + gamma_n*(f_1+k_n*X[2])+gamma_e*(f_2+k_e*X[3])
    )
    