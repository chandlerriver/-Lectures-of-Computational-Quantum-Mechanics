import numpy as np
from numpy import fft
np.random.seed(0)
import matplotlib.pyplot as plt
import copy

class waveevolve1D():
    def __init__(self, xmin, xmax, potential, psi0, 
                 dx=0.01, 
                 tmin=0,  tmax=1, dt=0.01, 
                 tol = 1e-6,
                 hbar=1, pi=1, m=1):

        self.xmin = xmin
        self.xmax = xmax
        self.dx = dx

        self.tmin = tmin
        self.tmax = tmax
        self.dt   = dt

        self.tol = tol
        self.potential = potential
        self.psi0 = psi0
        
        self.hbar = 1
        self.pi = 1
        self.m = 1

    #Crank-Nicolson Method
    def CN(self):
        T = np.arange(self.tmin, self.tmax, self.dt)
        X = np.arange(self.xmin, self.xmax, self.dx)
        psi = np.zeros(X.size, dtype="complex128")
        psi = self.psi0(X)/sum(abs(self.psi0(X))**2*self.dx)**0.5
        plt.title("prob")
        plt.ylabel("$|\psi(x)|^2$", rotation=90)
        plt.plot(X, abs(psi)**2)
        H= np.diag(potential(X)+1/self.dx**2, 0) + \
           np.diag(-1/2/self.dx**2*np.ones(X.size-1), +1) + \
           np.diag(-1/2/self.dx**2*np.ones(X.size-1), -1)
        M = np.diag(np.ones((X.size), dtype="complex128"), 0) + \
            complex(1j)/2/self.hbar*self.dt * H
        N = np.diag(np.ones((X.size), dtype="complex128"), 0) - \
            complex(1j)/2/self.hbar*self.dt * H
        # 从实际结果上看 np.dot(np.linalg.inv(N), M)
        # 还是 np.dot(np.linalg.inv(M), N)
        # 只影响波函数的复部的正负
        # CN 方法还需要学习
        U = np.dot(np.linalg.inv(N), M)
        for i in range(1, T.size):
            psi = np.dot(U, psi)
            plt.plot(X, abs(psi)**2)
            print("归一化结果：", sum(abs(psi)**2)*self.dx)
            plt.pause(0.3)
        return psi

    #Spectral Method
    def Spec(self):
        pass
            
A = 10
def potential(x, m=1, omega=1):
    return 0.5 * m * omega**2 * x**2

def psi0(x):
    return np.sqrt(np.exp(-x**2*40)/(np.sqrt(np.pi/10)/2) * 0.5)
    
a = waveevolve1D(xmin=-A, xmax=A, potential = potential, 
                 psi0 = psi0)
y = a.CN()
