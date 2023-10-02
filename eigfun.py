import numpy as np
import matplotlib.pyplot as plt

class eigfun1D():
    def __init__(self, xmin, xmax, potential, dx=0.01, eignum=1001):
        self.eignum = eignum
        self.xmin = xmin
        self.xmax = xmax
        self.dx = dx
        self.potential = potential

    #有限差分法
    def FDM(self, n=5):
        x = np.linspace(self.xmin, self.xmax, self.eignum)
        dx = x[1] - x[0]
        V= np.diag(potential(x[1:-1]), 0)
        L = np.diag(np.ones((self.eignum-2)), 0)/(dx**2) \
            - 0.5 * np.diag(np.ones((self.eignum-3)), +1)/(dx**2) \
            - 0.5 * np.diag(np.ones((self.eignum-3)), -1)/(dx**2)
        H = V + L

        w, v = np.linalg.eig(H)  
        idx_sorted = np.argsort(w) 
        eigE, eigV = w[idx_sorted], v[:, idx_sorted]

        psi = np.zeros(self.eignum)     
        psi[1:-1] = eigV[:, n]
        rho = np.zeros(self.eignum)
        rho[1:-1] = eigV[:, n] * eigV[:, n]

        plt.plot(x, rho/(sum(rho)*dx),label=r'$E_{%s}=%.3f\hbar \omega$'%(n, eigE[n]))
        plt.title(r'$E_{%s}=%.3f\hbar \omega$'%(n, eigE[n]))
        plt.xlabel(r'$x$')
        plt.pause(0.01)

        
        
A = 20
def potential(x, k=1):
    return 0.5 * k * x**2

a = eigfun1D(xmin=-A, xmax=A, potential = potential)
a.FDM(n = 10)
