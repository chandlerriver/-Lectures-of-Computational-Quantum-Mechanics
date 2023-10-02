import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
import copy

class eigfun1D():
    def __init__(self, xmin, xmax, potential, dx=0.01, \
                 eigmin=0.4, eigmax=2.6, eignum=1001, deig=0.001, \
                 psixmin=0, psixmax=0,
                 dpsixmin=1e-4, dpsixmax=-1e-4, 
                 tol = 1e-6,
                 hbar=1, pi=1, m=1):
        self.eignum = eignum
        self.eigmin = eigmin
        self.eigmax = eigmax
        self.deig = deig
        self.xmin = xmin
        self.xmax = xmax
        self.psixmin = psixmin
        self.psixmax = psixmax
        self.dpsixmin = dpsixmin
        self.dpsixmax = dpsixmax
        self.dx = dx
        self.tol = tol
        self.potential = potential
        self.hbar = 1
        self.pi = 1
        self.m = 1

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
        
        return eigE[n]


    #shooting method
    def Shooting(self, method="positive"):
        #method positive, negative, medium
        x = np.arange(self.xmin, self.xmax, self.dx)
        #psi_p 正向shooting
        #psi_n 反向shooting
        energy = np.arange(self.eigmin, self.eigmax, self.deig)
        
        if method == "positive":
            eigElist, psilist = [],[]
            psi_p = np.zeros(x.size)
            psi_p[0] = self.psixmin
            psi_p[1] = psi_p[0] + self.dx * self.dpsixmin            
            E = energy[0]
            k = np.zeros(x.size)
            k = 2 * self.m/self.hbar**2 * (E-self.potential(x))
            for i in range(2, x.size):
                psi_p[i] = (2*(1-5*self.dx**2/12*k[i-1])*psi_p[i-1]\
                    -(1+self.dx**2/12*k[i-2])*psi_p[i-2])/(1+self.dx**2*k[i]/12)
            flag = psi_p[-1] 
            for E in energy[1:]:
                psi_p = np.zeros(x.size)
                psi_p[0] = self.psixmin
                psi_p[1] = psi_p[0] + self.dx * self.dpsixmin
                k = np.zeros(x.size)
                k = 2 * self.m/self.hbar**2 * (E-self.potential(x))                
                for i in range(2, x.size):
                    psi_p[i] = (2*(1-5*self.dx**2/12*k[i-1])*psi_p[i-1]\
                              -(1+self.dx**2/12*k[i-2])*psi_p[i-2])/(1+self.dx**2*k[i]/12)
                if abs(psi_p[-1]-self.psixmax)<self.tol:    
                    eigElist.append(E)
                    psilist.append(psi_p)
                elif flag*psi_p[-1]<0:
                    eigElist.append(E)
                    psilist.append(psi_p)                      
                else:
                    pass
                flag = psi_p[-1]
            return eigElist

        elif method == "negative":
            eigElist, psilist = [],[]
            psi_n = np.zeros(x.size)
            psi_n[-1] = self.psixmax
            psi_n[-2] = psi_n[-1] - self.dx * self.dpsixmax           
            E = energy[0]
            k = np.zeros(x.size)
            k = 2 * self.m/self.hbar**2 * (E-self.potential(x))            
            for i in range(x.size-3, -1, -1):
                psi_n[i] = (2*(1-5*self.dx**2/12*k[i+1])*psi_n[i+1]\
                    -(1+self.dx**2/12*k[i+2])*psi_n[i+2])/(1+self.dx**2*k[i]/12)
            flag = psi_n[0]            
            for E in energy[1:]:
                psi_n = np.zeros(x.size)
                psi_n[-1] = self.psixmax
                psi_n[-2] = psi_n[-1] - self.dx * self.dpsixmax
                k = np.zeros(x.size)
                k = 2 * self.m/self.hbar**2 * (E-self.potential(x))                                
                for i in range(x.size-3, -1, -1):
                    psi_n[i] = (2*(1-5*self.dx**2/12*k[i+1])*psi_n[i+1]\
                              -(1+self.dx**2/12*k[i+2])*psi_n[i+2])/(1+self.dx**2*k[i]/12)
                if abs(psi_n[0]-self.psixmin)<self.tol:    
                    eigElist.append(E)
                    psilist.append(psi_n)                        
                elif flag*psi_n[0]<0:
                    eigElist.append(E)
                    psilist.append(psi_n)                      
                else:
                    pass                
                flag = psi_n[0]
            return eigElist

        elif method == "medium":
            midpoint = np.random.random(3) * (self.xmax - self.xmin) + self.xmin
            midpoint = np.array([int(round(i-self.xmin, int(np.log10(1/self.dx)))/self.dx) for i in midpoint])
            eigElist, psilist = [],[]
            psi_p = np.zeros(x.size)
            psi_p[0] = self.psixmin
            psi_p[1] = psi_p[0] + self.dx * self.dpsixmin      
            psi_n = np.zeros(x.size)
            psi_n[-1] = self.psixmax
            psi_n[-2] = psi_n[-1] - self.dx * self.dpsixmax           
            E = energy[0]
            k = np.zeros(x.size)
            k = 2 * self.m/self.hbar**2 * (E-self.potential(x))            
            for i in range(2, x.size):
                psi_p[i] = (2*(1-5*self.dx**2/12*k[i-1])*psi_p[i-1]\
                    -(1+self.dx**2/12*k[i-2])*psi_p[i-2])/(1+self.dx**2*k[i]/12)
            for i in range(x.size-3, -1, -1):
                psi_n[i] = (2*(1-5*self.dx**2/12*k[i+1])*psi_n[i+1]\
                    -(1+self.dx**2/12*k[i+2])*psi_n[i+2])/(1+self.dx**2*k[i]/12)
            flag = []
            for mid in midpoint:
                dN = (8*psi_n[mid+1]-psi_n[mid+2]-8*psi_n[mid-1]+psi_n[mid-2])/12*self.dx
                N  = psi_n[mid]
                dP = (8*psi_p[mid+1]-psi_p[mid+2]-8*psi_p[mid-1]+psi_p[mid-2])/12*self.dx
                P  = psi_p[mid]
                flag.append((dN/N)-(dP/P))
            flag = np.array(flag)
            for E in energy[1:]:
                psi_p = np.zeros(x.size)
                psi_p[0] = self.psixmin
                psi_p[1] = psi_p[0] + self.dx * self.dpsixmin      
                psi_n = np.zeros(x.size)
                psi_n[-1] = self.psixmax
                psi_n[-2] = psi_n[-1] - self.dx * self.dpsixmax           
                k = np.zeros(x.size)
                k = 2 * self.m/self.hbar**2 * (E-self.potential(x))            
                for i in range(2, x.size):
                    psi_p[i] = (2*(1-5*self.dx**2/12*k[i-1])*psi_p[i-1]\
                              -(1+self.dx**2/12*k[i-2])*psi_p[i-2])/(1+self.dx**2*k[i]/12)
                for i in range(x.size-3, -1, -1):
                    psi_n[i] = (2*(1-5*self.dx**2/12*k[i+1])*psi_n[i+1]\
                              -(1+self.dx**2/12*k[i+2])*psi_n[i+2])/(1+self.dx**2*k[i]/12)
                newflag = []
                for mid in midpoint:
                    dN = (8*psi_n[mid+1]-psi_n[mid+2]-8*psi_n[mid-1]+psi_n[mid-2])/12*self.dx
                    N  = psi_n[mid]
                    dP = (8*psi_p[mid+1]-psi_p[mid+2]-8*psi_p[mid-1]+psi_p[mid-2])/12*self.dx
                    P  = psi_p[mid]
                    newflag.append((dN/N)-(dP/P))
                newflag = np.array(newflag)
                for i in range(flag.size):
                    if flag[i]*newflag[i]<0:
                        eigElist.append(E)
                        psilist.append(psi_n)                      
                        break          
                    else:
                        pass
                flag = copy.deepcopy(newflag)
            return eigElist

    
        
A = 10
def potential(x, m=1, omega=1):
    return 0.5 * m * omega**2 * x**2

a = eigfun1D(xmin=-A, xmax=A, potential = potential)
a.FDM(n = 10)
##eigElist, psilist = a.Shooting(method = "positive")
##print(eigElist)
