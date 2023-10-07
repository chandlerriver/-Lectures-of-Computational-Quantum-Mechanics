import warnings
warnings.filterwarnings("error")
import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


#Gauss-Seidel
def GS(A, b, x):
    n = x.size
    X = np.zeros(n)
    for i in range(n):
        try:
            Sum = 0
            for j in range(n):
                Sum += A[i][j]*x[j]
            Sum -= b[i]
            Sum -= A[i][i]*x[i]
            X[i] = -1/A[i][i] * Sum
            print("X[", i, "]:", X[i])
        except Exception as e:
            print(e.__class__.__name__, "\n")
            print("i:", i, "j:", j)
    return X


#Relaxation coefficient
omega = 1.5


maxiter=30
A = np.random.randint(10, 20, (15, 15))
get_x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
b = np.dot(A, get_x.T)
tol = 1e-10

for omega in np.arange(0, 0.2, 0.1):

    x0 = np.zeros(15)

    
    for i in range(maxiter):
        x = (1-omega) * x0 + omega * GS(A, b, x0)
        try:
            if sum(abs(x-x0)**2) < tol:
                print(x)
                print("itertimes:", i, "the relaxation coefficient:", omega)
                break
            else:
                pass
        except Exception as e:
            print(e.__class__.__name__)
            print("x:", x)
            print("x0:", x0)
            
        x0 = x
        if i == maxiter - 1:
            print("not converged in itertimes:", i, "the relaxation coefficient:", omega)
