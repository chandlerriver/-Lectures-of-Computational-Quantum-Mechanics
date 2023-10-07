import warnings
warnings.filterwarnings("error")
import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import time


N = 800
R = 1000
A = np.random.uniform(-R, R, (N, N))
x_get = np.random.uniform(-R, R, (N, 1))
b = np.dot(A, x_get)

s = time.time()
result = np.linalg.solve(A, b)
e = time.time()
print("used time:", e-s)
print(sum((result-x_get)**2))

##errortest
a = 0
b = 10
try:
    a = b/a
except Exception as e:
    print(e.__class__.__name__)
