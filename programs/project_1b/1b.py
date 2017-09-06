import numpy as np, matplotlib.pyplot as plt, sys

n = int(sys.argv[1])
h = 1/(1+ float(n))
x = np.linspace(0, 1, n)
b = h**2*100*np.exp(-10*x) 

R = np.zeros((n, n+1))
R1 = R.copy()
R2 = R.copy()
R1[0] = R[0]
for i in range(n):  # setting up initial matrix
    if i < n-1:
        R[i, i+1] = -1
        R[i+1, i] = -1
    R[i, i] = 2
    R[i, -1] = b[i]


R1[0] = R[0]
for i in range(n-1):
    R1[i+1] = R[i+1] - R[i+1, i]*R1[i]/R1[i, i]  # forward
    
R2[-1] = R1[-1]    
for i in range(n)[-1:0:-1]:
    R2[i-1] = R1[i-1] - R1[i-1, i]*R2[i]/R2[i, i] # backward

result = R2[:, -1]
exact = 1-(1-np.exp(-10))*x - np.exp(-10*x)
plt.plot(x, result, label='numerical')
plt.plot(x, exact, label='exact')
plt.legend()
plt.title('N = %g' %n)
plt.show()
