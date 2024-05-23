import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

xdata = [0,
0.01,
0.02,
0.03,
0.05,
0.075,
0.1,
0.15,
0.2,
0.3,
0.4,
0.5,
0.6,
0.7,
0.8,
0.9,
1]
ydata = [37.17,
29.98,
26.21,
23.37,
19.15,
15.4,
12.65,
8.89,
6.46,
3.64,
2.16,
1.34,
0.85,
0.55,
0.37,
0.25,
0.17]

#Recast xdata and ydata into numpy arrays so we can use their handy features
xdata = np.asarray(xdata)
ydata = np.asarray(ydata)
plt.plot(xdata, ydata, 'o')

# Define the Gaussian function
def Magnet(x, A, B):
    y = B/((x+A)**3)
    return y

parameters, covariance = curve_fit(Magnet, xdata[4:-4], ydata[4:-4])

fit_A = parameters[0]
fit_B = parameters[1]
print(fit_A)
print(fit_B)

fit_y = Magnet(xdata, fit_A, fit_B)
plt.plot(xdata, ydata, 'o', label='data')
plt.plot(xdata, fit_y, '-', label='fit')
plt.legend()
plt.show()