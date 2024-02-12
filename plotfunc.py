#   Example of how to use scientific python to plot a function.
#   See https://www.scipy.org/ for all documentation and instruction
#   for installing scientific python. 
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
def f(x):         # Define your function here

    y = np.sin(x)/x
    return y
# Evenly spaced x'es
xl = -30
xh = 30
dx = (xh-xl)/200.
x = np.arange(xl,xh,dx)
y = f(x)    # the command to calculate the function values
plt.plot(x,y) # 'make' the plot (array of points connected with lines) 

plt.title('fu ')
plt.xlabel('x')
plt.ylabel('y')

plt.show()   # Show the plot (Without this you don't see anything on the screen)
