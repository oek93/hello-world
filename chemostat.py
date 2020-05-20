import numpy as np
import matplotlib.pyplot as plt
from math import sin
from scipy.integrate import odeint





# max growth rate of each genotype
a1 = 0.100 #AA
a2 = 0.125 #Aa
a3 = 0.150 #aa


# half velocity constant of the genotypes
H1 = 0.010 # 0.01  
H2=  0.125 #0.10  
H3 = 0.15 #0.15 

# death rate
d = 0.075

# max and min points for the nutrient concentration
E_max = 2
E_min = 1

# population density
A = [0.05]

# 'w' refers fluctuation rate 

for w in np.array(A):

    #w = 0.01
    #w = 0.05
    
    A0 = [0.1, 0.1, 0.1, 0.1]
    
    t_end = 5000
    t_step = 0.1
    

    
    
    def I(t):
        return E_min + 0.5*(1 + sin(w*t))*(E_max - E_min)
    
    def g1(E):
        return a1*E/(H1 + E)
    
    def g2(E):
        
        return a2*E/(H2+E)
    
    
    def g3(E):
        return a3*E/(H3 + E)
    
    for E in np.arange(0,10,0.1):
    
        def model(y, t):
            E, Y1, Y3, Y2 = y
            
                
        dEdt = I(t) - g1(E)*Y1 - g3(E)*Y3-g2(E)*Y2
        dY1dt = g1(E)*Y1 - d*Y1
        dY3dt = g3(E)*Y3 - d*Y3
        dY2dt = g2(E)*Y2 - d*Y2
                
            
                  
            return [dEdt, dY1dt, dY3dt, dY2dt]
    
    
    y0 = np.array(A0)
    time = np.arange(0, t_end, t_step)
    solution = odeint(model, y0, time)
    
    plt.plot(time, solution[:,0], linewidth=5,label='Env', alpha= 0.7)
    plt.plot(time, solution[:,1], linewidth=5, label='AA')
    plt.plot(time, solution[:,2], linewidth=5, label='aa', alpha = 0.5)
    plt.plot(time, solution[:,3], linewidth=5, label='Aa')
    
    
       
    
    plt.yscale('log')
    plt.ylim(0.001, 1000000)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,  ncol=6, mode="expand", borderaxespad=0.)
    plt.xlabel('Time (h)')
    plt.ylabel('Density')
    
    plt.show()
    
    
    
    
    
  
