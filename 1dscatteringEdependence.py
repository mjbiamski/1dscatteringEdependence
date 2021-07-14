#import libraries
import numpy as np 
import datetime
import matplotlib.pyplot as plt

#generate empty arrays
T_tot = []
T_sum = []
T_avg = []
n = 0 #define counter
realizations = 10
potentials = 50
E_max = 89

while n < realizations:
    start = datetime.datetime.now()
    E = 0.001
    g_array = [] #define empty random number array
    E_array = []
    kd_array = [] #define empty array for kd values
    MDtot_array = [] #define empty total transfer matrix array
    T_array = [] #define transmission coefficient array 
    g = 0  
    for i in range(potentials): 
        g = np.random.uniform(-0.01,0.01) #generate random number
        g_array.append(g) #populate g_array
    while E < E_max: 
        E_array.append(E) #populate E_array
        E += 0.001
    for i in range(len(E_array)):
        k = np.sqrt(E_array[i])
        kd_array.append(2*k) #populate kd_array
        MDtot = np.matrix([[1,0],[0,1]])
        for j in range(len(g_array)):
            D = 2
            d = 1
            k = D*np.sqrt(E_array[i])/d
            B = 2*D/d
            #define transfer matrix for one potential
            MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[j])*k), 
                          (-1j*B)/(2*k)], 
                       [(1j*B)/(2*k), 
                        (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[j])*k)]]) 
            MDtot = np.matmul(MDtot,MD0) #compute total transfer matrix
        MDtot_array.append(MDtot) #populate total transfer matrix array
    for i in range(len(E_array)):
        MD = MDtot_array[i]
        mdtot11 = MD[1,1]
        mdtot11star = np.conj(mdtot11)
        t1 = 1/(mdtot11*mdtot11star) #compute transmission coefficient
        T_array.append(t1) #populate transmission coefficient array
    T_tot.append(T_array)
    n += 1 #increment counter
    print(str(n) + " realization(s): " + str(datetime.datetime.now()-start) + 
          " seconds") #show time per realization

T_sum = [sum(i) for i in zip(*T_tot)] #sum transmission coefficients

#average transmission coefficients for realizations
for i in range(len(T_sum)):
    T_avg.append(T_sum[i]/realizations) 
    
print("-----------------------")

T_totB = []
T_sumB = []
T_avgB = []
n = 0
realizations = 10
potentials = 50
E_max = 89

while n < realizations:
    start = datetime.datetime.now()
    E = 0.001
    g_array = []
    E_array = []
    kd_array = []
    MDtot_array = []
    T_array = []
    g = 0  
    for i in range(potentials):
        g = np.random.uniform(-0.05,0.05)
        g_array.append(g)
    while E < E_max:
        E_array.append(E)
        E += 0.001
    for i in range(len(E_array)):
        k = np.sqrt(E_array[i])
        kd_array.append(2*k)
        MDtot = np.matrix([[1,0],[0,1]])
        for j in range(len(g_array)):
            D = 2
            d = 1
            k = D*np.sqrt(E_array[i])/d
            B = 2*D/d
            MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[j])*k), 
                          (-1j*B)/(2*k)], 
                       [(1j*B)/(2*k), 
                        (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[j])*k)]])
            MDtot = np.matmul(MDtot,MD0)
        MDtot_array.append(MDtot)
    for i in range(len(E_array)):
        MD = MDtot_array[i]
        mdtot11 = MD[1,1]
        mdtot11star = np.conj(mdtot11)
        t1 = 1/(mdtot11*mdtot11star)
        T_array.append(t1)
    T_totB.append(T_array)
    n += 1
    print(str(n) + " realization(s): " + str(datetime.datetime.now()-start) + " seconds")

T_sumB = [sum(i) for i in zip(*T_totB)]

for i in range(len(T_sumB)):
    T_avgB.append(T_sumB[i]/realizations)

print("-----------------------")

T_totC = []
T_sumC = []
T_avgC = []
n = 0
realizations = 10
potentials = 50
E_max = 89

while n < realizations:
    start = datetime.datetime.now()
    E = 0.001
    g_array = []
    E_array = []
    kd_array = []
    MDtot_array = []
    T_array = []
    g = 0  
    for i in range(potentials):
        g = np.random.uniform(-0.1,0.1)
        g_array.append(g)
    while E < E_max:
        E_array.append(E)
        E += 0.001
    for i in range(len(E_array)):
        k = np.sqrt(E_array[i])
        kd_array.append(2*k)
        MDtot = np.matrix([[1,0],[0,1]])
        for j in range(len(g_array)):
            D = 2
            d = 1
            k = D*np.sqrt(E_array[i])/d
            B = 2*D/d
            MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[j])*k), 
                          (-1j*B)/(2*k)], 
                       [(1j*B)/(2*k), 
                        (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[j])*k)]])
            MDtot = np.matmul(MDtot,MD0)
        MDtot_array.append(MDtot)
    for i in range(len(E_array)):
        MD = MDtot_array[i]
        mdtot11 = MD[1,1]
        mdtot11star = np.conj(mdtot11)
        t1 = 1/(mdtot11*mdtot11star)
        T_array.append(t1)
    T_totC.append(T_array)
    n += 1
    print(str(n) + " realization(s): " + str(datetime.datetime.now()-start) + " seconds")

T_sumC = [sum(i) for i in zip(*T_totC)]

for i in range(len(T_sumC)):
    T_avgC.append(T_sumC[i]/realizations)




plt.figure(figsize=(8,5),dpi=200)
plt.plot(kd_array,T_avgC,label='10%',linewidth=0.7)
plt.plot(kd_array,T_avgB,label='5%',linewidth=0.7)
plt.plot(kd_array,T_avg,label='1%',linewidth=0.7)
plt.xticks([0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi, 5*np.pi, 6*np.pi],
          [r'$0$', r'$\pi$', r'$2\pi$', r'$3\pi$', r'$4\pi$', r'$5\pi$', r'$6\pi$'])
plt.xlabel("$kd$")
plt.ylabel("$T_{N}$")
plt.legend(loc='upper left')
plt.ylim([0,1])
plt.show()