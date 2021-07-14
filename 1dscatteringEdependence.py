import numpy as np 
import datetime
import matplotlib.pyplot as plt
import warnings

#The program below is an updated version of the code I originally developed for my senior research project.  
#It is more efficient, utilizes classes, and allows users to input parameters for the calculation.

#See ElectronicTransmissionPeriodicSolid.pdf in https://mjbiamskiorg.wordpress.com/code/ for a detailed description my research as well as a derivation of the physics used below.

#This program calculates the probability of an electron passing through a 1D crystal with a disordered 
# lattice structure (disorder determined using random number). 
#The calculation is performed for a range of dimensionless energy values (defined here as "kd"). 
#The program allows for the user to input parameters where potentials represent the number of 
# lattices, disorder the non-uniformity in the lattice spacing, and realizations the number of calculations for 
# a given disorder (realizations are averaged for final caluclation).
#Once the calculations have been performed a single plot is created displaying the results.


warnings.filterwarnings("ignore", category=RuntimeWarning)
    
#average runtime for 10^3 potentials is 15 s per realization

#collect input data
disorder_array = list(map(float, input("Enter up to 5 different amounts of disorder (as decimal percentages, separated by a space): ").split(" ")))
pots = int(input("Enter number of potentials: "))
reals = int(input("Enter number of realizations for each set of data: "))


class TransmissionProb:
    def __init__(self, pevalue, potentials, realizations):
        self.pevalue = pevalue
        self.potentials = potentials
        self.realizations = realizations

        #define max energy value (dimensionless)
        E_max = 28*np.pi


        #generate empty arrays for transmission probability calculations
        T_tot = []
        T_sum = []
        self.T_avg = []

        n = 0 #define counter for realizations

        print("-----------------------")
        print("At {p} percent error ".format(p=self.pevalue*100))
        print("-----------------------")
        while n < self.realizations:
            start = datetime.datetime.now()
            E = 0.01
            g_array = [] #define random number array
            E_array = [] #define array for energy values
            self.kd_array = [] #define array for kd values
            MDtot_array = [] #define total transfer matrix array
            T_array = [] #define transmission coefficient array 
            g = 0  
            for i in range(self.potentials): 
                g = np.random.uniform(-1*self.pevalue,self.pevalue) #generate random number
                g_array.append(g) #populate g_array
            while E < E_max: 
                E_array.append(E) #populate E_array
                E += E_max/1000
            for i in range(len(E_array)):
                k = np.sqrt(E_array[i])
                self.kd_array.append(2*k) #populate kd_array
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
            self.T_avg.append(T_sum[i]/self.realizations) 


#create array for storing probabilities calculated with different disorder values
transprob_array = [[] for i in range(len(disorder_array))]

for i in range(len(disorder_array)):
    disordervalues = TransmissionProb(disorder_array[i], pots, reals)
    transprob_array[i].append(disordervalues.kd_array)
    transprob_array[i].append(disordervalues.T_avg)

print("")

#create plots 
plt.figure(figsize=(8,5),dpi=200)

for i in range(len(disorder_array)):
    plt.plot(transprob_array[i][0],transprob_array[i][1],label=str(str(disorder_array[i]*100)+"%"),linewidth=0.7)

plt.xticks([0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi, 5*np.pi, 6*np.pi],
          [r'$0$', r'$\pi$', r'$2\pi$', r'$3\pi$', r'$4\pi$', r'$5\pi$', r'$6\pi$'])
plt.xlabel("$kd$")
plt.ylabel("$T_{N}$")
plt.legend(loc='upper left')
plt.ylim([0,1])
plt.show()