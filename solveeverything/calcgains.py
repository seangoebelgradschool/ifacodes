#This code reads in the average ADU and variance data from stats.txt, compensates
#for the changing node capacitance, and then plots the resulting charge and avalanche
#gains.
#
#The functions to be called from the command line are:
# calcgains.use_ians_nodecap(), which uses #the node capacitances measured by Ian
#  Baker, or
# calcgains.minimize_cg_variance(), which tries to find a node capacitance which
#  reduces the dispersion of charge gain measurements, or
# calcgains.usemyownnodecap(), which allows the user to arbitrarily specify his/her
#   own node capacitances.
#The other functions (calcgains,etc) are called by these, and there generally isn't
#a need to call them directly from the command line.
#
#Near the top of this code, the user should specify the voltages at which the
#avalanche gain 'turns on'. In testing, thresh1=3.1 and thresh2=4.9 work well. These
#indicate where the avalanche gain fitting becomes an exponential function.
#
#Written 2018-12 by Sean Goebel, geekyrocketguy@gmail.com.

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pdb
from scipy import exp
import sys
import datetime

#these define where the transition between the functions for
#fitting to avalanche gains occurs
thresh1 = 3.1 #defines V for end of quadratic fit
thresh2 = 4.9 #defines V for start of exponential fit

def use_ians_nodecap():
    #This performs a exponential decay fit to the node capacitances mesaured by Ian Baker.
    #It then propagates that forward when calculating avalanche and charge gains.

    #Get node capacitance as a function of voltage.
    v_nc, nc = iannodecapacitance(showplot=False)

    #Calculate gains. Store the charge gains in cg.
    v_cg, cg = calcgains(v_nc, nc, showplot=True, savefig=False)

    
def minimize_cg_variance():
    #This varies the node capacitance function in order to minimize dispersion
    #in the charge gains.
    #It assumes that the node capacitance has the functional form of
    #y = a * exp( -x / tau ) + rho_0. It then performs a grid search to find which
    #values minimize the standard deviation of charge gains. Unfortunately, it is
    #not possible
    #to solve for a, tau, and rho_0 simultaneously, since the results are degenerate
    #and there are many well-fitting solutions. Therefore, I will asssume a value of
    #rho_0 (i.e. the node capacitance at high bias voltage) and solve for a and tau.
    #This is not a mathmatically rigorous optimization--I just do a grid search. 
    #Doing a finer grid spacing would cause the code to take longer but produce more
    #precise results.

    #Define the grid to search over here.
    a_grid = np.linspace(10,30, 21) #min, max, steps
    rho_0_grid = [30.6] #this could be replaced by an array if you assume a value
    # for one of the other variables
    tau_grid = np.linspace(1., 2., 21) #min, max, steps

    print("Assuming a node capacitance functional form of y = a * exp( -x / tau ) + rho_0 and searching for a, tau, and rho_0 that minimize standard deviation of the charge gains.")
    print("Due to degeneracy of the output, one of the three variables must be given, and then the other two can be optimized.")
    print("Searching for the best value of a in the range of " + str(np.min(a_grid)) + ' to ' + str(np.max(a_grid)) + " in " + str(len(a_grid)) + " steps.")
    print("Searching for the best value of tau in the range of " + str(np.min(tau_grid)) + ' to ' + str(np.max(tau_grid)) + " in " + str(len(tau_grid)) + " steps.")
    print("Searching for the best value of rho_0 in the range of " + str(np.min(rho_0_grid)) + ' to ' + str(np.max(rho_0_grid)) + " in " + str(len(rho_0_grid)) + " steps.")
    print('')
    
    v_nc = np.linspace(0,20,2e4+1) #generate voltages from 0 to 20V in 0.001V incs

    stddev_best = 123456 #initialize variable for later

    for x in range(len(a_grid)):
        a = a_grid[x]
        print(str(int(round(float(x)/len(a_grid)*100.))) + "% done.")
        
        for y in range(len(tau_grid)):
            tau = tau_grid[y]
            
            for z in range(len(rho_0_grid)):
                rho_0 = rho_0_grid[z]

                nc = expdec(v_nc, a, tau, rho_0) #calculate node capacitances
                
                #Calculate corresponding charge gains.
                v_cg, cg = calcgains(v_nc, nc, showplot=False, savefig=False, prompt=False)

                #if the stddev of the charge gains is lower than previously achieved,
                #store the values of a, tau, rho_0 that produced it.
                if np.std(v_cg, ddof=1) < stddev_best:
                    a_best, tau_best, rho_0_best = a, tau, rho_0
                    stddev_best = np.std(v_cg, ddof=1)

    #alert user if best value was on the edge of the grid space
    if ((a_best == np.min(a_grid)) or (a_best == np.min(a_grid))) & (len(a_grid) != 1):
        print("The best-fitting value of a was on the edge of the grid space.")
    if ((tau_best == np.min(tau_grid)) or (tau_best == np.min(tau_grid))) & (len(tau_grid) != 1):
        print("The best-fitting value of tau was on the edge of the grid space.")
    if ((rho_0_best == np.min(rho_0_grid)) or (rho_0_best == np.min(rho_0_grid))) & (len(rho_0_grid) != 1):
        print("The best-fitting value of rho_0 was on the edge of the grid space.")

    print('')
    print('A node capacitance of the form')
    print("y = " + str(a_best) + " exp(-x / " + str(tau_best) + ') + ' + str(rho_0_best))
    print('produces the least standard deviation in the charge gains')
    print('')

    #Now produce the plot showing the results.
    nc = expdec(v_nc, a_best, tau_best, rho_0_best) #calculate node capacitances
    v_cg, cg = calcgains(v_nc, nc, showplot=True, savefig=False, prompt=False)


def usemyownnodecap():
    #This demonstrates how you can manually input your own arbitrary values for
    #node capacitance. This bypasses any node cap curve fitting.
    #
    #The simplest thing to do is to define a node cap (in fF) at the same voltage
    #as each of your charge gains.
    #
    #This isn't strictly necessary, however, for charge gains at intermediate voltages.
    #The necessary thing is to input node
    #caps at the extreme voltages of your charge gains. In between these, calcgains()
    #will interpolate linearly for the intermediate node capacitances.

    #input the voltage of each node cap measurement
    v_nc = [0.5519, 2.0475, 6.025, 3.0405, 5.0262, 15.963] 

    #input the corresponding node capacitance
    nc = [53, 38.7, 31.4, 35.3, 32.1, 30.0]
    #^I've arbitrarily chosen these to demonstrate the code.^

    
    if len(nc) != len(v_nc):
        print('')
        print("Your arrays of node capacitances and voltages are not the same length! Cannot continue.")
        print('')
        return()

    #convert data type and sort so that v is in increasing order
    nc = np.array(nc)[np.argsort(v_nc)] 
    v_nc = np.sort(v_nc)
        
    calcgains(v_nc, nc, showplot=True, savefig=False, prompt=True)



    
    
def calcgains(v_nc, nc, showplot=True, savefig=False, prompt=True):
    #Given node capacitances and corresponding voltages, this reads in
    #the average ADU and variance measurements from stats.txt. It then
    #calculates the capacitance-corrected avalance and charge gains.
    #
    #showplot=True plots the resulting node cap, av gain, and charge gains.
    #savefig=True saves a .png file of that plot instead of displaying it.
    #prompt=True displays a plot of the uncorrected avalanche gains to see
    # if there are any aberrant points that need to be removed before
    # line fitting occurs.
    
    #read in mode and variance data
    v_g, avg, var = np.loadtxt('stats.txt', skiprows=1, unpack=True, comments='#')
    
    #calculate node capacitances at the voltage of each of the avg/var measurements
    nc_g = np.zeros(len(v_g)) #stores node capacitances at the right voltages
    for i in range(len(v_g)):
        if v_g[i] in v_nc: #check if the node capacitance has already been calced
            nc_g[i] = nc[np.where(v_nc == v_g[i])]
        #Perform a linear interpolation between the nearest node cap values.
        #This assumes that the node cap function is sampled finely enough
        #that a linear interpolation is accurate. This, as coded, sampled
        #it at 0.001V increments.
        else:
            loc1 = np.where(v_nc == np.max(np.array(v_nc)[v_nc < v_g[i]]))
            loc2 = np.where(v_nc == np.min(np.array(v_nc)[v_nc > v_g[i]]))
            
            nc_g[i] = ( nc[loc1] * (v_nc[loc2] - v_g[i]) + \
                        nc[loc2] * (v_g[i] - v_nc[loc1]) ) / \
                      (v_nc[loc2] -v_nc[loc1])

    #calculate capacitance-corrected avalanche gain with arbitrary scaling
    ag = var / avg * nc_g
    #normalize so that avgain=1 at 1.5V bias
    loc = np.where(abs(v_g - 1.5) == np.min(abs(v_g - 1.5) ))
    ag /= ag[loc]

    if prompt:
        plt.plot(v_g, ag, 'o')
        plt.yscale('log')
        plt.title('Should any points be expunged?')
        plt.xlabel("Volts")
        plt.ylabel("Avalanche aain")
        plt.ylim(0.8*np.min(ag), 1.1*np.max(ag))
        plt.tight_layout()
        plt.show()

        #freaking input and raw_input are different in python 2 and 3, so
        # I have to use this hacky technqiue.
        print("Would you like to continue fitting to this data? If a point is bad, you can comment out that line in stats.txt with a # and rerun this code. (y/n) ")
        cont = sys.stdin.readline()
        if cont.lower() != 'y\n':
            return 0, 0

    #calculate avalanche gain
    #add in a point at [0,1]
    ag = np.append([1], ag) #avalanche gain of 1 at 0v
    v_g = np.append([0], v_g)
    #The sigma keywords give the points a constant weighting in log space
    popt,pcov = curve_fit(combofit, v_g, ag, p0=[0.07, -0.2, 1.1, 2.8, 2.1], \
                          sigma=np.log(ag+1), absolute_sigma=True)
    
    #compute charge gain
    cg = avg / var * combofit(v_g[1:], popt[0], popt[1], popt[2], \
                              popt[3], popt[4])

    if showplot or savefig:
        fig = plt.figure(1, figsize=(7.5,10), dpi=100)
        plt.suptitle("Avalanche and Charge Gains Calculated Simultaneously", fontsize=18)

        v = np.linspace(0,20,2e4+1) #generate voltages from 0 to 20V in 0.001V incs
            
        #Node Capacitance plot
        ax1 = plt.subplot(311)
        v_ic,ic=np.transpose(np.loadtxt('iannodecapacitance.txt'))
        ic *= 1e15 #convert from F to fF
        ax1.plot(v_ic, ic, 'o', label="Ian Baker's measured values")
        ax1.plot(v_nc, nc,label='Node capacitance used in following panels')
                 #label=r'$C (V) = ' + str(popt[0])[:6] + \
                 #r'e^{-V/'+str(popt[1])[:6]+'} + '+ str(popt[2])[:6]+'$')
        plt.ylabel('Node Capacitance [fF]', fontsize=14)
        #plt.xlabel('Volts')
        plt.ylim(30,50)
        #ax1.xlim(0, 11)#1.1*np.max(v_g))
        plt.legend(fontsize=14)
        #plt.tick_params(axis='both', labelsize=14)

        #Avalanche Gain Plot
        ax2 = plt.subplot(312, sharex=ax1)
        ax2.plot(v_g[1:], ag[1:], 'o', label='Derived value')
        ax2.plot(v_g[0], ag[0], 'ro') #inserted value
        #Overplot fit to points
        loc1 = np.where(v < thresh1)
        loc2 = np.where((v >= thresh1) & (v <= thresh2))
        loc3 = np.where(v > thresh2)
        ax2.plot(v[loc1], combofit(v[loc1], \
                                   popt[0], popt[1], popt[2], \
                                   popt[3], popt[4]), 'r-', \
                 label=r'$A = '+str(popt[0])[:5]+'V^{\;2}+'+ \
                 str(popt[1])[:5] + 'V+' + str(popt[2])[:5]+'$'   )
        ax2.plot(v[loc2], combofit(v[loc2], \
                                   popt[0], popt[1], popt[2], \
                                   popt[3], popt[4]), 'c-', \
                 label="Transition region")
        ax2.plot(v[loc3], combofit(v[loc3], \
                                   popt[0], popt[1], popt[2], \
                                   popt[3], popt[4]), 'g-', \
                 label=r'$A = 2^{(V-'+str(popt[3])[:5]+')/'+str(popt[4])[:5]+'}$')

        plt.legend(loc=2, fontsize=14)
        plt.yscale('log')
        #plt.title("Capacitance Corrected Avalanche Gain", fontsize=14)
        #plt.xlabel("Bias Voltage")
        plt.ylabel("Avalanche Gain", fontsize=14)
        plt.ylim(0.8*np.min(ag), 1.5*np.max(ag))
        #plt.xlim(0, 1.1*np.max(v_g))
        #plt.tick_params(axis='both', labelsize=14)
        
        #Now charge gain plot
        ax3 = plt.subplot(313, sharex=ax1)
        ax3.plot(v_g[1:], cg, 'o', label="Derived value")
        #Now plot mean and stddev lines
        ax3.plot([np.min(v_g[1:]), np.max(v_g)], [np.mean(cg) , np.mean(cg)], \
             label='Mean = '+ str(np.mean(cg))[:5] + ' e-/ADU')
        ax3.plot([np.min(v_g[1:]), np.max(v_g)], \
                 [np.mean(cg) - np.std(cg, ddof=1), np.mean(cg) - np.std(cg, ddof=1)],\
                 'r--',
                 label=r'$\sigma=' + str(np.std(cg, ddof=1))[:5] + '$ e-/ADU')
        ax3.plot([np.min(v_g[1:]), np.max(v_g)], \
                 [np.mean(cg) + np.std(cg, ddof=1), np.mean(cg) + np.std(cg, ddof=1)],\
                 'r--')
    
        plt.legend(fontsize=14, loc=0, framealpha=0.7)
        #plt.title("Charge Gain", fontsize=14)
        
        plt.xlabel("Bias Voltage [v]", fontsize=14)
        plt.ylabel('Charge Gain [e-/ADU]', fontsize=14)
        #plt.tick_params(axis='both', labelsize=14)
        plt.xlim(0, 1.1*np.max(v_g))
    
        plt.subplots_adjust(top=0.94, bottom=0.06, left=0.10, right=0.99, hspace=0.1)

        if savefig:
            outtitle='Everything_at_Once_'+datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S") + '.png'
        
            plt.savefig(outtitle, dpi=150)
            print("Wrote", outtitle)
            plt.clf()
        else:
            plt.show()

    return v_g, cg

def iannodecapacitance(showplot=False):
    #Loads Ian's node capacitance measurements. Peforms an exponential
    # decline fit to them. Returns an array of volages (0-20V) and the
    # corresponding node capacitances. If plot=True is set, it plots
    # the fit and data points.
    
    #load Ian's measured node capacitance to use as initial guess for C(V)
    v_ic,ic=np.transpose(np.loadtxt('iannodecapacitance.txt'))
    ic *= 1e15 #convert from F to fF

    #get node capacitances by fitting to ian's data
    popt,pcov = curve_fit(expdec, v_ic, ic, p0=[26, 2, 30])

    print("The function that best fits Ian's node capacitance points is:")
    print("y = " + str(popt[0]) + " exp(-x / " + str(popt[1]) + ') + ' + str(popt[2]))
    
    v = np.linspace(0,20,2e4+1) #generate voltages from 0 to 20V in 0.001V incs
    nc = expdec(v, popt[0], popt[1], popt[2])

    if showplot:
        plt.plot(v_ic, ic, 'o', label="Ian Baker's Measured Values")
        plt.plot(v, nc, label=r'$C (V) = ' + str(popt[0])[:6] + \
                 r'e^{-V/'+str(popt[1])[:6]+'} + '+ str(popt[2])[:6]+'$')
        plt.ylabel('Node capacitance [fF]', fontsize=14)
        plt.xlabel('Volts', fontsize=14)
        #plt.ylim(30,50)
        plt.legend()
        plt.show()
          
    return v, nc

def expdec(x,a,tau, rho_0):
    return a*np.array(exp(-1. * np.array(x) / tau)) + rho_0
    
def combofit(x, a, b, c, d, e):
    #fits ax^2 + bx + c below thresh1
    #fits 2^((x-d)/e) above thresh2
    #computes linear interpolation between the two fcns between thresh1 and thresh2
    x = np.array(x)
    y = np.zeros(len(x))
    if np.min(x) <= thresh1:
        loc = np.where(x<=thresh1)
        y[loc] = a*(x[loc])**2 + b*(x[loc]) + c
    if np.max(x) >= thresh2:
        loc = np.where(x>=thresh2)
        y[loc] = 2**((x[loc]-d)/e)
    if np.max((x>thresh1) & (x<thresh2)) == True:
        #for points between the thresholds, compute weighted average
        # of the two functions
        loc = np.where((x>thresh1) & (x<thresh2))
        scaling = (x[loc] - thresh1) / (thresh2-thresh1)
        y[loc] =  (a*(x[loc])**2 + b*(x[loc]) + c) * (1-scaling) + \
                  (2**((x[loc]-d)/e)) * scaling
    return y
