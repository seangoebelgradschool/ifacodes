import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pdb

#updated to work with RR datasets

def main():
    volts = [3.5604, 3.5701, 3.5798, 3.5896, 3.5992, 3.6101, 3.6198]
    avg_adus = np.zeros(len(volts))
    
    #reset_freq = 200 #how often is detector reset?
    
    for j in range(len(volts)):
        voltage = volts[j]

        print "Working on dataset " + str(j+1) + " of " + str(len(volts)) + "."

        img = pyfits.getdata('/home/pizzabox/Desktop/seanstest_20180420/vg_20180517_prv='+str(voltage)+'.fits').astype(float)
        avg_adus[j] = np.median(img[1500:])

        if 1:
            meds = np.zeros(len(img))
            for i in range(len(img)):
                meds[i] = np.median(img[i])
            plt.plot(meds, 'o')
            plt.plot([0, len(meds)], [avg_adus[j], avg_adus[j]])
            plt.title(str(voltage)+ 'V, mean=' + str(avg_adus[j])[:7])
            plt.show()

    coeffs=np.polyfit(volts, avg_adus, 1)
    p = np.poly1d(coeffs)

    print "Volt gain is: ", 1.e6 / coeffs[0], 'uV/ADU'
    
    plt.plot(volts, avg_adus, 'o')
    plt.plot(volts, p(volts), '-', label=str(1.e6 / coeffs[0])[1:6]+' uV/ADU')
    plt.title("Volt Gain Mk14 SAPHIRA M06665-23 + Cryo Preamp + PB Rev3")
    plt.xlabel('PRV [V]')
    plt.ylabel('Mean ADUs')
    plt.legend(loc=1)
    plt.show()

    #pdb.set_trace()
