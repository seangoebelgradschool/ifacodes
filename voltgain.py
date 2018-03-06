import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pdb

def main():
    volts = [3.562, 3.572, 3.581, 3.591, 3.601, 3.612, 3.621]
    avg_adus = np.zeros(len(volts))
    
    reset_freq = 200 #how often is detector reset?
    
    for j in range(len(volts)):
        voltage = volts[j]
        meds = np.zeros(3998)
        print "Working on dataset " + str(j+1) + " of " + str(len(volts)) + "."
        for i in range(len(meds)):
            if i>1200: #exclude first 1000 frames
                if (i%reset_freq > 75): #not in first 75 frames after reset
                    filenum = str(i)
                    while len(filenum) <4: filenum = '0'+filenum
                    meds[i] = np.median(pyfits.getdata('/home/pizzabox/Desktop/SAPHIRA.26FEB2018/LABTEST-SAPHIRA/pbserver.lt/180302_Voltgain_PRV_'+str(voltage)+'V-02-'+filenum+'.fits'))

        meds = meds[np.where(meds != 0)] #remove empty elements
        avg_adus[j] = np.mean(meds)

        if 0:
            plt.plot(meds, 'o')
            plt.plot([0, len(meds)], [avg_adus[j], avg_adus[j]])
            plt.title(str(voltage)+ 'V, mean=' + str(avg_adus[j])[:7])
            plt.show()

    coeffs=np.polyfit(volts, avg_adus, 1)
    p = np.poly1d(coeffs)

    print "Volt gain is: ", 1.e6 / coeffs[0], 'uV/ADU'
    
    plt.plot(volts, avg_adus, 'o')
    plt.plot(volts, p(volts), '-', label=str(1.e6 / coeffs[0])[1:6]+' uV/ADU')
    plt.title("Volt Gain with Mk15 SAPHIRA M09105-27 + Cryo Preamp")
    plt.xlabel('PRV [V]')
    plt.ylabel('Mean ADUs')
    plt.legend(loc=1)
    plt.show()

    pdb.set_trace()
