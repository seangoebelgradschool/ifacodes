import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os
import pdb
import glob

#main() produces a map and histogram of charge gains for one image.
#doall() produces maps and histograms of charge gains for all
# fits files in the flux_var_files/ directory. See each code
# for its documentation. However, these aren't capacitance corrected.
#It writes the mode signal and variance to stats.txt files for use by
#3.calcgains.py.

def doall():
    files = np.array(glob.glob('flux_var_files/*illum*_avg.fits'))

    keep = []
    for i in range(len(files)):
        files[i] = files[i].replace('flux_var_files/', '')
        if 'unillum' not in files[i]:
            keep = np.append(keep, i)
    files = files[keep.astype(int)]
    
    print('The files found were:')
    for i in range(len(files)):
        print(files[i])
    print()
        
    prv=input('Please specify PRV\n')
    print("Hopefully the common voltages are specified in the filenames with c=[...]")

    text_file = open("stats.txt", "w")
    text_file.write("Bias_Voltage  Mode_of_CDS_Flux  Mode_of_Variance\n")
    text_file.close()
    
    for i in range(len(files)):
        loc1 = files[i].find('c=')+2
        loc2 = files[i].find('_avg')
        common = float(files[i][loc1:loc2])

        bias = float(prv) - common

        main(files[i], bias, save=True, savedatatotxt=True)


def main(infilename, bias, save=False, savedatatotxt=False):
#Displays a map and histogram of the gain = avg / var for a given bias voltage. Assumes that the files are in the flux_var_files/ directory. The figure can be saved by including save=True in the command, but this is optional and can be left out. The savedatatotxt keyword is used by the doall() command and saves the bias voltage, average adu, and variance to the stats.txt file.
#Syntax: gainmap.main(filename, bias_voltage).
# E.g. gainmap.main('illum_c=1.5435', 2.1, save=True)

    #Try to format the filename to what it should be
    infilename = infilename.replace('.fits', '')
    infilename = infilename.replace('_avg', '')
    infilename = infilename.replace('_var', '')
    infilename = infilename.replace('flux_var_files/', '')

    avg = pyfits.getdata('flux_var_files/'+infilename+'_avg.fits')
    var = pyfits.getdata('flux_var_files/'+infilename+'_var.fits')

    #What region of the detector should be used for the calculations?
    #cropzone = [32:-32, 64:-70] #this syntax doesn't work.
    cz1 = 32 #Distance from top of image
    cz2 = -32 #Distance from bottom of image
    cz3 = 64 #Distance from left of image
    cz4 = -70 #Distance from right of image

    #If the user wants to save flux and variance data to text file
    if savedatatotxt:
        #compute mode of flux and variance by creating a histogram

        plt.figure(1, figsize=(10, 4), dpi=100)
        plt.suptitle(infilename + ", Bias = " + str(bias)[:6] + ' V', size=14)

    ####FIRST DO THE FLUX    
        plt.subplot(121)
        plt.title('Flux')

        avgflat = avg[cz1:cz2, cz3:cz4].flatten() #crop to unmasked area
        avgflat = avgflat[np.isfinite(avgflat)]
        #try to figure out region over which to compute histogram
        mymed = np.median(avgflat)
        mystd = np.std(avgflat, ddof=1)
        #select values within 8 stddevs of the median
        loc = np.where((avgflat < mymed + 5 * mystd) &
                       (avgflat > mymed - 5 * mystd))
        #recompute median and std, rejecting the outliers
        mymed = np.median(avgflat[loc])
        mystd = np.std(avgflat[loc], ddof=1)
        mymin = 0 # mymed - 5*mystd
        mymax = mymed + 5*mystd

        
        junk = plt.hist(avgflat, bins=100, range=[mymin, mymax])
        loc = np.sum(junk[1] < 1) #exclude bad column
        bins = junk[1][loc: ]
        vals = junk[0][loc : ]
        loc = np.squeeze(np.where(vals == np.max(vals)))
        if np.size(loc) == 1:
            avgmode = np.mean(bins[loc : loc+2])
        else: #if the data is pathological and has two modes
            avgmode = np.mean(bins[loc[0] : loc[-1]+1])
        avgmode = str(int(round(avgmode)))
        plt.axvline(int(avgmode), label='Mode = ' + avgmode +' ADU', \
                    color='red', linestyle='dashed')
        plt.legend()

    ####NOW DO THE VARIANCE
        plt.subplot(122)
        plt.title('Variance')

        varflat = var[cz1:cz2, cz3:cz4].flatten() #crop to unmasked area
        varflat = varflat[np.isfinite(varflat)]
        #try to figure out region over which to compute histogram
        mymed = np.median(varflat)
        mystd = np.std(varflat, ddof=1)
        #select values within 8 stddevs of the median
        loc = np.where((varflat < mymed + 5 * mystd) &
                       (varflat > mymed - 5 * mystd))
        #recompute median and std, rejecting the outliers
        mymed = np.median(varflat[loc])
        mystd = np.std(varflat[loc], ddof=1)
        mymin = 0#mymed - 5*mystd
        mymax = mymed + 5*mystd

        
        junk = plt.hist(varflat, bins=100, range=[mymin, mymax])
        loc = np.sum(junk[1] < 1)
        bins = junk[1][loc: ]
        vals = junk[0][loc : ]
        loc = np.squeeze(np.where(vals == np.max(vals)))
        if np.size(loc) == 1:
            varmode = np.mean(bins[loc : loc+2])
        else: #if the data is pathological and has two modes
            varmode = np.mean(bins[loc[0] : loc[-1]+1])
        varmode = str(int(round(varmode)))
        plt.axvline(int(varmode), label='Mode = ' + varmode +' ADU^2', \
                    color='red', linestyle='dashed')
        plt.legend()

        plt.tight_layout()
        plt.subplots_adjust(top=0.86)
        plt.show()

        text_file = open("stats.txt", "a")
        text_file.write(str(bias)[:6] + ' ' + avgmode + ' ' + str(varmode) +'\n')
        text_file.close()
        
    gain = avg / var

    #Crop detector to unmasked region, select finite elements, and flatten it
    # into a 1D array
    gainflat=gain[cz1:cz2, cz3:cz4].flatten()
    gainflat = gainflat[np.isfinite(gainflat)]
    
    plt.figure(1, figsize=(10, 4), dpi=100)
    plt.suptitle(infilename + ", Bias = " + str(bias)[:6] + \
                 " V, Map of Mean / Variance", size=14)#, median=") + \
    #str(round(np.median(gainflat)*100)/100.), size=14)
            
    plt.subplot(121)
    #mymin = np.sort(gainflat)[int(round(0.04 * np.size(gainflat)))]
    #mymax = np.sort(gainflat)[int(round(0.96 * np.size(gainflat)))]
    mymin=0
    mymax=6
    #pdb.set_trace()
    #display image
    plt.imshow(gain, interpolation='none', vmin=mymin, vmax=mymax)
    plt.colorbar()#shrink=0.9)
    plt.title("Gain Map")

    plt.subplot(122)
    #mymin = np.sort(gainflat)[int(round(0.01 * np.size(gainflat)))]
    #mymax = np.sort(gainflat)[int(round(0.99 * np.size(gainflat)))]
    junk = plt.hist(gainflat, bins=100, range=[0,5])
    plt.title("Distribution of Gain Values in Unmasked Region")

    print()
    print("Statistics within masked region:")
    mymed = np.median(gainflat)
    print("Median:", mymed)
    mymean = np.mean(gainflat)
    print("Mean:", mymean)
    mymode = junk[1][np.where(junk[0] == np.max(junk[0]))][0]
    print("Mode:", mymode)
    print()

    plt.axvline(mymed, label='Median = ' + str(round(mymed*100)/100)+' e-/ADU', color='red')
    plt.axvline(mymode, label='Mode = ' + str(mymode)+' e-/ADU', color='yellow')
    plt.legend()

    
    plt.subplots_adjust(top=0.85, bottom=0.08, left=0.04, right=0.98)
    
    if save:
        outtitle='flux_var_files/gainmap_'+infilename+'.png'
        
        #plt.tight_layout()
        plt.savefig(outtitle, bbox_inches='tight')
        plt.clf()
        print("Wrote", outtitle)
    else:
        #plt.tight_layout()
        plt.show()
