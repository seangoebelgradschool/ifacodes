#started on 2018-02-10
#plots power spectrum of RFI, taking into account bonus clock cycles
#at ends of rows

import pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pdb
#from scipy.stats import mode
#import os
#import time
#from scipy import signal

#welch = True

def main(filename, controller='leach'):
    if (controller.lower() != 'leach') & (controller.lower() != 'pizzabox'):
        print "Please specify controller='leach' or controller='pizzabox'"
        return

    elif (controller.lower() == 'leach'):
        print "Using Leach controller timings."
        t_pix = 3.77e-6 #seconds per pixel
        t_endrow = 0.29e-6 #seconds at the end of each row
        t_clock = 10e-9 #duration of clock cycle
    else: #if pizzabox
        print "Using Pizza Box timings."
        t_pix = 1.0e-6 #seconds per pixel
        t_endrow = 210e-9 #seconds at the end of each row
        t_clock = 10e-9

    print "Reading image..."
    img = pyfits.getdata(filename)[20:320] #[700:]
    #img = np.random.normal(loc=1e4, scale=1000, size=(10, 256, 320))
    #img = np.random.uniform(low=100, high=200, size=(10, 256, 320))

    img_avg = np.median(img, 0)
    for z in range(np.shape(img)[0]):
        img[z] -= img_avg

    #inject a frequency to test recovery of it
    #img[:,:,32:64] = 1000

    print "Image read and median subtracted. Size is ", np.shape(img)
    print "Median, stddev:", np.median(img), np.std(img, ddof=1)

    if 0:
        for i in np.arange(0, len(img), 20):
            plt.imshow(img[i] , interpolation='none', vmin=-90, vmax=90)
            plt.colorbar()
            #plt.title("Example median-subtracted frame")
            plt.title(i)
            plt.show()

    row_clocks = np.shape(img)[2]/32*t_pix/t_clock + t_endrow/t_clock
    pixel_arr = np.zeros((np.shape(img)[0] , row_clocks * (np.shape(img)[1]) ))

    for z in range(np.shape(img)[0]):
        if z % (np.shape(img)[0] /100.) < np.shape(img)[0]/100.: #update progress
            print str(int(round(float(z) / np.shape(img)[0] * 100.)))+ "% complete."
        for y in range(np.shape(img)[1]):
            for x in range(np.shape(img)[2]/32):
                pixel_arr[z, 
                          y*row_clocks +     x*t_pix/t_clock : \
                          y*row_clocks + (x+1)*t_pix/t_clock] \
                    = np.median(img[z,y,x*32:(x+1)*32])

    print "Computing FFT."
    #pad array with 0s to make fft faster
#    print "before length", np.shape(pixel_arr)
    mylen = 2
    while mylen < np.shape(pixel_arr)[1]:
        mylen *= 2
    pad = mylen-np.shape(pixel_arr)[1]  

    pixel_arr = np.append(np.zeros((np.shape(pixel_arr)[0], pad/2)), pixel_arr, axis=1)
    pixel_arr = np.append(pixel_arr, np.zeros((np.shape(pixel_arr)[0], mylen - np.shape(pixel_arr)[1])), axis=1)
#    print "after length", np.shape(pixel_arr)

    fft = np.fft.fft(pixel_arr, axis=1) #/len(pixel_arr)
    xaxis = np.fft.fftfreq(np.shape(pixel_arr)[1], d=t_clock)
    powspec2d = (abs(fft))**2 #still 2d
    powspec = np.mean(powspec2d, axis=0) #make 1d

    cumsum = np.cumsum(powspec)
    
    #show image of power spectra over time
    loc = np.where((xaxis >= 0) & (xaxis <= 0.5*t_pix**(-1)))
    crop = powspec2d[:,np.min(loc):np.max(loc)]
    plt.imshow(crop/2, interpolation='none', norm=LogNorm(), \
                   vmin=np.min(powspec[30:np.max(loc)])/2, \
                   vmax=np.max(powspec[30:np.max(loc)])/2, aspect=4)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Frame Number')
    plt.xticks(np.linspace(0, np.shape(crop)[1], 11) , \
               np.round(np.linspace(np.min(xaxis[loc]), np.max(xaxis[loc]), 11)).astype('int'))
    plt.title('Power Spectra '+filename + ' ' +str(np.shape(img)))
    plt.colorbar(shrink=0.5)
    plt.tight_layout()
    plt.show()


    plt.figure(1, figsize=(12, 8), dpi=100)
    plt.subplot(211)
    plt.plot(xaxis, powspec/2, 'b')
    plt.xlim((row_clocks*np.shape(img)[1]*t_clock)**(-1) , 0.5*t_pix**(-1))
    plt.ylim(0, np.max(powspec[30:np.max(loc)])/2)
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.title('Power Spectrum Pizza Box 180219 1MHz Rate' +str(np.shape(img)))
    plt.title('Power Spectrum '+filename + ' ' +str(np.shape(img)))
    
    plt.ylabel('Power [arb. units]')
    plt.xlabel('Frequency [Hz]')

    plt.subplot(212)
    plt.plot(xaxis, cumsum)
    plt.xlim((row_clocks*np.shape(img)[1]*t_clock)**(-1) , 0.5*t_pix**(-1))
    plt.ylim(0, np.max(cumsum[30:len(powspec)/2]))
    plt.title('Cumulative sum of power')
    plt.ylabel('Cumulative Sum [arb. units]')
    plt.xlabel('Frequency [Hz]')
    #plt.yscale('log')
    #plt.xscale('log')

    plt.tight_layout() #prevent labels from going off the edge
    plt.show()
