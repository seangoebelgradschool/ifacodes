import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os
import pdb

#Syntax:
#calcfluxvar.main('mylistoffiles.txt', showframes=True, checkframedrops=True)
#
#See files.txt for an example of how to structure the input file list file.
#It should be structured as
#[PRV voltage] PRV
#[Common voltage] [filename (directory optional)]
#e.g. 10 ~/Desktop/myfile.fits
#Comments can be indicated with a #
#e.g. 3.5910 PRV
#If all files are in one directory, that directory can be specified after a 99.
#Otherwise it is assumed the directory info is in the filename.
#
#showframes=True will display the mean and variance frames for the user before
# saving them. Set =False for faster code operation.
#checkframedrops=True plots the mean of each frame in order to see whether
# resets are occurring where they should. Set checkframedrops=False for faster
# code operation.
#
#Written 2018-09 by Sean Goebel, geekyrocketguy@gmail.com.
    
def main(infilename, showframes=True, checkframedrops=True):
    #How many frames are in each ramp?
    reset_freq = 300

    #How many ramps are there?
    n_ramps = 100
    
    #Now, define the frame numbers in each ramp that should be used to produce the cds
    # pairs. Frames between cds_f1 and cds_f2 are averaged to reduce the read noise,
    # and that is subtracted from the average of frames cds_f3 through cds_f4.
    cds_f1 = 20
    cds_f2 = 120
    cds_f3 = 180
    cds_f4 = 280



    junk = np.loadtxt('files.txt', dtype='str')
    volts = junk[:,0].astype(float)
    filenames= junk[:,1]

    if 99 in volts:
        dir = filenames[np.squeeze(np.where(volts == 99))]
        filenames = np.delete(filenames, np.squeeze(np.where(volts == 99)))
        volts = np.delete(volts, np.squeeze(np.where(volts == 99)))
    else:
        dir=''

    #Print stuff out to the user to confirm it.
    lowercase_filenames = [v.lower() for v in filenames]
    if 'prv' in lowercase_filenames:
        prv=volts[lowercase_filenames.index('prv')]
        volts = np.delete(volts, lowercase_filenames.index('prv'))
        filenames = np.delete(filenames, lowercase_filenames.index('prv'))
    else:
        print('PRV not specified in'+infilename+'.')
        prv=input('Please specify PRV\n')
        print('Thank you. In the future, this can be avoided by adding to the text file "[value] PRV"')

    print('PRV =', prv)
    print()

    print("The data read was")
    for i in range(len(volts)):
        print( volts[i], dir+filenames[i])
    print()
        
    print("Frames", cds_f1,"through", cds_f2, "will be averaged and subtracted from the average of frames", cds_f3, "through", cds_f4)
    print()
#update
    for i in np.arange(1, 2):#range(len(filenames)):
        print("Working on file", i+1, "of", len(filenames))
        img = pyfits.getdata(dir+filenames[i])

        if not os.path.isdir('flux_var_files'):
            os.mkdir('flux_var_files')
            print('Created flux_var_files/ directory to store images.')

        imcube = np.zeros((n_ramps, np.shape(img)[1], np.shape(img)[2]))

        #populate the cube of cds pairs
        for j in range(n_ramps):
            imcube[j] = np.mean(img[j*reset_freq + cds_f3 : j*reset_freq + cds_f4], 0) - \
                        np.mean(img[j*reset_freq + cds_f1 : j*reset_freq + cds_f2], 0) 

        average = np.mean(imcube, 0)
        variance = np.var(imcube, axis=0, ddof=1)

        if checkframedrops: #plot median frame values and where resets should occur
            meds = np.zeros(len(img))
            for j in range(len(img)):
                meds[j] = np.median(img[j, 64:-64, 64:-64])
            plt.plot(meds, 'o')
            y0 = np.min(meds)
            y1 = np.max(meds)
            for j in range(n_ramps):
                plt.plot([reset_freq*j, reset_freq*j], [y0, y1], 'r-')
            plt.title('Median of Each Frame')
            plt.ylabel('ADUs')
            plt.xlabel('Frame number')
            plt.show()
        
        if showframes: #Display average and variance frames
            plt.figure(1, figsize=(10, 4), dpi=100)

            plt.suptitle(filenames[i])
            
            plt.subplot(121)
            #select min and max for scaling, presently 4% and 96%
            avgmin = np.sort(average.flatten())[int(round(0.04 * np.size(average)))]
            avgmax = np.sort(average.flatten())[int(round(0.96 * np.size(average)))]
            plt.imshow(average, interpolation='none', vmin=avgmin, vmax=avgmax)
            plt.title('Average')
            plt.colorbar(shrink=0.7)

            plt.subplot(122)
            #select min and max for scaling, presently 4% and 96%
            varmin = np.sort(variance.flatten())[int(round(0.01 * np.size(average)))]
            varmax = np.sort(variance.flatten())[int(round(0.99 * np.size(average)))]
            plt.imshow(variance, interpolation='none', vmin=varmin, vmax=varmax)
            plt.title('Variance')
            plt.colorbar(shrink=0.7)

            plt.tight_layout()
            plt.show()
            

        #save images
        if '/' in filenames[i]:
            #strip off directory information
            filenames[i] = filenames[i][filenames[i].rfind('/') : ]
        newfilename=filenames[i][:filenames[i].find('.fits')] #strip '.fits'
        pyfits.writeto('flux_var_files/'+newfilename+'_avg.fits', average, overwrite=True)
        print('Wrote flux_var_files/'+newfilename+'_avg.fits')
        pyfits.writeto('flux_var_files/'+newfilename+'_var.fits', variance, overwrite=True)
        print('Wrote flux_var_files/'+newfilename+'_var.fits')
