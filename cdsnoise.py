#!/usr/bin/env python
#Takes as input a filename and two frame numbers. Displays a CDS image
# and prints standard deviation of CDS pairs in the data cube.
#If in doubt about syntax, just type "python cdsnoise.py" or "cdsnoise.py" and it will print
# usage information.
#Updated by Sean on 2017/02/08 for use by Shane/Don.

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pdb
import sys
import time

#get arguments
args=sys.argv
if len(args) != 4: #If incorrect number of arguments were given
    print
    print "Usage: cdsnoise.py myfilename.fits frame1 frame2"
    print
    print "Displays a CDS image of the user-specified frames. Computes the standard deviation of that frame pair, and also the standard deviation of the image cube."
    print
    print "frame1 and frame2 should be two of the following: number, number, all, avg"
    print " If they are both numbers, the code forms a cds pair with those frames. "
    print " If one is a number and the other is avg, then the average frame from the cube is computed and subtracted from that number frame."
    print " If one is all and the other is avg, then the average frame is subtracted from all frames in the cube, and it loops to display the images to the user."
    print
    print "The standard deviation of the whole cube is also calculated using the selected method (average frame subtraction or temporally adjacent image pair subtraction)."
    print "In all cases, this code excludes the first 10 frames from stddev calculation to minimize settling effects."
    print
    #return

else:
    filename=str(args[1])
    img = pyfits.getdata(filename) #read in image

    args[2] = args[2].lower() #make lower case
    args[3] = args[3].lower()

####FIGURE OUT WHAT MODE THE USER WANTS
    if ('all' in args) & ('avg' in args):
        print
        print "Subtracting the average frame from all frames."
        f2 = np.median(img[10:], axis=0) #average frame
        looprange = range(10, len(img))
        title2 = 'avg' #for plot title

    elif ('all' in args): #'all' and a number given
        print
        print "Why do you want to subtract all frames from one frame?"
        print "Try again."
        print
        sys.exit(0)

    elif ('avg' in args): #'avg' and a number given
        print
        print "Subtracting the average frame from one frame."
        if args[2] == 'avg':
            f1 = int(args[3])
        else:
            f1 = int(args[2])

        f2 = np.median(img[10:], axis=0) #average frame
        looprange = range(f1,f1+1) #for for loop
        title2 = 'avg' #for plot title

    else: #two numbers given
        print
        print "Subtracting two user-specified frames."
        f1 = int(args[2])
        f2 = int(args[3])
        looprange = range(f1, f1+1) #for for loop
        title2 = str(f2) #for plot title


####DISPLAY CDS FRAME FOR THE USER TO VIEW.
    for f1 in looprange:
        if np.size(f2)==1: f2=img[f2] #f1 is a number, f2 is a 2d image.

        cds = img[f1] - f2 #form CDS pair with given frames
        #if np.median(cds) < 0: cds *= -1 #invert image if necessary

        stddev=np.std(cds[:,32:], ddof=1)
        print "Stddev of CDS pair formed from user-specified frames is:", stddev

        #Display CDS for user
        #Set image scaling in ADUs with these two lines.
        mymin=np.sort(cds.flatten())[0.01*np.size(cds)] #1st percentile
        mymax=np.sort(cds.flatten())[0.99*np.size(cds)] #99th percentile
        mymin = -15
        mymax = 15
        plt.imshow(cds, interpolation='none', vmin=mymin, vmax=mymax)
        #plt.imshow(cds, interpolation='none', norm=LogNorm(), vmin=2000, vmax=mymax)
        plt.title(filename + '\n CDS ' + str(f1) + "-" + title2 + ', stddev='+str(stddev)[:4])
        plt.colorbar()
        plt.show()


####NOW CALCULATE THE STANDARD DEVIATION IN THE IMAGE CUBE
    stddevs=np.array([]) #stores standard deviation of each pair
    if 'avg' in args:
        print "Calculating the standard deviation by subtracting the avg image from each frame."
        for i in range(10, np.shape(img)[0]): #excludes first 10 images, increments by 1
            stddevs=np.append(stddevs, np.std(img[i]-f2, ddof=1) )
    else:
        print "Calculating the standard deviation by subtracting adjacent images."
        for i in range(10, np.shape(img)[0], 2): #excludes first 10 images, increments by 2
            stddevs=np.append(stddevs, np.std(img[i+1]-img[i], ddof=1) )

    print "Average stddev is:", np.median(stddevs)
    print
    #time.sleep(5)
