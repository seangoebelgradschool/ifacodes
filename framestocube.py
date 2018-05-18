import astropy.io.fits as pyfits
import numpy as np
import glob
import pdb
import sys

def main(filestub='', outname='', first=-1, last=-1):
    if (filestub=='') | (outname==''):
        print "USAGE: framestocube.main('filestub', 'outputname', firstframe, lastframe')"
        print "The filestub should include all characters except the frame number and .fits."
        print "If the last two numbers are not included, all frames matching the stub will be processed into a cube."

    if ((first == -1) & (last != -1)) | \
       ((first != -1) & (last == -1)):
       print "Please specify first and last frames"

    files = glob.glob(filestub+'*.fits') #get list of filenames

    filenums = np.zeros(len(files))
    for i in range(len(files)):
        filenums[i] = framenum(files[i], filestub)

    print len(files), "files found."

    if len(files)==0: return #nothing found, therefore quit

    if (first==-1) & (last==-1): #if min and max not given, figure it out
        first= np.min(filenums)
        last = np.max(filenums)

    for i in range(len(files)):
        myframenum = framenum(files[i], filestub)
        if i%int(len(files)/100)==0:
            print str(int(round(float(i)/len(files)*100.)))+"% done."
            
        if (myframenum >= first) & (myframenum <= last):
            if 'cube' not in locals(): #first time, initialize cube
                im = pyfits.getdata(files[i])

                cube = np.zeros(np.append(int(last-first+1), np.shape(np.squeeze(im))))

            cube[int(myframenum-first), :,:] = pyfits.getdata(files[i])
            #print "file, layer", files[i], myframenum-first

    pyfits.writeto(outname, cube.astype('uint16'), overwrite=True)
    print "Wrote file", outname


def framenum(filename, filestub):
    num = filename.replace(filestub, '')
    num = num.replace('.fits', '')
    num = num.replace('-', '')
    num = int(num)
    return num
    
def batch(filestub2='', outname='', reset_freq=0):
    if filestub2=='' or outname=='' or reset_freq==0:
        print "Syntax: framestocube.batch(filestub, outname, reset_freq) \n"+\
            " e.g. framestocube.batch('blah-3-', 'mycubes.fits', 320)"
        return

    files = glob.glob(filestub2+'*.fits') 
    if len(files)%reset_freq != 0:
        print "Total number of frames is not divisible by reset frequency."
        print "N_frames:", len(files)
        print
        return
    
    for n in range(len(files)/reset_freq):
        outname2 = outname.replace('.fits', '-' + str(n) + '.fits')
        main(filestub=filestub2, outname=outname2, first=n*reset_freq, last=(n+1)*reset_freq-1)

     

