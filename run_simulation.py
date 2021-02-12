from astropy.nddata.utils import Cutout2D
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
from glob import glob
import os
import matplotlib.cm as cm
import seaborn as sea
from scipy.stats import norm

def make_model(gid, modelnum, xpos, ypos, field, flter):
    
#     print('Making mask...')
    
    if flter == 'F850LP':
        if field == 'GND':
            path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/North/850/'
            datafile = path+'GOODS-N_acsz_sci_sub.fits'
#            segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodsn_3dhst.v4.0.F160W_seg.fits'
        if field == 'GSD':
            path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/South/850/'
            datafile = path+'goodss_3dhst.v4.0.F850LP_orig_sci.fits'
#            segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodss_3dhst.v4.0.F160W_seg.fits'
            
    else:
    
        flternum = flter.split('F')[1].split('W')[0]

        if field =='GND':
            path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/North/'+flternum+'/'
            datafile = path+'goodsn_3dhst.v4.0.'+flter+'_orig_sci.fits'
#            segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodsn_3dhst.v4.0.F160W_seg.fits'
            
        if field =='GSD':
            path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/South/'+flternum+'/'
            datafile = path+'goodss_3dhst.v4.0.'+flter+'_orig_sci.fits'
#            segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodss_3dhst.v4.0.F160W_seg.fits'
            
    data = fits.open(datafile)[0].data
    cutout = Cutout2D(data, (xpos, ypos), (500, 500))
    
#    segmdata = fits.open(segmfile)[0].data
#    segmcutout = Cutout2D(segmdata, (xpos, ypos), (500, 500))
#
#    cutout.data[(segmcutout.data > 0)] = 0
    
    if not os.path.exists('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/masks/'):
        os.mkdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/masks/')
        
    if not os.path.exists('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/masks/'+str(gid)):
        os.mkdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/masks/'+str(gid))
    
    hdu = fits.PrimaryHDU(cutout.data, header = fits.open(datafile)[0].header)
    hdu.writeto('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/masks/'+str(gid)+'/'+str(gid)+'_'+flter+'_mask'+str(modelnum)+'.fits', overwrite = True)
    
Nsci850 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/GOODS-N_acsz_sci_sub.fits'
Nsci125 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodsn_3dhst.v4.0.F125W_orig_sci.fits'
Nsci160 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodsn_3dhst.v4.0.F160W_orig_sci.fits'
Ssci850 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodss_3dhst.v4.0.F850LP_orig_sci.fits'
Ssci125 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodss_3dhst.v4.0.F125W_orig_sci.fits'
Ssci160 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodss_3dhst.v4.0.F160W_orig_sci.fits'
mask_images = pd.DataFrame({'GND': [Nsci160,Nsci125,Nsci850], 'GSD': [Ssci160,Ssci125,Ssci850]}, index=['F160W', 'F125W', 'F850LP'])

def make_fits(gid, modelnum, flter, field, x, y, shape = (401,401)):

#    print(gid)

#    print('filter: {}'.format(flter))
#    print('field: {}'.format(type(field)))
#    print(mask_images.loc['F160W', 'GND'])
#    try:
#        print(mask_images.loc[flter,field])
#    except:
#        print(gid, modelnum, filter, field)
#    print(mask_images.loc[flter,'GSD'])
#    print(mask_images.loc[flter,field])

    file = mask_images.loc[flter,field]
    hdr = fits.open(file)[0].header
    
    make_model(gid, modelnum, x, y, field, flter)
    masked_file = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/masks/'+str(gid)+'/'+str(gid)+'_'+flter+'_mask'+str(modelnum)+'.fits'
    mask = fits.open(masked_file)[0].data
        
    cutout = Cutout2D(mask, (250, 250), shape, copy=True) #Make cutout

    modeldata = fits.open('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/model_fits/'+flter+'/'+field+'_'+str(gid)+'_'+flter+'_model'+str(modelnum)+'.fits')[0].data

    combarr = cutout.data + modeldata
    
    ###MASK###
    if field == 'GND':
        segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodsn_3dhst.v4.0.F160W_seg.fits'
    if field == 'GSD':
        segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodss_3dhst.v4.0.F160W_seg.fits'
        
    segmdata = fits.open(segmfile)[0].data
    segmcutout = Cutout2D(segmdata, (x, y), (400, 400))

    combarr[(segmcutout.data > 0)] = 0
    ###

    hdu = fits.PrimaryHDU(combarr, header = hdr)
    hdu.writeto('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/sim_fits/'+flter+'/'+field+'_'+str(gid)+'_'+flter+'_fakeim'+str(modelnum)+'.fits', overwrite = True)
    
    
    
    
    
    
Nsci850 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodsn_3dhst.v4.0.F850LP_psf.fits'
Nsci125 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodsn_3dhst.v4.0.F125W_psf.fits'
Nsci160 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodsn_3dhst.v4.0.F160W_psf.fits'
Ssci850 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodss_3dhst.v4.0.F850LP_psf.fits'
Ssci125 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodss_3dhst.v4.0.F125W_psf.fits'
Ssci160 = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/goodss_3dhst.v4.0.F160W_psf.fits'
psf_ims = pd.DataFrame({'GND': [Nsci160,Nsci125,Nsci850], 'GSD': [Ssci160,Ssci125,Ssci850]}, index=['F160W', 'F125W', 'F850LP'])





 #Write intro of .feedme file
def write_intro(flter, field, fitsfile, outputfile, is_model = False, use_mask = False):
    
    if is_model == True:
        model = 1
        minpos = 0
        maxpos = 400
    else:
        model = 0
        minpos = 10
        maxpos = 390
    
    outputfilefits = os.path.basename(outputfile).split('.feedme')[0]+'_output.fits'
    
    # zeropoints taken from https://iopscience.iop.org/article/10.1088/0067-0049/214/2/24/pdf
        
    if flter == 'F850LP':
        zeropoint_og = 24.871
        
    if flter == 'F160W':
        zeropoint_og = 25.946
        
    if flter == 'F125W':
        zeropoint_og = 26.230
            
    inputpsf = os.path.basename(psf_ims.loc[flter, field])
            
#    if use_mask == True:
#        make_mask(gids[0], field, flter)
#        fitsfile = 'masks/'+str(gids[0])+'_mask.fits'
#        xmin = 50
#        xmax = 250
#        ymin = 50
#        ymax = 250
    
    zeropoint=zeropoint_og
        
    
    text = '''================================================================================\n\
# IMAGE and GALFIT CONTROL PARAMETERS\n\
A) {0}           # Input data image (FITS file)\n\
B) {1}     					  # Output data image block\n\
C) none              					  # Sigma image name (made from data if blank or "none") \n\
D) {2}       # Input PSF image and (optional) diffusion kernel\n
E) none                                   # PSF fine sampling factor relative to data\n 
F) none              					  # Bad pixel mask (FITS image or ASCII coord list)\n\
G) none               					  # File with parameter constraints (ASCII file) \n\
H) {10} {11} {10} {11}       # Image region to fit (xmin xmax ymin ymax)\n\
I) 100    100        					  # Size of the convolution box (x y)\n\
J) {8}            					  # Magnitude photometric zeropoint \n\
K) 0.06  0.06      					  # Plate scale (dx dy)   [arcsec per pixel]\n\
O) regular           					  # Display type (regular, curses, both)\n\
P) {9}                 					  # Options: 0=normal run; 1,2=make model/imgblock & quit
\n'''.format(os.path.basename(fitsfile), outputfilefits, inputpsf, None, None, None, None, None, zeropoint, model, minpos, maxpos)
    
    return text


#Print out object info for .feedme file
def write_object(radius, axis_ratio, mag, sersic, PA, outputfile, is_model = False):
    
    if is_model == True:
        model = 0
    else:
        model = 1
    
    coord = '200 200'

    #Print text for .feedme file
    text='0) sersic            		    # Object type \n\
 1) {1}  {4} {4}                    # position x, y        [pixel]\n\
 3) {5}     {4}    		        # total magnitude    \n\
 4) {2}  {4}    		        #     R_e              [Pixels]\n\
 5) {6}   {4}                    # Sersic exponent (deVauc=4, expdisk=1)  \n\
 9) {3}   {4}                    # axis ratio (b/a)   \n\
10) {7}     {4}                    # position angle (PA)  [Degrees: Up=0, Left=90] \n\
\n'.format(None, coord, radius, axis_ratio, model, mag, sersic, PA)

    with open(outputfile, "a") as myfile:
        myfile.write(text)


#Write .feedme file for object and it's nearby objects
def write_feedme(flter, field, scifile, outputfile, rad = 2,  mag = 20, sersic = 1.5, PA = 45, is_model = False, axis_ratio = 0.5, use_mask = False,for_rewrite = None):

    text = write_intro(flter, field, scifile, outputfile, is_model = is_model, use_mask = use_mask)
    save = outputfile
    np.savetxt(save, np.array(text).reshape(1, ), fmt='%s')
            
    write_object(rad, axis_ratio, mag, sersic, PA, outputfile, is_model = is_model)

#For rewrite_galfit01
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 2

def get_output_data_lists(gid, flter, field, modelnum):
    
#    try:

    all_data = [flter, modelnum]

    galfit_file = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/outputs/data/'+flter+'/'+field+'_'+str(gid)+'_'+flter+'_'+'fakeim'+str(modelnum)+'.02'
    #find number of objects in file based on length of file
    number_of_objects=(file_len(galfit_file)-41)/13

    #read feedne file
    feedmepath = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/models/sim_feedmes/'+flter+'/'+field+'_'+str(gid)+'_'+flter+'_'+'fakeim'+str(modelnum)+'.feedme'
    with open(feedmepath) as feedmefile:
        feedme = feedmefile.readlines()

    #first object always on line 40
    object_lines = [40]
    comp = 40
    #make list of lines where data on objects begin
    for i in range(int(number_of_objects)-1):
        comp = comp+13
        object_lines.append(comp)

    with open(galfit_file, 'r') as file:
        data = file.readlines()

    #fix PA, position, and axis ratio for each object in file
    for index, i in enumerate(object_lines):

        mag = float(data[i+1].split()[1])
        rad = float(data[i+2].split()[1])
        sersic_index = float(data[i+3].split()[1])
        axisratio = float(data[i+7].split()[1])
        PA = float(data[i+8].split()[1])

        all_data.append(mag)
        all_data.append(rad)
        all_data.append(sersic_index)
        all_data.append(axisratio)
        all_data.append(PA)
            
#    except:
#        #make list of failed objects
#        if not os.path.exists('failed_run.txt'):
#            text = flter+'_fakeim'+str(modelnum)
#            np.savetxt('failed_run.txt', np.array(text).reshape(1, ), fmt='%s')
#
#        else:
#            with open('failed_run.txt', "a") as myfile:
#                myfile.write('\n'+flter+'_fakeim'+str(modelnum))
        
    return all_data

#For rewrite_galfit01
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 2

#rewrite galfit.01 file to keep position, PA and axis ratio fixed :)
#rewrite galfit.01 file to keep position, PA and axis ratio fixed :)
def rewrite_galfit01(galfitfile, feedmefile):

    feedmepath = feedmefile

    #find number of objects in file based on length of file
    number_of_objects=(file_len(galfitfile)-41)/13

    #first object always on line 40
    object_lines = [40]
    comp = 40
    #make list of lines where data on objects begin
    for i in range(int(number_of_objects)-1):
        comp = comp+13
        object_lines.append(comp)
        
    #read feedme file
    with open(feedmepath) as feedmefile:
        feedme = feedmefile.readlines()
            
    with open(galfitfile, 'r') as file:
        data = file.readlines()

    #fix PA, position, and axis ratio for each object in file
    for index, i in enumerate(object_lines):

        data[i] = data[i].replace(' 1 ', ' 0 ')
        data[i] = data[i].replace(' 1 ', ' 0 ')
        data[i+7] = data[i+7].replace(' 1 ', ' 0 ')
        data[i+8] = data[i+8].replace(' 1 ', ' 0 ')
        
        sersic_index = float(data[i+3].split()[1])
        
        if sersic_index > 8:
#            print('Replacing sersic index of object # '+str(index)+' with 8 and NOT fixing')
            sersic_index = '{:.4f}'.format(sersic_index)
            data[i+3] = data[i+3].replace(str(sersic_index), ' 8 ')

        with open(galfitfile, 'w') as file:
            file.writelines( data )

            
