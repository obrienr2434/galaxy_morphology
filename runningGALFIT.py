from astropy.io import ascii
import numpy as np
import pandas as pd
import subprocess
import os
import warnings
from astropy.nddata.utils import Cutout2D

from photutils import SExtractorBackground
from photutils import MADStdBackgroundRMS
from photutils import make_source_mask
from astropy.stats import SigmaClip

#for plotting
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
import seaborn as sea

class runGALFIT:

    def __init__(self, gid, field, flter, nearby_fitting_region, nearby_gids = None):

        if not (field == 'S' or field == 'N'):
            raise Exception('Field must be S or N.')
            
        if not type(flter) == str:
            raise Exception('Filter must be string.')

        self.gid = gid
        self.field = field
        self.flter = flter
        self.nearby_fitting_region = nearby_fitting_region

        if self.flter == '850':
            if self.field == 'N':
                self.fitsfile = 'GOODS-N_acsz_sci_sub.fits'
                self.fullfield = 'North'

            if self.field == 'S':
                self.fitsfile = 'goodss_3dhst.v4.0.F850LP_orig_sci.fits'
                self.fullfield = 'South'
                
        else:
            
            if self.field =='S':
                self.fitsfile = 'goodss_3dhst.v4.0.F'+flter+'W_orig_sci.fits'
                self.fullfield = 'South'

            if self.field =='N':
                self.fitsfile = 'goodsn_3dhst.v4.0.F'+flter+'W_orig_sci.fits'
                self.fullfield = 'North'

        ### CATALOG ###
        if (field == 'S') or (field == 'South'):
            catalog = '/Users/rosaliaobrien/research/website/catalogs/goodss_3dhst.v4.4.cat'
            
        if (field == 'N') or (field == 'North'):
            catalog = '/Users/rosaliaobrien/research/website/catalogs/goodsn_3dhst.v4.4.cat'

        self.cat = ascii.read(catalog)

        ### GET SKY ESTIMATE ###
        if not nearby_gids == None:
            xmin, xmax, ymin, ymax = self.get_boundaries(nearby_gids)
            path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/'
            skydata = fits.open(path+self.fitsfile)[0].data
            xmid = (xmin+xmax)/2
            ymid = (ymin+ymax)/2
            cutout = Cutout2D(skydata, (xmid, ymid), (400, 400), copy=True, mode='partial')
            masked_cutout = make_source_mask(cutout.data, snr = 1.5,
                            npixels=10, dilate_size=11)
            cutout.data[masked_cutout] = np.nan
            self.scidata = cutout.data
            sigma_clip = SigmaClip(sigma=3)
            bkg = SExtractorBackground(sigma_clip)
            rms = MADStdBackgroundRMS(sigma_clip)
            self.sky_value = bkg.calc_background(cutout.data)
            self.rms_value = rms.calc_background_rms(cutout.data)
        ########################


    #Get position for a SINGLE gid in format 'x y' for using in the .feedme file position
    def get_position(self, gid):
        '''
        Gets positions of objects in fitting region in pixel coord
        
        Parameters:
        -----------
        field : string
            GOODS_South or GOODS_North
        gid : float
            ID number of object
            
        Outputs
        -------
        coordinates of object in string 'x y' 
        '''
        
        #Return x and y coordinates (with a space) of single gid
        x = self.cat[self.cat['id']==gid]['x'].data[0]
        y = self.cat[self.cat['id']==gid]['y'].data[0]
        coord = str(x)+' '+str(y)
        return x, y, coord

    #Get positions for a LIST of gids so that we can determine the boundaries for the cutoutsize of the .feedme file
    def get_positions(self, gids):
        
        x=[]
        y=[]
        
        for gidt in gids:
            
            xt = self.cat[self.cat['id']==gidt]['x'].data[0]
            yt = self.cat[self.cat['id']==gidt]['y'].data[0]
            x.append(xt)
            y.append(yt)
            
        return x, y



     #Get boundaries for .feedme file (with 200pixel border at each edge)
    def get_boundaries(self, gids):
        
        x, y = self.get_positions(gids)

        # d = {'gid':gids, 'field':field, 'xpos':x, 'ypos':y}
        # df = pd.DataFrame(data=d)

        xmin = int(min(x))-200
        ymin = int(min(y))-200
        xmax = int(max(x))+200
        ymax = int(max(y))+200
        
        return xmin, xmax, ymin, ymax




    #Get objects nearby objects being fit to improve the fit
    def get_nearby_objects(self):
        '''
        Parameters
        ----------
        
        gid - float
            gid of object to get nearby objects of
            
        field - str
            field of object w/ gid
            
        nearby_fitting_region - str
            specifies number of nearby objects to grab based on a small, medium, or large radius from the object w/ gid
            
        Outputs
        -------
        list of nearby objects
        '''
            
        x = self.cat[self.cat['id']==self.gid]['x']
        y = self.cat[self.cat['id']==self.gid]['y']
        
        if self.nearby_fitting_region == 'masked':
            o = [self.gid]

        elif self.nearby_fitting_region != 'masked':

            if self.nearby_fitting_region == 'small':
                objects = self.cat[(x-50<self.cat['x']) & (self.cat['x']<x+50) & (y-50<self.cat['y']) & (self.cat['y']<y+50)]
            
            if self.nearby_fitting_region == 'medium':
                objects = self.cat[(x-70<self.cat['x']) & (self.cat['x']<x+70) & (y-70<self.cat['y']) & (self.cat['y']<y+70)]
                
            if self.nearby_fitting_region == 'large':
                objects = self.cat[(x-90<self.cat['x']) & (self.cat['x']<x+90) & (y-90<self.cat['y']) & (self.cat['y']<y+90)]
                
            if self.nearby_fitting_region == 'xlarge':
                objects = self.cat[(x-110<self.cat['x']) & (self.cat['x']<x+110) & (y-110<self.cat['y']) & (self.cat['y']<y+110)]
                
            if self.nearby_fitting_region == 'xxlarge':
                objects = self.cat[(x-140<cself.at['x']) & (self.cat['x']<x+140) & (y-140<self.cat['y']) & (self.cat['y']<y+140)]

            if self.flter == '160':
                fullfilter = 'F160W'
                
            if self.flter == '125':
                fullfilter = 'F125W'
                
            if self.flter == '850':
                fullfilter = 'F850LP'

            # print('nearby objects:', objects)

            # (25-2.5*np.log10(mag)).
            objects = objects[(25-2.5*np.log10(objects['f_'+fullfilter]) < 23) & (objects['flux_radius'] < 15)]

            # print(objects)
                
            o = objects['id'].data

            # print('nearby objects:', o)
            
            for i in o:
                if len(str(i)) > 5:
                    o = np.delete(o, np.where(o==i)[0][0])
            
        return o








    #Write intro of .feedme file
    def write_intro(self, gids, outputfile, sigma_image, use_mask = False):

        use_mask = True

        # print('Use mask?', use_mask, gid)

        # print(gid)
        
        outputfilefits = outputfile+'.fits'
        
        # if use_mask == False:
        xmin, xmax, ymin, ymax = self.get_boundaries(gids)
        
        # zeropoints taken from https://iopscience.iop.org/article/10.1088/0067-0049/214/2/24/pdf
            
        if self.flter == '850':
            if self.field == 'N':
                inputpsf = 'goodsn_3dhst.v4.0.F850LP_psf.fits'
                psffine = 'none'
                zeropoint_og = 24.871

            if self.field == 'S':
                zeropoint_og = 24.871
                inputpsf = 'goodss_3dhst.v4.0.F850LP_psf.fits'
                psffine = 'none'
                
        elif self.flter == '105':
            if self.field =='S':
                zeropoint_og = 23.972
                inputpsf = 'none'
                psffine = 'none'
                
        else:
            
            if self.flter == '160':
                zeropoint_og = 25.946
                
            if self.flter == '125':
                zeropoint_og = 26.230
                
            if self.flter == '606':
                zeropoint_og = 26.511
                
            if self.flter == '775':
                zeropoint_og = 25.671
            
            if self.field =='S':
                inputpsf = 'goodss_3dhst.v4.0.F'+self.flter+'W_psf.fits'
                psffine = 'goodss_3dhst.v4.0.F'+self.flter+'W_orig_wht.fits'

            if self.field =='N':
                inputpsf = 'goodsn_3dhst.v4.0.F'+self.flter+'W_psf.fits'
                psffine = 'goodsn_3dhst.v4.0.F'+self.flter+'W_orig_wht.fits'
                
        if use_mask == True:

            # print('Using mask!', gid)

            segmap = str(self.gid)+'_segm.fits'
            # if os.path.exists('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/segm_maps/'+segmap) == False:
            print('Creating segmentation map...')
            self.make_segm()

            # Move segmentation map to working directory
            os.rename('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/segm_maps/'+segmap, segmap)

        elif use_mask == False:

            # print('Not using mask!', gid)
            segmap = 'none'
        
        zeropoint=zeropoint_og
        
        if sigma_image == None:
            
            sigma_image = 'none'
            
        
        text = '''================================================================================\n\
    # IMAGE and GALFIT CONTROL PARAMETERS\n\
    A) {0}           # Input data image (FITS file)\n\
    B) {1}     					  # Output data image block\n\
    C) {9}              					  # Sigma image name (made from data if blank or "none") \n\
    D) {2}       # Input PSF image and (optional) diffusion kernel\n
    E) none                                   # PSF fine sampling factor relative to data\n 
    F) {10}              					  # Bad pixel mask (FITS image or ASCII coord list)\n\
    G) none               					  # File with parameter constraints (ASCII file) \n\
    H) {4} {5} {6} {7}       # Image region to fit (xmin xmax ymin ymax)\n\
    I) 200    200        					  # Size of the convolution box (x y)\n\
    J) {8}            					  # Magnitude photometric zeropoint \n\
    K) 0.06  0.06      					  # Plate scale (dx dy)   [arcsec per pixel]\n\
    O) regular           					  # Display type (regular, curses, both)\n\
    P) 0                 					  # Options: 0=normal run; 1,2=make model/imgblock & quit
    \n
    '''.format(self.fitsfile, outputfilefits, inputpsf, psffine, xmin, xmax, ymin, ymax, zeropoint, sigma_image, segmap)

        path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/'
        skydata = fits.open(path+self.fitsfile)[0].data
        xmid = (xmin+xmax)/2
        ymid = (ymin+ymax)/2
        cutout = Cutout2D(skydata, (xmid, ymid), (400, 400), copy=True, mode='partial')
        masked_cutout = make_source_mask(cutout.data, snr = 1.5,
                        npixels=10, dilate_size=11)
        cutout.data[masked_cutout] = np.nan
        self.scidata = cutout.data
        sigma_clip = SigmaClip(sigma=3)
        bkg = SExtractorBackground(sigma_clip)
        rms = MADStdBackgroundRMS(sigma_clip)
        self.sky_value = bkg.calc_background(cutout.data)
        self.rms_value = rms.calc_background_rms(cutout.data)

        ###PLOT MASKED IMAGE AND SEGMENTATION MAP###
        if use_mask == False: #Plot 1 subplot if not using segmentation map
            num_plots_x = 1
            num_plots_y = 1
        else:
            num_plots_x = 2
            num_plots_y = 1
        fig, ax = plt.subplots(num_plots_y, num_plots_x, figsize = (5,2))
        ax_og = ax #Make 'original' axis in case you splot segmentation map
        if use_mask == True: 
            ax = ax[0]
        im = ax.imshow(self.scidata, vmin = -100, vmax = 100)
        plt.colorbar(im, ax = ax)
        ax.set_title('gid: {}\nSky: {:.2f}\nRMS: {:.2f}'.format(self.gid, self.sky_value, self.rms_value), size=9)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if use_mask == False:
            plt.show()
        ########################

        if use_mask == True:
            ax = ax_og
            ### PLOT SEGMAP ###
            segmdata = fits.open(segmap)[0].data
            xmid = (xmin+xmax)/2
            ymid = (ymin+ymax)/2
            cutout = Cutout2D(segmdata, (xmid, ymid), (400, 400), copy=True, mode='partial')

            ax[1].imshow(cutout.data)
            ax[1].set_title('segmentation map'.format(self.gid, self.sky_value, self.rms_value), size=9)
            ax[1].set_xticklabels([])
            ax[1].set_yticklabels([])
            plt.show()
            ########################
            
        return text








    #Print out object info for .feedme file
    def write_object(self, gid, coord, outputfile, for_rewrite = None, mask = False, main_gid = True):
        
        #Estimate effective radius from catalog
        radius = self.cat[self.cat['id']==gid]['flux_radius'].data[0]
        
        #Estimate axis raio from catalog and round to 2 decimal places
        axis_ratio = '{:.2f}'.format(self.cat[self.cat['id']==gid]['b_image'].data[0]/self.cat[self.cat['id']==gid]['a_image'].data[0])
        
        # class_star = cat[cat['id']==gid]['class_star'].data[0]

        # column names and mag conversion taken from: 
        # http://cdsarc.u-strasbg.fr/viz-bin/ReadMe/J/ApJS/214/24?format=html&tex=true#sRM3.4
        if self.flter == '125':
            mag = (25-2.5*np.log10(self.cat[self.cat['id'] == gid]['f_F125W'])).item()
        if self.flter == '160':
            mag = (25-2.5*np.log10(self.cat[self.cat['id'] == gid]['f_F160W'])).item()
        if self.flter == '850':
            mag = (25-2.5*np.log10(self.cat[self.cat['id'] == gid]['f_F850LP'])).item()

        if main_gid == True:
                #Print text for .feedme file
                text='#ID: {0} \n\
            0) sersic            		    # Object type \n\
             1) {1}  1 1                    # position x, y        [pixel]\n\
             3) {4}     1    		        # total magnitude    \n\
             4) {2}  1    		        #     R_e              [Pixels]\n\
             5) 2.00   1                    # Sersic exponent (deVauc=4, expdisk=1)  \n\
             9) {3}   1                    # axis ratio (b/a)   \n\
            10) 0     1                    # position angle (PA)  [Degrees: Up=0, Left=90] \n\
            \n'.format(gid, coord, radius, axis_ratio, mag)

        if main_gid == False:
                #Print text for .feedme file
                text='#ID: {0} \n\
            0) sersic                       # Object type \n\
             1) {1}  1 1                    # position x, y        [pixel]\n\
             3) {4}     1                   # total magnitude    \n\
             4) {2}  1                  #     R_e              [Pixels]\n\
             5) 2.00  1                    # Sersic exponent (deVauc=4, expdisk=1)  \n\
             9) {3}   1                    # axis ratio (b/a)   \n\
            10) 0     1                    # position angle (PA)  [Degrees: Up=0, Left=90] \n\
            \n'.format(gid, coord, radius, axis_ratio, mag)

        with open(outputfile+'.feedme', "a") as myfile:
            myfile.write(text)




    #Write .feedme file for object and it's nearby objects
    def write_feedme(self, outputfile, sigma_image, use_mask = False):
        
        if use_mask == True:
            
            objects = [self.gid]
        
        else:
            nearby_objects = self.get_nearby_objects()
            objects = nearby_objects
            self.nearby_objects = nearby_objects
            
        fields = [self.field]*len(objects)

        print('Writing intro...')
        text = self.write_intro(objects, outputfile, sigma_image = sigma_image, use_mask=use_mask)
        save = outputfile+'.feedme'
        np.savetxt(save, np.array(text).reshape(1, ), fmt='%s')

        for i in objects:
            x, y, coord = self.get_position(i)
            
            # if use_mask == True:
            #     coord = str(150)+' '+str(150)
            
            if i == self.gid:
                main_gid = True
            if i != self.gid:
                main_gid = False
                
            self.write_object(i, coord, outputfile, mask = use_mask, main_gid = main_gid)            

            print('Writing object '+str(i)+'...')

        return objects








    def organize_files(self,delete_segm = True):

        if delete_segm == True:
            path = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+str(self.gid)+'_segm.fits'

            if os.path.exists(path) == True:
                print('Deleting segmentation map...')
                os.remove(path)

        if delete_segm == False:
            os.rename(segmap,'/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/segm_maps/'+segmap)

        if not os.path.exists('outputs'):
            os.mkdir('outputs')
            
        if not os.path.exists('outputs/fits'):
            os.mkdir('outputs/fits')
            
        if not os.path.exists('outputs/'+str(self.gid)):
            os.mkdir('outputs/'+str(self.gid))

        feedme = str(self.gid)+'.feedme'
        fits = str(self.gid)+'.fits'
        galfitrun = 'galfit.01'
        galfitrun2 = 'galfit.02'
        galfitrun2c = 'pre-2ndcomponent_fit.txt'
        segmap = str(self.gid)+'_segm.fits'

        if os.path.exists(fits):
            os.rename(str(self.gid)+'.fits', 'outputs/fits/'+str(self.gid)+'.fits')
            
        if os.path.exists(feedme):
            os.rename(str(self.gid)+'.feedme', 'outputs/'+str(self.gid)+'/'+str(self.gid)+'.feedme')
            
        if os.path.exists(galfitrun):
            os.rename('galfit.01', 'outputs/'+str(self.gid)+'/galfit.01')

        if os.path.exists(galfitrun2):
            os.rename('galfit.02', 'outputs/'+str(self.gid)+'/galfit.02')

        if os.path.exists(galfitrun2c):
            os.rename(galfitrun2c, 'outputs/'+str(self.gid)+'/'+galfitrun2c)





    #For rewrite_galfit01
    def file_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 2




    #rewrite galfit.01 file to keep position, PA and axis ratio fixed :)
    def rewrite_galfit01(self, use_galfit02 = None):
        
        if use_galfit02 == True:
            filename = 'galfit.02'
            
        else:
            filename = 'galfit.01'

        feedmepath = str(self.gid)+'.feedme'

        #find number of objects in file based on length of file
        number_of_objects=(self.file_len(filename)-41)/13

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

        #make list of gids based on readme file
        gidlist_from_files = []
        for index, line in enumerate(feedme):
            try:
                if line.split()[0] == '#ID:':
                    gidlist_from_files.append(int(line.split()[1]))

            except:
                pass

        #set fitting_object_index as the index from gidlist_from_files
        for k, j in enumerate(gidlist_from_files):
            if j == self.gid:
                fitting_object_index = k
                
        with open(filename, 'r') as file:
            data = file.readlines()

        #fix PA, position, and axis ratio for each object in file
        for index, i in enumerate(object_lines):

            # print('------------'+str(index)+'------------')
            # print(data[i])
            # print(data[i+1])
            # print(data[i+2])
            # print(data[i+3])
            # print(data[i+4])
            # print(data[i+5])
            # print(data[i+6])
            # print(data[i+7])
            # print(data[i+8])
            
            # print(data[i+3].split())
            sersic_index = float(data[i+3].split()[1])

            data[i] = data[i].replace(' 1 ', ' 0 ')
            data[i] = data[i].replace(' 1 ', ' 0 ')
            data[i+7] = data[i+7].replace(' 1 ', ' 0 ')
            data[i+8] = data[i+8].replace(' 1 ', ' 0 ')

            if (sersic_index > 8) | (sersic_index < 1):
            
                if fitting_object_index == object_lines.index(i):
                    # print('Replacing sersic index of object # '+str(index)+' with 8 and NOT fixing')
                    sersic_index = '{:.4f}'.format(sersic_index)
                    data[i+3] = data[i+3].replace(str(sersic_index), ' 8 ')
                
                # If the sersic index for a specfic object is too large, get the true values from the catalog
                # and fix this object to these values
                else:

                    calc_radius = float(data[i+2].split()[1])
                    calc_ar = float(data[i+7].split()[1])
                    calc_mag = float(data[i+1].split()[1])

                    gid_specific = int(gidlist_from_files[index])
                    
                    #Estimate effective radius from catalog
                    cat_radius = self.cat[self.cat['id']==gid_specific]['flux_radius'].data[0]
                    
                    #Estimate axis raio from catalog and round to 2 decimal places
                    cat_axis_ratio = '{:.2f}'.format(self.cat[self.cat['id']==gid_specific]['b_image'].data[0]/self.cat[self.cat['id']==gid_specific]['a_image'].data[0])

                    # column names and mag conversion taken from: 
                    # http://cdsarc.u-strasbg.fr/viz-bin/ReadMe/J/ApJS/214/24?format=html&tex=true#sRM3.4
                    if self.flter == '125':
                        cat_mag = (25-2.5*np.log10(self.cat[self.cat['id'] == gid_specific]['f_F125W'])).item()
                    if self.flter == '160':
                        cat_mag = (25-2.5*np.log10(self.cat[self.cat['id'] == gid_specific]['f_F160W'])).item()
                    if self.flter == '850':
                        cat_mag = (25-2.5*np.log10(self.cat[self.cat['id'] == gid_specific]['f_F850LP'])).item()

                    # print('Replacing sersic index of object # '+str(index)+' with 8 and fixing')
                    sersic_index = '{:.4f}'.format(sersic_index)
                    data[i+3] = data[i+3].replace(str(sersic_index), ' 8 ')
                    # data[i+3] = data[i+3].replace(' 1 ', ' 0 ')

                    data[i+2] = data[i+2].replace(str(calc_radius), ' '+str(cat_radius)+' ')
                    data[i+2] = data[i+2].replace(' 1 ', ' 0 ')

                    data[i+7] = data[i+7].replace(str(calc_ar), ' '+str(cat_axis_ratio)+' ')
                    # data[i+7] = data[i+7].replace(' 1 ', ' 0 ')

                    data[i+1] = data[i+1].replace(str(calc_mag), ' '+str(cat_mag)+' ')
                    data[i+1] = data[i+1].replace(' 1 ', ' 0 ')

                    # Let position and PA be free sinzce you are fixing everything else
                    # data[i] = data[i].replace(' 0 ', ' 1 ')
                    # data[i] = data[i].replace(' 0 ', ' 1 ')
                    # data[i+7] = data[i].replace(' 0 ', ' 1 ')
                    # data[i+8] = data[i+8].replace(' 0 ', ' 1 ')

            with open(filename, 'w') as file:
                file.writelines( data )



    def check_units(self, file_path=None, update=False):

        # print('\nChecking fits headers...')
        
        if file_path == None:

            # if self.flter == '850':
                
            #     if self.field == 'N':
            #         fitsfile = 'GOODS-N_acsz_sci_sub.fits'
            #         inputpsf = 'goodsn_3dhst.v4.0.F850LP_psf.fits'
            #         fullfield = 'North'

            #     if self.field == 'S':
            #         fitsfile = 'goodss_3dhst.v4.0.F850LP_orig_sci.fits'
            #         inputpsf = 'goodss_3dhst.v4.0.F850LP_psf.fits'
            #         fullfield = 'South'
                    
            # elif self.flter == '105':
            #     if self.field =='S':
            #         fitsfile = 'goodss-F105W-astrodrizzle-v4.4_drz_sci.fits'
            #         inputpsf = 'none'
            #         psffine = 'none'
            #         fullfield = 'South'

            # else:
            #     if self.field =='S':
            #         fitsfile = 'goodss_3dhst.v4.0.F'+flter+'W_orig_sci.fits'
            #         inputpsf = 'goodss_3dhst.v4.0.F'+flter+'W_psf.fits'
            # #         psffine = 'goodss_3dhst.v4.0.F'+flter+'W_orig_wht.fits'
            #         fullfield = 'South'

            #     if self.field =='N':
            #         fitsfile = 'goodsn_3dhst.v4.0.F'+flter+'W_orig_sci.fits'
            #         inputpsf = 'goodsn_3dhst.v4.0.F'+flter+'W_psf.fits'
            # #         psffine = 'goodsn_3dhst.v4.0.F'+flter+'W_orig_wht.fits'
            #         fullfield = 'North'

            os.chdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter)
            
            files = [self.fitsfile]
            
        else:
            
            fitsfile = file_path
            
            files = [fitsfile]
                
        for file in files:

            data, header = fits.getdata(file, header=True)

            try:
                print(file, 'UNITS: ', header['IMAGEUN'])
                changed = True
            except KeyError:
                print('Units of '+file+' have not been changed to counts')
                changed = False
                continue
                
        if update == False:
            print('You are choosing to not update the units.')
                
        if update == True:

            for file in files:
                data, header = fits.getdata(file, header=True)

                if changed == False:
                    exptime = header['EXPTIME']
                    gain = header['CCDGAIN']

                    data = data*exptime/gain #convert from e/s to counts

                    header['IMAGEUN'] = 'counts'

                    if files.index(file) == 0:
                        savedfile = field+flter+'_sci.fits'
                        fits.writeto(savedfile, data, header)
                        print('Units of '+str(file)+' changed to counts (data * EXPTIME) / GAIN')

                    os.rename(savedfile, file)

                else:
                    print('Image already changed to counts')
                
                

    def run_galfit(self, write_files=True, update=True, sigma_image = None, use_mask=False, delete_segm = True):

        self.check_units(update=update)
        
        #get rid of galfit.number files that might already exist
        numbers = ['01', '02', '03', '04', '05']
        for number in numbers:
            if os.path.exists('galfit.'+number):
                trash = '/usr/local/bin/trash'
                subprocess.call([trash, 'galfit.'+number])
            
        #Field must be 'S' or 'N'
        if not (self.field == 'S' or self.field == 'N'):
            raise Exception('Field must be S or N.')
            
        if not type(self.flter) == str:
            raise Exception('Filter must be string.')
        
        #change to directory where you want to save all outputs
        os.chdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/')
        #assign galfit excecutable
        galfit = '/Users/rosaliaobrien/Bin/galfit'

        if write_files == True:
            #write each feedme file
            nearby_objects = self.write_feedme(str(self.gid), sigma_image = sigma_image, use_mask = use_mask)
                
        else:
            os.rename('outputs/'+str(self.gid)+'/'+str(self.gid)+'.feedme', str(self.gid)+'.feedme')
            
        if os.path.exists('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/failed_first_run.txt'):
            
            try:
            
                # get rid of gid if it's in failed_first_run.txt
                rewrite_failed_first_run = pd.read_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/failed_first_run.txt', 
                            header=None)
                rewrite_list1 = rewrite_failed_first_run[0].tolist()

                if self.gid in rewrite_list1:
                    print('Removing '+str(self.gid)+' from failed_first_run.txt')
                    rewrite_list1.remove(self.gid)

                new_failed_first_run = pd.DataFrame(rewrite_list1)
                new_failed_first_run.to_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/failed_first_run.txt', 
                            header=None, index=False)
                
            except:
                pass

        if os.path.exists('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/failed_second_run.txt'):
            
            try:
            
                # get rid of gid if it's in failed_second_run.txt
                rewrite_failed_second_run = pd.read_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/failed_second_run.txt', 
                            header=None)
                rewrite_list2 = rewrite_failed_second_run[0].tolist()

                if gid in rewrite_list2:
                    print('Removing '+str(self.gid)+' from failed_second_run.txt')
                    rewrite_list2.remove(self.gid)

                new_failed_second_run = pd.DataFrame(rewrite_list2)
                new_failed_second_run.to_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/failed_second_run.txt', 
                            header=None, index=False)
                
                
            except:
                pass

        print('Running GALFIT...')
        #run galfit for each feedme file!
        subprocess.run([galfit+' -skyped {:.2f} -skyrms {:.2f} '.format(self.sky_value, self.rms_value)+str(self.gid)+'.feedme'],shell = True)
        
        #if galfit.01 wasnt produced, GALFIT failed first run
        if not os.path.exists('galfit.01'):
            print('!!!OBJECT FAILED FIRST RUN!!!')
            
            #make list of failed objects
            if not os.path.exists('failed_first_run.txt'):
                text = str(self.gid)
                np.savetxt('failed_first_run.txt', np.array(text).reshape(1, ), fmt='%s')

            else:
                with open('failed_first_run.txt', "a") as myfile:
                    myfile.write('\n'+str(self.gid))

        #if galfit didnt fail first run and also didnt do a two component model, continue
        if os.path.exists('galfit.01'):
            
            print('Rewriting galfit.01...')
            self.rewrite_galfit01()
            print('Rerunning GALFIT...')
            subprocess.run([galfit+' -skyped {:.2f} - skyrms {:.2f} '.format(self.sky_value, self.rms_value)+'galfit.01'], shell = True)

            
            if not os.path.exists('galfit.02'):
                print('!!!OBJECT FAILED SECOND RUN!!!')
                
                #make list of failed objects
                if not os.path.exists('failed_second_run.txt'):
                    text = str(self.gid)
                    np.savetxt('failed_second_run.txt', np.array(text).reshape(1, ), fmt='%s')

                else:
                    with open('failed_second_run.txt', "a") as myfile:
                        myfile.write('\n'+str(self.gid))

        self.organize_files(delete_segm = delete_segm)
        
        print('DONE!')
        print('\n')
                
                
                
    def normalize_colors(self, sci, model, resid):
        num_add = abs(min(sci.flatten()))
        denom = max(sci.flatten())
        norm_arr = (sci+num_add)/denom
        norm_model = (model+num_add)/denom
        norm_resid = (resid+num_add)/denom
        return norm_arr, norm_model, norm_resid

    def plot(self, norm_colors = True, show_object = False, show_plot=False, 
        save_fig=True, used_mask=False, saveas_png = False):
            
        if self.flter == '160':
            fullfilter = 'F160W'
            
        if self.flter == '125':
            fullfilter = 'F125W'
            
        if self.flter == '850':
            fullfilter = 'F850LP'
        
        os.chdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+self.fullfield+'/'+self.flter+'/')
            
        # print(str(gid)+'...')

        file = 'outputs/fits/'+str(self.gid)+'.fits'
    #         savefig = 'outputs/pngs/output_gid'+str(gid)+'.png'

        if not os.path.exists('outputs/pngs'):
                os.mkdir('outputs/pngs')
            
        image = fits.getdata(file, ext=1)
        model = fits.getdata(file, ext=2)
        residual = fits.getdata(file, ext=3)

        fig = plt.figure(figsize=(20, 7))
        # plt.rcParams.update({'font.size': 35})

        gs = gridspec.GridSpec(1, 3, wspace=0.01, hspace=0)

        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1], sharey = ax1)
        ax3 = plt.subplot(gs[0,2], sharey = ax1)

        ax1.set_title('G'+self.field+'D-'+str(self.gid)+' ('+fullfilter+')', size = 32)
                          
        ax2.set_title('GALFIT Model', size = 35)
        ax3.set_title('Residual', size = 35)

        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])

        if norm_colors == False:
            #get min and max values for images
            # maskimg2 = np.ravel(image)
            # maskimg2 = maskimg2[maskimg2 < 0.5] 
            vmin = np.percentile(image,4)
            vmax = np.percentile(image,96)

        if norm_colors == True:
            image, model, residual = self.normalize_colors(image, model, residual)
            vmin = np.percentile(image,4)
            vmax = np.percentile(image,96)

        ax1.imshow(image, vmin=vmin, vmax=vmax, cmap='Greys_r', origin='lower', aspect="auto")
        ax2.imshow(model, vmin=vmin, vmax=vmax, cmap='Greys_r', origin='lower', aspect="auto")
        ax3.imshow(residual, vmin=vmin, vmax=vmax, cmap='Greys_r', origin='lower', aspect="auto")
        
        if show_object == True:

            x, y, coord = self.get_position(self.gid)
            nearby_objects = self.get_nearby_objects()
            xmin, xmax, ymin, ymax = self.get_boundaries(nearby_objects)

            ax1.scatter(x-xmin, y-ymin, marker='+', color='red', s=300)
            ax2.scatter(x-xmin, y-ymin, marker='+', color='red', s=300)

            ax1.set_xlim(left=x-xmin-200, right=x-xmin+200)
            ax2.set_xlim(left=x-xmin-200, right=x-xmin+200)
            ax3.set_xlim(left=x-xmin-200, right=x-xmin+200)
            ax1.set_ylim(bottom=y-ymin-200, top=y-ymin+200)
            ax2.set_ylim(bottom=y-ymin-200, top=y-ymin+200)
            ax3.set_ylim(bottom=y-ymin-200, top=y-ymin+200)
                
        if save_fig == True:
            savefig = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/plots/residuals/pdfs/output_gid'+str(self.gid)+'_'+self.flter+'.pdf'
            plt.savefig(savefig, bbox_inches = 'tight')
            
            if saveas_png == True:
                savefig = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/plots/residuals/pngs/output_gid'+str(self.gid)+'_'+self.flter+'.png'
                plt.savefig(savefig, bbox_inches = 'tight')
            
        if show_plot == True:
            plt.show()
            
        if show_plot == False:
            plt.close()

    def update_csv(file, update_gid, column = None, cell_value = None):
        data = pd.read_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+file)
        for index, gid in zip(data.index, data['gid']):
            if gid == update_gid:
                data.at[index, column] = cell_value
                
        data.to_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+file, index=False)
        
        return data

    def make_temp_df(df, tempgids):
        dflist = []
        for index, gid in zip(df.index, df['gid']):
            if gid in tempgids:
                dflist.append(index)
                
        return df.loc[dflist]

    def make_segm(self, ignore = [00000]):

        nearby_objects = self.get_nearby_objects()
        
        # print('Making mask...')
        
        # if flter == '850':
        #     if field == 'N':
        #         # datafile = 'GOODS-N_acsz_sci_sub.fits'
        #         segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodsn_3dhst.v4.0.F160W_seg.fits'
        #     if field == 'S':
        #         # datafile = 'goodss_3dhst.v4.0.F850LP_orig_sci.fits'
        #         segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/goodss_3dhst.v4.0.F160W_seg.fits'
                
        # else:
            
        if self.field =='S':
            # datafile = 'goodss_3dhst.v4.0.F'+flter+'W_orig_sci.fits'
            segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/GOODS-S/goodss_3dhst.v4.0.F160W_seg.fits'

        if self.field =='N':
            # datafile = 'goodsn_3dhst.v4.0.F'+flter+'W_orig_sci.fits'
            segmfile = '/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/og_images/GOODS-N/goodsn_3dhst.v4.0.F160W_seg.fits'
                
                
        # data = fits.open(datafile)[0].data
        # x,y,coord = get_position(field, gid)
        # cutout = Cutout2D(data, (x, y), (300, 300))
        
        segmdata = fits.open(segmfile)[0].data
        # segmcutout = Cutout2D(segmdata, (x, y), (300, 300))
        segmdata[segmdata == self.gid] = 0

        if self.nearby_fitting_region != 'masked':
            for obj in nearby_objects:
                if not obj in ignore:
                    segmdata[segmdata == obj] = 0
        
        # cutout.data[(segmcutout.data != gid) & (segmcutout.data > 0)] = 0
        if not os.path.exists('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/segm_maps/'):
            os.mkdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/segm_maps/')
        
        hdu = fits.PrimaryHDU(segmdata, header = fits.open(segmfile)[0].header)
        hdu.writeto('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/segm_maps/'+str(self.gid)+'_segm.fits', overwrite = True)