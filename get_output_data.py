import os
import pandas as pd

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 2

def get_output_data_lists(fullfield, flter, gidlist, fitting_regions):

    print(fullfield, flter)

    all_data = []
            
    for gid, fitting_reg in zip(gidlist, fitting_regions):
        
        # print(str(gid)+'...')
        
        os.chdir('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+fullfield+'/'+flter+'/outputs/'+str(gid))

        failed_second_run = pd.read_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_GALFIT/'+fullfield+'/'+flter+'/failed_second_run.txt', header=None)
        
        try:
            galfit_file = 'galfit.02'
            #find number of objects in file based on length of file
            number_of_objects=(file_len(galfit_file)-41)/13
        except:
            print('Not using galfit.02')
        
            try:
                #if object failed second run, use galfit.01 instead
                for index, fgid in failed_second_run.iterrows():
                    if gid == fgid[0]:
                        print('Using galfit.01 for id'+str(gid))
                        galfit_file = 'galfit.01'
                        number_of_objects=(file_len(galfit_file)-41)/13
            except:
                pass
        
        #read feedne file
        with open(str(gid)+'.feedme') as feedmefile:
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
            if j == gid:
                fitting_object_index = k

        #first object always on line 40
        object_lines = [40]
        comp = 40
        #make list of lines where data on objects begin
        for i in range(int(number_of_objects)-1):
            comp = comp+13
            object_lines.append(comp)

        data_list = [gid]
        
        with open(galfit_file, 'r') as file:
            data = file.readlines()

        #fix PA, position, and axis ratio for each object in file
        for index, i in enumerate(object_lines):

            if fitting_object_index == object_lines.index(i):
                
                mag = float(data[i+1].split()[1])
                rad = float(data[i+2].split()[1])
                sersic_index = float(data[i+3].split()[1])
                axisratio = float(data[i+7].split()[1])
                PA = float(data[i+8].split()[1])
                
                data_list.append(mag)
                data_list.append(rad)
                data_list.append(sersic_index)
                data_list.append(axisratio)
                data_list.append(PA)

        all_data.append(data_list)
        
    return all_data

def get_output_data(fullfield, gidlist, fitting_regions, diff_fitting_regions = False):
    
    # print('Getting output data...')
    
    if diff_fitting_regions == True:
        if not type(fitting_regions) == list:
            raise Exception('fitting_region must be a list of fitting regions [F160W, F125W, 850LP]')

        data160 = get_output_data_lists(fullfield, '160', gidlist, fitting_regions[0])
        data125 = get_output_data_lists(fullfield, '125', gidlist, fitting_regions[1])
        data850 = get_output_data_lists(fullfield, '850', gidlist, fitting_regions[2])
        
    else:
        data160 = get_output_data_lists(fullfield, '160', gidlist, fitting_regions)
        data125 = get_output_data_lists(fullfield, '125', gidlist, fitting_regions)
        data850 = get_output_data_lists(fullfield, '850', gidlist, fitting_regions)
    
    column_names = ['gid', 'mag', 'radius', 'sersic_index', 'axis_ratio', 'PA']
    
    df160 = pd.DataFrame(data=data160, columns = column_names)
    df125 = pd.DataFrame(data=data125, columns = column_names)
    df850 = pd.DataFrame(data=data850, columns = column_names)
    
    df160.to_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_galfit/'+fullfield+'/F160W_output_data.csv', index=False)
    df125.to_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_galfit/'+fullfield+'/F125W_output_data.csv', index=False)
    df850.to_csv('/Users/rosaliaobrien/research/GALFIT_CLEAR/running_galfit/'+fullfield+'/F850LP_output_data.csv', index=False)
    
    # print('csv files saved to running_galfit/')
    
    return df160, df125, df850


