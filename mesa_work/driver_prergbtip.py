from __future__ import print_function
import re
import numpy as np
import os
import sys
from astropy.io import ascii
import glob
from time import sleep
import h5py
import pandas as pd

def set_gyre_inlist(inlist_path, outlist_path, inputFileName, summaryFileName, freqMinRadial, freqMaxRadial, freqMinNonRadial, freqMaxNonRadial):
    # reads in the template inlist and writes out a new inlist with the 
    # parameters set appropriately
    try:
        inlist = open(inlist_path,'r')
        outlist = open(outlist_path,'w')
    except:
        sleep(2.0)
        inlist = open(inlist_path,'r')
        sleep(2.0)
        outlist = open(outlist_path,'w')


    for line in inlist.read().split('\n'):

        first = line.split()
        if len(first)>0:
            if first[0] != '!':

                if re.search("file = 'spb.mesa'",line):
                    line = "\tfile = '%s'  !set by driver.py" % inputFileName

                if re.search("summary\_file = 'summary.txt'",line):
                    line = "\tsummary_file = '%s'  !set by driver.py" % summaryFileName

                if re.search("freq\_min\_radial = 100",line):
                    line = "\tfreq_min={:0.5f}  !set by driver.py".format(freqMinRadial)

                if re.search("freq\_max\_radial = 1800",line):
                    line = "\tfreq_max={:0.5f}  !set by driver.py".format(freqMaxRadial)

                if re.search("freq\_min\_nonradial = 100",line):
                    line = "\tfreq_min={:0.5f}  !set by driver.py".format(freqMinNonRadial)

                if re.search("freq\_max\_nonradial = 1800",line):
                    line = "\tfreq_max={:0.5f}  !set by driver.py".format(freqMaxNonRadial)

        print(line,file=outlist)

    outlist.close()
    inlist.close()


class readTable:
    '''
    A parent class to be wrapped by other class, in order to read in such as mesa history file.
    These files have very similar structures, typically header+table.

    '''
    def __init__(self, filepath, verbose=True):
        '''
        Args:
            filepath: the path of .history file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            track: a structure numpy array containing the evol track
            colnames: a tuple containing column names

        '''  
        self.filepath = filepath
        if verbose: 
            print('Processing :', self.filepath)
        return
    
    def readFile(self, filepath, headerNameLine=1, headerDataLine=2, tableHeaderLine=6):
        '''
        Reads in a file.
        '''

        with open(self.filepath) as f:
            content = [line.split() for line in f]
        header = {content[headerNameLine-1][i]:content[headerDataLine-1][i] for i in range(len(content[headerNameLine-1]))}
        table = np.genfromtxt(self.filepath, skip_header=tableHeaderLine-1, names=True)
        colnames = table.dtype.names

        return header, table, colnames



class history(readTable):
    '''

    A class to read mesa history files, store the data within, and offer useful routines (?).

    '''
    
    def __init__(self, filepath, verbose=True, ifReadProfileIndex=False):
        '''
        Args:
            filepath: the path of .history file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            track: a structure numpy array containing the evol track
            colnames: a tuple containing column names
            profileIndex: a structured array containing the map between model_number and profile_number

        '''
        super().__init__(filepath, verbose)
        self.header, self.track, self.colnames = self.readFile(filepath, headerNameLine=2, headerDataLine=3, tableHeaderLine=6)
        
        if ifReadProfileIndex:
            self.profileIndex = self.read_profile_index()
        return
    
    def read_profile_index(self):
        '''
        Reads in the profile.index file
        '''
        filepath = self.filepath.split('.history')[0] + 'profile.index'
        profileIndex = np.genfromtxt(filepath, skip_header=1, names=('model_number', 'priority', 'profile_number'))
        return profileIndex


class profile(readTable):
    '''

    A class to read mesa history files, store the data within, and offer useful routines.

    '''

    
    def __init__(self, filepath, verbose=True):
        '''
        Args:
            filepath: the path of *profile*.data file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            profile: a structure numpy array containing the structure profile
            colnames: a tuple containing column names

        '''
        super().__init__(filepath, verbose)
        self.header, self.profile, self.colnames = self.readFile(filepath, headerNameLine=2, headerDataLine=3, tableHeaderLine=6)

        return


class sums(readTable):
    '''

    A class to read gyre mode summary file, store the data within, and offer useful routines.

    '''

    
    def __init__(self, filepath, verbose=True):
        '''
        Args:
            filepath: the path of .sums file.
            verbose: whether print info to device. default: True.

        Attributes:
            header: a dictionary containing history file header
            modeSummary: a structure numpy array containing the summary table
            colnames: a tuple containing column names

        '''
        super().__init__(filepath, verbose)
        self.header, self.modeSummary, self.colnames = self.readFile(filepath, headerNameLine=3, headerDataLine=4, tableHeaderLine=6)

        return

def set_mesa_inlist(index, inlist_path, outlist_path, 
                    mass, Xinit, Yinit, Zinit, amlt, 
                    fov_shell, fov0_shell, 
                    fov_core, fov0_core, 
                    ifsetfinalmodel):
    # reads in the template inlist and writes out a new inlist with the 
    # parameters set appropriately
    
    try:
        inlist = open(inlist_path,'r')
        outlist = open(outlist_path,'w')
    except:
        sleep(2.0)
        inlist = open(inlist_path,'r')
        sleep(2.0)
        outlist = open(outlist_path,'w')


    for line in inlist.read().split('\n'):

        first = line.split()
        if len(first)>0:
            if first[0] != '!':

                if re.search('initial\_mass',line):
                    line = "\tinitial_mass = %g  !set by driver.py" % mass

                if re.search('initial\_z =',line):
                    line = "\tinitial_z = %g  !set by driver.py" % Zinit
                
                if re.search('Zbase =',line):
                    line = "\tZbase =%g  !set by driver.py" % Zinit

                if re.search('initial\_y',line):
                    line = "\tinitial_y = %g  !set by driver.py" % Yinit

                if re.search('mixing\_length\_alpha',line):
                    line = "\tmixing_length_alpha = %g  !set by driver.py" % amlt

                if fov_shell >0:
                    if re.search('overshoot\_scheme\(1\)',line):
                        line = "\tovershoot_scheme(1) = '%s' !set by driver.py " % "exponential" 
                    if re.search('overshoot\_zone\_type\(1\)',line):
                        line = "\tovershoot_zone_type(1) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_zone\_loc\(1\)',line):
                        line = "\tovershoot_zone_loc(1) = '%s' !set by driver.py " % "shell" 
                    if re.search('overshoot\_bdy\_loc\(1\)',line):
                        line = "\tovershoot_bdy_loc(1) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_f\(1\)',line):
                        line = "\tovershoot_f(1) = %g  !set by driver.py" %fov_shell
                    if re.search('overshoot\_f0\(1\)',line):
                        line = "\tovershoot_f0(1) = %g !set by driver.py" %fov0_shell

                else:
                    if re.search('overshoot\_scheme\(1\)',line):
                        line = "\t!overshoot_scheme(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_type\(1\)',line):
                        line = "\t!overshoot_zone_type(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_loc\(1\)',line):
                        line = "\t!overshoot_zone_loc(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_bdy\_loc\(1\)',line):
                        line = "\t!overshoot_bdy_loc(1) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_f\(1\)',line):
                        line = "\t!overshoot_f(1) = %g  !set by driver.py" % 0.
                    if re.search('overshoot\_f0\(1\)',line):
                        line = "\t!overshoot_f0(1) = %g !set by driver.py" % 0.               

                if fov_core >0:
                    if re.search('overshoot\_scheme\(2\)',line):
                        line = "\tovershoot_scheme(2) = '%s' !set by driver.py " % "exponential" 
                    if re.search('overshoot\_zone\_type\(2\)',line):
                        line = "\tovershoot_zone_type(2) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_zone\_loc\(2\)',line):
                        line = "\tovershoot_zone_loc(2) = '%s' !set by driver.py " % "core" 
                    if re.search('overshoot\_bdy\_loc\(2\)',line):
                        line = "\tovershoot_bdy_loc(2) = '%s' !set by driver.py " % "any" 
                    if re.search('overshoot\_f\(2\)',line):
                        line = "\tovershoot_f(2) = %g  !set by driver.py" %fov_core
                    if re.search('overshoot\_f0\(2\)',line):
                        line = "\tovershoot_f0(2) = %g !set by driver.py" %fov0_core

                else:
                    if re.search('overshoot\_scheme\(2\)',line):
                        line = "\t!overshoot_scheme(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_type\(2\)',line):
                        line = "\t!overshoot_zone_type(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_zone\_loc\(2\)',line):
                        line = "\t!overshoot_zone_loc(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_bdy\_loc\(2\)',line):
                        line = "\t!overshoot_bdy_loc(2) = '%s' !set by driver.py " % "" 
                    if re.search('overshoot\_f\(2\)',line):
                        line = "\t!overshoot_f(2) = %g  !set by driver.py" % 0.
                    if re.search('overshoot\_f0\(2\)',line):
                        line = "\t!overshoot_f0(2) = %g !set by driver.py" % 0.   

                # smass = 'm{:03.0f}'.format(mass*100)
                header = 'index{:06.0f}'.format(index)

                if re.search('star\_history\_name',line):
                    output_history_name = header + '.history'
                    line = "\tstar_history_name = '%s'  !set by driver.py" % output_history_name

                if re.search('profile\_data\_prefix',line):
                    output_profile = header + 'profile'
                    line = "\tprofile_data_prefix = '%s' !set by driver.py " % output_profile 
                
                if re.search('profiles\_index\_name =', line):
                    output_index = header + 'profile.index'
                    line = "\tprofiles_index_name = '%s' !set by driver.py " % output_index 

                if ifsetfinalmodel:
                    if re.search('save\_model\_filename =', line):
                        output_final_model = header + 'final.mod'
                        line = "\tsave_model_filename = '%s' !set by driver.py " % output_final_model 

        print(line,file=outlist)

    outlist.close()
    inlist.close()

# read in input params
tracks = ascii.read('coarse_grid_input_params.txt', delimiter=',')

# get from input which track to run for this simulation
index_this_run = int(sys.argv[1])
idx = (tracks['index'] >= index_this_run) & (tracks['index'] < (index_this_run+1))

for track in tracks[idx]:

    # get input params for this simulation
    index, mass = track['index'], track['star_mass']
    Xinit, Yinit, Zinit, amlt = track['Xinit'], track['Yinit'], track['Zinit'], track['amlt']
    fov_shell, fov0_shell, fov_core, fov0_core = track['fov_shell'], track['fov0_shell'], track['fov_core'], track['fov0_core']
    mh = np.log10(Zinit/Xinit) - np.log10(0.0134/0.7381) # asplund 09 current solar abundance scale
    
    output_history_name = 'index{:06.0f}'.format(index)+'.history'
    output_final_model_name = 'index{:06.0f}final.mod'.format(index)
    output_profile_index_name = 'index{:06.0f}profile.index'.format(index)


    # if os.path.exists('../pre_rgb_tip/finalmodels/'+output_final_model_name): 
    #     print('File exists, skipping ', output_final_model_name)
    #     continue
    # else:
    #     print('Now calculating ', output_final_model_name)
    print('Now calculating ', output_final_model_name)

    # # # Step 1: modify inlist and run MESA.

    set_mesa_inlist(index, 'inlist_template', 'inlist', 
                    mass, Xinit, Yinit, Zinit, amlt, 
                    fov_shell, fov0_shell, fov_core, fov0_core, True)
    #os.system('\\rm -r LOGS; \\rm -r png; \\rm -r photos')

    print('------ MESA start ------')
    os.system('./rn1 > mesa_output_index{:06.0f}.txt'.format(index))
    print('------ MESA done ------')

    # if os.path.exists(output_final_model_name):
    #     os.system('mv {:s} ../pre_rgb_tip/finalmodels/'.format(output_final_model_name))

    #     # # List all files in the home directory
    #     # files = glob.glob(os.path.expanduser("photos/x*"))
    #     # if len(files) != 0:
    #     #     # Sort by modification time (mtime) descending
    #     #     latest_binary_photo = sorted(files, key=lambda t: -os.stat(t).st_mtime)[0]
    #     #     os.system('mv {:s} ../finalmodels/{:s}.x'.format(latest_binary_photo, output_final_model_name))
    # else:
    #     os.system('touch ../pre_rgb_tip/finalmodels/{:s}'.format(output_final_model_name))


    # # # Step 2: collect FGONGs and run GYRE. Only run if history and index file exists.
    if (os.path.exists('LOGS/'+output_history_name) & os.path.exists('LOGS/'+output_profile_index_name)):

        # calculate oscillation modes

        h = history('LOGS/'+output_history_name, ifReadProfileIndex=True)
        track = h.track # 'delta_nu', 'nu_max'
        profileIndex = h.profileIndex # 'model_number', 'priority', 'profile_number'

        fgongPaths = [f for f in os.listdir('LOGS/') if (f.endswith('.FGONG') & f.startswith('index{:06.0f}'.format(index)))]

        for fgongPath in fgongPaths:
            # if os.path.exists(summaryFileName): continue

            profileNo = int(fgongPath.split('profile')[-1].split('.data')[0])
            modelNo = profileIndex['model_number'][profileIndex['profile_number'] == profileNo][0]
            delta_nu = track['delta_nu'][track['model_number'] == modelNo][0]
            nu_max = track['nu_max'][track['model_number'] == modelNo][0]
            delta_Pg = track['delta_Pg'][track['model_number'] == modelNo][0]
            Teff = 10.0**track['log_Teff'][track['model_number'] == modelNo][0]
            
            # correct nu_max with the metallicity dependence
            p = np.array([-0.00718002, -0.05279003,  0.9892853 ])
            nu_max = nu_max * np.polyval(p, mh)
            
            k, b = 0.9638, -1.7145
            width = np.exp(k*np.log(nu_max) + b)
            freqMinRadial = nu_max - width*4
            freqMaxRadial = nu_max + width*4
            freqMinNonRadial = nu_max - width*4
            freqMaxNonRadial = nu_max + width*4
            # freqMin = nu_max - 10.*delta_nu
            # freqMax = nu_max + 10.*delta_nu
            if freqMinRadial <= 0.: freqMinRadial = 0.0001
            if freqMinNonRadial <= 0.: freqMinNonRadial = 0.0001


            inputFileName = 'LOGS/'+fgongPath
            summaryFileName = 'LOGS/'+fgongPath+'.sum'
            if (Teff<5500) & ((delta_nu/(nu_max**2.0 * delta_Pg * 1e-6))>3.5):
                inlist_path = 'gyre_template_pi.in'
            else:
                inlist_path = 'gyre_template_mixed.in'
            outlist_path = 'gyre.in'
            set_gyre_inlist(inlist_path, outlist_path, inputFileName, summaryFileName, freqMinRadial, freqMaxRadial, freqMinNonRadial, freqMaxNonRadial)

            # print('------ GYRE start ------')
            os.system('$GYRE_DIR/bin/gyre gyre.in >> mesa_output_index{:06.0f}.txt'.format(index))
            # print('------ GYRE done ------')



    # # # Step 3: create a .h5 file to store history and frequencies. Only run if history file exists.
    if os.path.exists('LOGS/'+output_history_name):

        # os.system('mv LOGS/{:s} ../pre_rgb_tip/history/'.format(output_history_name))
        # os.system('mv LOGS/{:s} ../pre_rgb_tip/history/'.format(output_profile_index_name))
        # if (not os.path.exists('../pre_rgb_tip/freqs/index{:06.0f}/'.format(index))): os.mkdir('../pre_rgb_tip/freqs/index{:06.0f}/'.format(index))
        # os.system('mv LOGS/*.sum ../pre_rgb_tip/freqs/index{:06.0f}/'.format(index))


        # # os.system('mv LOGS/*.FGONG ../freqs/index{:06.0f}/'.format(index))
        # os.system('rm LOGS/*.FGONG')
        # os.system('rm LOGS/*profile*.data')


        # # read in models
        if os.path.exists('LOGS/'+output_profile_index_name):
            h = history('LOGS/'+output_history_name, ifReadProfileIndex=True)
            ifprofile = True 
        else:
            h = history('LOGS/'+output_history_name, ifReadProfileIndex=False)
            ifprofile = False 


        track = h.track 
        _, idx = np.unique(track['model_number'], return_index=True)
        track = track[idx]

        table = pd.DataFrame(track)

        if ifprofile:
            table = table.merge(pd.DataFrame(h.profileIndex), on='model_number', how='left')

        # # append grid initial parameters as new columns
        table['index'], table['Xinit'], table['Yinit'], table['Zinit'], table['amlt'] = index, Xinit, Yinit, Zinit, amlt
        table['fov_shell'], table['fov0_shell'], table['fov_core'], table['fov0_core'] = fov_shell, fov0_shell, fov_core, fov0_core


        # # append phase as a new column
        phase = np.zeros(len(table))
        turnoffidx = table['center_h1'] < 1.e-7
        msidx = (10.0**(table['log_Lnuc']-table['log_L'])>0.99) & (~turnoffidx)
        pmsidx = (10.0**(table['log_Lnuc']-table['log_L'])<=0.99) & (~turnoffidx)
        sgidx = (turnoffidx) & (table['nu_max']>=300)
        rgbidx = (turnoffidx) & (table['nu_max']<300)
        hebidx = (turnoffidx) & ((table['center_he4']+table['center_he3']) <0.95)
        tfidx = (turnoffidx) & (~sgidx) & (~rgbidx) & (~hebidx)

        phase[pmsidx] = -1
        phase[msidx] = 0
        phase[sgidx] = 1
        phase[rgbidx] = 2
        phase[hebidx] = 3
        phase[tfidx] = -2
        table['phase'] = phase

        # # append log properties
        table['luminosity'] = 10.0**table['log_L']
        table['radius'] = 10.0**table['log_R']
        table['Teff'] = 10.0**table['log_Teff']

        # # append seismic scaling quantities
        Dnu_sun, numax_sun, Teff_sun = 135.1, 3090., 5777.
        table['delta_nu_scaling'] = table['star_mass']**0.5 * table['radius']**-1.5 * Dnu_sun
        table['numax_scaling'] = table['star_mass'] * table['radius']**-2.0 * (table['Teff']/Teff_sun)**-0.5 * numax_sun
        
        # # append surface quantities
        Zsun, Xsun = 0.0134, 0.7381 # 0.0134, 0.7381, a09 # 0.0169, 0.7345, gs98
        table['FeH'] = np.log10((1-table['surface_h1']-table['surface_he4']-table['surface_he3'])/table['surface_h1']) - np.log10(Zsun/Xsun)


        # # # assign a prior
        # age = np.array(table['star_age'])
        # prior = np.concatenate([[0],(age[2:]-age[:-2])/2.,[0]])
        # prior[~np.isfinite(prior)] = 0. 
        # prior = prior/np.sum(prior)
        # table['prior'] = prior


        sumPaths = ['LOGS/'+f for f in os.listdir('LOGS/') if f.endswith('.FGONG.sum')]
        sumDirs = np.array([f.split('/index')[0]+'/' for f in sumPaths])
        sumNames = np.array([f.split('/')[-1] for f in sumPaths])

        if ifprofile:
            # # set a seismic flag
            table['flag_seismo'] = np.array(np.isfinite(table['profile_number']), dtype=int)

            # # read in radial mode frequencies
            seismicCols = ['l', 'n_p', 'n_g', 'n_pg', 'E_p', 'E_g', 'E_norm', 'freq']
            # seismicCols = ['l', 'n_p', 'n_g', 'n_pg', 'E_norm', 'freq']
            seismicData = [[] for i in range(len(seismicCols))]

            for imod, mod in table.loc[:,:].iterrows():
                profileIndex = mod['profile_number']
                if (not np.isfinite(profileIndex)):
                    for i in range(len(seismicData)):
                        seismicData[i].append(np.nan)
                else:
                    sumFile = 'index{:06.0f}profile{:0.0f}.data.FGONG.sum'.format(index, profileIndex)
                    if len(sumDirs[sumNames==sumFile])==0:
                        for i in range(len(seismicData)):
                            seismicData[i].append(np.nan)
                        table.loc[imod, 'flag_seismo'] = 0
                        table.loc[imod, 'profile_number'] = np.nan
                    else:
                        sumPath = sumDirs[sumNames==sumFile][0] + sumFile
                        s = sums(sumPath, verbose=False)
                        for i in range(len(seismicData)-1):
                            seismicData[i].append(np.array([s.modeSummary[seismicCols[i]]]))
                        seismicData[-1].append(np.array([s.modeSummary['Refreq']]))
        else:
            # # set a seismic flag
            table['flag_seismo'] = np.zeros(len(table), dtype=int)

        # #  write out the table
        with h5py.File('index{:06.0f}.history.h5'.format(index), 'w') as h5f:
            for col in table.columns:
                h5f.create_dataset(col, data=table[col])
            if ifprofile:
                for imod, mod in table.iterrows():
                    profileIndex = mod['profile_number']
                    if (not np.isfinite(profileIndex)): continue
                    for i in range(len(seismicCols)):
                        h5f.create_dataset('profile{:0.0f}/{:s}'.format(profileIndex, seismicCols[i]), data=seismicData[i][imod])

