import os
import glob
import pickle
import pandas
import shutil
import argparse
import numpy as np
from astroquery import sha
from astropy import coordinates as coord
from astropy import units as u
from collections import Counter

from tools import query_exo_archive

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process (default = all)"
    parser.add_argument("-t", "--target", help=help_, type=str, default="all")

    help_ = "Directory containing data (default = DATA/)"
    parser.add_argument("-d", "--datadir", help=help_, default="DATA/", type=str)

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    # directory to save data
    DATADIR = args.datadir

    data = query_exo_archive()
    
    # only deal with transiting exoplanets
    pmask = (data['pl_orbper']<120) & (data['tran_flag']==1)
    data = data[pmask]
    data = data.reset_index(drop=True)

    #initialize new column for image counts
    data['spitzer count'] = 0.0

    # keys to save from spitzer archive query
    # https://irsa.ipac.caltech.edu/onlinehelp/heritage/api.html#output
    keys = ['modedisplayname','wavelength','exposuretime','scet','filesize', 'accessUrl']

    # loop through targets
    for i in range(data.shape[0]):

        # for resuming with a %paste
        if data.at[i,'spitzer count'] != 0:
            continue 

        # ignore KOIs
        if data.iloc[i].hostname.split('-')[0] == 'KOI':
            data.at[i,'spitzer count'] = -1
            continue
        
        if args.target not in data.iloc[i].hostname:
            continue

        pname = data.iloc[i].pl_name
        star = ''.join(pname.split(' ')[:-1])
        print("Querying Spitzer Archive for: ",pname)
        #if pname
        # query archive
        try:
            query = sha.query(ra=data.iloc[i].ra, dec=data.iloc[i].dec, size=0.01, verbosity=3, dataset=1)
            print('query finished {}'.format(len(query)))

            # alloc arrays
            alldata = {'aor':[]}
            for k in keys:
                alldata[k] = []

            # store info for each spitzer db object
            for j in range(len(query)):

                # only accept valid images
                if query[j]['primaryfield'] == 1 and query[j]['hasAccess'].strip() == 'true' and query[j]['filetype'].strip() == 'Image':
                    if query[j]['wavelength'].strip() in ['IRAC 3.6um','IRAC 4.5um','IRAC 5.8um','IRAC 8.0um']:
                        # store information
                        for k in keys:
                            alldata[k].append( query[j][k] )
                        alldata['aor'].append( query[j]['reqkey'] )

            # create directory to save data            
            if not os.path.exists(DATADIR+star):
                os.mkdir(DATADIR+star)
            
            pickle.dump(alldata, open(DATADIR+star+"/filelist.pkl",'wb'))
            data.at[i,'spitzer count'] = len(alldata[k])

        except Exception as e:
            print("something failed",e)
            data.at[i,'spitzer count'] = -1