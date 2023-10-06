import os
import json
import pickle
import argparse
from astroquery import sha

def parse_args():
    parser = argparse.ArgumentParser()
    
    help_ = "Meta data, json file with {target:[aor..]} "
    parser.add_argument("-i", "--input", help=help_, type=str)

    help_ = "Directory containing data (default = DATA/)"
    parser.add_argument("-d", "--datadir", help=help_, default="DATA/", type=str)

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    # meta data 
    data = json.load(open(args.input,"r"))

    keys = ['modedisplayname','wavelength','exposuretime','scet','filesize', 'accessUrl']

    targets = list(data.keys())

    for t in targets:
        
        # alloc arrays
        alldata = {'aor':[]}
        for k in keys:
            alldata[k] = []

        for i in range(len(data[t])):

            try:
                print("Querying {} @ AOR:{}".format(t,data[t][i])) 
                query = sha.query(reqkey=data[t][i], verbosity=3, dataset=1)
                print(" {}".format(len(query)))
                # store info for each spitzer db object
                for j in range(len(query)):

                    # only accept valid images
                    if query[j]['primaryfield'] == 1 and query[j]['hasAccess'].strip() == 'true' and query[j]['filetype'].strip() == 'Image':

                        if query[j]['wavelength'].strip() in ['IRAC 3.6um','IRAC 4.5um','IRAC 5.8um','IRAC 8.0um']:

                            # no duplicates
                            if query[j]['accessUrl'] not in alldata['accessUrl']:
                                alldata['aor'].append(data[t][i])
                                # store information
                                for k in keys:
                                    alldata[k].append( query[j][k] )

                savedir = os.path.join(args.datadir, ''.join(t.split(' ')))
                if not os.path.exists(savedir):
                    os.mkdir(savedir)

                pickle.dump(alldata, open(savedir+"/filelist.pkl",'wb'))
            except:
                print("pid failed:",data[t][i])
    
        print("images found:",len(alldata['wavelength']))