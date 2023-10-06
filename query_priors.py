import os
import glob
import json
import argparse

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

    translation = {
        "pl_orbper":"period",
        "pl_orbpererr1":"period_uperr",
        "pl_orbpererr2":"period_lowerr",

        "pl_orbsmax":"sma",
        "pl_orbsmaxerr1":"sma_lowerr",
        "pl_orbsmaxerr2":"sma_uperr",
        
        "pl_ratdor":"ars",
        "pl_ratdorerr1":"ars_lowerr",
        "pl_ratdorerr2":"ars_upper",

        "pl_orbeccen":"ecc",
        "pl_orbeccenerr1":"ecc_uperr",
        "pl_orbeccenerr2":"ecc_lowerr",
        
        "pl_orbincl":"inc",
        "pl_orbinclerr1":"inc_uperr",
        "pl_orbinclerr2":"inc_lowerr",

        "pl_bmassj":"mass",
        "pl_bmassjerr1":"mass_uperr",
        "pl_bmassjerr2":"mass_lowerr",

        "pl_radj":"rp",
        "pl_radjerr1":"rp_uperr",
        "pl_radjerr2":"rp_lowerr",
        
        "pl_tranmid":"t0",
        "pl_tranmiderr1":"t0_uperr",
        "pl_tranmiderr2":"t0_lowerr",

        "pl_orblper":"omega",
        "pl_orblpererr1":"omega_uperr",
        "pl_orblpererr2":"omega_lowerr",

        "st_teff":"T*",
        "st_tefferr1":"T*_uperr",
        "st_tefferr2":"T*_lowerr",

        "st_mass":"M*",
        "st_masserr1":"M*_uperr",
        "st_masserr2":"M*_lowerr",
        
        "st_rad":"R*",
        "st_raderr1":"R*_uperr",
        "st_raderr2":"R*_lowerr",
        
        "st_met":"FEH*",
        "st_meterr1":"FEH*_uperr",
        "st_meterr2":"FEH*_lowerr",
        
        "st_logg":"LOGG*",
        "st_loggerr1":"LOGG*_uperr",
        "st_loggerr2":"LOGG*_lowerr",

        "ra":"RA", # TODO check if this is saving
        "dec":"DEC",

        "pl_refname":"P.Reference",
        "st_refname":"S.Reference"
    }

    fixed_planet = {
        "period_units":"[days]",
        "sma_units":"[AU]",
        "ecc_units":"",
        "inc_units":"[degree]",
        "omega_units":"[degree]",
        "mass_units":"[Jupiter mass]",
        "rp_units":"[Jupiter radius]",
        "t0_units":"[Julian Day]",
    }

    fixed_star = {
        "T*_units":"[K]",
        "M*_units":"[Solar mass]",
        "R*_units":"[Solar radius]",
        "FEH*_units":"[Fe/H]",
        "LOGG*_units":"log10[cm.s-2]",
    }

    api_keys = [
        "hostname", "pl_name", "pl_letter"
    ]

    for k in translation.keys():
        api_keys.append(k)

    if args.target == "all":
        # query nasa exoplanet archive
        data = query_exo_archive(api_keys)
    else:
        data = query_exo_archive(api_keys,target=args.target)
        
    # same naming format
    data['dirname'] = [''.join(data.iloc[i].pl_name.split(' ')[:-1]) for i in range(len(data))]

    # get planet names from directories
    dirs = glob.glob("{}/*/".format(args.datadir))

    for i in range(len(dirs)):

        # process specific target
        if args.target != 'all':
            if ''.join(args.target.lower().split(' ')) not in dirs[i].lower():
                continue

        name = dirs[i].split('/')[1]

        mask = data['dirname'] == name
        prior= {}
        for k in fixed_star:
            prior[k] = fixed_star[k]

        for k in translation:
            if "st" in k:
                prior[translation[k]] = data[mask].iloc[0][k]
        
        # for each planet
        prior['planets'] = []
        for j in range(len(data[mask])):
            p = data[mask].iloc[j].pl_letter.strip()
            prior['planets'].append(p)
            prior[p] = {}

            # keys for units
            for k in fixed_planet:
                prior[p][k] = fixed_planet[k]

            # parameters from NASA exoplanet archive
            for k in translation:
                if "pl" in k:
                    prior[p][translation[k]] = data[mask].iloc[j][k]

        with open(dirs[i]+"prior.json","w") as outfile:
            json.dump(prior, outfile, indent=4)
