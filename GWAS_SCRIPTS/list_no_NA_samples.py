#!/usr/bin/env python3

#Script to extract individuals for a specific phenotype (in plink/regenie input) without NA, for filtering purposes

import argparse

arg_parser = argparse.ArgumentParser()

arg_parser.add_argument('-p', '--pheno_file', type=str, required=True, help='The phenotype file')
arg_parser.add_argument('-pn', '--pheno', type=str, required=True, help='The phenotype name to scan for and count NAs in')
arg_parser.add_argument('-c', '--covar_file', type=str, required=True, help='The covariate file to scan for and count NAs in') #TODO: Assumes covar file always present, fix later
arg_parser.add_argument('-cn', '--covars', type=str, required=True, help='The covariate names to check for missingness, comma separated, will be stripped of spaces') #TODO: Assumes covar file always present, fix later
arg_parser.add_argument('-o', '--out_file', type=str, required=True, help='The list of individuals without NA for the given phenotype')

args = arg_parser.parse_args()

pheno_path  = args.pheno_file
pheno_name  = args.pheno
covar_path  = args.covar_file
covar_names = args.covars
output_path = args.out_file

no_na_inds = {}

#Check phenotype first
with open(pheno_path, 'r') as data_file:
    header = data_file.readline()

    #Find index of phenotype of interest
    headers = header.strip().split("\t")
    pheno_index = headers.index(pheno_name)

    for row in data_file:
        fields = row.split("\t")
        fam_id = fields[0]
        ind_id = fields[1]
        pheno = fields[pheno_index]
        if (pheno != "NA"):
            no_na_inds[fam_id + '_'  + ind_id] = [True, fam_id, ind_id]

data_file.close()

covars = covar_names.replace(' ', '').split(",")
#Check covars
with open(covar_path, 'r') as covar_file:
    header = covar_file.readline()
    headers = header.strip().split("\t")
    for covar in covars:
        covar_index = headers.index(covar) 
        for row in covar_file:
            fields = row.strip().split("\t")
            fam_id = fields[0]
            ind_id = fields[1]
            covar_val = fields[covar_index]
            if (covar_val == "NA" and ((fam_id + "_" + ind_id) in no_na_inds)):
                no_na_inds[fam_id + "_" + ind_id] = [False, fam_id, ind_id]
        covar_file.seek(1)


with open(output_path, 'w') as output_file:
    for key in no_na_inds:
        if(no_na_inds[key][0] == True):
            output_file.write(no_na_inds[key][1] + "\t" + no_na_inds[key][2] + "\n")

