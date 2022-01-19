"""
Author: Samantha Piekos
Date: 1/12/22
Last Edited Date: 1/19/22
Name: get_genetic_variant_clinical_significance.py
Description: "This queries Biomedical Data Commons using the python API to get
the clinical significance of each genetic variant from an input file that
contains a list of genetic variants of interest. It counts and prints the
number of genetic variants in the list that belong to each category of clinical
significance. Then it writes this info to a tab delimited file containing two 
columns: rsID and clinical significance"

python3 get_genetic_variant_clinical_significance.py file_import file_export

@filename_import  filepath to the import file of the genetic variant list
@filename_export  filepath of where to write the genetic variants and thier
                  associated clinical significance

"""

# Import Data Commons
import datacommons as dc

# Import pandas
import pandas as pd

# Specify if import file of SNP list contains a header
HEADER = 0  # change to None if no header in import file

def query_data_commons(list_gv):
  # query Biomedical Data Commons to retrieve the clinical signficance of each
  # genetic variant
  dict_clin_sig = {}
  for snp in list_gv:
    try:
      temp = dc.get_property_values(['bio/' + snp], 'clinicalSignificance')
      if len(temp['bio/' + snp]) == 0: 
        dict_clin_sig[snp] = ['Unknown']
      else:
        dict_clin_sig[snp] = temp['bio/' + snp]
    except:
      dict_clin_sig[snp] = ['Unknown']
  return dict_clin_sig


def reformat_clin_sig(list_clin_sig):
  # reduce clinical significance types to 13 categories
  dict_rename = {'ClinSigAffects': 'Associated',
                 'ClinSigAssociation': 'Associated',
                 'ClinSigAssociationNotFound': 'Not Associated',
                 'ClinSigBenign': 'Benign',
                 'ClinSigBenignLikelyBenign': 'Likely Benign',
                 'ClinSigConflictingPathogenicity': 'Conflicting Reports',
                 'ClinSigDrugResponse': 'Drug Response',
                 'ClinSigHistocompatability': 'Histocompatibility',
                 'ClinSigLikelyBenign': 'Likely Benign',
                 'ClinSigLikelyPathogenic': 'Likely Pathogenic',
                 'ClinSigProtective': 'Protective',
                 'ClinSigRiskFactor': 'Risk Factor',
                 'ClinSigUncertain': 'Unknown',
                 'ClinSigUnknown': 'Unknown',
                 'ClinSigUntested': 'Unknown',
                 'Unknown': 'Unknown'}
  set_renamed = set()
  for item in list_clin_sig:
    set_renamed.add(dict_rename[item])
  # remove unknown if another report of clinical significance exists
  if len(set_renamed) > 1 and 'Unknown' in set_renamed:
    set_renamed.remove('Unknown')
  return list(set_renamed)


def remove_zero_values(d):
  # remove any key, value pairs from dictionary if value is 0
  d_new = {}
  for k, v in d.items():
    if v != 0:
      d_new[k] = v
  return d_new


def print_dict(d):
  for k, v in d.items():
    print('# ' + k + ': ' + str(v))


def reformat_instances(dict_clin_sig):
  # count the number of instances of each clinical significance that's
  # associated with genetic variants in the input list
  dict_counts = {'Associated': 0,
                 'Benign': 0,
                 'Conflicting Reports': 0,
                 'Drug Responce': 0,
                 'Histocompatibility': 0,
                 'Likely Benign': 0,
                 'Likely Pathogenic': 0,
                 'Not Associated': 0,
                 'Other': 0,
                 'Pathogenic': 0,
                 'Protective': 0,
                 'Risk Factor': 0,
                 'Unknown': 0}
  dict_clin_sig_formatted = {}
  for snp, clin_sig in dict_clin_sig.items():
      clin_sig = reformat_clin_sig(clin_sig)
      if len(clin_sig) > 1:
        dict_counts['Conflicting Reports'] += 1
        clin_sig = ['Conflicting Reports']
      else:
        dict_counts[clin_sig[0]] += 1
      dict_clin_sig_formatted[snp] = clin_sig[0]
  dict_counts = remove_zero_values(dict_counts)
  print('# Number of SNPs of Each Clinical Significance Category:')
  print_dict(dict_counts)
  return dict_clin_sig_formatted


def print_dict_to_file(filename_export, d):
  # write results of final formatted dictionary to export file
  w = open(filename_export, 'w')
  w.write('rsID\tclinical_significance\n')
  for k, v in d.items():
    w.write(k + '\t' + v + '\n')


def identify_clinical_significance(filename_import, filename_export):
  df = pd.read_csv(filename_import, header=HEADER)
  list_gv = list(df.iloc[:, 0])
  dict_clin_sig = query_data_commons(list_gv)
  dict_clin_sig_formatted = reformat_instances(dict_clin_sig)
  # write clinical significance of genetic variants to export file
  print_dict_to_file(filename_export, dict_clin_sig_formatted)
    

def main():
  import sys

  filename_import = sys.argv[1]
  filename_export = sys.argv[2]

  identify_clinical_significance(filename_import, filename_export)


if __name__ =='__main__':
  main()
