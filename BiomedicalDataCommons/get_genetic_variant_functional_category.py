"""
Author: Samantha Piekos
Date: 1/19/22
Name: get_genetic_variant_functional_category.py
Description: "This queries Biomedical Data Commons using the python API to get
the functional category of each genetic variant from an input file that
contains a list of genetic variants of interest. It counts and prints the
number of genetic variants in the list that belong to each functional category
Then it writes this info to a tab delimited file containing two columns:
rsID and functional category"

python3 get_genetic_variant_functional_category.py file_import file_export

@filename_import  filepath to the import file of the genetic variant list
@filename_export  filepath of where to write the genetic variants and thier
                  associated functional category

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
      temp = dc.get_property_values(['bio/' + snp], 'functionalCategory')
      if len(temp['bio/' + snp]) == 0: 
        dict_func_cat[snp] = ['Unknown']
      else:
        dict_clin_sig[snp] = temp['bio/' + snp]
    except:
      dict_clin_sig[snp] = ['Unknown']
  return dict_clin_sig


def clean_entries(s):
  # remove unknown if another report of clinical significance exists
  if len(s) > 1 and 'Unknown' in s:
    s.remove('Unknown')
  # remove Near Gene 3' report if 3' UTR report also exists
  if ("3' UTR" in s) and ("Near Gene 3'" in s):
    s.remove("Near Gene 3'")
  # remove Near Gene 5' report if 5' UTR report also exists
  if ("5' UTR" in s) and ("Near Gene 5'" in s):
    s.remove("Near Gene 5'")
  # remove Intron report if ncRNA report also exists
  if ('ncRNA' in s) and ('Intron' in s):
    s.remove('Intron')
  return(s)

def reformat_category(list_category):
  # rename functional category types
  dict_rename = {'GeneticVariantFunctionalCategoryUTR3': "3' UTR",
                 'GeneticVariantFunctionalCategoryUTR5': "5' UTR",
                 'GeneticVariantFunctionalCategoryCodingSynon': 'Coding Synonomous',
                 'GeneticVariantFunctionalCategoryFrameshift': 'Frameshift',
                 'GeneticVariantFunctionalCategoryIntron': 'Intron',
                 'GeneticVariantFunctionalCategoryMissense': 'Missense',
                 'GeneticVariantFunctionalCategoryNearGene3': "Near Gene 3'",
                 'GeneticVariantFunctionalCategoryNearGene5': "Near Gene 5'",
                 'GeneticVariantFunctionalCategorySplice3': "Splice 3'",
                 'GeneticVariantFunctionalCategorySplice5': "Splice 5'",
                 'GeneticVariantFunctionalCategoryStopLoss': 'Stop Loss',
                 'GeneticVariantFunctionalCategoryUnknown': 'Unknown',
                 'GeneticVariantFunctionalCategoryCDSIndel': 'CDS Indel',
                 'GeneticVariantFunctionalCategoryCDSReference': 'CDS Reference',
                 'GeneticVariantFunctionalCategoryncRNA': 'ncRNA',
                 'Unknown': 'Unknown'}
  set_renamed = set()
  for item in list_category:
    set_renamed.add(dict_rename[item])
  # clean up functional categories associated with a genetic variant
  set_renamed = clean_entries(set_renamed)
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


def reformat_instances(dict_func_cat):
  # count the number of instances of each clinical significance that's
  # associated with genetic variants in the input list
  dict_counts = {"3' UTR": 0,
                 "5' UTR": 0,
                 'Coding Synonomous': 0,
                 'Multiple Reports': 0,
                 'Frameshift': 0,
                 'Intron': 0,
                 'Missense': 0,
                 "Near Gene 3'": 0,
                 "Near Gene 5'": 0,
                 "Splice 3'": 0,
                 "Splice 5'": 0,
                 'Stop Loss': 0,
                 'Unknown': 0,
                 'CDS Indel': 0,
                 'CDS Refernce': 0,
                 'ncRNA': 0}
  dict_func_cat_formatted = {}
  for snp, func_cat in dict_func_cat.items():
    func_cat = reformat_category(func_cat)
    # report as new multiple reports category if multiple different
    # functional categories are associated with a given genetic variant
    if len(func_cat) > 1:
      dict_counts['Multiple Reports'] += 1
    else:
      dict_counts[func_cat[0]] += 1
    dict_func_cat_formatted[snp] = func_cat
  dict_counts = remove_zero_values(dict_counts)
  print('# Number of SNPs of Each Functional Category:')
  print_dict(dict_counts)
  return dict_func_cat_formatted


def print_dict_to_file(filename_export, d):
  # write results of final formatted dictionary to export file
  w = open(filename_export, 'w')
  w.write('rsID\tfunctional_category\n')
  for k, v in d.items():
    w.write(k + '\t' + (';').join(v) + '\n')


def identify_clinical_significance(filename_import, filename_export):
  df = pd.read_csv(filename_import, header=HEADER)
  list_gv = list(df.iloc[:, 0])
  dict_func_cat = query_data_commons(list_gv)
  dict_func_cat_formatted = reformat_instances(dict_func_cat)
  # write functional category of genetic variants to export file
  print_dict_to_file(filename_export, dict_func_cat_formatted)
    

def main():
  import sys

  filename_import = sys.argv[1]
  filename_export = sys.argv[2]

  identify_clinical_significance(filename_import, filename_export)


if __name__ =='__main__':
  main()
