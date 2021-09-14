API_KEY = <my_api_key>

def set_up_environment(API_KEY):
  # Import Data Commons
  import datacommons as dc

  # Import other required libraries
  import matplotlib.pyplot as plt
  import matplotlib.patches as mpatches
  import pandas as pd
  import requests
  import json
  from google.colab import files
  from google.colab import drive
  uploaded = files.upload()
  from pywaffle.waffle import Waffle
  import io
  import time
  import multiprocessing
  # Mount the Drive
  drive.mount('/content/drive', force_remount=True)

  # REPLACE THIS with the path to your key.
  #key_path = '/content/drive/My Drive/DataCommons/secret.json'
  dc.set_api_key(API_KEY)
  # Read the key in and provide it to the Data Commons API
  #with open(key_path, 'r') as f:
  #  secrets = json.load(f)
  #  dc.set_api_key(secrets['dc_api_key'])

# The API server root URL
root = 'https://api.datacommons.org'

def read_candidate_snps():
  filename = 'Type_1_Diabetes_GM12878_all_trios.txt'
  df = pd.read_csv(io.BytesIO(uploaded[filename]), sep='\t')
  list_snps = list(df["SNP_id"])
  list_genes = list(df["gene_id"])
  count = 0
  dict_snp_genes = {}
  while count < len(list_snps):
    if list_snps[count] not in list(dict_snp_genes.keys()):
      dict_snp_genes[list_snps[count]] = set(list_genes[count])
    else:
      dict_snp_genes[list_snps[count]].add(list_genes[count])
    count += 1
  return(dict_snp_genes)

def query_data_commons(dict_snp_genes):
  dict_tissue_count = {'Skin Not Sun Exposed Suprapubic': 0, 'Whole Blood':0, 'Pancreas':0, 'Thyroid':0, 'None': 0}
  for snp in dict_snp_genes:
    try:
      temp = dc.get_property_values(['bio/' + snp], 'referenceSNPClusterID', out=False)['bio/' + snp]
      list_genes = set()
      for item in temp:
        tissue = dc.get_property_values([item], 'associatedTissue')[item][0]
        gene = dc.get_property_values([item], 'geneSymbol')[item][0].split('_')[1]
        if tissue in dict_tissue_count and gene in list(dict_snp_genes[snp]):
          dict_tissue_count[tissue] += 1
          list_genes.add(gene)
      dict_tissue_count['None'] += len(dict_snp_genes[snp]) - len(list(genes))
    except:
      dict_tissue_count['None'] += len(dict_snp_genes[snp])
return(dict_tissue_count)

 

def main():
  set_up_environment(API_KEY)
  start_time=time.time()
  dict_snp_genes=read_candidate_snps()
  dict_tissue_count = query_data_commons(dict_snp_genes)
  print("The number of CPUs used: " + str(multiprocessing.cpu_count()))
  print("The amount of time to process: " + str(time.time()-start_time))
   
if __name__ == '__main__':
    main()
    
