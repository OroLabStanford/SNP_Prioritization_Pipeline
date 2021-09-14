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

def query_data_commons(snp, dict_genetic_variant):
	try:
		temp = dc.get_property_values(["bio/"+snp], 'functionalCategory')
		for key in list(temp.keys()):
			dict_genetic_variant[snp] = temp[key]
	except:
		dict_genetic_variant[snp] = ['GeneticVariantFunctionalCategoryUnknown']
	return(dict_genetic_variant)

def get_variant_functional_category(df):
  """
  Given a file contaning genetic variant names, compute the functional category 
  for each genetic variant.
  @waf is a boolean indicating weither or not to display the waffle graph. 
  Default : True
  @hist is a boolean indicating weither or not to display the histogram
  Default : False
  @Return : A dataframe with the generic variant DCIDs and their associated 
  functional category
  """
  dict_genetic_variant = {}
  for index, row in df.iterrows():
    dict_genetic_variant = query_data_commons(row["SNPs"], dict_genetic_variant)
  for key in list(dict_genetic_variant.keys()):
    if 'GeneticVariantFunctionalCategoryCodingSynon' in dict_genetic_variant[key]:
      dict_genetic_variant[key].remove('GeneticVariantFunctionalCategoryCodingSynon')
      dict_genetic_variant[key].append('GeneticVariantFunctionalCategoryCodingSynonomous')
    if len(dict_genetic_variant[key]) > 1:
      if 'GeneticVariantFunctionalCategoryUnknown' in dict_genetic_variant[key]:
        dict_genetic_variant[key].remove('GeneticVariantFunctionalCategoryUnknown')
      if ('GeneticVariantFunctionalCategoryNearGene5' in dict_genetic_variant[key]) & ('GeneticVariantFunctionalCategoryUTR5' in dict_genetic_variant[key]) : 
          #If both FCs "NearGene5" and "UTR5" exists, we keep only "UTR5"
          dict_genetic_variant[key].remove('GeneticVariantFunctionalCategoryNearGene5')
      if ('GeneticVariantFunctionalCategoryNearGene3' in dict_genetic_variant[key])&('GeneticVariantFunctionalCategoryUTR3' in dict_genetic_variant[key]) : 
          #If both FCs "ncRNA" and "Intron" exists, we keep only "ncRNA"
          dict_genetic_variant[key].remove('GeneticVariantFunctionalCategoryNearGene3')
      if ('GeneticVariantFunctionalCategoryncRNA' in dict_genetic_variant[key])&('GeneticVariantFunctionalCategoryIntron' in dict_genetic_variant[key]) : 
          #If both FCs "NearGene3" and "UTR3" exists, we keep only "UTR3"
          dict_genetic_variant[key].remove('GeneticVariantFunctionalCategoryIntron')
    if len(dict_genetic_variant[key]) > 1:
      dict_genetic_variant[key] = ['GeneticVariantFunctionalCategoryConflicting']
    dict_genetic_variant[key] = dict_genetic_variant[key][0].replace('GeneticVariantFunctionalCategory',"")
  return(dict_genetic_variant)

def count_instances(d):
  '''
  counts = {"3'UTR":0, "5'UTR":0, 'CodingSynonomous':0, 'Frameshift':0, 'Intron':0,\
            'Missense':0, "NearGene3":0, "NearGene5":0, "Nonsense":0, "Splice3":0,\
            "Splice5":0, "StopLoss":0, "Unknown":0, "cdsIndel":0, "cdsReference":0,\
            "ncRNA":0, "Conflicting":0}
  for key, value in d.items():
    if value in counts:
      counts[value] += 1
    else:
      print("Warning: Functional category type " + value + " not counted!")
  '''
  counts = {}
  for key, value in d.items():
    if value not in counts:
      counts[value] = 1
    else:
      counts[value] +=1
  # rename keys
  dict_new_keys = {"3'UTR": "3' UTR", "5'UTR": "5' UTR", "CodingSynonomous": "Synonomous Codon", 
                   "NearGene3": "Near Gene 3'", "NearGene5": "Near Gene 5'", 
                   "Splice3": "Splice Acceptor Variant", "Splice5": "Splice Donor Variant",
                   "StopLoss": "Stop Loss", "cdsIndel": "Coding Reading Frame Synonomous Indel",
                   "cdsReference": "Allele on the Contig"}
  for key in list(counts.keys()):
    if key in dict_new_keys:
      counts[dict_new_keys[key]] = counts.pop(key)
  '''
  counts["5' UTR"] = counts.pop("5'UTR")
  counts["Synonomous Codon"] = counts.pop("CodingSynonomous")
  counts["Near Gene 3'"] = counts.pop("NearGene3")
  counts["Near Gene 5'"] = counts.pop("NearGene5")
  counts["Splice Acceptor Variant"] = counts.pop("Splice3")
  counts["Splice Donor Variant"] = counts.pop("Splice5")
  counts["Stop Loss"] = counts.pop("StopLoss")
  counts["Coding Reading Frame Synonomous Indel"] = counts.pop("cdsIndel")
  counts["Allele on the Contig"] = counts.pop("cdsReference")
  '''
  return(counts)
  

def main():
	set_up_environment(API_KEY)
	start_time=time.time()
	dict_genetic_variant=get_variant_functional_category(pd.read_csv(io.BytesIO(uploaded['Type_1_Diabetes_GM12878_cohesin_SNPs.txt']), names=["SNPs"]))
	dict_counts = count_instances(dict_genetic_variant)
	print("The number of CPUs used: " + str(multiprocessing.cpu_count()))
	print("The amount of time to process: " + str(time.time()-start_time))
    
if __name__ == '__main__':
    main()
