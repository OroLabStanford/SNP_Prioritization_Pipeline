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

def query_data_commons(list_gv):
	dict_clin_sig = {}
	for snp in list_gv:
	  try:
	    temp = dc.get_property_values(['bio/' + snp], 'clinicalSignificance')
	    dict_clin_sig[snp] = temp['bio/' + snp]
	  except:
	    dict_clin_sig[snp] = ['Unknown']
	return(dict_clin_sig)

def count_instances(dict_clin_sig):
	# Count the number of snps in each clinicalSignificance Category
	rename_dict = {'ClinSigBenign': 'Benign', 'ClinSigUncertain': 'Unknown'}
	counts_dict = {'Unknown': 0, 'Conflicting': 0}
	for snp, clinSig in dict_snps_clin_sig.items():
	  if len(clinSig) > 1:
	    counts_dict['Conflicting'] += 1
	    continue
	  clinSig = clinSig[0]
	  if clinSig == 'Unknown':
	    counts_dict['Unknown'] += 1
	  else:
	    clinSig = rename_dict[clinSig]
	    if clinSig in counts_dict:
	      counts_dict[clinSig] += 1
	    else:
	      counts_dict[clinSig] = 1
	return(counts_dict)

  

def main():
	set_up_environment(API_KEY)
	# start clock
	start_time = time.time()

	filename = 'Type_1_Diabetes_GM12878_candidate_SNPs.txt'
	df = pd.read_csv(io.BytesIO(uploaded[filename]), header=None)
	list_gv = list(df.iloc[:, 0])
	dict_clin_sig = query_data_commons(list_gv)
	dict_counts = count_instances(dict_clin_sig)
	print("The number of CPUs used: " + str(multiprocessing.cpu_count()))
	print("The amount of time to process: " + str(time.time()-start_time))

    
if __name__ == '__main__':
    main()
