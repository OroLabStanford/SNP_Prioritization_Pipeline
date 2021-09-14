"""
Author:	        Pranav Bhardwaj
Date: 		5/20/2019
Title: 		all_steps.py
Edited By:      Samantha Piekos
"""

# Import Modules
import os
import numpy as np
import pandas as pd
import sys
import time

# File open
def read_file(data_file, file_type):
    '''The read_file function reads the given file to a pandas 
    dataframe with desired headers'''

    if file_type == "HiChIP":
        data = []
        with open(data_file, 'r') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split()
                data.append(line) 
        
        headers = ["bin1_chr", "bin1_start", "bin1_stop", "bin2_chr", "bin2_start", "bin2_stop", "counts", "FDR", "ID"]
        df = pd.DataFrame(data, columns = headers) 

        cast_cols = ["bin1_start", "bin1_stop","bin2_start", "bin2_stop"]
        cast_dict = {}
        for col in cast_cols:
            cast_dict[col] = int
        
        return df.astype(cast_dict)

    elif file_type == "ChIPseq":
        data = []
        with open(data_file, 'r') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split()
                data.append(line[:4]) 

        headers = ["region_chr", "region_start", "region_stop", "region_id"]
        df = pd.DataFrame(data, columns = headers) 

        cast_cols = ["region_start", "region_stop"]
        cast_dict = {}
        for col in cast_cols:
            cast_dict[col] = int
        
        return df.astype(cast_dict)

    elif file_type == "SNP":
        data = []
        with open(data_file, 'r') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split()
                data.append(line) 

        headers = ["SNP_chr", "SNP_start", "SNP_stop", "SNP_id", "SNP_MAF"]
        df = pd.DataFrame(data[1:], columns = headers) 

        cast_cols = ["SNP_start", "SNP_stop"]
        cast_dict = {}
        for col in cast_cols:
            cast_dict[col] = int
        
        return df.astype(cast_dict)

    elif file_type == "Gene":
        data = []
        with open(data_file, 'r') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split()
                data.append(line) 

        headers = ["gene_chr", "gene_start", "gene_stop", "gene_id", "pm"]
        df = pd.DataFrame(data, columns = headers) 

        cast_cols = ["gene_start", "gene_stop"]
        cast_dict = {}
        for col in cast_cols:
            cast_dict[col] = int
        
        return df.astype(cast_dict)

# Step 1 - Anchor Loops
def anchor_loops(HiChIP_df, ChIPseq_df):
    '''Anchors peaks in a bin. Returns a dataframe with these bin - peak pairs.
    Has a column region_bin with the bin that the peak falls in'''

    output_df = pd.DataFrame(columns = ["bin1_chr", "bin1_start", "bin1_stop", "bin2_chr", 
                                        "bin2_start", "bin2_stop", "counts", "FDR", "ID",
                                        "region_chr", "region_start", "region_stop", "region_id",
                                        "region_bin"])
    counter = 0

    # Iterate over ChIPseq peaks
    for i in range(ChIPseq_df.shape[0]):

        # Extract peak information
        peak_chr = ChIPseq_df['region_chr'][i]
        peak_start = ChIPseq_df['region_start'][i]
        peak_stop = ChIPseq_df['region_stop'][i]

        # Find bin1s where peak falls
        bin1_df = HiChIP_df[ HiChIP_df['bin1_chr'] == peak_chr]
        bin1_df = bin1_df[ peak_start <= bin1_df['bin1_stop'] ]
        bin1_df = bin1_df[ bin1_df['bin1_start'] <= peak_stop ]

        # Write row to output dataframe
        for j in range(bin1_df.shape[0]):
            index = bin1_df.index[j]
            HiChIP_row = list(HiChIP_df.iloc[index])
            peak_row = list(ChIPseq_df.iloc[i])
            output_df.loc[counter] = HiChIP_row + peak_row + [1]
            counter += 1

        # Find bin2s where peak falls
        bin2_df = HiChIP_df[ HiChIP_df['bin2_chr'] == peak_chr]
        bin2_df = bin2_df[ peak_start <= bin2_df['bin2_stop'] ]
        bin2_df = bin2_df[ bin2_df['bin2_start'] <= peak_stop ]

        # Write row to output dataframe
        for j in range(bin2_df.shape[0]):
            index = bin2_df.index[j]
            HiChIP_row = list(HiChIP_df.iloc[index])
            peak_row = list(ChIPseq_df.iloc[i])
            output_df.loc[counter] = HiChIP_row + peak_row + [2]
            counter += 1

    print('Finished Step 1')
    print("Number of seconds passed: {}".format(time.time()-start_time))
    print('')
    return output_df

# Step 2 - Match SNPs
def match_snps(step1_df, SNP_df):
    '''Matches SNPs falling partially in anchored peaks'''

    output_df = pd.DataFrame(columns = ["bin1_chr", "bin1_start", "bin1_stop", "bin2_chr", 
                                        "bin2_start", "bin2_stop", "counts", "FDR", "ID",
                                        "region_chr", "region_start", "region_stop", "region_id",
                                        "region_bin",
                                        "SNP_chr", "SNP_start", "SNP_stop", "SNP_id", "SNP_MAF"])
    counter = 0

    # Iterate over SNPs
    for i in range(SNP_df.shape[0]):

        # Extract SNP information 
        snp_chr = SNP_df['SNP_chr'][i]
        snp_start = SNP_df['SNP_start'][i]
        snp_stop = SNP_df['SNP_stop'][i]

        # Find peaks where SNP falls
        peak_df = step1_df[ step1_df['region_chr'] == snp_chr ]
        peak_df = peak_df[ snp_start <= peak_df['region_stop'] ]
        peak_df = peak_df[ peak_df['region_start'] <= snp_stop ]

        # Write row to output dataframe
        for j in range(peak_df.shape[0]):
            index = peak_df.index[j]
            step1_row = list(step1_df.iloc[index])
            SNP_row = list(SNP_df.iloc[i])
            output_df.loc[counter] = step1_row + SNP_row
            counter += 1

    print('Finished Step 2')
    print("Number of seconds passed: {}".format(time.time()-start_time))
    print('')
    return output_df

# Step 3 - Match Genes
def match_genes(step2_df, gene_df, wings):
    '''Matches genes falling partially in bin opposite of matched SNPs'''

    # Initialize headers of output file
    headers = ["bin1_chr", "bin1_start", "bin1_stop", "bin2_chr", 
                "bin2_start", "bin2_stop", "counts", "FDR", "ID",
                "region_chr", "region_start", "region_stop", "region_id",
                "region_bin",
                "SNP_chr", "SNP_start", "SNP_stop", "SNP_id", "SNP_MAF",
                "gene_chr", "gene_start", "gene_stop", "gene_id", 'pm']

    output_df = pd.DataFrame(columns = headers)
    counter = 0

    # Check for genes in bin 2 if SNP falls in bin 1
    bin2_df = step2_df[ step2_df['region_bin'] == 1 ]
    # Check for genes in bin 1 if SNP falls in bin 2
    bin1_df = step2_df[ step2_df['region_bin'] == 2 ]

    # Iterate over Genes
    for i in range(gene_df.shape[0]):

        # Extract Gene information 
        gene_chr = gene_df['gene_chr'][i]
        gene_start = gene_df['gene_start'][i] - wings
        gene_stop = gene_df['gene_stop'][i] + wings

        # Find bin1s where gene falls opposite of SNP
        bin1_iter_df = bin1_df[ bin1_df['bin1_chr'] == gene_chr ]
        bin1_iter_df = bin1_iter_df[ gene_start <= bin1_iter_df['bin1_stop'] ]
        bin1_iter_df = bin1_iter_df[ bin1_iter_df['bin1_start'] <= gene_stop ]

        # Write row to output dataframe
        for j in range(bin1_iter_df.shape[0]):
            index = bin1_iter_df.index[j]
            step2_row = list(step2_df.iloc[index])
            gene_row = list(gene_df.iloc[i])
            output_df.loc[counter] = step2_row + gene_row
            counter += 1

        # Find bin2s where gene falls opposite of SNP
        bin2_iter_df = bin2_df[ bin2_df['bin2_chr'] == gene_chr ]
        bin2_iter_df = bin2_iter_df[ gene_start <= bin2_iter_df['bin2_stop'] ]
        bin2_iter_df = bin2_iter_df[ bin2_iter_df['bin2_start'] <= gene_stop ]

        # Write row to output dataframe
        for j in range(bin2_iter_df.shape[0]):
            index = bin2_iter_df.index[j]
            step2_row = list(step2_df.iloc[index])
            gene_row = list(gene_df.iloc[i])
            output_df.loc[counter] = step2_row + gene_row
            counter += 1
            
    print('Finished Step 3')
    print("Number of seconds passed: {}".format(time.time()-start_time))
    print('')
    return output_df

def post_process(step3_df, output_file):
    '''Returns useful information on matched values, saves output to file'''

    # Save output to file
    step3_df.to_csv(output_file, sep='\t', index = False)

    # Return unique elements in trios
    unique_peaks = set(step3_df['region_id'])
    unique_snps = set(step3_df['SNP_id'])
    unique_genes = set(step3_df['gene_id'])
    unique_trios = set()
    for i in range(step3_df.shape[0]):
        trio = (step3_df['region_id'][i], step3_df['SNP_id'][i], step3_df['gene_id'][i])
        unique_trios.add(trio)

    return unique_peaks, unique_snps, unique_genes, unique_trios

if __name__ == '__main__':

    # Start program timer
    start_time = time.time()
    print("Job initiated")
    print("Number of seconds passed: {}".format(time.time()-start_time))
    
    # Usage message
    if (len(sys.argv) < 6):
        print("Usage: ")
        print("        python all_steps.py <HiChIP_file> <ChIPseq_file> <SNP_file> <Gene_file> <Output_file> <Optional: Gene_Wings>")
        sys.exit()

    # Ensure all information carried along in pandas dataframes
    pd.set_option('display.max_colwidth', -1)

    # Read inputs
    HiChIP_df = read_file(sys.argv[1], "HiChIP")
    ChIPseq_df = read_file(sys.argv[2], "ChIPseq")
    SNP_df = read_file(sys.argv[3], "SNP")
    Gene_df = read_file(sys.argv[4], "Gene")

    # Perform step 1 of pipeline
    step1_df = anchor_loops(HiChIP_df, ChIPseq_df)

    # Perform step 2 of pipeline
    step2_df = match_snps(step1_df, SNP_df)

    # Perform step 3 of pipeline 
    if len(sys.argv) >= 7:
        step3_df = match_genes(step2_df, Gene_df, int(sys.argv[6]))
    else:
        step3_df = match_genes(step2_df, Gene_df, 5000)

    # Extract results from data, save step3_df to output file
    unique_peaks, unique_snps, unique_genes, unique_trios = post_process(step3_df, sys.argv[5])

    # Print summary
    print("Summary of Run:")
    print("HiChIP File:                {}".format(sys.argv[1]))
    print("ChIPseq File:               {}".format(sys.argv[2]))
    print("SNP File:                   {}".format(sys.argv[3]))
    print("Gene File:                  {}".format(sys.argv[4]))
    print("Output File:                {}".format(sys.argv[5]))
    print("Unique Peaks identified:    {}".format(len(unique_peaks)))
    print("Unique SNPs identified:     {}".format(len(unique_snps)))
    print("Unique Genes identified:    {}".format(len(unique_genes)))
    print("Unique Trios identified:    {}".format(len(unique_trios)))
    print("Time for execution:         {} seconds".format(time.time() - start_time))
    print("Unique Peaks:")
    print(unique_peaks)
    print("Unique SNPs:")
    print(unique_snps)
    print("Unique Genes:")
    print(unique_genes)
