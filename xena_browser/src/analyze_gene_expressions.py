"""
Author:
Clinical Genomics Analysis Branch,
National Cancer Center of Korea

Update:
Last Updated: Aug 17, 2017

Purpose:
The purpose of this program is to create a .csv file that
contains the following information based on gene expression folder
dwonloaded from UCSC Xena Browser
"""

import Tkinter, tkFileDialog
import gzip
import os
import numpy as np
import errno
from make_dataset import *
import constants

root = Tkinter.Tk()
root.withdraw()

print("\n\n===Clinical Genomics Analysis Branch, National National Cancer of Korea. All rights reserved. 2017===\n")
print("The purpose of this program is to calculate average gene expression of tumor samples and normal samples\n" +
      "given a zipped gene expression file from UCSC Xena Browser.\n\n")
print("---> Please select the zipped gene expression file downloaded from UCSC Xena Browser <---\n\n")
zipped_gene_expression_file_path = tkFileDialog.askopenfilename()
root.update()

if zipped_gene_expression_file_path is None:
    exit()

# Process zipped gene expression file
with gzip.open(zipped_gene_expression_file_path, 'r') as f:
    gene_expression = GeneExpression(zipped_file=f)
    # Calculate average gene expression value ratio. Avg(tumor) and avg(normal) for each gene
    avg_gene_exp_ratio = {} # key = gene, value = avg(tumor) / avg(normal)
    tumor_samples_gene_expressions = {} # key = gene, value = []
    normal_samples_gene_expressions = {} # key = gene, value = []
    genes_calculated = 0
    for gene_name in gene_expression.get_genes_list():
        if genes_calculated % 200 == 0:
            progress_percentage = 100*(float(genes_calculated)/float(len(gene_expression.get_genes_list())))
            sys.stdout.write("\rCalculating average gene expression ratios between tumor and normal samples %d%%" % progress_percentage)
            sys.stdout.flush()
        genes_calculated = genes_calculated + 1
        tumor_samples_gene_expressions[gene_name] = []
        normal_samples_gene_expressions[gene_name] = []
        for sample_id in gene_expression.get_samples_list():
            sample_gene_expression_value = gene_expression.get_gene_exp_matrix_dataframe()[sample_id][gene_name]
            if gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_TUMOR:
                tumor_samples_gene_expressions[gene_name].append(sample_gene_expression_value)
            elif gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_NORMAL:
                normal_samples_gene_expressions[gene_name].append(sample_gene_expression_value)
        if float(np.mean(np.array(normal_samples_gene_expressions[gene_name]))) == 0.0:
            avg_gene_exp_ratio[gene_name] = "n/a"
        else:
            avg_gene_exp_ratio[gene_name] = float(np.mean(np.array(tumor_samples_gene_expressions[gene_name]))) / \
                                            float(np.mean(np.array(normal_samples_gene_expressions[gene_name])))
    print("\nFinished calculating average gene expression value ratio")
    # Write .csv files
    print("Started writing data to .csv files")
    summary_file = os.getcwd() + "/results/analysis.csv"
    raw_gene_exp_file = os.getcwd() + "/results/gene_expressions.csv"
    # Create the files
    if not os.path.exists(os.path.dirname(summary_file)):
        try:
            os.makedirs(os.path.dirname(summary_file))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    if not os.path.exists(os.path.dirname(raw_gene_exp_file)):
        try:
            os.makedirs(os.path.dirname(raw_gene_exp_file))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    # Write to the files
    print("Writing data to analysis.csv file")
    with open(summary_file, "w") as csv_save_file:
        # Header
        csv_save_file.write("gene_name,tumor_samples_count,normal_samples_count,avg(tumor),avg(normal),avg(tumor)/avg(normal),avg(tumor)=0_count,avg(normal)=0_count\n")
        # Gene expression data and avg ratio
        genes_written = 0
        for gene_name in gene_expression.get_genes_list():
            if genes_calculated % 200 == 0:
                progress_percentage = 100 * (float(genes_calculated) / float(len(gene_expression.get_genes_list())))
                sys.stdout.write("\r%d%%" % progress_percentage)
                sys.stdout.flush()
            csv_save_file.write(str(gene_name))
            csv_save_file.write("," + str(len(tumor_samples_gene_expressions[gene_name])))
            csv_save_file.write("," + str(len(normal_samples_gene_expressions[gene_name])))
            csv_save_file.write("," + str(np.mean(np.array(tumor_samples_gene_expressions[gene_name]))))
            csv_save_file.write("," + str(np.mean(np.array(normal_samples_gene_expressions[gene_name]))))
            csv_save_file.write("," + str(avg_gene_exp_ratio[gene_name]))
            zero_gene_exp_values_tumor = [a for a in tumor_samples_gene_expressions[gene_name] if a == 0]
            zero_gene_exp_values_normal = [a for a in normal_samples_gene_expressions[gene_name] if a == 0]
            csv_save_file.write("," + str(len(zero_gene_exp_values_tumor)))
            csv_save_file.write("," + str(len(zero_gene_exp_values_normal)) + "\n")
    print("Finished writing data to analysis.csv file")
    print("Writing data to gene_expressions.csv file")
    with open(raw_gene_exp_file, "w") as csv_save_file:
        # Header - 1 (sample tissue type; tumor or normal)
        csv_save_file.write("Tissue type (T=tumor;N=normal)")
        for sample_id in gene_expression.get_samples_list():
            if gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_TUMOR:
                csv_save_file.write(",T")
        for sample_id in gene_expression.get_samples_list():
            if gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_NORMAL:
                csv_save_file.write(",N")
        csv_save_file.write("\n")
        # Header - 2 (sample ids)
        csv_save_file.write("gene_name")
        for sample_id in gene_expression.get_samples_list():
            if gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_TUMOR or \
                gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_NORMAL:
                csv_save_file.write("," + str(sample_id))
        csv_save_file.write("\n")
        # Gene expression data and avg ratio
        genes_written = 0
        for gene_name in gene_expression.get_genes_list():
            if genes_calculated % 200 == 0:
                progress_percentage = 100 * (float(genes_calculated) / float(len(gene_expression.get_genes_list())))
                sys.stdout.write("\r%d%%" % progress_percentage)
                sys.stdout.flush()
            csv_save_file.write(str(gene_name))
            for sample_id in gene_expression.get_samples_list():
                csv_save_file.write("," + str(gene_expression.get_gene_exp_matrix_dataframe()[sample_id][gene_name]))
            csv_save_file.write("\n")
    print("Finished writing data to .csv files")