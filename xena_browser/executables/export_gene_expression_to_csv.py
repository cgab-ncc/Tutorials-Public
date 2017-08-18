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

import sys
sys.path.append("../")
import Tkinter, tkFileDialog
import gzip
import numpy as np
from make_dataset import *

root = Tkinter.Tk()
root.withdraw()

zipped_gene_expression_file_path = tkFileDialog.askopenfilename()
if zipped_gene_expression_file_path is None:
    exit()

# Process zipped gene expression file
with gzip.open(zipped_gene_expression_file_path, 'r') as f:
    gene_expression = GeneExpression(zipped_file=f)
    # Calculate average gene expression value ratio. Avg(tumor) and avg(normal) for each gene
    print("Started calculating average gene expression value ratio")
    avg_gene_exp_ratio = {} # key = gene, value = avg(tumor) / avg(normal)
    tumor_samples_gene_expressions = {} # key = gene, value = []
    normal_samples_gene_expressions = {} # key = gene, value = []
    for gene_name in gene_expression.get_genes_list():
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
    print("Finished calculating average gene expression value ratio")
    # Write the gene expression to .csv file
    print("Started writing data to .csv file")
    csv_save_file = tkFileDialog.asksaveasfile(mode='w', defaultextension=".csv")
    if csv_save_file is None:
        exit()
    # Header - 1 (sample tissue type; tumor or normal)
    csv_save_file.write("Tissue type (T=tumor,N=normal)")
    for sample_id in gene_expression.get_samples_list():
        if gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_TUMOR:
            csv_save_file.write(",T")
        elif gene_expression.get_sample_tissue_types_dict()[sample_id] == constants.SAMPLE_TYPE_NORMAL:
            csv_save_file.write(",N")
    csv_save_file.write("\n")
    # Header - 2 (sample ids)
    csv_save_file.write("gene_name")
    for sample_id in gene_expression.get_samples_list():
        csv_save_file.write("," + str(sample_id))
    csv_save_file.write(",avg(tumor),avg(normal),avg(tumor)/avg(normal)\n")
    # Gene expression data and avg ratio
    for gene_name in gene_expression.get_genes_list():
        csv_save_file.write(str(gene_name))
        for sample_id in gene_expression.get_samples_list():
            csv_save_file.write("," + str(gene_expression.get_gene_exp_matrix_dataframe()[sample_id][gene_name]))
        csv_save_file.write("," + str(np.mean(np.array(tumor_samples_gene_expressions[gene_name]))))
        csv_save_file.write("," + str(np.mean(np.array(normal_samples_gene_expressions[gene_name]))))
        csv_save_file.write("," + str(avg_gene_exp_ratio[gene_name]) + "\n")
    csv_save_file.close()
    print("Finished writing data to .csv file")