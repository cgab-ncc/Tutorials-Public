import pandas as pd
import numpy as np
import constants


class GeneExpression:

    def __init__(self, zipped_file, verbose=False):
        """
        REQUIRES:	zipped_file = string (full path to the zipped HiSeqV2 file)

        """
        if verbose:
            print("## Started initializing GeneExpression ##")
        # Raw data
        self.__verbose = verbose
        self.__gene_exp_matrix_df = None
        # Processed data
        self.__sample_tissue_types_dict = {} # key = sample id, value = tissue type
        self.__sample_cancer_stages_dict = {} # key = sample id, value = cancer stage
        # Construct data
        self.__read_file(zipped_file=zipped_file)
        self.__compute_sample_tissue_types_dict()
        if verbose:
            print("## Finished initializing GeneExpression ##\n")

    """
    PRIVATE functions
    """
    def __read_file(self, zipped_file):
        """
        This function reads the zipped gene expression file
        REQUIRES:   zipped_file = string (should include the path to the zipped file)
        MODIFIES:   self.__gene_exp_matrix_df
        EFFECTS:    populates self.participants
        """
        if self.__verbose:
            print("Started reading zipped gene expression file")
        self.__gene_exp_matrix_df = pd.read_csv(zipped_file, sep="\t", header=0, index_col=0)
        if self.__verbose:
            print("Finished reading zipped gene expression file")

    def __compute_sample_tissue_types_dict(self):
        """
        This function populates self.__sample_tissue_types_dict based on self.__gene_exp_matrix_df

        REQUIRES:   self.__gene_exp_matrix_df is populated
        MODIFIES:   self.__sample_tissue_types_dict
        EFFECTS:    populates self.__sample_tissue_types_dict where
                    key = sample_id
                    value = sample type where
                    "1" = "primary solid tumor"
                    "6" = "metastatic"
                    "11" = "solid tissue normal"
        """
        if self.__verbose:
            print("Started computing sample types dictionary")
        for column in self.__gene_exp_matrix_df.columns:
            id = str(column)
            id_components = id.split("-")
            sample_type = int(id_components[-1])
            self.__sample_tissue_types_dict[id] = str(sample_type)
        if self.__verbose:
            print("Finished computing sample types dictionary")

    """
    GET functions
    """
    def get_gene_exp_matrix_dataframe(self):
        return self.__gene_exp_matrix_df

    def get_genes_list(self):
        return self.__gene_exp_matrix_df.iloc[:,0].index.values.tolist()

    def get_samples_list(self):
        return self.__gene_exp_matrix_df.columns.values.tolist()

    def get_sample_tissue_types_dict(self):
        return self.__sample_tissue_types_dict

    """
    PRINT functions
    """