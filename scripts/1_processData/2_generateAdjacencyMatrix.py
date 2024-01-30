#####################################
# Title: generateAdjacencyMatrix.py
# Author: Vivek Sriram 
# Last revised: 03/07/23
# -----------------------------------
# Function: Create the adjacency matrix corresponding to a network where nodes
#             represent diseases and edges represent shared SNPs from the input
#             PheWAS summary data. This script is applied one-at-a-time to the
#             male and female PheWAS data.


# Import statements
import numpy as np
import pandas as pd
import os
import sklearn
import sklearn.metrics

##############################################
# 0. GET OUR INPUT DATA AND DISTANCE ARRAY
# Read in our raw data
rawDataDirectory = "/Users/viveksrm/Desktop/bernabeuData/bernabeu_female"
# Specify the p-value threshold to count a SNP as associated with a phenotype
pvalueThreshold = 1e-4

# CREATE OUR ASSOCIATION MATRIX
# The rows of the matrix will be the phenotypes, and the columns will be the SNPs
# Initialize a dictionary
phenotypeToSNPDictionary = {}
# Initialize a list of unique phenotypes and SNPs in the input data
uniquePhenotypes = []
uniqueSNPs = []

# Iterate over the input data directory
i = 1
numFiles = len(os.listdir(rawDataDirectory))
for filename in os.listdir(rawDataDirectory):
    print("Processing file " + filename + ", " + str(i) + " out of " + str(numFiles))
    # Read in the current file into a pandas dataframe and get the phenotype name
    currentFile = pd.read_csv(rawDataDirectory + "/" + filename, sep = "\t")
    phenotypeName = os.path.splitext(filename)[0]
    
    # Get all associated SNPs and their corresponding p-values for the current phenotype
    snpsForThisPhenotype = currentFile['snpid'].tolist()
    pValuesForThisPhenotype = currentFile['pval'].tolist()
    # If the current phenotype has associated SNPs...
    if(len(snpsForThisPhenotype) > 0):
        uniquePhenotypes.append(phenotypeName)
        phenotypeToSNPDictionary[str(phenotypeName)] = []
    
        # Iterate over all associated SNPs and check that they pass our max p-value threshold
        rowIndex = 0
        for elem in snpsForThisPhenotype:
            if(pValuesForThisPhenotype[rowIndex] <= pvalueThreshold):
            	phenotypeToSNPDictionary[str(phenotypeName)].append(elem)
            	print("Adding SNP " + str(elem) + " with p-value " + str(pValuesForThisPhenotype[rowIndex]))
            	if(elem not in uniqueSNPs):
                    uniqueSNPs.append(elem)
            rowIndex = rowIndex+1
    i = i+1

# Now we have a dictionary of phenotypes and their corresponding SNPs
# We also have a list of unique phenotypes and SNPs
# Initialize a dataframe that is the right dimensions: number of rows = number of phenotypes, number of columns = number of SNPs
# Label all the rows and columns according to our uniquePhenotypes and uniqueSNPs lists
associationMatrix = pd.DataFrame(np.zeros((len(uniquePhenotypes), len(uniqueSNPs))),
                     columns=uniqueSNPs, 
                     index=uniquePhenotypes)

for phen in phenotypeToSNPDictionary:
    currentListOfSnps = phenotypeToSNPDictionary[phen]
    for snp in currentListOfSnps:
        associationMatrix[snp][phen] = 1.0

# Drop any SNPs from our association matrix that have no associations with phenotypes in the matrix
droplist = [i for i in associationMatrix if associationMatrix[i].sum()<=1]
associationMatrix.drop(droplist,axis=1,inplace=True)

# Use cosine similarity to generate a weighted similarity matrix (W_w) from the association matrix
W_w = sklearn.metrics.pairwise.cosine_similarity(associationMatrix)
W_w_pandas = pd.DataFrame(W_w, index=associationMatrix.index, columns=associationMatrix.index)

# Get rid of all rows and columns in the similarity matrix 
W_w_pandas_noZeroes = W_w_pandas.loc[(W_w_pandas!=0).any(axis=1)]
W_w_pandas_noZeroes = W_w_pandas_noZeroes.loc[:, (W_w_pandas_noZeroes != 0).any(axis=0)]

# Convert the weighted similarity matrix (W_w) to a weighted adjacency matrix (A_w)
A_w = 1-W_w_pandas_noZeroes
A_w.to_csv("/Users/viveksrm/Desktop/bernabeuData/femaleBlock_assocMat_1neg4.csv")
