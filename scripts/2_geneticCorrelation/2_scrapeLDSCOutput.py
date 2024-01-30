#####################################
# Title: scrapeLDSCOutput.py
# Author: Vivek Sriram
# Last revised: 04/08/2023
# ------------------------------------------------
# Function: Script to collect genetic correlation (rg) values
#           from LDSR output log files. This output is used
#           to define edge weights for a genetic correlation
#           disease-disease network.

# Import statements
import os

# Get the directory of rg results
in_rgResultsPath = "/Users/viveksrm/Desktop/rgResults/male/"
rgResults = os.fsencode(in_rgResultsPath)

# Establish the output file of pairwise genetic correlations
outputFilePath = "/Users/viveksrm/Desktop/rgResults/gcDDN_allEdges_male.tsv"
# Open the output file
currentOutputFile = open(outputFilePath, "a")

# Iterate over all .log files in the directory of rg results
i = 0
for file in os.listdir(rgResults):
    filename = os.fsdecode(file)
    if filename.endswith(".log"):
        print("Processing " + filename)
        # Get the first phenotype ID and the second phenotype ID from the filename (formated as x_x_x_ID1_x_x_x_ID2_x_x.log)
        firstPhenotypeIndex = 3
        secondPhenotypeIndex = 7
        firstPhenotype = filename.split("_")[firstPhenotypeIndex]
        secondPhenotype = filename.split("_")[secondPhenotypeIndex]

        # Get correlations for phenotype pairs that have two different phenotypes
        if(firstPhenotype != secondPhenotype):
            # Open the file and get the correlation info for the phenotype pair
            fullFilePath = in_rgResultsPath + filename
            with open(fullFilePath, "r") as currentRgFile:
                x = currentRgFile.readlines()
                x = [elem.strip() for elem in x[x.index('Summary of Genetic Correlation Results\n') + 1:] if elem != '\n' ]
                header = x[0].split()
                values = x[1].split()
                currentRgFile.close()

                # If this correlation isn't NA...
                phenotype1Index = 0
                phenotype2Index = 1
                rgIndex = 2
                if(values[rgIndex] != "NA"):
                    # If we haven't created the output file yet, then we need to add the header
                    if(i == 0):
                        currentOutputFile.write("\t".join(header))
                        currentOutputFile.write("\n")

                    # Clean up values for p1 (index 0) and p2 (index 1)
                    values[phenotype1Index] = values[phenotype1Index].replace('munged_male_bernabeu_blocks/munged_cc_', '')
                    values[phenotype1Index] = values[phenotype1Index].replace('.sumstats.gz', '')
                    values[phenotype2Index] = values[phenotype2Index].replace('munged_male_bernabeu_blocks/munged_cc_', '')
                    values[phenotype2Index] = values[phenotype2Index].replace('.sumstats.gz', '')

                    # Write the correlation info to our output file
                    currentOutputFile.write("\t".join(values))
                    currentOutputFile.write("\n")
                    i = i+1
currentOutputFile.close()
