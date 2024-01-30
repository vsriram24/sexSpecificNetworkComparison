#####################################
# Title: splitDataBySex.py
# Author: Vivek Sriram
# Last revised: 03/10/2023
# ------------------------------------------------
# Function: The raw UKBB data from Bernabeu et al. includes both male-
#            and female-specific association data within the same file
#            for each phenotype. This script serves to split relevant
#            data into separate files by sex.

# Import statements
import os
import argparse

# Input and output directory paths
inDirectoryPath = "/Users/viveksrm/Desktop/bernabeuData/"
outputDirectory_female = "/Users/viveksrm/Desktop/bernabeu_female/"
outputDirectory_male = "/Users/viveksrm/Desktop/bernabeu_male/"

# Specify the delimiter used in the input files
columnDelimiter = "\t"

# Iterate over all files in the input directory
for filename in os.listdir(inDirectoryPath):
    print("Working on " + str(filename))
    # Open new files in the two filter directories that correspond to the filtered versions of the current file
    currentMaleOutput = open(outputDirectory_male + "male_" + filename, 'a')
    currentFemaleOutput = open(outputDirectory_female + "female_" + filename, 'a')

    # Set the header lines for the filtered files
    headerLine = "snpid" + columnDelimiter + "chr" + columnDelimiter + "effect" + columnDelimiter + "ref" + columnDelimiter + "beta" + columnDelimiter + "SE" + columnDelimiter + "pval" + "\n"
    currentMaleOutput.write(headerLine)
    currentFemaleOutput.write(headerLine)

    # Read through the lines of the current file
    f = open(inDirectoryPath + filename, 'r')
    lines = f.readlines()
    lineCounter = 0

    # Get the elements of each line (and make sure to remove the EOL character in the last column)
    for line in lines[1:]:
        lineAsArray = line.split(columnDelimiter)
        lineAsArray[len(lineAsArray)-1] = lineAsArray[len(lineAsArray)-1].rstrip()

        snpInformation = str(lineAsArray[0]) + columnDelimiter + str(lineAsArray[3]) + columnDelimiter + str(lineAsArray[1]) + columnDelimiter + str(lineAsArray[2])
        maleOutputLine = snpInformation + columnDelimiter + str(lineAsArray[9]) + columnDelimiter + str(lineAsArray[11]) + columnDelimiter + str(lineAsArray[13]) + "\n"
        femaleOutputLine = snpInformation + columnDelimiter + str(lineAsArray[8]) + columnDelimiter + str(lineAsArray[10]) + columnDelimiter + str(lineAsArray[12]) + "\n"
                
        currentMaleOutput.write(maleOutputLine)
        currentFemaleOutput.write(femaleOutputLine)

        lineCounter = lineCounter +1

    f.close()
    currentMaleOutput.close()
    currentFemaleOutput.close()
    
print("Finished filtering data")

