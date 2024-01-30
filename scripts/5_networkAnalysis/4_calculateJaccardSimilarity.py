"""
Title: calculateJaccardSimilarity.py
Author: Vivek Sriram
Date: 04/28/23

-----------------------------------------------------
Function: Calculate Jaccard Similarity (statistic to gauge similarity of two sets)
            between the male- and female-specific edge maps. Also, determine the
            types of edges (connecting same or different disease categories) for
            both graphs.
"""

# Import packages
import os
import pandas as pd
import csv

# Input DDNs (edge map is in TSV format)
in_edgelist_f = pd.read_csv('/Users/viveksrm/Desktop/female_edgeMap.tsv', index_col=None, sep = "\t")
in_edgelist_m = pd.read_csv('/Users/viveksrm/Desktop/male_edgeMap.tsv', index_col=None, sep = "\t")

# Get descriptions (including disease categories) for all phenotypes
in_traitDescriptions = pd.read_csv('/Users/viveksrm/Desktop/traitDescription_diseaseCategories.csv')

# Initialize a dictionary of disease-category values
diseaseCategoryMap = {}

##########################################################
edges_f = set()
edges_m = set()

# Function to convert edge maps to sets of edges
def getEdgesFromEdgeList(edgeList):
    outputEdgeSet = set()
    for index, row in edgeList.iterrows():
        source = row['Source']
        target = row['Target']
        if (source < target):
            currentEdge = row['Source'] + " " + row['Target']
        else:
            currentEdge = row['Target'] + " " + row['Source']
        outputEdgeSet.add(currentEdge)

    return outputEdgeSet

edges_f = getEdgesFromEdgeList(in_edgelist_f)
edges_m = getEdgesFromEdgeList(in_edgelist_m)

# Get intersection and union of edge sets for the two sexes
edgeIntersection = edges_m.intersection(edges_f)
edgeUnion = edges_m.union(edges_f)

# Calculate Jaccard similarity
print("Intersection count is: ")
print(len(edgeIntersection))
print()
print("Union count is: ")
print(len(edgeUnion))
print()
print("Jaccard similarity is:")
print(len(edgeIntersection)/len(edgeUnion))


#####################################################################
# Get disease categories for each phenotype
for index, row in in_traitDescriptions.iterrows():
    diseaseCategoryMap[row['Trait']] = row['Disease Category']

# Function to calculate number of edges from same and different disease categories
def getCommonAndDifferentCategoryEdgeCounts(edgeList):
    sharedCategoryCount = 0
    differentCategoryCount = 0
    for index, row in edgeList.iterrows():
        firstCategory = diseaseCategoryMap[row['Source']]
        secondCategory = diseaseCategoryMap[row['Target']]

        if(firstCategory == secondCategory):
            sharedCategoryCount = sharedCategoryCount + 1
        else:
            differentCategoryCount = differentCategoryCount + 1

    return([sharedCategoryCount, differentCategoryCount])

# Get number of edges with shared and different disease categories
femaleCategoryCounts = getCommonAndDifferentCategoryEdgeCounts(in_edgelist_f)
print()
print("Female same phenotype category count = " + str(femaleCategoryCounts[0]))
print("Female different phenotype category count = " + str(femaleCategoryCounts[1]))

maleCategoryCounts = getCommonAndDifferentCategoryEdgeCounts(in_edgelist_m)
print()
print("Male same phenotype category count = " + str(maleCategoryCounts[0]))
print("Male different phenotype category count = " + str(maleCategoryCounts[1]))

