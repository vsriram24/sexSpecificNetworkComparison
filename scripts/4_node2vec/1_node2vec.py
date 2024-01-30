#####################################
# Title: node2vec.py
# Author: Vivek Sriram
# Last revised: 04/28/23
# -----------------------------------------------------------------------
# Function: Node2vec is a machine learning framework to learn feature
#           embeddings for nodes in a graph. It is derived from
#           the word2vec method in NLP for deriving vector
#           representations of words in a sentence. For node2vec, paths
#           of nodes are represented as sentences, and the context of
#           each node's connections is represented in a vector.
#
#           Here, we apply node2vec to our two graphs to generate 
#           a female-specific and male-specific vector representation of
#           each disease in our dataset. We then apply a method known as
#           "temporal embedding matching" to align our vector representations
#           across the sexes to the same axes. Finally, we use t-SNE dimensionality
#           reduction to visualize our aligned embeddings. These embeddings
#           can be compared across the sexes to evaluate differences in
#           disease connectivity patterns. 

# Import packages
import os
import pandas as pd
import numpy as np
import networkx as nx
from node2vec import Node2Vec
import csv
from scipy.spatial import distance
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from sklearn.cluster import KMeans
from umap import UMAP
import plotly.express as px
from sklearn.manifold import TSNE

# Parameters for node2vec random walks
returnParam = 1.0 # Options: 0.5, 1, 2, 4
inOutParam = 1.0
outFilePath = '/cubic/projects/kimlab/projects/viveksrm_DDNComp/embeddingDistances.csv'
outCategoryPath = '/cubic/projects/kimlab/projects/viveksrm_DDNComp/categoryEmbeddingDistances.csv'
numRuns = 100

# Parameters for node2vec random walks
dim = 16 # Options: {16, 32, 64, 128}
walkLength = 20 # diameter for female ssDDN is 4, diameter for male ssDDN is 6
n_walks = 10000

# Parameters for node2vec embedding calculations from random walks
slides = 10
numEpochs = 20

# Input DDNs (pandas format)
#A. ssDDN (1e-4, female) - 104 nodes, 676 edges
in_ssDDN_pandas_female = pd.read_csv('/cubic/projects/kimlab/projects/viveksrm_DDNComp/runningEmbeddings/femaleBlock_assocMat_neg4.csv', index_col = 0)
#B. ssDDN (1e-4, male) - 104 nodes, 598 edges
in_ssDDN_pandas_male = pd.read_csv('/cubic/projects/kimlab/projects/viveksrm_DDNComp/runningEmbeddings/maleBlock_assocMat_neg4.csv', index_col = 0)

distanceDict = {}

"""
# Temporal Embedding Matching
:param emb_to_align: embedding vectors to be align
:param emb_base: base embedding vectors
:return: aligned_embeddings of emb_to_align

Example of embedding's data structure
    embs = {
        a: [0,0,1,....,0],
        b: [1,0,0,....,0],
        ...
    }
"""
def align_two_embs(emb_to_align, emb_base, common_keys=None):
    if not common_keys:
        common_keys = list(set(emb_to_align.keys()).intersection(set(emb_base.keys())))

    A = np.array([emb_to_align[key] for key in common_keys]).T
    B = np.array([emb_base[key] for key in common_keys]).T
    M = B.dot(A.T)
    u, sigma, v_t = np.linalg.svd(M)
    rotation_matrix = u.dot(v_t)
    aligned_embedding = {k: rotation_matrix.dot(v) for k, v in emb_to_align.items()}

    return aligned_embedding



###############################################################################################
## PREPARE DISEASE-DISEASE NETWORKS
print("Processing imported male and female ssDDNs...")
print()

# Get the set of diseases that are common between the two adjacency matrices
femaleRowNames = in_ssDDN_pandas_female.index
maleRowNames = in_ssDDN_pandas_male.index
commonDiseases = list(set(femaleRowNames) & set(maleRowNames))

# Filter down the female adjacency matrix to only the common diseases
in_ssDDN_pandas_female = in_ssDDN_pandas_female[in_ssDDN_pandas_female.columns.intersection(commonDiseases)]
femaleDiseasesToDrop = list(set(femaleRowNames) - set(commonDiseases))
for elem in femaleDiseasesToDrop:
    in_ssDDN_pandas_female = in_ssDDN_pandas_female.drop(elem)
in_ssDDN_pandas_female[in_ssDDN_pandas_female < 0] = 0

# Filter down the male adjacency matrix to only the common diseases
in_ssDDN_pandas_male = in_ssDDN_pandas_male[in_ssDDN_pandas_male.columns.intersection(commonDiseases)]
maleDiseasesToDrop = list(set(maleRowNames) - set(commonDiseases))
for elem in maleDiseasesToDrop:
    in_ssDDN_pandas_male = in_ssDDN_pandas_male.drop(elem)
in_ssDDN_pandas_male[in_ssDDN_pandas_male < 0] = 0

# Convert pandas to numpy
ssDDN_numpy_female = in_ssDDN_pandas_female.to_numpy()
ssDDN_numpy_male = in_ssDDN_pandas_male.to_numpy()

# Convert numpy to networkx
G_ssDDN_female = nx.from_numpy_array(ssDDN_numpy_female)
G_ssDDN_male = nx.from_numpy_array(ssDDN_numpy_male)


for i in range(1, numRuns):
    print("On iteration " + str(i) + " :")
    # Set random seed (for random walks and for embedding calculations from random walks)
    randomSeed = i

    ###############################################################################################
    ## SIMULATE RANDOM BIASED WALKS
    print("Simulating random walks for female ssDDN...")
    G_f_n2v = Node2Vec(G_ssDDN_female, # Networkx graph
                        dimensions = dim, # Embedding dimensions (default: 128)
                        walk_length = walkLength, # Number of nodes in each walk (default: 80)
                        num_walks = n_walks, # Number of walks per node (default: 10)
                        p = returnParam, # Return hyperparameter to promote BFS (default: 1)
                        q = inOutParam, # In-out parameter to promote DFS (default: 1)
                        seed = randomSeed) # seed for the random number generator (defualt: None)
    print("Simulating random walks for male ssDDN...")
    G_m_n2v = Node2Vec(G_ssDDN_male, # Networkx graph
                        dimensions = dim, # Embedding dimensions (default: 128)
                        walk_length = walkLength, # Number of nodes in each walk (default: 80)
                        num_walks = n_walks, # Number of walks per node (default: 10)
                        p = returnParam, # Return hyperparameter to promote BFS (default: 1)
                        q = inOutParam, # In-out parameter to promote DFS (default: 1)
                        seed = randomSeed) # seed for the random number generator (defualt: None)
    print()


    ###############################################################################################
    ## LEARN NODE EMBEDDINGS FROM GENERATED RANDOM WALKS
    print("Learning node embeddings for female ssDDN...")
    model_f = G_f_n2v.fit(min_count = 1, # ignore all words with total absolute freuqency lower than this; minimum count of words to consider when training the model
                        window = slides, # max distance between the target node and its neighbors
                        epochs = numEpochs, # Number of iterations (epochs) over the corpus
                        seed = randomSeed) # seed for the random number generator

    print("Learning node embeddings for male ssDDN...")
    model_m = G_m_n2v.fit(min_count = 1, # ignore all words with total absolute freuqency lower than this; minimum count of words to consider when training the model
                        window = slides, # max distance between the target node and its neighbors
                        epochs = numEpochs, # Number of iterations (epochs) over the corpus
                        seed = randomSeed) # seed for the random number generator
    print()

    embedding_matrix_f = model_f.wv.vectors
    embedding_matrix_m = model_m.wv.vectors


    ###############################################################################################
    ## PERFORM TEMPORAL MATCHING
    # Prepare embedding dictionaries
    print("Generating embedding dictionaries for temporal matching...")
    female_embedding_dict = {}
    for i in range(0, len(in_ssDDN_pandas_female.index)):
         female_embedding_dict[in_ssDDN_pandas_female.index[i]] = embedding_matrix_f[i]

    male_embedding_dict = {}
    for i in range(0, len(in_ssDDN_pandas_male.index)):
         male_embedding_dict[in_ssDDN_pandas_male.index[i]] = embedding_matrix_m[i]
    print()


    print("Calculating new alignment for the secondary embedding...")
    # Use the male embedding as a base and align the female embedding to the male's
    alignedEmbeddingForSecondaryDDN = align_two_embs(emb_to_align = female_embedding_dict, 
                                                    emb_base = male_embedding_dict, 
                                                    common_keys = commonDiseases)
    print()

    print("Calculating distances between the embeddings for the two ssDDNs...")
    for key in alignedEmbeddingForSecondaryDDN:
        femaleVector = alignedEmbeddingForSecondaryDDN[key]
        maleVector = male_embedding_dict[key]
        dist = distance.cosine(femaleVector, maleVector)

        if key not in distanceDict.keys():
            distanceDict[key] = [dist]
        else:
            distanceDict[key].append(dist)

    print("Saved results!")
    print()


traitDescriptions = pd.read_csv('/cubic/projects/kimlab/projects/viveksrm_DDNComp/traitDescription_diseaseCategories.csv')
diseaseCategoryMap = {}
for index, row in traitDescriptions.iterrows():
    disease = row['Trait']
    category = row['Disease Category']
    diseaseCategoryMap[disease] = category

category_outputDistanceDict = {}


print("Calculating average distances across runs...")
# Calculate average and std deviation of distances for each key
outputDistanceDict = {}

for key in distanceDict.keys():
    currentListOfDistances = distanceDict[key]
    currentDiseaseCategory = diseaseCategoryMap[key]

    for elem in currentListOfDistances:
        if currentDiseaseCategory in category_outputDistanceDict:
            category_outputDistanceDict[currentDiseaseCategory].append(elem)
        else:
            category_outputDistanceDict[currentDiseaseCategory] = [elem]

    avgDist = np.mean(currentListOfDistances)

    stdev = np.std(currentListOfDistances)
    outputDistanceDict[key] = [avgDist, stdev]


print("Saving to output file " + outFilePath)
with open(outFilePath, 'w') as outFile:
    outFile.write("disease" + "\t" + "meanDist" + "\t" + "stdDevDist" + "\n")
    for key in outputDistanceDict.keys():
        outFile.write(str(key) + "\t" + str(outputDistanceDict[key][0]) + "\t" + str(outputDistanceDict[key][1]) + "\n")
    outFile.close()


for elem in category_outputDistanceDict:
    category_outputDistanceDict[elem] = [np.mean(category_outputDistanceDict[elem]), np.std(category_outputDistanceDict[elem])]

print("Saving to output file " + outCategoryPath)
with open(outCategoryPath, 'w') as outCategory:
    outCategory.write("diseaseCategory" + "\t" + "meanDist" + "\t" + "stdDevDist" + "\n")
    for key in category_outputDistanceDict.keys():
        outCategory.write(str(key) + "\t" + str(category_outputDistanceDict[key][0]) + "\t" + str(category_outputDistanceDict[key][1]) + "\n")
    outCategory.close()
print("All done!")



###############################################################################################
## VISUALIZE EMBEDDINGS

# Convert matrices to pandas dataframes
embedding_matrix_m_pd = pd.DataFrame(embedding_matrix_m)
embedding_matrix_f_pd = pd.DataFrame(embedding_matrix_f)
aligned_matrix_f_pd = pd.DataFrame(alignedEmbeddingForSecondaryDDN)

# Perform t-SNE to visualize the embeddings in 2D
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results_maleEmbeddingOriginal = tsne.fit_transform(embedding_matrix_m_pd)
tsne_results_femaleEmbeddingOriginal = tsne.fit_transform(embedding_matrix_f_pd)
tsne_results_femaleEmbeddingAligned = tsne.fit_transform(aligned_matrix_f_pd)

print('t-SNE done!')

# Save t-SNE results to the pandas dataframes
embedding_matrix_m_pd['tsne-2d-one'] = tsne_results_maleEmbeddingOriginal[:,0]
embedding_matrix_m_pd['tsne-2d-two'] = tsne_results_maleEmbeddingOriginal[:,1]

embedding_matrix_f_pd['tsne-2d-one'] = tsne_results_femaleEmbeddingOriginal[:,0]
embedding_matrix_f_pd['tsne-2d-two'] = tsne_results_femaleEmbeddingOriginal[:,1]

aligned_matrix_f_pd['tsne-2d-one'] = tsne_results_femaleEmbeddingAligned[:,0]
aligned_matrix_f_pd['tsne-2d-two'] = tsne_results_femaleEmbeddingAligned[:,1]
    
# Plot t-SNE using seaborn
plt.figure(figsize=(16,4))
ax1 = plt.subplot(1, 3, 1)
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-one",
    data=embedding_matrix_m_pd,
    legend="full",
    alpha=0.3,
    ax=ax1
)
    
ax2 = plt.subplot(1, 3, 2)
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    data=embedding_matrix_f_pd,
    legend="full",
    alpha=0.3,
    ax=ax2
)
    
ax3 = plt.subplot(1, 3, 3)
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    data=aligned_matrix_f_pd,
    legend="full",
    alpha=0.3,
    ax=ax3
)