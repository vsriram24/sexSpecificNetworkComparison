# The interplay of sex and genotype in disease associations: a comprehensive network analysis in the UK Biobank
Authored by Vivek Sriram, PhD, MA

## Introduction
This project is the final first-author paper I completed for my PhD in Biomedical Informatics and Computational Genomics at the University of Pennsylvania. My dissertation was broadly focused on the network-based integration of genomic and phenotypic data from electronic health record (EHR)-linked biobanks for the identification of genetic contributors to disease comorbidities. 

In this project, I developed sex-stratified disease-association networks from phenome-wide association study (PheWAS) data from the UK Biobank for the identification of genotype-by-sex (GxS) effects in disease comorbidities. Networks were generated in the two following ways: (1) An edge was constructed between two phenotypes if they shared common associated genetic variants (single nucleotide polymorphisms / SNPs) according to the PheWAS. (2) An edge was constructed between two phenotypes if a statistically significant genetic correlation was identified through the application of linkage disequilibrium score regression (LDSR).

After generating my two sex-stratified networks, I applied network comparison methods to identify GxS effects in the cross-phenotype associations. This code sample includes scripts for the following three network comparison methods: (1) direct comparisons of edge sets between the two graphs through Jaccard Similarity, (2) comparison of clustering behavior through K-means/Louvain clustering after dimensionality reduction, and (3) comparison of network embeddings through the application of the NLP-inspired node2vec representation learning algorithm and temporal embedding matching.

The corresponding manuscript for this work is currently under review with _Genome Medicine_. The attached files provide an overview of the project as well as a sample of some of the data processing and analysis steps that were involved. Since this work is still being reviewed for publication, not all data and scripts for the project have been shared.

Key elements presented in this code sample include:
- Generating sex-stratified networks representing genetic associations between diseases from large-scale EHR-linked biobank data
- Directly comparing the two sex-stratified networks by their edge sets and clustering behavior
- Applying graph representation learning using node2vec, and aligning generated node embeddings using temporal embedding matching to compare networks to each other in a common embedding space.


## Directory Organization
Directories are organized as follows:

1. **Background (background)**
This folder includes background information to provide context for the research project and analysis.
- projectAbstract.docx: An abstract for the manuscript summarizing this research project. This document provides a summary of the biological context of the work and highlights a few of the final translational takeaways.
- methodOverview.png: An overview figure highlighting different methods of network comparison that were applied in this research project.

2. **Sample Data (sampleData)**
This folder includes sample input and output data to give context for the data used in and output by the analyses scripts.
- inData:
	- femaleSummaryData/: PheWAS summary data for the female UKBB population. Each file corresponds to a distinct phenotype, and each row of each file corresponds to association information between a SNP and the current phenotype.
	- maleSummaryData/: PheWAS summary data for the male UKBB population.
	- sampleCounts.csv: Case/Control counts for the male and female populations for each disease.
	- traitDescription_diseaseCategories.csv: Phenotype descriptions and disease categories for each phenotype.
	
- outData:
	- embeddingDistances.csv: Distances between the embeddings for each disease between the male- and female-specific versions of the data.
	- embeddingVisualizations.pdf: Knitted version of visualizeEmbeddingsInGgplot.Rmd, depicting differences in graph embedding patterns between the sexes
	- female_assocMat.csv: Association matrix for the female-specific disease network, including normalized distances for each pair of diseases in the network. A value of 0 indicates no connection between the phenotypes, while a value of 1 indicates the strongest connection in the network.
	- female_edgeMap.tsv: Edgemap for the female-specific disease network. Each row includes a source node ID, a target node ID, and a list of shared SNPs between the two diseases.
	- male_assocMat.csv: Association matrix for the male-specific disease network.
	- male_edgeMap.tsv: Edgemap for the male-specific disease network.
	- rgResults/: Log file outputs of Linkage Disequilibrium Score Regression.

3. **Scripts (scripts)**
This folder includes samples of the code used for data analysis. Scripts are broken down into five stages: (1) Data Processing, (2) Genetic Correlations, (3) Interactive Visualization, (4) Node2Vec, and (5) Network Analysis. 
- 1_processData/: Convert source PheWAS data into female- and male-specific disease networks.
- 2_geneticCorrelation/: Apply Linkdage Disequilibrium Score Regression (LDSR) to identify genetic correlations between phenotypes.
- 3_interactiveVis/: Prepare network files for interactive visualization in the Gephi network visualization API.
- 4_node2vec/: Generate node embeddings from sex-stratified disease networks using node2vec.
- 5_networkAnalysis/: Evaluate sex-stratified networks against one another with edge set, clustering, and node embedding comparisons.

More detail regarding the individual functions of each script in each directory can be found in the **Code Overview** section.


## Code Overview

### 1. Process Data
- 1_splitDataBySex.py: This Python script takes the original source data, a directory of files where each file corresponds to SNP-association data for a single phenotype for both sexes, and splits each file into two (one file for each sex per phenotype).

- 2_generateAdjacencyMatrix.py: This Python script takes PheWAS data for a single sex and generates a corresponding adjacency matrix to represent the network of shared genetic associations across diseases. Edge weights are determined according to the number of SNPs that are shared between each pair of traits, and then cosine similarity is applied to normalize edge weights to between 0 and 1.


### 2. Genetic Correlation
- 1_runLDSC.sh: This shell script facilitates the application of LDSR through batch processing to calculate genetic correlations between all pairs of phenotypes.

- 2_scrapeLDSCOutput.py: This Python script takes the standardized output of LDSR and collects pairwise genetic correlation values into a network of edges defined by the correlations between all pairs of phenotypes.


### 3. Interactive Visualization
- 1_makeNodeAndEdgeMaps.py: This Python script takes PheWAS data for a single sex and converts the information into a node map and edge map that can be plugged into the Gephi user interface for interactive network visualization.


### 4. Node2Vec
- 1_node2vec.py: This Python script applies the NLP-inspired graph representation learning method "node2vec" to define vectors for each phenotype based upon the topology of the input graph adjacency matrix. The script also includes methodology for temporal embedding matching, allowing us to map generated vectors across the two graphs to one another and compare how diseases differ from one another across the sexes in terms of their co-associations with other phenotypes.


### 5. Network Analysis
- 1_compareEmbeddingDistances.ipynb: This iPython notebook includes exploratory data analysis (EDA) to evaluate clustering behavior and embedding distances across the male- and female-specific networks.

- 2_visualizeDDNEmbeddings.ipynb: This iPython notebook includes EDA to visualize clustering behavior and embedding patterns across the two networks.

- 3_visualizeEmbeddingsInGgplot.Rmd: This RMarkdown document replicates the visualizations completed within visualizeDDNEmbeddings.ipynb using ggplot2.

- 4_calculateJaccardSimilarity.py: This Python script compares edge sets for the male- and female-specific graphs to provide insight into similarity of the types of edges that appear across the two networks.


## Software Dependencies
Scripts were run using Python 3.11.3 and R version 4.0.2. 
