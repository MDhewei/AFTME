# AFTME

AFTME is an alignment-free method for the automatic mapping of 3D motifs to different FGs of a specific ligand through two-dimensional clustering. This repository contains the source codes, datasets and results for large-scale identification of FG-based protein-ligand binding motifs.

Any questions or requests, please contact us: hwkobe.1027@gmail.com or yangliupoplar@gmail.com

## AFTME source codes and usage

### Requirements:
- AFTME is written in Python, Python>=3.6 is needed. 

- Python Packages: Biopython, pandas, seaborn, Matplotlib, SciPy, NumPy

#### [1. DistMat.py](https://github.com/MDhewei/AFTME/blob/master/DistMat.py): Source codes to create distance matrix for each of ligands from given liganfinfo and nonredundant pdbID.

#### Required arguments of the program:
* -i/--inputfile:
the 2 paths in order:   
> the absolute path of an excel file recording all 233 ligands information. 

> the absolute path of a folder containing txt files recording ligand nonredundant pdbID for each of 233 ligands.

* -o/--outfile:  
the absolute path of a folder to store distance matrix excel files created for each of 233 ligands.

#### Example to run the Dismat.py:
```
python DistMat.py -i "Ligand_information.xlsx" "NonredundantLigID_folder" -o "DistanceMatrix_folder"
```

#### [2. AFTME.py](https://github.com/MDhewei/AFTME/blob/master/AFTME.py): Source codes to extract FG-binding motifs based on ligand information and distance matrix.

#### there are 2 functions: motifGen and assign_score. The second function takes the output fils of the 1st function as input.  
#### - MotifGen: obtain the clustering heatmaps, corresponding clustering results for each of all 233 ligands and the first stage scores for each of the ligands comprising 2 or more function groups. 

##### Required arguments of motifGen:
* -i/--inputfile:
the 2 paths in order:  

> the absolute path of an excel file recording all 233 ligands information  

> the absolute path of a folder placing distance matrix excel files created for each of 233 ligands

* -o/--outfile:
the 3 paths in order:  

> the absolute path of a folder for saving clustering 233 heatmaps of the 233 ligands 

> the absolute path of a folder for saving clustering result texts recording atoms of function group and binding motifs 

> the absolute path of a folder for saving excel file recording the first stage scores of 181 ligands with 2 or more funtion groups out of 233 ligands.   

##### Example to run the motifGen:
```
python AFTME.py motifGen -i "Ligand_information.xlsx" "DistanceMatrix_folder" -o "heatmapfolder" "statisfolder" "scorefolder"
```
#### - assign_score: obtain assignment between artificially defined function groups and clustered function groups of each of 233 ligands, and all ligands scores based on the first stage scores. 

##### Required arguments of assign_score:  
* -i/--inputfile:
the 3 paths in order:  

> the absolute path of an excel file recording inherent atom composition for each manualy predefined function group of each of 233 ligands.

> the absolute path of a folder placing txt files recording clustered atoms of function group and binding motif for each of 233 ligands, which is genetated from motifGen.

> the absolute path of a folder placing excel files recording 1st stage scores, which is generated from motifGen.(note: the folder is also used for saving the second stage scores of outfile)

* -o/--outfile:    
the absolute path of a folder for saving an excel file recording assignment between artificially defined function groups and clustered function groups of each of 233 ligands.  

##### Example to run the assign_score:
```
python AFTME.py assign_score -i "atom_233.xlsx" "statisfolder" "scorefolder" -o "assignfolder"
```

## Datasets and large-scale-analysis results

[1. Motif_all](https://github.com/MDhewei/AFTME/tree/master/%20Motif_all): This folder contains the heatmaps, statistics results and figures for all the FG-binding motifs.

- **Heatmaps**: Heatmaps for all the ligands generated with AFTME 

- **motif_analyse**: Bar plots, pie plots and venn plots to visualize different FG-binding motifs

- **statistics**: Dict files recording the motif informations for all the ligands

[2. Datasets](https://github.com/MDhewei/AFTME/tree/master/Datasets): This folder contains the datasets used for large-scale FG-motif identification.

- **NonredundantLigID**: The .txt files containing all the PDB codes used for motif extraction for each ligand

- [Ligand_information.xlsx](https://github.com/MDhewei/AFTME/blob/master/Datasets/Ligand_information.xlsx): The excel file recording the information: name, functional groups, 2d figures etc. for all the ligands involved in our analysis.

- [PDBCodes.xlsx](https://github.com/MDhewei/AFTME/blob/master/Datasets/PDBCodes.xlsx): The excel file recording all the PDB codes used for analysis

[3. DistanceMatrix](https://github.com/MDhewei/AFTME/tree/master/DistanceMatrix): This folder contains the excel files recording all the distance matrix generated for all the ligands.

[4. Systematic_analysis](https://github.com/MDhewei/AFTME/tree/master/Systematic_analysis): This folder contains the files and source codes to reproduce the results for systematic analysis.

- [Cluster_file](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/Cluster_file): Dict file recording the motifs in each cluster.

- [FG_binding_motif_all.csv](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/FG_binding_motif_all.csv): Table recroding the information: ID, atom properties, aa compositions, etc. for all the identified motifs.

- [motif_assign.xlsx](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/motif_assign.xlsx): Table recording the clustering assignments for each FG in the analysis.

- [Motif_clustering.ipynb](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/Motif_clustering.ipynb): Source codes for the motif clustering: feature matrix construction, "elbow" plot, t-SNE plot, etc. and characterization (bar and pie plots).

- [ClusteringAnalysis.ipynb](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/ClusteringAnalysis.ipynb): Source codes for motif analysis within each motif cluster: FG-motif overlaps, distribution of different binding modes, affinity comparison, etc. 

















