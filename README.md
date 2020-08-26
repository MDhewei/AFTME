# AFTME

AFTME is an alignment-free method for the automatic mapping of 3D motifs to different FGs of a specific ligand through two-dimensional clustering.

## 

[1.Motif_all](https://github.com/MDhewei/AFTME/tree/master/%20Motif_all): This folder contains the heatmaps, statistics results and figures for all the FG-binding motifs.

- **Heatmaps**: Heatmaps for all the ligands generated with AFTME 

- **motif_analyse**: Bar plots, pie plots and venn plots to visualize different FG-binding motifs

- **statistics**: Dict files recording the motif informations for all the ligands

[2.Datasets](https://github.com/MDhewei/AFTME/tree/master/Datasets): This folder contains the datasets used for large-scale FG-motif identification.

- **NonredundantLigID**: The .txt files containing all the PDB codes used for motif extraction for each ligand

- [Ligand_information.xlsx](https://github.com/MDhewei/AFTME/blob/master/Datasets/Ligand_information.xlsx): The excel file recording the information: name, functional groups, 2d figures etc. for all the ligands involved in our analysis.

- [PDBCodes.xlsx](https://github.com/MDhewei/AFTME/blob/master/Datasets/PDBCodes.xlsx): The excel file recording all the PDB codes used for analysis

[3.DistanceMatrix](https://github.com/MDhewei/AFTME/tree/master/DistanceMatrix): This folder contains the excel files recording all the distance matrix generated for all the ligands.

[4.Systematic_analysis](https://github.com/MDhewei/AFTME/tree/master/Systematic_analysis): This folder contains the files and source codes to reproduce the results for systematic analysis.

- [Cluster_file](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/Cluster_file): Dict file recording the motifs in each cluster.

- [FG_binding_motif_all.csv](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/FG_binding_motif_all.csv): Table recroding the information: ID, atom properties, aa compositions, etc. for all the identified motifs.

- [motif_assign.xlsx](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/motif_assign.xlsx): Table recording the clustering assignments for each FG in the analysis.

- [Motif_clustering.ipynb](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/Motif_clustering.ipynb): Source codes for the motif clustering: feature matrix construction, "elbow" plot, t-SNE plot, etc. and characterization (bar and pie plots).

- [ClusteringAnalysis.ipynb](https://github.com/MDhewei/AFTME/blob/master/Systematic_analysis/ClusteringAnalysis.ipynb): Source codes for motif analysis within each motif cluster: FG-motif overlaps, distribution of different binding modes, affinity comparison, etc. 















