from filepath import *
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import fcluster
import json


def exMotif(lig_cls, G):
    '''this function is to obtain the first stage scoring, ligand atoms clustering, and functional atoms clustering
       o input: (1) lig_cls: the clustering result generated directly from the function sns.clustermap()
                (2) G: the number of function groups of the ligand , type: int
       o output: return three dicts
                 the first one is the first stage scores comprising the clustering quality for each function group of ligand
                 the second one is the clustering number with its clustered ligand atoms
                 the last one is the clustering number which is mapped into ligand clustering numbers, and its clustered
                 function atoms
    '''
    # get the reordered data(i.e., the data after clustering)
    df_cls = lig_cls.data2d

    # assign the reordered column indices to the single clusters ranging from 1 to G
    cluster_colind = fcluster(lig_cls.dendrogram_col.linkage, t=G, criterion="maxclust")
    col_cls = [cluster_colind[c] for c in lig_cls.dendrogram_col.reordered_ind]

    # assign the reordered row indices to the single clusters ranging from 1 to G
    cluster_rowind = fcluster(lig_cls.dendrogram_row.linkage, t=G, criterion="maxclust")
    row_cls = [cluster_rowind[r] for r in lig_cls.dendrogram_row.reordered_ind]

    # obtain ligand function groups and corresponding functional atoms
    if sorted(col_cls) == col_cls and sorted(row_cls) == row_cls:
        # obtain the index of the first atom in each cluster
        y_axis = [row_cls.index(g) for g in range(1, G+1)]
        y_axis.append(len(row_cls))
        x_axis = [col_cls.index(g) for g in range(1, G+1)]
        x_axis.append(len(col_cls))

        # obtain function group atoms for each function group of the ligand
        ligatom_order = {"funcG{i}".format(i=g): list(df_cls.columns)[x_axis[g-1]:x_axis[g]] for g in range(1, G+1)}
        # score and map the clusters of fas into function groups
        Mat = []
        for i in range(G):
            for j in range(G):
                mat = df_cls.iloc[y_axis[j]:y_axis[j+1], x_axis[i]:x_axis[i+1]]
                # calculate the relative proportion of grids with the value no more than 0.4 in each square
                val = mat[mat <= 0.4].count().sum() / float(mat.shape[0] * mat.shape[1])
                Mat.append(val)
        val_mat = np.array(Mat).reshape(G, G).T
        val_max = val_mat.max(axis=1)
        pro_cls = val_mat.argmax(axis=1)
        # save the clustering quality as scores
        score = dict(zip([pp+1 for pp in pro_cls], val_max))

        if sorted(pro_cls) == list(range(G)):
                pro2lig_cls = {i: pro_cls[i-1]+1 for i in range(1, G+1)}
                proatom_order = {"cluster{i}".format(i=pro2lig_cls[g]):
                                     list(df_cls.index)[y_axis[g-1]:y_axis[g]] for g in range(1, G+1)}
                return score, ligatom_order, proatom_order
        # # filter out the bad clustering data
        # else:
        #     print("different rows have same column")


def motifGen(LigInfo = ligInfo, DistFolder=distFolder, HeatFolder=heatFolder, StatisFolder=statisFolder,
             ScoreFolder = scoreFolder):
    '''this function is to assess the clustering heatmap, corresponding clustering result for each of all ligands
       and the first stage scores for each of the ligands comprising 2 or more function groups
       o input: (1) ligInfo: the location of file recording all 233 ligands information
                (2) DistFolder: the location of file recording all 233 ligands information
                (3) HeatFolder: the location of folder where clustering heatmap of the ligand to be saved
                (4) StatisFolder: the location of the folder where derived clustering results txt
                                  recording atoms of function group and motif binding to be written
                (5) ScoreFolder: the location of folder placing excel files to be recording all 233 ligands
                                 clustering quality for the first stage
       o output: plot all ligands clustering heatmaps,
                 write all ligands' clustering results into their corresponding single .txt files,
                 and write the first stage scores of function group for ligands comprising 2 or more function groups
                 into .xlsx file
                 the clustering results of each ligand include:
                 a) ligand function group: clustered atoms for each function group of the ligand
                 b) protein binding motif: clustered function atoms for each motif
    '''
    df_info = pd.read_excel(LigInfo, index_col=0)
    d_score = {"ligands": [], "scores": []}
    for ligand in df_info.index:
        # read the distance matrix of the ligand
        df2 = pd.read_excel(os.path.join(DistFolder, ligand+"_dismat.xlsx"),
                            index_col=0)
        # obtain the the number of ligand function group
        g_num = df_info.loc[ligand]["Number"]
        # two-way hierarchical clustering on the distance matrix
        clsmap = sns.clustermap(df2, standard_scale=0, method="ward", metric="euclidean", cmap="bwr")
        # write the heatmap into .png file
        plt.savefig(HeatFolder + ligand + "_cls.png", dpi=300, bbox_inches='tight')
        try:
            # get score, ligand function groups, protein binding motifs
            s, lig_order, atom_order = exMotif(clsmap, g_num)
            # score ligands comprising 2 or more function groups
            if g_num >= 2:
                d_score["ligands"].append(ligand)
                d_score["scores"].append(s)
            Stat_dic = {ligand: lig_order, "pro": atom_order}
            data = json.dumps(Stat_dic)
            # write the ligand function groups and corresponding protein binding motifs into statistic folder
            with open(os.path.join(StatisFolder, ligand), "w") as fi:
                fi.write(data)

        except TypeError:
            print("%s has got bad cluster result" % ligand)

    # write the fist stage scoring for the clustering quality of the ligand into .xlsx
    pd.DataFrame(d_score).to_excel(os.path.join(ScoreFolder, "score1st.xlsx"), index=False)




def dic_num2fg(fuc, cls):
    '''this function is to map the clustered function group of the ligand atoms into the manually defined function group
       of the ligand
       o input: (1) fuc: manually defined function group atoms
                (2) cls: clustered function group atoms
       o output: return a dict storing the defined function group name with its mapped clustered number of function atoms
                 and the dict values include 2 section:
                 1) the manually defined function group names
                 2) a tuple scoring the match degrees between defined atoms and clustered atoms for each function
                    group, which is used for second scoring
    '''
    dic_lig = {}
    for k1, v1 in cls.items():
        s1 = set(v1)
        for k2, v2 in fuc.items():
            s2 = set(v2.split(","))
            # proportion for intersection of clustered atoms and defined atoms relative to clustered atoms
            a1 = len(s1 & s2) / len(s1)
            # proportion for intersection of clustered atoms and defined atoms relative to defined atoms
            a2 = len(s1 & s2) / len(s2)
            # if two proportions both are more than 0.5, the two sets can be matched together
            if a1 >= 0.5 and a2 >= 0.5:
                dic_lig["cluster"+str(k1[-1])] = [k2, (a1, a2)]
    return dic_lig


def assign_score(AtomsLigand=atomsLigand, StatisFolder=statisFolder, ScoreFolder=scoreFolder,
                 AssignFolder=assignFolder):
    '''this function is to obtain all ligands' function group names assignment to the clustered function groups
       and all ligands score combined with the first stage scores
       o input: (1) AtomsLigand: the location of excel file recording inherent atom composition for each
                    predefined function group of each ligand
                (2) StatisFolder: the location of the folder placing files recording clustered atoms of function group
                                  and motif binding
                (3) ScoreFolder: the location of folder placing excel files recording scores for ligands
                (4) AssignFolder: the location of folder placing excel file recording assignment of functional group
                                  names and cluster numbers
       o output: write the assignment of functional group names and cluster numbers for
                 each of all ligands into .xlsx file, and write the second scores combined with the first scores
                 for each of all ligands into .xlsx file.
                 return 2 dict recording scores and assignments for all ligands.
    '''
    atom_lib = pd.read_excel(AtomsLigand, index_col=0)
    df_score1st = pd.read_excel(os.path.join(ScoreFolder, "score1st.xlsx"), index_col=0)
    score2nd = {"ligands": [], "scores": []}
    assign = {"ligands": [], "assign": []}

    for ligand in atom_lib.index:
        assign["ligands"].append(ligand)
        score2nd["ligands"].append(ligand)
        # assign the ligands in df_score1st comprising 2 or more function groups, and score
        if ligand in df_score1st.index:
            value = []
            ass = {}
            fuc_atoms = eval(atom_lib.loc[ligand]["atoms"])
            with open(os.path.join(StatisFolder, ligand), "r") as f:
                clsted_atoms = json.load(f)[ligand]
            result = dic_num2fg(fuc_atoms, clsted_atoms)
            # make sure every cluster number is assigned to each atom
            if len(result) == len(clsted_atoms):
                score1 = eval(df_score1st.loc[ligand]["scores"])
                for g in sorted(score1.keys()):
                    # calculate the 2nd score combined with the 1st score
                    s2 = score1[g] * 0.5 + sum(result["cluster%d" % g][1]) * 0.25
                    value.append((result["cluster%d" % g][0], s2))
                    ass["cluster%d" % g] = result["cluster%d" % g][0]
                score2nd["scores"].append(value)
                assign["assign"].append(ass)
            else:
                print(ligand+" get wrong")

        # add ligands in which the function groups count=1, whose quality score=1
        else:
            fname = list(eval(atom_lib.loc[ligand]["atoms"]).keys())[0]
            assign["assign"].append({"cluster1": fname})
            score2nd[ligand].append([(fname, 1.0)])

    # write the assignment and score into corresponding excel files
    pd.DataFrame(assign).to_excel(os.path.join(AssignFolder, "assign.xlsx"), index_col=0)
    pd.DataFrame(score2nd).to_excel(os.path.join(ScoreFolder, "score2nd.xlsx"), index_col=0)


