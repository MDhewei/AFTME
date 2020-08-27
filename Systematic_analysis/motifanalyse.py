from __future__ import division
from filepath import *
import os, itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2_unweighted, venn3_unweighted
from scipy import stats
import json

#  set  figure font and resolution.
params = {"font.sans-serif": 'Helvetica', "font.family": "sans-serif", "figure.dpi": 300}
plt.rcParams.update(params)

atomtypes = ['hydrophilic', 'acceptor', 'donor', 'hydrophobic', 'aromatic', 'neutral']
AminoAcids = ['ALA', 'GLY', 'PRO', 'ILE', 'LEU', 'VAL', 'PHE', 'TYR', 'TRP', 'ASN', 'CYS',
              'MET', 'THR', 'SER', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS']


def GetMotif(moti):
    '''this function is to extract the pdb code, amino acids, atom properties from motif
       o input: (1) moti: a list of atoms labeled as "pdbID_AminoAcid+AminoAcidID_atom"
       o output: return 3 lists.
                 the first one is a list of pdbcode occurring in the given motif
                 the second one is a list of amino acids occurring in the given motif
                 the third one is a list of atom properties occurring in the given motif
    '''
    Protein = []
    AA = []
    Prop = []
    for atom in moti:
        Protein.append(atom[:4])
        aa = atom[5:atom.rindex("_")][:3]
        at = atom[atom.rindex("_") + 1:]
        if aa[:2] in ["CA", "FE", "MN", "ZN", "MG", "CU"]:
            aa = " " + aa[:2]
            at = aa[1:3]
        AA.append(aa)
        if aa in atom_type.keys() and at in atom_type[aa].keys():
            pp = atom_type[aa][at]
            Prop.append(pp)
    # filter out the pdbcode that appears little in the given motif
    Prop_in = [pro for pro in set(Protein) if Protein.count(pro) >= 5]

    return Protein, AA, Prop


def PlotBarChart(ax, moti_AA):
    '''this function is to make a horizontal bar plot displaying the amino acids distribution of one motif
       o input: (1) ax: a single Axes object
                (2) moti_AA:  a list of amino acids occurring in the given motif
       o output: return a bar plot object
    '''
    labels = []
    values = []
    y_cor = []
    y = 0.1
    for i in ['ALA', 'GLY', 'PRO', 'ILE', 'LEU', 'VAL', 'PHE', 'TYR', 'TRP', 'ASN',
              'CYS', 'MET', 'THR', 'SER', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS']:
        values.append(moti_AA.count(i) / len(moti_AA))
        labels.append(i)
        y_cor.append(y)
        y += 0.3

    ax.barh(y_cor, values, height=0.1, tick_label=labels)
    ax.set_xlabel('Proportion of different amino acid')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis="y", labelsize=8)


def PlotPieChart(ax, props):
    '''this function is to make a pie chart displaying the proportion of atom properties in one motif
       o input: (1) ax: a single Axes object
                (2) props: a list of atom properties occurring in the given motif
       o output: return a pie chart object
    '''
    labels = []
    sizes = []
    for i in list(set(props)):
        labels.append(prop_dic[i])
        sizes.append(props.count(i))
    colors = ['gold', 'yellowgreen', 'lightcoral', 'grey', 'silver', 'white', 'lightskyblue', 'yellow']
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=45, colors=colors[0:len(sizes)])


def motifplot(G, proatom_order, lig_name, MotifFolder):
    '''this function is to make a picture combined bar plot with pie chart characterizing the motif biochemical properties,
       and make a venn plot to show overlaps of protein pockets containing different motifs for one ligand
       o input: (1) G: the number of function groups of the ligand , type: int
                (2) proatom_order: a dict of different binding motifs for a ligand
                (3) lig_name: the ligand name
                (4) MotifFolder: the location of a folder to place plots that display motif features
       o output: save the pictures into motifFolder folder.
                 return 2 pictures, one is a picture comprising bar plot and pie chart,
                 and the other one is a venn plot
                 note: if the number of function groups for one ligand is more than 3, there will be only one picture.
    '''
    fig, axes = plt.subplots(G, 2, figsize=(8, 3 * G), dpi=300)
    fig1, axes1 = plt.subplots(1, 1)
    atomorder = sorted(proatom_order.items())
    cluster_set = []
    labels = []
    # plot bar and pie chart displaying motif amino acids distribution and atom properties distribution respectively
    for cluster in atomorder:
        A, aminoAcids, props = GetMotif(cluster[1])
        PlotBarChart(axes[atomorder.index(cluster), 0], aminoAcids)
        PlotPieChart(axes[atomorder.index(cluster), 1], props)
        cluster_set.append(set(A))
        labels.append(cluster[0])
    fig.savefig(os.path.join(MotifFolder, "%s_analysis.png" % lig_name), dpi=300, bbox_inches="tight")

    # make a venn plot to show overlaps of protein pockets containing 2 motifs for one ligand
    if G == 2:
        venn2_unweighted(cluster_set, labels, ax=axes1)
        fig1.savefig(os.path.join(MotifFolder, "%s_provenn.png" % lig_name), dpi=300)

    # make a venn plot to show overlaps of protein pockets containing 3 motifs for one ligand
    elif G == 3:
        venn3_unweighted(cluster_set, labels, ax=axes1)
        fig1.savefig(os.path.join(MotifFolder, "%s_provenn.png" % lig_name), dpi=300)


# read the function group and ligands that contain the function group
df_fg = pd.read_excel(FG_inLigand, index_col=0)
# read the ligand and its function groups with matched clustering number
df_assign = pd.read_excel(os.path.join(assignFolder, "assign.xlsx"), index_col=0)


def set_bar(axes, labels, label_num, y_label):
    '''this function is to set the plot appearance
       o input: (1) axes: a single Axes object
                (2) labels: a list
                (3) label_num: int
                (4) y_label: string
    '''
    # remove right side and top side of frame
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    # manually set Axes by size and position
    axes.set_xticks([0.52 + 1.5 * i for i in range(label_num)])
    axes.set_xticklabels(labels, fontsize=15)
    axes.set_ylabel(y_label, size=20)


def freq_list(name):
    '''this function is to calculate the frequencies of amino acids occurring in the given motif, and the frequencies
       of atom properties occurring in the given motif
       o input: (1)name : string in scheme with: "FunctionGroup_ligand_clusterID", e.g., "Ribose_GDU_cluster2"
       o output: return 2 list
                 the first one is a list of frequencies of amino acids
                 the second one is a list of frequencies of atom properties
    '''
    # extract the ligand name and the cluster ID given the function group
    lig, cls = name[name.rindex("_")-3:name.rindex("_")], name[-8:]
    # obtain the motif that binding ligand's function group specified
    data = json.loads(open(os.path.join(statisFolder, lig)).read())["pro"][cls]
    # extract the motif amino acids, atom properties and pdb code
    pro, aa, pp = GetMotif(data)
    # calculate the frequencies of amino acids
    aa_freq = [aa.count(AA) / len(aa) for AA in AminoAcids]
    #  calculate the frequencies of atom properties
    pp_freq = [pp.count(num) / len(pp) for num in range(1, 7)]
    return aa_freq, pp_freq


def func_cls(func):
    '''this function is to obtain all motifs binding the given function group
       o input: (1) string, denotes the function group
       o output: return a list of motifs which each is named with
                 "FunctionGroup_ligand_clusterID", e.g., "Ribose_GDU_cluster2"
    '''
    lst = []
    data = eval(df_fg.loc[func]["ligand"])
    for i in data:
        for j in data[i]:
            lst.append(func + "_" + i + "_" + j)
    return lst


def fg_reuse(func, fuc_list, labels_lig):
    '''this function is to visualize the reusability of different motifs binding the given function group
       o input: (1) func: string, denotes the function group
                (2) func_list: a list of motif name that binds the function group
                (3) labels_lig: a list of ligand names
       o output: return 3 plot given the function group
                 the fist one is to display the consistency of 20 amino acids distribution for given motifs
                 the second one is to dispaly the consistency of 6 atop properties distribution for given motifs
                 the third one is to  display the consistency of integration atom and amino acid distribution for all
                 motifs binding the given function groups
    '''
    loc1 = np.arange(0.2, 29, step=1.5)
    # calculate the frequencies of the given motifs
    h1_aa, h1_prop = freq_list(fuc_list[0])
    h2_aa, h2_prop = freq_list(fuc_list[1])
    h3_aa, h3_prop = freq_list(fuc_list[2])
    # create a grouped bar chart, and each group displays same amino acids proportions in the different motifs.
    fig1, axes1 = plt.subplots(1, 1, figsize=(10, 6))
    ax1 = axes1.bar(loc1, h1_aa, width=0.1, color='lightcoral', ec="k", lw=0.5)
    ax2 = axes1.bar(loc1+0.2, h2_aa, width=0.1, color='lightskyblue', ec="k", lw=0.5)
    ax3 = axes1.bar(loc1+0.4, h3_aa, width=0.1, color='gold', ec="k", lw=0.5)
    axes1.legend([ax1, ax2, ax3], labels_lig, fontsize='small', loc="best")
    set_bar(axes1, AminoAcids, 20, "Amino Acids")

    loc2 = np.arange(0.2, 8, step=1.5)
    # create a grouped bar chart, and each group displays same atom properties proportions in the different motifs.
    fig2, axes2 = plt.subplots(figsize=(10, 6))
    ax4 = axes2.bar(loc2, h1_prop, width=0.1, color='lightcoral', ec="k", lw=0.5)
    ax5 = axes2.bar(loc2+0.2, h2_prop, width=0.1, color='lightskyblue', ec="k", lw=0.5)
    ax6 = axes2.bar(loc2+0.4, h3_prop, width=0.1, color='gold', ec="k", lw=0.5)
    axes2.legend([ax4, ax5, ax6], labels_lig, fontsize='small', loc="best")
    set_bar(axes2, atomtypes, 6, "Atom Properties")

    # obtain all motifs binding the function group
    fuc_lig_cls = func_cls(func)
    # calculate the pearson correlation pairwise among all motifs according to amino acid+ atom property
    p_fuc = []
    for i, j in itertools.combinations(fuc_lig_cls, 2):
        p_fuc.append(stats.pearsonr(freq_list(i)[0] + freq_list(i)[1], freq_list(j)[0] + freq_list(j)[1])[0])

    # calculate the pearson correlation pairwise from randomly selected motifs in the fg-contained ligands
    p_fuc_rand = []
    fuc_ligcls_rand = []
    fuc_ligands = [i[i.rindex("_")-3:i.rindex("_")] for i in fuc_lig_cls]
    for lig in fuc_ligands:
        k = eval(df_assign.loc[lig]["assign"]).keys()
        cls = np.random.choice(list(k))
        fuc_ligcls_rand.append(func + "_" + lig + "_" + cls)
    for m, n in itertools.combinations(fuc_ligcls_rand, 2):
        p_fuc_rand.append(stats.pearsonr(freq_list(m)[0] + freq_list(m)[1], freq_list(n)[0] + freq_list(n)[1])[0])

    # create a boxplot displaying distributions of fg-bind motif pairs pearson correlation
    # vs random motif pars pearson correlation
    fig3, axes3 = plt.subplots()
    axes3.boxplot(np.array(p_fuc),  boxprops=dict(facecolor="#659435"), medianprops=dict(lw=2), showfliers=False,
                  patch_artist=True, positions=[0.25])
    axes3.boxplot(np.array(p_fuc_rand), boxprops=dict(facecolor="#954535"), medianprops=dict(lw=2), showfliers=False,
                  patch_artist=True, positions=[0.75])
    axes3.set_xticks([0.25, 0.75])
    axes3.set_xticklabels(["%s pair" % func, "Random pair"], size="large")
    axes3.set_ylabel("Pearson correlation", size="large")
    axes3.spines['top'].set_visible(False)
    axes3.spines['right'].set_visible(False)
    axes3.hlines(y=1.05, xmin=0.25, xmax=0.75, linewidth=1)
    # Calculate the T-test for the two distributions
    p_value = stats.ttest_ind(np.array(p_fuc), np.array(p_fuc_rand))[1]
    axes3.text(0.33, 1.07, s="%.3e" % p_value, fontsize=14)

    # save the 3 plots with p-value
    fig1.savefig(os.path.join(atpFolder, func, "AA.png"), dpi=300, bbox_inches="tight")
    fig2.savefig(os.path.join(atpFolder, func, "ATOM.png"), dpi=300, bbox_inches="tight")
    fig3.savefig(os.path.join(atpFolder, func, "PCC.png"), dpi=300, bbox_inches="tight")


def PPC_overall():
    '''this function is to quantify Reusability for overall function groups by pearson correlation,
       and generate a box plot comparing the pearson values distribution between motifs binding same FG and motifs
       randomly selected from specified ligands
       o output: return a boxplot labeled with p-value of T-test for the two distributions
    '''
    # choose the fg-binding motifs from different function groups

    # calculate the pearson-values pairwise from the FGs which count == 1
    p_1 = []
    func_1 = [i for i in df_fg.index if df_fg.loc[i]["count"] == 1]
    comb_len1 = []
    for k in func_1:
        comb_len1 += func_cls(k)
    for ki, kj in itertools.combinations(comb_len1, 2):
        p_1.append(stats.pearsonr(freq_list(ki)[0]+freq_list(ki)[1], freq_list(kj)[0]+freq_list(kj)[1])[0])
    # calculate the pearson-values pairwise from the FGs which count > 1
    p_m2 = []
    func_m2 = [i for i in df_fg.index if df_fg.loc[i]["count"] != 1]
    for i, j in itertools.combinations(func_m2, 2):
        l1 = func_cls(i)
        l2 = func_cls(j)
        for m in l1:
            for n in l2:
                p_m2.append(stats.pearsonr(freq_list(m)[0]+freq_list(m)[1], freq_list(n)[0]+freq_list(n)[1])[0])
    # calculate the pearson-values pairwise from the count == 1 and count >1
    p_1_m2 = []
    comb_len2 = []
    for i in func_m2:
        comb_len2 += func_cls(i)
    for fc1 in comb_len1:
        for fc2 in comb_len2:
            p_1_m2.append(stats.pearsonr(freq_list(fc1)[0]+freq_list(fc1)[1], freq_list(fc2)[0]+freq_list(fc2)[1])[0])
    p_from_diff = p_1 + p_m2 + p_1_m2


    # choose the fg-binding motifs from same function groups
    p_from_same = []
    for i in func_m2:
        for j, k in itertools.combinations(func_cls(i), 2):
            p_from_same.append(stats.pearsonr(freq_list(j)[0]+freq_list(j)[1], freq_list(k)[0]+freq_list(k)[1])[0])

    # obtain all pairs
    p_all = p_from_same + p_from_diff

    sames = np.array(p_from_same)  # 3228
    all = np.array(p_all)  #141778
    np.random.shuffle(all)
    # randomly select pairs from all pairs with equal size of same pairs
    all_random = np.random.choice(all, size=3228)

    # create a boxplot displaying distributions of same fg-bind motif pairs pearson correlation
    # vs random motif pars pearson correlation
    fig, axes = plt.subplots()
    axes.boxplot(sames, patch_artist=True, boxprops=dict(facecolor="#659435"), medianprops=dict(lw=2),
                 showfliers=False, positions=[0.1], widths=0.45)
    axes.boxplot(all_random, patch_artist=True, boxprops=dict(facecolor="#954535"), medianprops=dict(lw=2),
                 showfliers=False, positions=[0.6], widths=0.45)
    axes.set_xticks([0.1, 0.6])
    axes.set_xticklabels(["Same FG pair", "Random FG pair"], size=18)
    axes.set_ylabel("Pearson correlation", size=18)
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.hlines(y=1.02, xmin=0.1, xmax=0.6, linewidth=1, colors="k")
    # Calculate the T-test for the two distributions
    p_value = stats.ttest_ind(sames, all_random)[1]
    axes.text(0.1, 1.04, s="p=%.2e" % p_value, fontsize=20)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    # save the boxplot with p-value
    fig.savefig(os.path.join(aftDir, "aftme", "PPC_overall.png"), dpi=300, bbox_inches="tight")


def func_cal(func):
    '''this function is to calculate unique FG reusability and the p-value
       o input: (1) string, denotes the function group
       o output: return 2 values
                 the first one is average pearson correlations
                 the seond one is average p-values
    '''
    # obtain all motifs binding the function group
    fuc_lig_cls = func_cls(func)
    # calculate the pairwise pearson correlation with corresponding t-test
    p_fuc = []
    for i, j in itertools.combinations(fuc_lig_cls, 2):
        p_fuc.append(stats.pearsonr(freq_list(i)[0] + freq_list(i)[1], freq_list(j)[0] + freq_list(j)[1]))
    # calculate the average pearson correlations, average p-values
    pp = np.array(p_fuc).mean(axis=0)
    return "%.3f" % pp[0], "%.3e" % pp[1]


############################## generation of files and plots using above functions#####################

def featplot(ligInfo, statisFolder, motifAnalyse):
    '''this function is to plot the bar, pie chart for ligands comprising 2 or more function groups, plot the venn diagrams for the
       ligands comprising 2 or 3 function groups, save pictures into the motifAnalyse folder
       o input: (1) ligInfo: an excel file recording all 233 ligands information
                (2) statisFolder: a folder placing txt files recording clustered atoms of function group and binding motif
                    for each of 233 ligands
                (3) motifAnalyse: a folder for saving plots that display motif features
       o output: save the pictures into motifAnalyse folder.
                 for one ligand: return 2 pictures, one is a picture comprising bar plot and pie chart,
                 and the other one is a venn plot.
                 note: if the number of function groups for one ligand is more than 3, there will be no corresponding venn
                       plot.
    '''

    df = pd.read_excel(ligInfo, index_col=0)
    for i in df.index:
        g_num = df.loc[i]["Functional Group Number"]
        if g_num != 1:
            with open(os.path.join(statisFolder, i), "r") as f:
                pro_atom = json.loads(f.readlines()[0])["pro"]
                motifplot(g_num, pro_atom, i, motifAnalyse)



# Reusability for Adenine, Ribose and Triphosphate of ATP
# save pictures into the picForATP folder
aden_name = ["Adenine_ATP_cluster3", "Adenine_ADP_cluster3", "Adenine_AMP_cluster1"]
ribo_name = ["Ribose_ATP_cluster2", "Ribose_ADP_cluster2", "Ribose_AMP_cluster3"]
tripp_name = ["Triphosphate_ATP_cluster1", "Triphosphate_GTP_cluster1", "Triphosphate_UTP_cluster1"]
fg_reuse("Adenine", aden_name, ["ATP", "ADP", "AMP"])
fg_reuse("Ribose", ribo_name, ["ATP", "ADP", "AMP"])
fg_reuse("Triphosphate", tripp_name, ["ATP", "GTP", "UTP"])

# quantify Reusability for overall function groups by pearson correlation
# save picture into the Current working directory, named in "PCC_overall.png"
PPC_overall()

# obtain the average pearson correlations and p-values for all fgs which count >=3
# given the file "FG_3more.xlsx"from fg_reuse folder, save the files into the fg_reuse folder.
fg = pd.read_excel(os.path.join(fg3mFolder, "FG_3more.xlsx"), index_col=0)
for i in fg.index:
    x1, x2 = func_cal(i)
    fg.loc[i, "Conservation"] = json.loads(x1)
    fg.loc[i, "P-value"] = json.loads(x2)
fg.to_excel(os.path.join(fg3mFolder, "FG_3more_result.xlsx"), index_col=0)

