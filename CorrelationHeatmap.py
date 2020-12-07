import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
import scipy.stats as sp
import numpy as np
#%% Filtering out data from speadsheet, setting up vars
CorrTitles = ["Day 2 Correlation","Day 4 Correlation","Day 10 293T Correlation","Day 10 A549 Correlation", "Day 4 Chip Correlation"]
PTitls = ["Day 2 P-value","Day 4 P-value","Day 10 293T P-value","Day 10 A549 P-value","Day 4 ChIP P-value"]
Location = ['/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay2CorrelationData/Day2MasterCorrelationRawData.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay2CorrelatedWithRNADay4/Day2vsDay4RNAMasterCorrelationRawData.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay10Corr/ChipDay10Corr293T.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay10Corr/ChipDay10A549.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay4/Day4MasterCorrelationRawData.csv']
for file in range(5):  
    df = pd.read_csv(Location[file])
    df = df.iloc[:, 0:11]
    DKK1_RNA = ((df.loc[df["Gene"]=="DKK1_RNA"]))
    DKK1_RNA = DKK1_RNA.to_numpy()[0,2:]
    IGF_RNA = ((df.loc[df["Gene"]=="IGFBPL1_RNA"]))
    IGF_RNA = IGF_RNA.to_numpy()[0,2:]
    MLNR_RNA = ((df.loc[df["Gene"]=="MLNR_RNA"]))
    MLNR_RNA = MLNR_RNA.to_numpy()[0,2:]
    DKK1_CHIP = ((df.loc[df["Gene"]=="DKK1_CHIP"]))
    DKK1_CHIP = DKK1_CHIP.to_numpy()[:,2:]
    IGF_CHIP = ((df.loc[df["Gene"]=="IGFBPL1_CHIP"]))
    IGF_CHIP = IGF_CHIP.to_numpy()[:,2:]
    MLNR_CHIP = ((df.loc[df["Gene"]=="MLNR_CHIP"]))
    MLNR_CHIP = MLNR_CHIP.to_numpy()[:,2:]
    
    def swap(arr):
        t = [[],[],[],[]]
        temp0 = arr[0]
        temp1 = arr[1]
        temp2 = arr[2]
        temp3 = arr[3]
        t[0] = temp2
        t[1] = temp3
        t[2] = temp0
        t[3] = temp1
        return t
    
    DKK1_CHIP=swap(DKK1_CHIP)
    IGF_CHIP=swap(IGF_CHIP)
    MLNR_CHIP=swap(MLNR_CHIP)
    
    DKK1_CORR = np.zeros(4)
    DKK1_PVAL = np.zeros(4)
    IGF_CORR = np.zeros(4)
    IGF_PVAL = np.zeros(4)
    MLNR_CORR = np.zeros(4)
    MLNR_PVAL = np.zeros(4)
    #%% Running spearman correlations
    for i in range (0,4):
        rho,pval=sp.mstats.spearmanr(np.array(DKK1_RNA,dtype=float),np.array(DKK1_CHIP[i],dtype=float),nan_policy="omit")
        DKK1_CORR[i] = rho
        DKK1_PVAL[i] = pval
    for i in range (0,4):
        rho,pval=sp.mstats.spearmanr(np.array(IGF_RNA,dtype=float),np.array(IGF_CHIP[i],dtype=float),nan_policy="omit")
        IGF_CORR[i] = rho
        IGF_PVAL[i] = pval
    for i in range (0,4):
        rho,pval=sp.mstats.spearmanr(np.array(MLNR_RNA,dtype=float),np.array(MLNR_CHIP[i],dtype=float),nan_policy="omit")
        MLNR_CORR[i] = rho
        MLNR_PVAL[i] = pval
    
    corr = np.vstack((DKK1_CORR,IGF_CORR,MLNR_CORR))
    pval = np.vstack((DKK1_PVAL,IGF_PVAL,MLNR_PVAL))
    corr = corr.T
    pval = pval.T
    #%% Graphing Section. Uncomment lines 80 and 81 and comment lines 66 and 67 for PVAL heatmap
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(10,10))
    sns.set(font_scale=2.0)
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220,10,sep=20,as_cmap=True)
    x_axis_labels = ["DKK1","IGFBPL1","MLNR"]# labels for x-axis
    y_axis_labels = ["H3K9me3", "H3K9ac","H4K16ac","H3K27ac"]
    ax.tick_params(axis='x', which='major', pad=10, rotation = 0,labelsize=30)
    ax.tick_params(axis='y', labelsize=25)
    plt.yticks(va="center")
    # Draw the heatmap with the mask and correct aspect ratio
    #ax.set_title(CorrTitles[file], fontsize = 30,pad=20)
    #sns.heatmap(corr, cmap=cmap, center=0,annot = True, square = True, linewidths=.5,vmin=-1, vmax=1,xticklabels=x_axis_labels, yticklabels=y_axis_labels)
    #f.savefig(CorrTitles[file]+".tiff")
    ax.set_title(PTitls[file], fontsize = 30,pad=20)
    sns.heatmap(pval, cmap=cmap, center=0,annot = True, square = True, linewidths=.5,vmin=-0, vmax=1,xticklabels=x_axis_labels, yticklabels=y_axis_labels)
    f.savefig(PTitls[file]+".tiff")
