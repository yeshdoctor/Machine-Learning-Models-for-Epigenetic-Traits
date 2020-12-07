import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
import scipy.stats as sp
import numpy as np

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
    
#%% Filtering out data from speadsheet, setting up vars
    
Titles = ["Day 2","Day 4","Day 10 293T","Day 10 A549","Day 4 ChIP"]
Location= ['/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay2CorrelationData/Day2MasterCorrelationRawData.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay2CorrelatedWithRNADay4/Day2vsDay4RNAMasterCorrelationRawData.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay10Corr/ChipDay10Corr293T.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay10Corr/ChipDay10A549.csv','/Users/yeshdoctor/Desktop/Duke/LAB/ChipDay4/Day4MasterCorrelationRawData.csv']
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
    DKK1_CHIP=swap(DKK1_CHIP)
    IGF_CHIP=swap(IGF_CHIP)
    MLNR_CHIP=swap(MLNR_CHIP)
    DKK1 = np.array(np.append(DKK1_RNA,DKK1_CHIP).reshape(5,9).T[1:],dtype=float)
    IGF = np.array(np.append(IGF_RNA,IGF_CHIP).reshape(5,9).T[1:],dtype=float)
    MLNR = np.array(np.append(MLNR_RNA,MLNR_CHIP).reshape(5,9).T[1:],dtype=float)
    data = [DKK1,IGF,MLNR]
    labels = ["DKK1","IGFBPL1","MLNR"]
    for i in range(3):
        plt.figure(i)
        f, ax = plt.subplots(figsize=(20,15))
        sns.set(font_scale=2.0)
        cmap = sns.diverging_palette(220,10,sep=20,as_cmap=True)
        # Generate a custom diverging colormap
        y_axis_labels = ["dCas9","KRAB","HDAC8","HDAC8-mut.","Sirt6","Sirt6-mut.","SetDB2","SetDB2-mut."]# labels for x-axis
        x_axis_labels = ["RNA","H3K9me3", "H3K9ac","H4K16ac","H3K27ac"]
        ax.tick_params(axis='x', which='major', pad=10, rotation = 0,labelsize=30)
        ax.tick_params(axis='y', labelsize=25)
        plt.yticks(va="center")
        # Draw the heatmap with the mask and correct aspect ratio
        ax.set_title("Log2 " + Titles[file] + " " + labels[i], fontsize = 40,pad=20)
        sns.heatmap(np.log2(data[i]), cmap=cmap,center=0,annot = True, square = True, linewidths=.5,vmin=-1.7, vmax=1.7,xticklabels=x_axis_labels, yticklabels=y_axis_labels)
        f.savefig("Log2 " + Titles[file] + " " + labels[i]+".tiff",bbox_inches = "tight")