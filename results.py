# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 10:20:19 2019

@author: Anik
"""

import json
from inspect import getsourcefile
import os
import pandas as pd

path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))
essentialProteinData=[]
nonEssentialProteinData=[]


def accuracy():
    
    path2 = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))+'\\Input\\essential_non_essential_data.xls'
    xl = pd.ExcelFile(path2)
    # Load a sheet into a DataFrame by name: df1
    df1 = xl.parse('essential_proteins',header=None) 
    df2 = xl.parse('non_essential_proteins',header=None)
    #print(df1)
    for i in df1.itertuples(index=False):
        essentialProteinData.append(i[0])
    for j in df2.itertuples(index=False):
        nonEssentialProteinData.append(j[0])
    #print(essentialProteinData)
    #print(nonEssentialProteinData)
    Y1=essentialProtein(1)
    nEP1=nonEssentialProtein(1)
    Y2=essentialProtein(2)
    nEP2=nonEssentialProtein(2)
    Y3=essentialProtein(3)
    nEP3=nonEssentialProtein(3)
    preDict=open(path+'\\Output\\details.txt','a+')
    first = second = third = 0
    preDict.write('\n\nFor Top 100 '+'\n\n')
    for i in range(100):
        if Y1[i] in essentialProteinData:
            first = first + 1
            #print(Y1[i])
        if Y2[i] in essentialProteinData:
            second = second + 1
        if Y3[i] in essentialProteinData:
            third = third + 1
    print(first ,'-->',second ,'-->',third)
    
    preDict.write("For TH1 "+str(first))
    preDict.write("\n")
    preDict.write("For TH2 "+str(second))
    preDict.write("\n")
    preDict.write("For TH3 "+str(third))
    preDict.write("\n")
    
    preDict.write('\n\nFor Top 200 '+'\n\n')
    first = len(set(Y1).intersection(essentialProteinData))
    second = len(set(Y2).intersection(essentialProteinData))
    third = len(set(Y3).intersection(essentialProteinData))
    
    preDict.write("For TH1 "+str(first))
    preDict.write("\n")
    preDict.write("For TH2 "+str(second))
    preDict.write("\n")
    preDict.write("For TH3 "+str(third))
    preDict.write("\n")
    
    preDict.close()
    
   
    sixstatcal(Y1,nEP1,1)
    sixstatcal(Y2,nEP2,2)
    sixstatcal(Y3,nEP3,3)

def essentialProtein(k):
    path2 = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))+'\\Output\\essentialProtein'+str(k)+'.txt'
    eProtein=open(path2)
    eProteinSet=json.load(eProtein)
    eProtein.close()
    preDict=open(path+'\\Output\\details.txt','a+')
    preDict.write('\n************************\nEssential Protein Count For k = '+str(k)+' is:  '+str(len(eProteinSet))+'\n************************\n')
    preDict.close()
    return eProteinSet

def nonEssentialProtein(k):
    path3 = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))+'\\Output\\nonEssentialProtein'+str(k)+'.txt'
    eProtein=open(path3)
    eProteinSet=json.load(eProtein)
    preDict=open(path+'\\Output\\details.txt','a+')
    preDict.write('\n************************\nNon Essential Protein Count For k = '+str(k)+' is:  '+str(len(eProteinSet))+'\n************************\n')
    preDict.close()
    print('Non Essential Protein',len(eProteinSet))
    eProtein.close()
    return eProteinSet
 
def sixstatcal(EP,nEP,k):
    #print(nEP)
    TP = FP = TN = FN = 0
    for i in range(len(EP)):
        if EP[i] in essentialProteinData:
            TP = TP + 1
        elif EP[i] in nonEssentialProteinData:
            FP = FP + 1
    for i in range(len(nEP)):
        if nEP[i] in nonEssentialProteinData:
            TN = TN + 1
        elif nEP[i] in essentialProteinData:
            FN = FN + 1
    Sensitivity = TP / (TP + FN)
    Specificity = TN / (FP + TN)
    NPV = TN / (FN + TN)
    PPV = TP / (TP + FP)
    ACC = (TP + TN)/(len(EP) + len(nEP))
    F = (2 * Sensitivity * PPV) / (Sensitivity + PPV)
    preDict=open(path+'\\Output\\details.txt','a+')
    preDict.write('\n\nFor K = '+str(k))
    preDict.write("\n\n")
    preDict.write("Sensitivity = "+str(Sensitivity))
    preDict.write("\n")
    preDict.write("Specificity = "+str(Specificity))
    preDict.write("\n")
    preDict.write("NPV = "+str(NPV))
    preDict.write("\n")
    preDict.write("PPV = "+str(PPV))
    preDict.write("\n")
    preDict.write("ACC = "+str(ACC))
    preDict.write("\n")
    preDict.write("F - measure = "+str(F))
    preDict.write("\n")
    preDict.write("TP = "+str(TP))
    preDict.write("\n")
    preDict.write("FP = "+str(FP))
    preDict.write("\n")
    preDict.write("TN = "+str(TN))
    preDict.write("\n")
    preDict.write("FN = "+str(FN))
    preDict.write("\n")
    preDict.close()