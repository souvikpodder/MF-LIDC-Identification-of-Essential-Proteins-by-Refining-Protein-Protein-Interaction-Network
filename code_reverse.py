# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:45:41 2019

@author: souvi
"""

import json
import statistics
from inspect import getsourcefile
import os
import pandas as pd
import results as rs
path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))

data='data'       # data is the data set file name without extension
          # graph is the variable to store the network usin dictionary

output=path+'\\Output'
input=path+'\\Input'

def Intersection(lst1, lst2): 
    return set(lst1).intersection(lst2)  

def Union(lst1, lst2): 
    final_list = lst1 + lst2 
    t=set(final_list)
    return t

def graphBuild(resultPath,dataPath):
    edgeF=open(input+'\\'+dataPath+'.txt')    
    lvl=[]
    flag=0
    i=0
    graph={}
    while True:
        a=edgeF.readline()
        if a == '':
            break
        b=a.split()
        i+=1
        ##print (b)
        if flag == 0:
            lvl.append(b[1])
            graph[b[0]]=lvl.copy()
            lvl.clear()
            lvl.append(b[0])  
            graph[b[1]]=lvl.copy()
            flag=1
        else:
            if b[0] in graph:
                graph[b[0]].append(b[1])
            else:
                lvl.append(b[1])
                graph[b[0]]=lvl.copy()
                lvl.clear()
            if b[1] in graph:
                graph[b[1]].append(b[0])
            else:
                lvl.append(b[0])
                graph[b[1]]=lvl.copy()
        lvl.clear()
    edgeF.close()
    file3=open(output+'\\details.txt','w+')
    file3.write('At begin :\n'+'node: '+ str(len(graph))+'  Edge: '+ str(i)+'\n')
    file3.close()
    file2=open(output+'\\'+resultPath+'.txt','w+')
    file2.write(json.dumps(graph,indent=4, separators=(",",":")))
    file2.close()
    return graph

def graphBuildList(resultPath,dataPath):
    edgeF=open(output+'\\'+dataPath+'.txt')
    edgeList=json.load(edgeF)    
    lvl=[]
    lvl.clear()
    graph={}
    flag=0
    ##print('Edge in reduced graph: ',edgeList)
    for b in edgeList:
        ###print (b)
        if flag == 0:
            lvl.append(b[1])
            graph[b[0]]=lvl.copy()
            lvl.clear()
            lvl.append(b[0])  
            graph[b[1]]=lvl.copy()
            flag=1
        else:
            if b[0] in graph:
                graph[b[0]].append(b[1])
            else:
                lvl.append(b[1])
                graph[b[0]]=lvl.copy()
                lvl.clear()
            if b[1] in graph:
                graph[b[1]].append(b[0])
            else:
                lvl.append(b[0])
                graph[b[1]]=lvl.copy()
        lvl.clear()
    file2=open(output+'\\'+resultPath+'.txt','w+')
    file2.write(json.dumps(graph,indent=4, separators=(",",":")))
    file2.close()
    return graph
    
#%%
def uniqueProtein(filePath):
    uniqueProteinF=open(input+'\\'+filePath+'.txt')
    s1=set()
    while True:
        a=uniqueProteinF.readline()
        if a == '':
            break
        b=a.split()
        s1.add(b[0])
        s1.add(b[1])
    ##print(s1)
    uniqueProteinList=list(s1)
    uniqueProteinList.sort()
    #print('Unique protein list: ',uniqueProteinList)
    uniqueProteinF.close()
    return uniqueProteinList

def uniqueProteinL(filePath):
    uniqueProteinF=open(output+'\\'+filePath+'.txt')
    uniqueProteinEdgeList=json.load(uniqueProteinF)
    s1=set()
    for b in uniqueProteinEdgeList:
        s1.add(b[0])
        s1.add(b[1])
    ##print(s1)
    uniqueProteinList=list(s1)
    uniqueProteinList.sort()
    #print('Unique protein list: ',uniqueProteinList)
    uniqueProteinF.close()
    return uniqueProteinList
#%%
#   Non Essential Protein
def essentialNonEssential(k):
    reducedProteinF=open(output+'\\LIDC'+str(k)+'.txt')
    reducedProteinList=json.load(reducedProteinF)
    #print(reducedProteinList)
    totalProtein=len(reducedProteinList)
    print(totalProtein)
    essentialProteinNo=int(totalProtein*0.20)
    essentialProteinList=[]
    nonEssentialProteinList=[]
    count=0
    for i in reducedProteinList:
        if count<=essentialProteinNo :
            essentialProteinList.append(i[0])
        else:
            nonEssentialProteinList.append(i[0])
        count+=1
    #print(s1)
    
    essentialProteinF=open(output+'\\essentialProtein'+str(k)+'.txt','w+')
    essentialProteinF.write(json.dumps(essentialProteinList,indent=4 ,separators=(",",":")))
    essentialProteinF.close()
    
    nonEssentialProteinF=open(output+'\\nonEssentialProtein'+str(k)+'.txt','w+')
    nonEssentialProteinF.write(json.dumps(nonEssentialProteinList,indent=4 ,separators=(",",":")))
    nonEssentialProteinF.close()
    
#%%
#%%
def nodeWeight(data):
    graph=graphBuild('mainGraph',data)
    ###print(graph)
    proteinList=uniqueProtein(data)        
    totalProtein=len(proteinList)
    ##print("totalProtein: ",totalProtein)
    ###print(proteinList)
    proteinDegList=[]
    for p in proteinList:
        levList=graph[p]
        deg=0
        for child in levList:
            childList=graph[child]
            deg+=len(childList)
        proteinDegList.append(deg/totalProtein)
    ##print(proteinDegList)
    degF=open(output+'\\nodeDeg'+'.txt','w+')
    degF.write(json.dumps(proteinDegList,indent=4, separators=(",",":")))
    degF.close()    
    alpha = statistics.mean(proteinDegList)
    SD = statistics.stdev(proteinDegList)
    for k in range(1,4):
        print('started k = ',k)
        nodeDelete(alpha,SD,proteinList.copy(),proteinDegList.copy(),graph,k)
        print('finished k = ',k)

def nodeDelete(alpha,SD,proteinList,proteinDegList,graph,k):
    
    THk=alpha+k*SD*(1-(1/(1+(SD**2))))    
    print('For k =',k,'Th = ',THk)
    #print(uniqueProteinList)
    #print(proteinDegList)
    graphTemp=graph.copy()
    for node,val in list(zip(proteinList,proteinDegList)):
        #print('Node :',node,'val: ',val)
        if val < THk:
            graphTemp.pop(node)
            #print('R Node :',node,'val: ',val)
            for parent in graphTemp:
                child=set(graphTemp[parent])
                if node in child:
                    graphTemp[parent].remove(node)
            proteinList.remove(node)
    proteinListSet=set(proteinList)
    edge=open(input+'\\data.txt')
    reducedEdgeList=[]
    while True:
        a=edge.readline()
        if a == '':
            break
        i=a.split()
        if i[0] in proteinListSet and i[1] in proteinListSet:
            reducedEdgeList.append(i)
    edge.close()
    
    df=pd.DataFrame(reducedEdgeList)
    df.to_csv(output+'\\edgesAfterNodeReduction'+str(k)+'.csv',sep=',',header=None,index=None)
    reducedEdgeF=open(output+'\\edgesAfterNodeReduction'+str(k)+'.txt','w+')
    reducedEdgeF.write(json.dumps(reducedEdgeList,indent=4,separators=(",",":")))
    reducedEdgeF.close()
    
    
    
    reducedGraphF=open(output+'\\reducedNodeGraph'+str(k)+'.txt','w+')
    reducedGraphF.write(json.dumps(graphTemp,indent=4,separators=(",",":")))
    reducedGraphF.close()
       
    file3=open(output+'\\details.txt','a+')
    file3.write('For k = '+ str(k) +'\n After node reduce :\n'+'node: '+ str(len(proteinList))+'  Edge: '+ str(len(reducedEdgeList))+'\n')
    file3.close()
    edgeWeight(k,graphTemp)
    
#%%
def edgeWeight(k,graph):
    #print(data)    
    j=0
    #print(data)
    edgeF=open(output+'\\edgesAfterNodeReduction'+str(k)+'.txt')
    edgeList=json.load(edgeF)    
    #print(edgeList)
    edgeWeightVal=[]
    ##print(edgeList[1][0])
    for j in range(len(edgeList)):
        lst1 = graph.get(edgeList[j][0])
        ##print(lst1)
        lst2 = graph.get(edgeList[j][1])
        ##print(lst2)
        w=len(Intersection(lst1,lst2))/len(Union(lst1,lst2))
        ##print(w)
        edgeWeightVal.append(w)
        #lst1=lst2
        j = j + 1
    edgeF.close()
    #print('Edge weight: ',edgeWeightVal)
    #print('Edge: ',edgeList)
    alpha = statistics.mean(edgeWeightVal)
    SD = statistics.stdev(edgeWeightVal)
    ##print(alpha)
    ##print(SD)
    
    TH1=edgeWeightVal
    T1=edgeList
    THk=alpha+k*SD*(1-(1/(1+(SD**2))))
    #print('Threshhold for k = ',k,': ',THk)
    for ele,val in list(zip(edgeList,edgeWeightVal)):
        ##print(ele," :",val)
        if val < THk:
            ##print(ele," :")
           # #print(val)
            T1.remove(ele)
            TH1.remove(val)
    print('Edges after reduction: ',len(T1))
    
    
    
    file2=open(output+'\\reducedEdge'+str(k)+'.txt','w+')
    file2.write(json.dumps(T1,indent=4, separators=(",",":")))
    file2.close()
    
    
    proteinList=uniqueProteinL('reducedEdge'+str(k))
    
    reducedProteinF=open(output+'\\reducedNodeAfterEdgeReduce'+str(k)+'.txt','w+')
    reducedProteinF.write(json.dumps(proteinList,indent=4,separators=(",",":")))
    reducedProteinF.close()
    
    file3=open(output+'\\details.txt','a+')
    file3.write('For k = '+ str(k) +'\n After edge reduce :\n'+'node: '+ str(len(proteinList))+'  Edge: '+ str(len(T1))+'\n')
    file3.close()
    df=pd.DataFrame(T1)
    df.to_csv(output+'\\edgesAfterEdgeReduction'+str(k)+'.csv',sep=',',header=None,index=None)
    LIDCDriver(k,graphBuildList('reducedEdgeGraph'+str(k),'reducedEdge'+str(k)))
    essentialNonEssential(k)
#%%
def LID(k,graph):    
    node=set()
    LIDList={}
    for key in graph:
        v=graph[key]
        s1=set(v)
        edge=0
        for i in v:
            v1=graph[i]
            edge+=len(set(v1).intersection(s1))
            node.update(set(v1).intersection(s1))
            #j loop end
        #i loop end
        ##print(node)
        edge=edge/2
        if(len(node)>0):
            w=edge/len(node)
        else:
            w=0
        node.clear()
        LIDList[key]=w
    LIDListSorted=sorted(LIDList.items(), key=lambda t: t[1],reverse=True)
    LID=open(output+'\\LID'+str(k)+'.txt','w+')
    LID.write(json.dumps(LIDList,indent=4 ,separators=(",",":")))
    LID.close()
    #print(LIDList)
    return(LIDListSorted)
#%%
def IDC(k,graph):
    path2 = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))+'\\Input\\pComplex.xlsx'
    reducedProteinF=open(output+'\\reducedNodeAfterEdgeReduce'+str(k)+'.txt')
    uniqueProteinSet=json.load(reducedProteinF)
    xl = pd.ExcelFile(path2)
#    #print(xl.sheet_names)
    s=[] 
    IDCList={}
    #temp
    s1=[]
    # Load a sheet into a DataFrame by name: df1
    df1 = xl.parse('proteinComplex',header=None)
#    #print(df1)
    #print(type(uniqueProteinSet))
    for row in df1.itertuples(index=False):
            s1.clear()
            for x in row:        
                if x==x and x in uniqueProteinSet:
                    s1.append(x)
            s.append(s1)
            for node in s1:
    #            #print(node)
                s2=set(graph[node])
                IDC=len(s2.intersection(set(s1)))
                if node in IDCList:
                    IDCList[node]+=IDC
                else:
                    IDCList[node]=IDC
    #print(IDCList)
    IDC=open(output+'\\IDC'+str(k)+'.txt','w+')
    IDC.write(json.dumps(IDCList,indent=4 ,separators=(",",":")))
    IDC.close()
    return IDCList

#%%
#LIDC
    
def LIDC(LIDList,IDCList,k):
    rank=1
    N=len(LIDList)
    LIDCList={}
    for tuple in LIDList:
#        #print(type(tuple))
        LID=tuple[1]
        #print(tuple[0])
        if(tuple[0] in IDCList):
            IDC=IDCList[tuple[0]]
        else:
            IDC=0
        #print('LID: ',LID,'IDC: ',IDC,'Rank: ',rank,'N: ',rank/N)
        LIDC=(LID*(1-(rank/N)))+(IDC*(rank/N))
        #print('LIDC: ',LIDC)
        LIDCList[tuple[0]]=LIDC
        rank+=1
    #print(LIDCList)
    LIDCListSorted=sorted(LIDCList.items(), key=lambda t: t[1],reverse=True)
    LIDC=open(output+'\\LIDC'+str(k)+'.txt','w+')
    LIDC.write(json.dumps(LIDCListSorted,indent=4 ,separators=(",",":")))
    LIDC.close()
    print('k= ',k,' finished')

def LIDCDriver(k,graph):
    IDCList=IDC(k,graph)
    LIDList=LID(k,graph)
    LIDC(LIDList,IDCList,k)
 


#%%

def mainDriver():
    nodeWeight(data)
mainDriver()
rs.accuracy()