# -*- coding: utf-8 -*-
"""
Created on Fri May  3 21:59:12 2019

@author: souvi
"""
import json

from inspect import getsourcefile
import os
import pandas as pd
path = os.path.dirname(os.path.abspath(getsourcefile(lambda:0)))

data='data'       # data is the data set file name without extension
          # graph is the variable to store the network usin dictionary

output=path+'\\Output'
input=path+'\\Input'
def uniqueProtein(filePath):
    print('in')
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
    #uniqueProteinList.sort()
    print('Unique protein list: ',uniqueProteinList)
    uniqueProteinF.close()
    return uniqueProteinList
L=uniqueProtein('data')
df=pd.DataFrame(L)
df.to_csv(output+'\\uniqueProtein'+'.csv',sep=',',header=None,index=None)
file2=open(output+'\\'+'uniqueProtein'+'.txt','w+')
file2.write(json.dumps(L,indent=4, separators=("\n",":")))
file2.close()