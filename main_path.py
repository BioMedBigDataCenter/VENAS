import pandas as pd
import numpy as np
import networkx as nx
import cdlib
import matplotlib.pyplot as plt
from networkx.algorithms import community
import itertools
from networkx import edge_betweenness_centrality as betweenness
from operator import itemgetter
from cdlib import algorithms, viz, NodeClustering

def impcsv(csvfile):
    t  = pd.read_csv(csvfile,sep=",")
    return t

def generateGraph(pdt):
    G=nx.from_pandas_edgelist(pdt,source="Source",target="Target")
    return G

def getCluster(G):
    coms = cdlib.algorithms.louvain(G,resolution=2.0,randomize=True)
    return coms

'''
def getClusterNewman(G):
    coms = nx.community.girvan_newman(G)
    return coms
'''


def process(G,coms,mode="big"):
    import math
    n_com = len(coms.communities)
    #find key node in each cluster
    deg =dict(G.degree())
    keyNode={}
    keyNodeList=[]
    comLengthCount={}
    for com in coms.communities:
        score=0
        keyNum=0
        for node in com:
            if score < deg[node]:
                score = deg[node]
                keyNum = node
        keyNode[keyNum]=score
        comLengthCount[keyNum]=len(com)
        keyNodeList.append(keyNum)
    #if key node link to each other?
    #filter keyNode use threshold or cluster count

    G_smaller=G.copy()
    removeNode=[]
    if mode=="big":
        for node in G.node():
            if deg[node] ==1:
                removeNode.append(node)
    for node in removeNode:
        G_smaller.remove_node(node)

    filterNode={}
    if mode=="small":
        threshold = math.log(len(G.node()))*2
    else:
        threshold = len(G_smaller.node())/float(len(coms.communities))
    '''
    for node in keyNode.keys():
        if keyNode[node] > threshold:
            filterNode[node]=keyNode[node]
    '''
    for node in keyNode.keys():
        if comLengthCount[node]>threshold:
            filterNode[node]=comLengthCount[node]
    #for node in keyNode.keys():
    #    if deg[node]>threshold:
    #        filterNode[node]=comLengthCount[node]

    print(len(filterNode.keys()))
    #check dijkstra_path
    mainPath={}
    for nodestart in filterNode.keys():
        for nodeend in filterNode.keys():
            if nodestart != nodeend:
                thiskey = (nodestart,nodeend)
                reversekey = (nodeend,nodestart)
                if thiskey not in mainPath.keys() and reversekey not in mainPath.keys():
                    mainPath[thiskey]=None

    return_Path={}
    for thiskey in mainPath.keys():
        thispath=nx.dijkstra_path(G_smaller,thiskey[0],thiskey[1])
        if len(thispath)>2:
            flag = 0
            for node in thispath:
                if node != thiskey[0] and node != thiskey[1]:
                    if node  in  filterNode.keys():
                        flag +=1
            if flag==0:
                return_Path[thiskey]=thispath
        else:
            return_Path[thiskey]=thispath
    
    f=open("edgeTable.csv","w")
    f.write("Source,Target\n")
    #print(return_Path.keys())
    for eachPath in return_Path.keys():
        l = len(return_Path[eachPath])
        i=0
        j=1
        item = return_Path[eachPath]
        #print(item)
        while (j<l):
                f.write(str(item[i])+","+str(item[j])+"\n")
                i=i+1
                j=j+1
    f.close()
    return return_Path,filterNode,keyNodeList


#def createGephiNodeCSV(coms,rPath,keyNodeList,filterNode,G):
#    i=0
#    prDict={}
#    originDict={}
#    for node in list(G.node):
#        prDict[node]=-3
#    for com in coms.communities:
#        for node in com:
#            prDict[node]=i
#            originDict[node]=i
#        i=i+1
#    transNode=[]
#    for p in rPath.keys():
#        if len(rPath[p]) >2:
#            for node in rPath[p]:
#                if (node!=rPath[p][0] and node!=rPath[p][-1]):
#                    thisCom=prDict[node]
#                    if len(coms.communities[thisCom])>1:
#                        if node not in keyNodeList:
#                            transNode.append(node)
#                            #print(node)
#                            prDict[node]=i
#                            i=i+1
#    for tnode in transNode:
#        thisCom = originDict[tnode]
#        for node in coms.communities[thisCom]:
#            if node != tnode:
#                toKnodePath=nx.dijkstra_path(G,node,keyNodeList[thisCom])
#                if tnode  in toKnodePath:
#                    prDict[node]=prDict[tnode]
#    #process node not in filterNode but in keynode
#    #assign the sam com num direct link to it
#    small_node=[]
#    for node in keyNodeList:
#        if node not in filterNode.keys():
#            if node not in transNode:
#                small_node.append(node)
#                print(node)
#                
#    small_node_color={}          
#    for node in small_node:
#        flag=0
#        for nbr in G[node]:
#            if nbr in filterNode.keys():
#                small_node_color[node]=originDict[nbr]
#                #print(node,nbr)
#                flag=1
#                break
#        if flag==0:
#            small_node_color[node]=-1
#    for node in small_node_color.keys():
#        if small_node_color[node]==-1:
#            min_len=len(G.node)
#            nearest_fnode=0
#            for fNode in filterNode.keys():
#                toFnodePath = nx.dijkstra_path(G,node,fNode)
#                if (len(toFnodePath)<min_len):
#                    min_len=len(toFnodePath)
#                    nearest_fnode=fNode
#            small_node_color[node]=originDict[nearest_fnode]
#    for node in small_node_color.keys():
#        thisCom=coms.communities[originDict[node]]
#        for thisComNode in thisCom:
#            prDict[thisComNode]=small_node_color[node]
#
#    mainPathNode=set()
#    for eachPath in rPath.keys():
#        item = rPath[eachPath]
#        for node in item:
#            mainPathNode.add(node)
#    valueDict={}
#    #valueDict.keys()=list(G.node)
#    for node in list(G.node):
#        if node not in mainPathNode:
#            valueDict[node]=0
#        else:
#            valueDict[node]=1000
#    
#    f=open("nodeTable.csv","w")
#    f.write("Id,Label,ClusterId,Value\n")
#    for node in prDict.keys():
#        f.write(str(node)+","+str(node)+","+str(prDict[node])+","+str(valueDict[node])+"\n")
#    f.close()

def createGephiNodeCSV(coms,rPath,keyNodeList,filterNode,G):
    i=0
    prDict={}
    originDict={}
    for com in coms.communities:
        for node in com:
            prDict[node]=i
            originDict[node]=i
        i=i+1
    #transNode=[]
    transNode=set()
    for p in rPath.keys():
        if len(rPath[p]) >2:
            for node in rPath[p]:
                if (node!=rPath[p][0] and node!=rPath[p][-1]):
                    #if  prDict[node] < len(coms.communities):
                    if node not in transNode and node not in filterNode.keys():
                        #if len(coms.communities[prDict[node]])>1:
                            transNode.add(node)
                            #print(node)
                            prDict[node]=i
                            i=i+1
    for tnode in transNode:
        thisCom = originDict[tnode]
        for node in coms.communities[thisCom]:
            if node != tnode and node not in transNode and prDict[node] < len(coms.communities):
                toKnodePath=nx.dijkstra_path(G,node,keyNodeList[thisCom])
                if tnode  in toKnodePath:
                    minPath=len(nx.dijkstra_path(G,node,tnode))
                    nearestTNode=tnode
                    for ttnode in toKnodePath:
                        if ttnode in transNode and ttnode != tnode:
                            thisPath = len(nx.dijkstra_path(G,node,ttnode))
                            if thisPath < minPath:
                                minPath=thisPath
                                nearestTNode = ttnode
                    prDict[node]=prDict[nearestTNode]
    #process node not in filterNode but in keynode
    #assign the sam com num direct link to it
    small_node=[]
    for node in keyNodeList:
        if node not in filterNode.keys():
            if node not in transNode:
                small_node.append(node)
                print(node)
                
    small_node_color={}          
    for node in small_node:
        flag=0
        for nbr in G[node]:
            if nbr in filterNode.keys():
                small_node_color[node]=originDict[nbr]
                #print(node,nbr)
                flag=1
                break
        if flag==0:
            small_node_color[node]=-1
    for node in small_node_color.keys():
        if small_node_color[node]==-1:
            min_len=len(G.node)
            nearest_fnode=0
            for fNode in filterNode.keys():
                toFnodePath = nx.dijkstra_path(G,node,fNode)
                if (len(toFnodePath)<min_len):
                    min_len=len(toFnodePath)
                    nearest_fnode=fNode
            small_node_color[node]=originDict[nearest_fnode]
    for node in small_node_color.keys():
        thisCom=coms.communities[originDict[node]]
        for thisComNode in thisCom:
            prDict[thisComNode]=small_node_color[node]

    mainPathNode=set()
    for eachPath in rPath.keys():
        item = rPath[eachPath]
        for node in item:
            mainPathNode.add(node)
    valueDict={}
    #valueDict.keys()=list(G.node)
    for node in prDict.keys():
        if node not in mainPathNode:
            valueDict[node]=0
        else:
            valueDict[node]=1000
    
    f=open("nodeTable.csv","w")
    f.write("Id,Label,ClusterId,Value\n")
    for node in prDict.keys():
        f.write(str(node)+","+str(node)+","+str(prDict[node])+","+str(valueDict[node])+"\n")
    f.close()
