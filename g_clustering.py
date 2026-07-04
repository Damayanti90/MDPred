# -*- coding: utf-8 -*-




drug=[]
protein=[]
go=[]
pathway=[]
symb="`"
fs=open(r"reduced_merged_embeddings_hgt.txt","r")
for line in fs:
    x=line.split(symb)
    if x[0].startswith("DB"):
       drug.append(x[0])
    elif x[0].startswith("GO"):
       go.append(x[0])
    elif x[0].startswith("R-HSA"):
       pathway.append(x[0])
    else:
       protein.append(x[0])
fs.close()
print(len(drug),len(protein),len(go),len(pathway))

count=0
c=0
import re
path_dic={}
fs=open(r"ReactomePathways.gmt",'r')
for line in fs:
    count+=1
    x=re.split("\t",line)
    if x[1].strip() in pathway:
       c+=1
       genes=[a.strip() for a in x[2:]]
       path_dic[x[1]]=genes
       if c<=10:
          print(count,len(x),len(genes))
          print(line)
fs.close()
print(count,c,len(path_dic.keys()))

def jaccard(path1_genes, path2_genes):
    

    
    set1 = set(path1_genes)
    set2 = set(path2_genes)

    
    inter = set1.intersection(set2)
    inter = len(inter)

    
    union = set1.union(set2)
    union = len(union)

    
    if union == 0:
        
        return 0.0
    else:
        jaccard= inter / union
        return jaccard

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering
import networkx as nx

def hier_clus(go_terms: list, sim_mat: np.ndarray, threshold: float = 0.5):
   

    
    dist_mat = 1 - sim_mat
    np.fill_diagonal(dist_mat, 0)

    if np.any(dist_mat < 0):
        dist_mat[dist_mat < 0] = 0

    
    thresh = 1 - threshold

    
    clustering = AgglomerativeClustering(
        n_clusters=None,
        metric='precomputed',
        linkage='average',
        distance_threshold=thresh
    )


    
    clustering.fit(dist_mat)

    
    #return clustering.n_clusters_
    return clustering.labels_
    #return len(np.unique(cluster_labels))



import pandas as pd

def sim_mat(path_lst):

    sim_matrix = pd.DataFrame(index=path_lst, columns=path_lst, dtype=float)
    for term1 in path_lst:
        for term2 in path_lst:
            score = jaccard(path_dic[term1], path_dic[term2])
            sim_matrix.loc[term1, term2] = score

    
    return sim_matrix

go_coords={}
fp=open(r"2d_pathway_SBERT_embeddings.txt","r")
for line in fp:
    x=line.split("`")
    go_coords[x[0]]=[round(float(x[1]),4),round(float(x[2]),4)]
fp.close()
print(len(go_coords.keys()))
#for i in go_coords.keys():
    #print(go_coords[i])

def inertia(X, labels):
    
    ulabels = np.unique(labels)
    inertia = 0
    for label in ulabels:
        if label == -1:
            continue
        clus_pts = X[labels == label]
        if clus_pts.size > 0:
            centroid = np.mean(clus_pts, axis=0)
            inertia += np.sum((clus_pts - centroid)**2)
    return inertia

import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator

symb="`"
fv=open(r"drugs_threshold_cluster_optimal_pathway.txt","w")
ft=open(r"drugs_threshold_cluster_all_pathway.txt","w")
fs=open(r"predicted_pathway_per_multi_drug.txt","r")
count=0
doubtful=0
for line in fs:
    count+=1   
    print("Drug no",count)
    print(line.strip())
    x=line.split("|")    
    temp=x[1].split(",")
    go_terms=[aa.strip() for aa in temp if aa.strip() in path_dic.keys()]
    print("no of pathways involved",len(go_terms))
    data=sim_mat(go_terms)
    #print(data)
    a=data.to_numpy()
    #print(a)
    #similarity_threshold=0.5
    flag=0
    s=str(x[0]).strip()+"`"
    t=str(x[0]).strip()+"`"
    thresholds = np.linspace(0.1, 0.5, 50)
    inertias = []
    num_clus = []
    for sim_thresh in thresholds:
        #print("checking for threshold=",similarity_threshold)
        labels = hier_clus(go_terms, a, threshold=sim_thresh)
        n_clus = len(np.unique(labels))
        num_clus.append(n_clus)
        #go_terms=np.array(go_terms)
        temp=[]
        for g in go_terms:
            temp.append(go_coords[g])
        temp=np.array(temp)
        inert = inertia(temp, labels)
        s+="threshold="+str(sim_thresh)+",WCSS="+str(inert)+", #clusters="+str(num_clus)+symb
        inertias.append(inert)
    s=s.strip()
    ft.write(s)
    ft.write("\n")
    print(s)
    kn = KneeLocator(thresholds, inertias, curve='convex', direction='decreasing', S=1.0)
    opt_thresh= kn.elbow    
    if opt_thresh is not None:
       idx = list(thresholds).index(opt_thresh)
       opt_clus = num_clus[idx]

       print("optimal threshold: ",opt_thresh)
       print("optimal clusters: ",opt_clus)
       t+="threshold="+str(opt_thresh)+", WCSS="+str(inertias[idx])+",#clusters="+str(opt_clus)+symb
    else:
       print("elbow not found.")
       t+="threshold=NA,WCSS=NA,#clusters=NA`"
    t=t.strip()
    print(t)
    fv.write(t)
    fv.write("\n")
    
fs.close()
fv.close()
ft.close()


from pygosemsim import similarity
from pygosemsim import graph


G = graph.from_obo(r"go-basic.obo")

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering
import networkx as nx


def sim_mat_GO(go_terms):
    sim_matrix = pd.DataFrame(index=go_terms, columns=go_terms, dtype=float)

    
    for term1 in go_terms:
        for term2 in go_terms:
            score = similarity.wang(G, term1, term2)
            sim_matrix.loc[term1, term2] = score

    
    return sim_matrix

go_coords={}
fp=open(r"2d_GO_SBERT_embeddings.txt","r")
for line in fp:
    x=line.split("`")
    go_coords[x[0]]=[round(float(x[1]),4),round(float(x[2]),4)]
fp.close()
print(len(go_coords.keys()))
#for i in go_coords.keys():
    #print(go_coords[i])





import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator
symb="`"
fv=open(r"drugs_threshold_cluster_optimal_mf.txt","w")
ft=open(r"drugs_threshold_cluster_all_mf.txt","w")
fs=open(r"predicted_mf_per_multi_drug.txt","r")
count=0
doubtful=0
for line in fs:
    count+=1
    print("Drug no",count)
    print(line.strip())
    x=line.split("|")    
    temp=x[1].split(",")
    go_terms=[aa.strip() for aa in temp]
    print("no of GO terms involved",len(go_terms))

    data=sim_mat_GO(go_terms)
    #print(data)
    a=data.to_numpy()
    #print(a)
    #similarity_threshold=0.5
    flag=0
    s=str(x[0]).strip()+"`"
    t=str(x[0]).strip()+"`"
    thresholds = np.linspace(0.1, 0.5, 50)
    inertias = []
    num_clus = []
    for sim_thresh in thresholds:
        #print("checking for threshold=",similarity_threshold)
        labels = hier_clus(go_terms, a, threshold=sim_thresh)
        n_clus = len(np.unique(labels))
        num_clus.append(n_clus)
        #go_terms=np.array(go_terms)
        temp=[]
        for g in go_terms:
            temp.append(go_coords[g])
        temp=np.array(temp)
        inert = inertia(temp, labels)
        s+="threshold="+str(sim_thresh)+",WCSS="+str(inert)+", #clusters="+str(num_clus)+symb
        inertias.append(inert)
    s=s.strip()
    ft.write(s)
    ft.write("\n")
    print(s)
    kn = KneeLocator(thresholds, inertias, curve='convex', direction='decreasing', S=1.0)
    opt_thresh= kn.elbow

    
    if opt_thresh is not None:
       idx = list(thresholds).index(opt_thresh)
       opt_clus = num_clus[idx]

       print("optimal threshold: ",opt_thresh)
       print("optimal clusters: ",opt_clus)
       t+="threshold="+str(opt_thresh)+", WCSS="+str(inertias[idx])+",#clusters="+str(opt_clus)+symb
    else:
       print("elbow not found.")
       t+="threshold=NA,WCSS=NA,#clusters=NA`"
    t=t.strip()
    print(t)
    fv.write(t)
    fv.write("\n")
    
    
fs.close()
fv.close()
ft.close()

import numpy as np
import matplotlib.pyplot as plt
from kneed import KneeLocator
symb="`"
fv=open(r"drugs_threshold_cluster_optimal_bp.txt","w")
ft=open(r"drugs_threshold_cluster_all_bp.txt","w")
fs=open(r"predicted_bp_per_multi_drug.txt","r")
count=0
doubtful=0
for line in fs:
    count+=1    
    print("Drug no",count)
    print(line.strip())
    x=line.split("|")
    
    temp=x[1].split(",")
    go_terms=[aa.strip() for aa in temp]
    print("no of GO terms involved",len(go_terms))
    print(go_terms)
    data=sim_mat_GO(go_terms)
    #print(data)
    a=data.to_numpy()
    #print(a)
    #similarity_threshold=0.5
    flag=0
    s=str(x[0]).strip()+"`"
    t=str(x[0]).strip()+"`"
    thresholds = np.linspace(0.1, 0.5, 50)
    inertias = []
    num_clus = []
    for sim_thresh in thresholds:
        #print("checking for threshold=",sim_thresh)
        labels = hier_clus(go_terms, a, threshold=sim_thresh)
        n_clus = len(np.unique(labels))
        num_clus.append(n_clus)
        #go_terms=np.array(go_terms)
        temp=[]
        for g in go_terms:
            temp.append(go_coords[g])
        temp=np.array(temp)
        inert = inertia(temp, labels)
        s+="threshold="+str(sim_thresh)+",WCSS="+str(inert)+", #clusters="+str(num_clus)+symb
        inertias.append(inert)
    s=s.strip()
    ft.write(s)
    ft.write("\n")
    print(s)
    kn = KneeLocator(thresholds, inertias, curve='convex', direction='decreasing', S=1.0)
    opt_thresh= kn.elbow   
    if opt_thresh is not None:
       idx = list(thresholds).index(opt_thresh)
       opt_clus = num_clus[idx]

       print("optimal threshold: ",opt_thresh)
       print("optimal clusters: ",opt_clus)
       t+="threshold="+str(opt_thresh)+", WCSS="+str(inertias[idx])+",#clusters="+str(opt_clus)+symb
    else:
       print("elbow not found.")
       t+="threshold=NA,WCSS=NA,#clusters=NA`"
    t=t.strip()
    print(t)
    fv.write(t)
    fv.write("\n")
    
    
fs.close()
fv.close()
ft.close()




drug_clus={}

count=0
symb="`"
fv=open(r"drugs_threshold_cluster_optimal_pathway.txt","r")
min_wcss=99999
max_wcss=0
min_thres=1
max_thres=0
min_clus=99999
max_clus=0
for line in fv:
    count+=1
    print(line.strip())
    x=line.split(symb)
    y=x[1].split(",")
    if y[0].strip().endswith("NA"):
       continue
    clus=float(y[2].split("=")[1])
    wcss=float(y[1].split("=")[1])
    thres=float(y[0].split("=")[1])
    if wcss>max_wcss:
       max_wcss=wcss
    if wcss<min_wcss:
       min_wcss=wcss
    if thres>max_thres:
       max_thres=thres
    if thres<min_thres:
       min_thres=thres
    if clus>max_clus:
       max_clus=clus
    if clus<min_clus:
       min_clus=clus
    #if count<=10:

    print(y[2])
    print(clus)
    drug_clus[x[0]]=[clus]
fv.close()
print(count,len(drug_clus.keys()))
print("threshold",min_thres,max_thres)
print("wcss",min_wcss,max_wcss)
print("#clusters",min_clus,max_clus)

count=0
symb="`"
fv=open(r"drugs_threshold_cluster_optimal_mf.txt","r")
min_wcss=99999
max_wcss=0
min_thres=1
max_thres=0
min_clus=99999
max_clus=0
for line in fv:
    count+=1
    print(line.strip())
    x=line.split(symb)
    y=x[1].split(",")
    if y[0].strip().endswith("NA"):
       continue
    clus=float(y[2].split("=")[1])
    wcss=float(y[1].split("=")[1])
    thres=float(y[0].split("=")[1])
    if wcss>max_wcss:
       max_wcss=wcss
    if wcss<min_wcss:
       min_wcss=wcss
    if thres>max_thres:
       max_thres=thres
    if thres<min_thres:
       min_thres=thres
    if clus>max_clus:
       max_clus=clus
    if clus<min_clus:
       min_clus=clus
    #if count<=10:

    print(y[2])
    print(clus)
    #drug_clus[x[0]]=[clus]

    temp=drug_clus[x[0]]
    temp.append(clus)
    drug_clus[x[0]]=temp
fv.close()
print(count,len(drug_clus.keys()))
print("threshold",min_thres,max_thres)
print("wcss",min_wcss,max_wcss)
print("#clusters",min_clus,max_clus)

count=0
symb="`"
fv=open(r"drugs_threshold_cluster_optimal_bp.txt","r")
min_wcss=99999
max_wcss=0
min_thres=1
max_thres=0
min_clus=99999
max_clus=0
for line in fv:
    count+=1
    print(line.strip())
    x=line.split(symb)
    y=x[1].split(",")
    if y[0].strip().endswith("NA"):
       continue
    clus=float(y[2].split("=")[1])
    wcss=float(y[1].split("=")[1])
    thres=float(y[0].split("=")[1])
    if wcss>max_wcss:
       max_wcss=wcss
    if wcss<min_wcss:
       min_wcss=wcss
    if thres>max_thres:
       max_thres=thres
    if thres<min_thres:
       min_thres=thres
    if clus>max_clus:
       max_clus=clus
    if clus<min_clus:
       min_clus=clus
    #if count<=10:

    print(y[2])
    print(clus)
    #drug_clus[x[0]]=[clus]

    temp=drug_clus[x[0]]
    temp.append(clus)
    drug_clus[x[0]]=temp
fv.close()
print(count,len(drug_clus.keys()))
print("threshold",min_thres,max_thres)
print("wcss",min_wcss,max_wcss)
print("#clusters",min_clus,max_clus)

import re
id_name={}
fs=open(r"structure_links.tsv","r")
for line in fs:
    x=re.split('\t',line)
    id_name[x[0]]=x[1]
fs.close()
print(len(id_name.keys()))

for i in drug_clus.keys():
    print(i,drug_clus[i])

mult_drug=0
symb="`"
fs=open(r"predicted_multifunctional_drugs.txt","w")
for d in drug_clus.keys():
    ttt=len(drug_clus[d])
    if ttt!=3:
       continue
    [a,b,c]=drug_clus[d]
    #print(d,a,b,c)
    if a>1 and b>1 and c>1:
       mult_drug+=1
       s=str(d)+symb+str(id_name[d])+symb+str(c)+symb+str(b)+symb+str(a)+symb
       s=s.strip()
       #print(d,id_name[d],a,b,c)
       print(s)
       fs.write(s)
       fs.write("\n")
fs.close()
print(mult_drug)
