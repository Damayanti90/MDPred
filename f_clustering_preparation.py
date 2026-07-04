# -*- coding: utf-8 -*-





threshold=0.7
no_of_trials=50
drug_go_dic={}
drugs_covered=[]
count=0
fs=open(r"result_drug_go_prediction_data_XGBoost.txt","r")
for line in fs:
    count+=1
    x=line.split("`")
    s=x[0]+"`"+x[1]+"`"
    s=s.strip()
    drugs_covered.append(x[0])
    #if s in drug_go:
    drug_go_dic[s]=[float(x[4].strip())]
fs.close()
print(len(drug_go_dic.keys()))
drugs_covered=list(set(drugs_covered))
print(len(drugs_covered))


count=0
fs=open(r"result_drug_go_prediction_data_lightGBM.txt","r")
for line in fs:
    count+=1
    x=line.split("`")
    s=x[0]+"`"+x[1]+"`"
    s=s.strip()
    drugs_covered.append(x[0])
    #if s in drug_go:
    temp=drug_go_dic[s]
    temp.append(float(x[4].strip()))
    drug_go_dic[s]=temp
fs.close()
print(len(drug_go_dic.keys()))
drugs_covered=list(set(drugs_covered))
print(len(drugs_covered))





c=0
d=0
fs=open(r"predicted_drug_go_associations.txt","w")
for line in drug_go_dic.keys():
    [a,b]=drug_go_dic[line]
    if round(a,4)>=threshold and round(b,4)>=threshold:
       c+=1
       s=line.strip()+str(a)+"`"+str(b)+"`"
       s=s.strip()
       if c<=10:
          print(s)
       fs.write(s)
       fs.write("\n")
fs.close()
print(c,c/len(drug_go_dic.keys()))
print(d)

drug_go_dic={}

count=0
fs=open(r"result_drug_pathway_prediction_data_XGBoost.txt","r")
for line in fs:
    count+=1
    x=line.split("`")
    s=x[0]+"`"+x[1]+"`"
    s=s.strip()
    #if s in drug_go:
    drug_go_dic[s]=[float(x[4].strip())]
fs.close()
print(len(drug_go_dic.keys()))


count=0
fs=open(r"result_drug_pathway_prediction_data_lightGBM.txt","r")
for line in fs:
    count+=1
    x=line.split("`")
    s=x[0]+"`"+x[1]+"`"
    s=s.strip()
    #if s in drug_go:
    temp=drug_go_dic[s]
    temp.append(float(x[4].strip()))
    drug_go_dic[s]=temp
fs.close()
print(len(drug_go_dic.keys()))



c=0
d=0
fs=open(r"predicted_drug_pathway_associations.txt","w")
for line in drug_go_dic.keys():
    [a,b]=drug_go_dic[line]
    if round(a,4)>=threshold and round(b,4)>=threshold:
       c+=1
       s=line.strip()+str(a)+"`"+str(b)+"`"
       s=s.strip()
       if c<=10:
          print(s)
       fs.write(s)
       fs.write("\n")
fs.close()
print(c,c/len(drug_go_dic.keys()))
print(d)

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
fs.close()
print(len(drug),len(go),len(pathway))

count=0
drug_go=[]
symb="`"
fs=open(r"predicted_drug_go_associations.txt","r")
for line in fs:
    count+=1
    x=line.split(symb)
    s=x[0]+symb+x[1]+symb
    s=s.strip()
    drug_go.append(s)
fs.close()
print(count)

count=0
par_dic={}
fs=open(r"go_terms_with_parents-all.txt","r")
for line in fs:
    count+=1
    x=line.split("|")
    if x[0].strip() in go:
       y=x[1].split(",")
       if y[-1].strip() in par_dic.keys():
          temp=par_dic[y[-1].strip()]
          temp.append(x[0])
          par_dic[y[-1].strip()]=temp
       else:
          par_dic[y[-1].strip()]=[x[0].strip()]
    if count<=10:
       print(line)
fs.close()
print(count,len(par_dic.keys()))

bp_go=[]
mf_go=[]
cc_go=[]
for parent in par_dic.keys():
    print(parent, len(par_dic[parent]))

    if parent=="GO:0008150:0":

      bp_go=par_dic[parent]
    elif parent=='GO:0003674:0':

      mf_go=par_dic[parent]
    else:

      cc_go=par_dic[parent]

bp_drug={}
mf_drug={}
cc_drug={}
count=0
for line in drug_go:
  count+=1
  if count%100000==0:
     print(count,len(bp_drug.keys()))
  x=line.split(symb)
  if x[1] in bp_go:
    if x[0] in bp_drug.keys():
       temp=bp_drug[x[0]]
       temp.append(x[1])
       bp_drug[x[0]]=temp
    else:
       bp_drug[x[0]]=[x[1]]
  if x[1] in mf_go:
    if x[0] in mf_drug.keys():
       temp=mf_drug[x[0]]
       temp.append(x[1])
       mf_drug[x[0]]=temp
    else:
       mf_drug[x[0]]=[x[1]]
  if x[1] in cc_go:
    if x[0] in cc_drug.keys():
       temp=cc_drug[x[0]]
       temp.append(x[1])
       cc_drug[x[0]]=temp
    else:
       cc_drug[x[0]]=[x[1]]
print(len(bp_drug.keys()),len(mf_drug.keys()),len(cc_drug.keys()))

print(len(bp_drug.keys()))
print(len(mf_drug.keys()))
print(len(cc_drug.keys()))

count=0
symb="`"
drug_pathway=[]
fs=open(r"predicted_drug_pathway_associations.txt","r")
for line in fs:
    count+=1
    x=line.split(symb)
    s=x[0]+symb+x[1]+symb
    s=s.strip()
    drug_pathway.append(s)
fs.close()
print(count,len(drug_pathway))

pathway_drug={}
count=0
for line in drug_pathway:
  count+=1
  if count%100000==0:
     print(count,len(pathway_drug.keys()))
  x=line.split(symb)
  if x[1] in pathway:
    if x[0] in pathway_drug.keys():
       temp=pathway_drug[x[0]]
       temp.append(x[1])
       pathway_drug[x[0]]=temp
    else:
       pathway_drug[x[0]]=[x[1]]
print(len(pathway_drug.keys()))

multi_drug_go=[]
for p in bp_drug.keys():
    if p not in mf_drug.keys():
       continue
    if len(bp_drug[p])>=2 and len(mf_drug[p])>=2 and len(pathway_drug[p])>=2:
       multi_drug_go.append(p)
print(len(multi_drug_go))

involved=[]
fs=open(r"predicted_bp_per_multi_drug.txt","w")
for d in multi_drug_go:
    temp=bp_drug[d]
    s=str(d)+"|"
    for p in temp:
        s+=str(p)+","
    s=s.strip()
    s=s[:-1]
    involved.append(len(temp))
    #print(len(temp),s)
    fs.write(s)
    fs.write("\n")
fs.close()
import numpy as np
involved=np.array(involved)
print(np.max(involved),np.min(involved),np.mean(involved))

involved=[]
fs=open(r"predicted_mf_per_multi_drug.txt","w")
for d in multi_drug_go:
    temp=mf_drug[d]
    s=str(d)+"|"
    for p in temp:
        s+=str(p)+","
    s=s.strip()
    s=s[:-1]
    involved.append(len(temp))
    #print(len(temp),s)
    fs.write(s)
    fs.write("\n")
fs.close()
import numpy as np
involved=np.array(involved)
print(np.max(involved),np.min(involved),np.mean(involved))





multi_drug_pathway=multi_drug_go

involved=[]
fs=open(r"predicted_pathway_per_multi_drug.txt","w")
for d in multi_drug_pathway:
    temp=pathway_drug[d]
    s=str(d)+"|"
    for p in temp:
        s+=str(p)+","
    s=s.strip()
    s=s[:-1]
    involved.append(len(temp))
    #print(len(temp),s)
    fs.write(s)
    fs.write("\n")
fs.close()
import numpy as np
involved=np.array(involved)
print(np.max(involved),np.min(involved),np.mean(involved))






import numpy as np

import pandas as pd
import pickle
from sklearn.preprocessing import Normalizer

df = pd.read_csv(r"pathways-SBERT-embeddings.txt", sep="`", header=None)
#df = pd.read_csv(r"consolidated_fet.txt", sep="`", header=None)
x = df.iloc[:, 1:-1].values


print(df.head().to_string())


from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

x= sc.fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(0.99)
#pca=PCA(n_components=64)


unpickled_df = pca.fit_transform(x)
explained_variance = pca.explained_variance_ratio_


print(explained_variance)

print(unpickled_df)
print(type(unpickled_df))
print(unpickled_df.shape)


cumulative_variance = np.sum(explained_variance)
print("explained variance:", cumulative_variance)


mother=[]
x=df.values.tolist()
y=unpickled_df.tolist()
print(x[0])
print (len(x),len(x[0]))
print(y[0])
print(len(y),len(y[0]))
for i,j in zip(x,y):
    son=[]
    son.append(i[0])

    for k in j:
        son.append(k)
    mother.append(son)
print(mother[0])
print(len(mother),len(mother[0]))
dh=pd.DataFrame(mother)
print(dh)


symb="`"
count=0
#fp=open(r"reduced_consolidated_fet_PCA_1.txt","w")
fp=open(r"reduced_pathways-SBERT-embeddings_PCA.txt","w")



for i in mother:
    count+=1
    st=""
    for j in i:
        st=st+str(j)+symb
    if count<=10:
       print(st.strip())
    fp.write(st.strip())
    fp.write("\n")
    #count+=1
fp.close()
print(count)

count=0
ids=[]
feat_lst=[]
symb="`"
#fs=open(r"reduced_consolidated_fet_PCA_1.txt","r")

fs=open(r"reduced_pathways-SBERT-embeddings_PCA.txt","r")
for line in fs:
    count+=1
    if count<=10:
       print(line)
    x=line.split(symb)
    #if count<=10:
    ids.append(x[0])
    feat_lst.append(x[1:-1])
fs.close()
print(len(ids),len(feat_lst))

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import optuna
from optuna.pruners import MedianPruner, PatientPruner
import numpy as np


input=len(feat_lst[0])
data=np.array(feat_lst).astype(np.float32)


scale = MinMaxScaler()
data = scale.fit_transform(data)
X_train, X_val = train_test_split(data, test_size=0.1, random_state=42)


train_data= TensorDataset(torch.from_numpy(X_train))
val_data= TensorDataset(torch.from_numpy(X_val))


class AE(nn.Module):
    def __init__(self, orig, inter, latent):
        super(AE, self).__init__()


        self.encoder = nn.Sequential(
            nn.Linear(orig, inter),
            nn.ReLU(),
            nn.Linear(inter, inter // 2),
            nn.ReLU(),
            nn.Linear(inter // 2, latent)
        )


        self.decoder = nn.Sequential(
            nn.Linear(latent, inter// 2),
            nn.ReLU(),
            nn.Linear(inter // 2, inter),
            nn.ReLU(),
            nn.Linear(inter, orig),
            nn.Sigmoid()
        )

    def forward(self, x):
        z = self.encoder(x)
        recon= self.decoder(z)
        return recon, z

    def encode(self, x):
        return self.encoder(x)


def train_func(model, train_load, opt):
    model.train()
    tot_loss= 0
    for bx, in train_load:
        opt.zero_grad()
        recon, _ = model(bx)
        loss = F.mse_loss(recon, bx, reduction='mean')
        loss.backward()
        opt.step()
        tot_loss += loss.item()
    return tot_loss/ len(train_load)

def val_func(model, val_load):
    model.eval()
    tot_loss = 0
    with torch.no_grad():
        for bx, in val_load:
            recon, _ = model(bx)
            loss= F.mse_loss(recon, bx, reduction='mean')
            tot_loss += loss.item()
    return tot_loss / len(val_load)

def objective(trial):
    latent= trial.suggest_categorical('latent_dim', [32])
    inter= trial.suggest_int('intermediate_dim', 128, 512)
    lr = trial.suggest_float('learning_rate', 1e-5, 1e-2, log=True)
    batch= trial.suggest_categorical('batch_size', [32, 64, 128, 256])
    epoch = trial.suggest_int('epochs', 30, 100)

    train_load = DataLoader(train_data, batch_size=batch, shuffle=True)
    val_load= DataLoader(val_data, batch_size=batch)

    model = AE(orig=input, inter=inter, latent=latent)
    opt = optim.Adam(model.parameters(), lr=lr)

    for ep in range(epoch):
        trainl_loss = train_func(model, train_load, opt)
        val_loss = val_func(model, val_load)

        trial.report(val_loss, ep)
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    return val_loss

if __name__ == "__main__":
    pr= PatientPruner(wrapped_pruner=MedianPruner(n_startup_trials=5, n_warmup_steps=5), patience=5)
    study = optuna.create_study(direction='minimize', pruner=pr)
    study.optimize(objective, n_trials=no_of_trials, show_progress_bar=True)

    best_par = study.best_params
    print("\nBest hyp", best_par)


    fmodel= AE(input, best_par['intermediate_dim'], best_par['latent_dim'])
    fopt = optim.Adam(fmodel.parameters(), lr=best_par['learning_rate'])

    full_load = DataLoader(TensorDataset(torch.from_numpy(data)),
                             batch_size=best_par['batch_size'], shuffle=True)

    for epoch in range(best_par['epochs']):
        train_func(fmodel, full_load, fopt)


    fmodel.eval()
    with torch.no_grad():
        full_data = torch.from_numpy(data)
        compr = fmodel.encode(full_data).numpy()
        print("\nCompressed Shape:", compr.shape)


count=0
symb="`"
fet_lst=[]
fs=open(r"2d_pathway_SBERT_embeddings.txt","w")
for id,line in zip(ids,compr):
    count+=1
    #if count<=10:
       #print(len(line),line[:5],line[2])
    s=""
    s+=str(id)+symb
    t=""
    for i in range(len(line)):
        s+=str(line[i])+symb
        t+=str(line[i])+symb
    s=s.strip()
    t=t.strip()
    fet_lst.append(t)
    fs.write(s)
    fs.write("\n")
    if count<=10:
       print(s)
fet_lst=list(set(fet_lst))
print(count,len(fet_lst))
fs.close()




#fs=open(r"trial_info_AE_1.txt","w")
fs=open(r"trial_info_2d_pathway_SBERT_embeddings.txt","w")
for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    #print("lr",trial.params['learning_rate'],"batch size", trial.params['batch_size'],"encoding dimension",trial.params["latent_dim"],"hidden_dim",trial.params["intermediate_dim"],"epochs",trial.params['epochs'])
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['learning_rate'])+"`"
    s+="batch_size="+str(trial.params['batch_size'])+"`"
    s+="encoding_dim="+str(trial.params['latent_dim'])+"`"
    s+="hidden_dim="+str(trial.params['intermediate_dim'])+"`"
    s+="epochs="+str(trial.params['epochs'])+"`"
    #s+="beta="+str(trial.params['beta'])+"`"
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()





import numpy as np

import pandas as pd
import pickle
from sklearn.preprocessing import Normalizer

df = pd.read_csv(r"GO-SBERT-embeddings.txt", sep="`", header=None)
#df = pd.read_csv(r"consolidated_fet.txt", sep="`", header=None)
x = df.iloc[:, 1:-1].values


print(df.head().to_string())


from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

x= sc.fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(0.99)
#pca=PCA(n_components=64)


unpickled_df = pca.fit_transform(x)
explained_variance = pca.explained_variance_ratio_


print(explained_variance)

print(unpickled_df)
print(type(unpickled_df))
print(unpickled_df.shape)


cumulative_variance = np.sum(explained_variance)
print("explained variance:", cumulative_variance)


mother=[]
x=df.values.tolist()
y=unpickled_df.tolist()
print(x[0])
print (len(x),len(x[0]))
print(y[0])
print(len(y),len(y[0]))
for i,j in zip(x,y):
    son=[]
    son.append(i[0])

    for k in j:
        son.append(k)
    mother.append(son)
print(mother[0])
print(len(mother),len(mother[0]))
dh=pd.DataFrame(mother)
print(dh)


symb="`"
count=0
#fp=open(r"reduced_consolidated_fet_PCA_1.txt","w")
fp=open(r"reduced_GO-SBERT-embeddings_PCA.txt","w")



for i in mother:
    count+=1
    st=""
    for j in i:
        st=st+str(j)+symb
    if count<=10:
       print(st.strip())
    fp.write(st.strip())
    fp.write("\n")
    #count+=1
fp.close()
print(count)

count=0
ids=[]
feat_lst=[]
symb="`"
#fs=open(r"reduced_consolidated_fet_PCA_1.txt","r")

fs=open(r"reduced_GO-SBERT-embeddings_PCA.txt","r")
for line in fs:
    count+=1
    if count<=10:
       print(line)
    x=line.split(symb)
    #if count<=10:
    ids.append(x[0])
    feat_lst.append(x[1:-1])
fs.close()
print(len(ids),len(feat_lst))

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import optuna
from optuna.pruners import MedianPruner, PatientPruner
import numpy as np


input=len(feat_lst[0])
data=np.array(feat_lst).astype(np.float32)


scale = MinMaxScaler()
data = scale.fit_transform(data)
X_train, X_val = train_test_split(data, test_size=0.1, random_state=42)


train_data= TensorDataset(torch.from_numpy(X_train))
val_data= TensorDataset(torch.from_numpy(X_val))


class AE(nn.Module):
    def __init__(self, orig, inter, latent):
        super(AE, self).__init__()


        self.encoder = nn.Sequential(
            nn.Linear(orig, inter),
            nn.ReLU(),
            nn.Linear(inter, inter // 2),
            nn.ReLU(),
            nn.Linear(inter // 2, latent)
        )


        self.decoder = nn.Sequential(
            nn.Linear(latent, inter// 2),
            nn.ReLU(),
            nn.Linear(inter // 2, inter),
            nn.ReLU(),
            nn.Linear(inter, orig),
            nn.Sigmoid()
        )

    def forward(self, x):
        z = self.encoder(x)
        recon= self.decoder(z)
        return recon, z

    def encode(self, x):
        return self.encoder(x)


def train_func(model, train_load, opt):
    model.train()
    tot_loss= 0
    for bx, in train_load:
        opt.zero_grad()
        recon, _ = model(bx)
        loss = F.mse_loss(recon, bx, reduction='mean')
        loss.backward()
        opt.step()
        tot_loss += loss.item()
    return tot_loss/ len(train_load)

def val_func(model, val_load):
    model.eval()
    tot_loss = 0
    with torch.no_grad():
        for bx, in val_load:
            recon, _ = model(bx)
            loss= F.mse_loss(recon, bx, reduction='mean')
            tot_loss += loss.item()
    return tot_loss / len(val_load)

def objective(trial):
    latent= trial.suggest_categorical('latent_dim', [32])
    inter= trial.suggest_int('intermediate_dim', 128, 512)
    lr = trial.suggest_float('learning_rate', 1e-5, 1e-2, log=True)
    batch= trial.suggest_categorical('batch_size', [32, 64, 128, 256])
    epoch = trial.suggest_int('epochs', 30, 100)

    train_load = DataLoader(train_data, batch_size=batch, shuffle=True)
    val_load= DataLoader(val_data, batch_size=batch)

    model = AE(orig=input, inter=inter, latent=latent)
    opt = optim.Adam(model.parameters(), lr=lr)

    for ep in range(epoch):
        trainl_loss = train_func(model, train_load, opt)
        val_loss = val_func(model, val_load)

        trial.report(val_loss, ep)
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    return val_loss

if __name__ == "__main__":
    pr= PatientPruner(wrapped_pruner=MedianPruner(n_startup_trials=5, n_warmup_steps=5), patience=5)
    study = optuna.create_study(direction='minimize', pruner=pr)
    study.optimize(objective, n_trials=no_of_trials, show_progress_bar=True)

    best_par = study.best_params
    print("\nBest hyp", best_par)


    fmodel= AE(input, best_par['intermediate_dim'], best_par['latent_dim'])
    fopt = optim.Adam(fmodel.parameters(), lr=best_par['learning_rate'])

    full_load = DataLoader(TensorDataset(torch.from_numpy(data)),
                             batch_size=best_par['batch_size'], shuffle=True)

    for epoch in range(best_par['epochs']):
        train_func(fmodel, full_load, fopt)


    fmodel.eval()
    with torch.no_grad():
        full_data = torch.from_numpy(data)
        compr = fmodel.encode(full_data).numpy()
        print("\nCompressed Shape:", compr.shape)


count=0
symb="`"
fet_lst=[]
fs=open(r"2d_GO_SBERT_embeddings.txt","w")
for id,line in zip(ids,compr):
    count+=1
    #if count<=10:
       #print(len(line),line[:5],line[2])
    s=""
    s+=str(id)+symb
    t=""
    for i in range(len(line)):
        s+=str(line[i])+symb
        t+=str(line[i])+symb
    s=s.strip()
    t=t.strip()
    fet_lst.append(t)
    fs.write(s)
    fs.write("\n")
    if count<=10:
       print(s)
fet_lst=list(set(fet_lst))
print(count,len(fet_lst))
fs.close()




#fs=open(r"trial_info_AE_1.txt","w")
fs=open(r"trial_info_2d_GO_SBERT_embeddings.txt","w")
for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    #print("lr",trial.params['learning_rate'],"batch size", trial.params['batch_size'],"encoding dimension",trial.params["latent_dim"],"hidden_dim",trial.params["intermediate_dim"],"epochs",trial.params['epochs'])
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['learning_rate'])+"`"
    s+="batch_size="+str(trial.params['batch_size'])+"`"
    s+="encoding_dim="+str(trial.params['latent_dim'])+"`"
    s+="hidden_dim="+str(trial.params['intermediate_dim'])+"`"
    s+="epochs="+str(trial.params['epochs'])+"`"
    #s+="beta="+str(trial.params['beta'])+"`"
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()


