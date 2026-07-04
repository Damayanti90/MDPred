

no_of_iterations=50
no_of_splits=10
no_of_trials=50



count=0
ft=open(r"merged_embeddings.txt","w")
fs=open(r"modified_total_chembert_embeddings_approved_drugs.txt","r")
for line in fs:
    count+=1
    ft.write(line.strip())
    ft.write("\n")

fs.close()

fs=open(r"GO-SBERT-embeddings.txt","r")
for line in fs:
    count+=1
    ft.write(line.strip())
    ft.write("\n")
fs.close()


fs=open(r"pathways-SBERT-embeddings.txt","r")
for line in fs:
    count+=1
    ft.write(line.strip())
    ft.write("\n")
fs.close()

ft.close()

import numpy as np

import pandas as pd
import pickle
from sklearn.preprocessing import Normalizer

df = pd.read_csv(r"merged_embeddings.txt", sep="`", header=None)

x = df.iloc[:, 1:-1].values


print(df.head().to_string())


from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

x= sc.fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(0.99)



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
fp=open(r"reduced_merged_embeddings_PCA.txt","w")



for i in mother:
    count+=1
    st=""
    for j in i:
        st=st+str(j)+symb
    if count<=10:
       print(st.strip())
    fp.write(st.strip())
    fp.write("\n")
    
fp.close()
print(count)

count=0
ids=[]
feat_lst=[]
symb="`"


fs=open(r"reduced_merged_embeddings_PCA.txt","r")
for line in fs:
    count+=1
    if count<=10:
       print(line)
    x=line.split(symb)
    
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
fs=open(r"reduced_merged_embeddings_AE.txt","w")
for id,line in zip(ids,compr):
    count+=1
    
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





fs=open(r"trial_info_reduced_merged_embeddings_AE.txt","w")
for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['learning_rate'])+"`"
    s+="batch_size="+str(trial.params['batch_size'])+"`"
    s+="encoding_dim="+str(trial.params['latent_dim'])+"`"
    s+="hidden_dim="+str(trial.params['intermediate_dim'])+"`"
    s+="epochs="+str(trial.params['epochs'])+"`"
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()



fet_lst=[]
id_lst=[]
count=0
l=[]
fs=open(r"reduced_merged_embeddings_AE.txt","r")
for line in fs:
    count+=1
    x=line.split("`")
    l.append(len(x))  
    id_lst.append(x[0])
    fet_lst.append(x[1:-1])
    if count<=10:
       print(line)
fs.close()
print(count)
import numpy as np
fet_lst=np.array(fet_lst)
fet_lst=fet_lst.astype(np.float32)
fet_lst=np.expand_dims(fet_lst, axis=0)
l=list(set(l))
emb_len=l[0]-2
print(len(id_lst),len(fet_lst[0]),l,emb_len)


import torch
import torch.nn as nn
import torch.nn.functional as F
import optuna

input=torch.from_numpy(fet_lst)

class MHSA(nn.Module):
    def __init__(self, emb=emb_len, num_heads=12, dropout=0.1):
        super().__init__()
        self.emb_dim = emb
        self.num_heads = num_heads
        self.head_dim = emb // num_heads
        self.q_proj = nn.Linear(emb, emb)
        self.k_proj = nn.Linear(emb, emb)
        self.v_proj = nn.Linear(emb, emb)
        self.dropout = nn.Dropout(dropout)
        self.out_proj = nn.Linear(emb, emb)

    def forward(self, x):
        batch, seq_len, _ = x.shape
        q = self.q_proj(x).view(batch, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(batch, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(batch, seq_len, self.num_heads, self.head_dim).transpose(1, 2)
        with torch.backends.cuda.sdp_kernel(enable_flash=True, enable_math=True, enable_mem_efficient=True):
            att_out = F.scaled_dot_product_attention(q, k, v, dropout_p=self.dropout.p if self.training else 0.0)
        att_out = att_out.transpose(1, 2).contiguous().view(batch, seq_len, self.emb_dim)
        return self.out_proj(att_out)

def objective(trial):
    num_heads = trial.suggest_categorical("num_heads", [4, 8, 16])
    dropout = trial.suggest_float("dropout", 0.0, 0.3)
    lr = trial.suggest_float("lr", 1e-5, 1e-3, log=True)
    model = MHSA(num_heads=num_heads, dropout=dropout)
    opt = torch.optim.AdamW(model.parameters(), lr=lr)
    targ=input.clone()
    return train_func(model,opt,targ)

def train_func(model,opt,targ):
    model.train()
    opt.zero_grad()
    out = model(input)
    loss = F.mse_loss(out, targ)
    loss.backward()
    opt.step()
    return loss.item()

study = optuna.create_study(direction="minimize")
study.optimize(objective, n_trials=no_of_trials)
print("Best Hyperparameters: ",study.best_params)
best_par = study.best_params
fmodel = MHSA(num_heads=best_par['num_heads'], dropout=best_par['dropout'])
opt = torch.optim.AdamW(fmodel.parameters(), lr=best_par['lr'])
train_func(fmodel,opt,input)
fmodel.eval()
with torch.no_grad():
    fout = fmodel(input)
    femb = fout.squeeze(0)






l=[]
ft=open(r"reduced_merged_embeddings_MHSA.txt","w")
for i in range(len(femb)):

    temp=femb[i]
    
    temp=temp.numpy()
    temp=[str(float(str(t))) for t in temp]
    s=id_lst[i]+"`"+"`".join(temp)+"`"
    if i<=10:
       print(s)
    ft.write(s)
    ft.write("\n")
    fgh=s.split("`")
    l.append(len(fgh))
ft.close()
l=list(set(l))
print(l)




fs=open(r"trial_info_reduced_merged_embeddings_MHSA.txt","w")

for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")  
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['lr'])+"`"
    s+="num_heads="+str(trial.params['num_heads'])+"`"
    s+="dropout="+str(trial.params['dropout'])+"`"
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()

import numpy as np

import pandas as pd
import pickle
from sklearn.preprocessing import Normalizer

df = pd.read_csv(r"reduced_merged_embeddings_MHSA.txt", sep="`", header=None)

x = df.iloc[:, 1:-1].values


print(df.head().to_string())


from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

x= sc.fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(0.99)



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

fp=open(r"reduced_merged_embeddings_MHSA_PCA.txt","w")



for i in mother:
    count+=1
    st=""
    for j in i:
        st=st+str(j)+symb
    if count<=10:
       print(st.strip())
    fp.write(st.strip())
    fp.write("\n")
    
fp.close()
print(count)

count=0
ids=[]
feat_lst=[]
symb="`"

fs=open(r"reduced_merged_embeddings_MHSA_PCA.txt","r")
for line in fs:
    count+=1
    if count<=10:
       print(line)
    x=line.split(symb)
    
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
fs=open(r"reduced_merged_embeddings_MHSA_AE.txt","w")
for id,line in zip(ids,compr):
    count+=1
    
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





fs=open(r"trial_info_reduced_merged_embeddings_MHSA_AE.txt","w")

for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['learning_rate'])+"`"
    s+="batch_size="+str(trial.params['batch_size'])+"`"
    s+="encoding_dim="+str(trial.params['latent_dim'])+"`"
    s+="hidden_dim="+str(trial.params['intermediate_dim'])+"`"
    s+="epochs="+str(trial.params['epochs'])+"`"
    
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()







fet_lst=[]
id_lst=[]
count=0
l=[]
fs=open(r"reduced_merged_embeddings_AE.txt","r")
for line in fs:
    count+=1
    x=line.split("`")
    l.append(len(x)) 
    id_lst.append(x[0])
    fet_lst.append(x[1:-1])
    if count<=10:
       print(line)
fs.close()
print(count)
import numpy as np
fet_lst=np.array(fet_lst)
fet_lst=fet_lst.astype(np.float32)
fet_lst=np.expand_dims(fet_lst, axis=0)
l=list(set(l))
emb_len=l[0]-2
print(len(id_lst),len(fet_lst[0]),l,emb_len)


import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import optuna
import math

input = torch.from_numpy(fet_lst).float()

class Luong(nn.Module):
    def __init__(self, emb, mlen=count):
        super().__init__()
        self.emb = emb
        self.W = nn.Linear(emb, emb, bias=False)
        self.pos_emb = nn.Parameter(torch.zeros(1, mlen, emb))
        nn.init.trunc_normal_(self.pos_emb, std=0.02)
        self.lnorm = nn.LayerNorm(emb)
        self.dropout = nn.Dropout(0.1)

    def forward(self, x):
        batch, seq, _ = x.shape
        x = x + self.pos_emb[:, :seq, :]
        wa = self.W(x)
        scale = math.sqrt(self.emb)
        score = torch.matmul(wa, x.transpose(1, 2)) / scale
        att_wt = F.softmax(score, dim=-1)
        att_wt = self.dropout(att_wt)
        context_vec = torch.matmul(att_wt, x)
        out = self.lnorm(x + context_vec)
        return out


def objective(trial):
    lr = trial.suggest_float("lr", 1e-5, 1e-3, log=True)
    model = Luong(emb=emb_len)
    opt = torch.optim.AdamW(model.parameters(), lr=lr)
    targ= input.clone()
    return train_func(model,opt,targ)

def train_func(model,opt,targ):
    model.train()
    for _ in range(100):
        opt.zero_grad()
        out = model(input)
        loss = F.mse_loss(out, targ)
        loss.backward()
        opt.step()
    return loss.item()


study = optuna.create_study(direction="minimize")
study.optimize(objective, n_trials=no_of_trials)
print("Best Hyperparameters: ",study.best_params)
best_par = study.best_params
fmodel = Luong(emb=emb_len)
opt = torch.optim.AdamW(fmodel.parameters(), lr=best_par['lr'])
train_func(fmodel,opt,input)
fmodel.eval()
with torch.no_grad():
    fout = fmodel(input)
    femb = fout.squeeze(0).numpy()
    
    
    
l=[]

ft=open(r"reduced_merged_embeddings_luong.txt","w")
for i in range(len(femb)):

    temp=femb[i]
    temp=[str(float(str(t))) for t in temp]
    s=id_lst[i]+"`"+"`".join(temp)+"`"
    if i<=10:
       print(s)
    ft.write(s)
    ft.write("\n")
    fgh=s.split("`")
    l.append(len(fgh))
ft.close()
l=list(set(l))
print(l)



fs=open(r"trial_info_reduced_merged_embeddings_luong.txt","w")
for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['lr'])+"`"
    
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()

import numpy as np

import pandas as pd
import pickle
from sklearn.preprocessing import Normalizer

df = pd.read_csv(r"reduced_merged_embeddings_luong.txt", sep="`", header=None)

x = df.iloc[:, 1:-1].values


print(df.head().to_string())


from sklearn.preprocessing import StandardScaler
sc = StandardScaler()

x= sc.fit_transform(x)


from sklearn.decomposition import PCA
pca = PCA(0.99)



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

fp=open(r"reduced_merged_embeddings_luong_PCA.txt","w")



for i in mother:
    count+=1
    st=""
    for j in i:
        st=st+str(j)+symb
    if count<=10:
       print(st.strip())
    fp.write(st.strip())
    fp.write("\n")
    
fp.close()
print(count)

count=0
ids=[]
feat_lst=[]
symb="`"


fs=open(r"reduced_merged_embeddings_luong_PCA.txt","r")
for line in fs:
    count+=1
    if count<=10:
       print(line)
    x=line.split(symb)
    
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
    print("\nBest hyperparameters: ", best_par)


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
fs=open(r"reduced_merged_embeddings_luong_AE.txt","w")
for id,line in zip(ids,compr):
    count+=1
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





fs=open(r"trial_info_reduced_merged_embeddings_luong_AE.txt","w")

for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['learning_rate'])+"`"
    s+="batch_size="+str(trial.params['batch_size'])+"`"
    s+="encoding_dim="+str(trial.params['latent_dim'])+"`"
    s+="hidden_dim="+str(trial.params['intermediate_dim'])+"`"
    s+="epochs="+str(trial.params['epochs'])+"`"
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()

count=0
ft=open(r"merged_embeddings_for_HGT.txt","w")
node_dic={}
fs=open(r"reduced_merged_embeddings_MHSA_AE.txt","r")
for line in fs:
    count+=1
    x=line.split("`",1)
    node_dic[x[0]]=x[1].strip()
fs.close()

fs=open(r"reduced_merged_embeddings_luong_AE.txt","r")
for line in fs:
    count+=1
    x=line.split("`",1)
    temp=node_dic[x[0]]
    temp+=x[1].strip()
    node_dic[x[0]]=temp
fs.close()


l=[]
count=0
for p in node_dic.keys():
    count+=1

    s=p+"`"+node_dic[p]
    s=s.strip()
    if count<=10:
       print(s)
    t=s.split("`")
    l.append(len(t))
    ft.write(s)
    ft.write("\n")


ft.close()
l=list(set(l))
print(count,l)


