

no_of_iterations=50
no_of_splits=10
no_of_trials=50

count=0
import numpy as np
import torch
ids_drug=[]
feat_lst_drug=[]
ids_protein=[]
feat_lst_protein=[]
ids_go=[]
feat_lst_go=[]
ids_path=[]
feat_lst_path=[]
symb="`"
count_dr=0
count_go=0
count_path=0
fs=open(r"merged_embeddings_for_HGT.txt","r")
for line in fs:
    count+=1
    if count<=10:
       print(line)
    x=line.split(symb)
    if x[0].startswith("DB"):
       ids_drug.append(x[0])
       feat_lst_drug.append(x[1:-1])
       count_dr+=1
       
    elif x[0].startswith("GO"):
       ids_go.append(x[0])
       feat_lst_go.append(x[1:-1])
       count_go+=1
    elif x[0].startswith("R-HSA"):
       ids_path.append(x[0])
       feat_lst_path.append(x[1:-1])
       count_go+=1
fs.close()
print(len(ids_path),len(feat_lst_path))
data=np.array(feat_lst_drug).astype(np.float32)
feat_tensor_drug=torch.from_numpy(data)
data=np.array(feat_lst_protein).astype(np.float32)
feat_tensor_protein=torch.from_numpy(data)
data=np.array(feat_lst_go).astype(np.float32)
feat_tensor_go=torch.from_numpy(data)
data=np.array(feat_lst_path).astype(np.float32)
feat_tensor_path=torch.from_numpy(data)
print(len(feat_tensor_drug),len(feat_tensor_drug[0]))
print(len(feat_tensor_go),len(feat_tensor_go[0]))
print(len(feat_tensor_path),len(feat_tensor_path[0]))

drug_go=[]
drug_pathway=[]
drug_drug=[]
drug_protein=[]
protein_protein=[]
protein_go=[]
protein_pathway=[]
symb="`"
fs=open(r"positive_interaction_data.txt","r")
for line in fs:
    x=line.split(symb)
    t=x[1]+symb+x[0]+symb
    t=t.strip()
    if x[0] in ids_drug and x[1] in ids_go:
       drug_go.append(line.strip())
    elif x[0] in ids_drug and x[1] in ids_path:
       drug_pathway.append(line.strip())
    elif x[0] in ids_drug and x[1] in ids_drug:
       drug_drug.append(line.strip())
    elif x[0] in ids_protein and x[1] in ids_protein:
       protein_protein.append(line.strip())
    elif x[0] in ids_drug and x[1] in ids_protein:
       drug_protein.append(line.strip())
    elif x[0] in ids_protein and x[1] in ids_go:
       protein_go.append(line.strip())
    elif x[0] in ids_protein and x[1] in ids_path:
       protein_pathway.append(line.strip())

    elif x[1] in ids_drug and x[0] in ids_go:
       drug_go.append(t.strip())
    elif x[1] in ids_drug and x[0] in ids_path:
       drug_pathway.append(t.strip())
    elif x[1] in ids_drug and x[0] in ids_protein:
       drug_protein.append(t.strip())
    elif x[1] in ids_protein and x[0] in ids_go:
       protein_go.append(t.strip())
    elif x[1] in ids_protein and x[0] in ids_path:
       protein_pathway.append(t.strip())
fs.close()
print(len(drug_drug),len(drug_go),len(drug_pathway),len(drug_protein),len(protein_protein),len(protein_go),len(protein_pathway))

coord_a=[]
coord_b=[]
for line in drug_drug:
    x=line.split(symb)
    coord_a.append(ids_drug.index(x[0]))
    coord_b.append(ids_drug.index(x[1]))
drug_drug_tensor=[coord_a,coord_b]
drug_drug_tensor=np.array(drug_drug_tensor)
drug_drug_tensor=torch.from_numpy(drug_drug_tensor)
print(len(drug_drug_tensor),len(drug_drug_tensor[0]))



coord_a=[]
coord_b=[]
for line in drug_go:
    x=line.split(symb)
    coord_a.append(ids_drug.index(x[0]))
    coord_b.append(ids_go.index(x[1]))
drug_go_tensor=[coord_a,coord_b]
drug_go_tensor=np.array(drug_go_tensor)
drug_go_tensor=torch.from_numpy(drug_go_tensor)
print(len(drug_go_tensor),len(drug_go_tensor[0]))


coord_a=[]
coord_b=[]
for line in drug_pathway:
    x=line.split(symb)
    coord_a.append(ids_drug.index(x[0]))
    coord_b.append(ids_path.index(x[1]))
drug_pathway_tensor=[coord_a,coord_b]
drug_pathway_tensor=np.array(drug_pathway_tensor)
drug_pathway_tensor=torch.from_numpy(drug_pathway_tensor)
print(len(drug_pathway_tensor),len(drug_pathway_tensor[0]))



neg_drug_go=[]
neg_drug_protein=[]
neg_drug_pathway=[]
neg_drug_drug=[]
neg_protein_protein=[]
neg_protein_go=[]
neg_protein_pathway=[]
symb="`"
fs=open(r"hard_negative_samples.txt","r")
for line in fs:
    x=line.split(symb)
    t=x[1]+symb+x[0]+symb
    t=t.strip()
    if x[0] in ids_drug and x[1] in ids_go:
       neg_drug_go.append(line.strip())
    elif x[0] in ids_drug and x[1] in ids_path:
       neg_drug_pathway.append(line.strip())
    elif x[0] in ids_drug and x[1] in ids_drug:
       neg_drug_drug.append(line.strip())
    elif x[0] in ids_protein and x[1] in ids_protein:
       neg_protein_protein.append(line.strip())
    elif x[0] in ids_drug and x[1] in ids_protein:
       neg_drug_protein.append(line.strip())
    elif x[0] in ids_protein and x[1] in ids_go:
       neg_protein_go.append(line.strip())
    elif x[0] in ids_protein and x[1] in ids_path:
       neg_protein_pathway.append(line.strip())

    elif x[1] in ids_drug and x[0] in ids_go:
       neg_drug_go.append(t.strip())
    elif x[1] in ids_drug and x[0] in ids_path:
       neg_drug_pathway.append(t.strip())
    elif x[1] in ids_drug and x[0] in ids_protein:
       neg_drug_protein.append(t.strip())
    elif x[1] in ids_protein and x[0] in ids_go:
       neg_protein_go.append(t.strip())
    elif x[1] in ids_protein and x[0] in ids_path:
       neg_protein_pathway.append(t.strip())
fs.close()
print(len(neg_drug_drug),len(neg_drug_go),len(neg_drug_pathway),len(neg_drug_protein),len(neg_protein_protein),len(neg_protein_go),len(neg_protein_pathway))

coord_a=[]
coord_b=[]
for line in neg_drug_drug:
    x=line.split(symb)
    coord_a.append(ids_drug.index(x[0]))
    coord_b.append(ids_drug.index(x[1]))
neg_drug_drug_tensor=[coord_a,coord_b]
neg_drug_drug_tensor=np.array(neg_drug_drug_tensor)
neg_drug_drug_tensor=torch.from_numpy(neg_drug_drug_tensor)
print(len(neg_drug_drug_tensor),len(neg_drug_drug_tensor[0]))



coord_a=[]
coord_b=[]
for line in neg_drug_go:
    x=line.split(symb)
    coord_a.append(ids_drug.index(x[0]))
    coord_b.append(ids_go.index(x[1]))
neg_drug_go_tensor=[coord_a,coord_b]
neg_drug_go_tensor=np.array(neg_drug_go_tensor)
neg_drug_go_tensor=torch.from_numpy(neg_drug_go_tensor)
print(len(neg_drug_go_tensor),len(neg_drug_go_tensor[0]))


coord_a=[]
coord_b=[]
for line in neg_drug_pathway:
    x=line.split(symb)
    coord_a.append(ids_drug.index(x[0]))
    coord_b.append(ids_path.index(x[1]))
neg_drug_pathway_tensor=[coord_a,coord_b]
neg_drug_pathway_tensor=np.array(neg_drug_pathway_tensor)
neg_drug_pathway_tensor=torch.from_numpy(neg_drug_pathway_tensor)
print(len(neg_drug_pathway_tensor),len(neg_drug_pathway_tensor[0]))



from torch_geometric.data import HeteroData

data = HeteroData()
data['drug'].x =feat_tensor_drug
data['go'].x = feat_tensor_go
data['pathway'].x = feat_tensor_path


data['drug', 'assoc', 'drug'].edge_index = drug_drug_tensor
data['drug', 'relates', 'go'].edge_index = drug_go_tensor
data['drug', 'impacts', 'pathway'].edge_index = drug_pathway_tensor
data['drug', 'no_interaction', 'drug'].edge_index = neg_drug_drug_tensor
data['drug', 'no_interaction', 'go'].edge_index = neg_drug_go_tensor
data['drug', 'no_interaction', 'pathway'].edge_index = neg_drug_pathway_tensor


import torch_geometric.transforms as T
transform = T.ToUndirected()
data_undirected = transform(data)
data=data_undirected

print(data)
print(f"Number of nodes: {data.num_nodes}")
print(f"Number of edges: {data.num_edges}")
print(f"Has isolated nodes: {data.has_isolated_nodes()}")
print(f"Has self loops: {data.has_self_loops()}")
print(f"Is undirected: {data.is_undirected()}")


import torch
import torch.nn.functional as F
from torch_geometric.nn import HGTConv
import optuna
from torch_geometric.nn import Linear


class HGT(torch.nn.Module):
    def __init__(self, metadata, hidden_channel, out_channel, tot_heads, tot_layers):
        super().__init__()
        self.pre= torch.nn.ModuleDict()
        for node_type in metadata[0]:
            self.pre[node_type] = Linear(-1, hidden_channel)

        self.block = torch.nn.ModuleList()
        for _ in range(tot_layers):
            self.block.append(HGTConv(hidden_channel, hidden_channel, metadata, heads=tot_heads))

        self.post = Linear(hidden_channel, out_channel)

    def forward(self, x_dict, edge_index_dict):
        temp = {}
        for node_type, x in x_dict.items():
             a= self.pre[node_type]
             b = a(x)
             temp[node_type] = b.relu()
        x_dict=temp
        for bl in self.block:
            x_dict = bl(x_dict, edge_index_dict)
            temp_1={}
            for node_type,x in x_dict.items():
                temp_1[node_type]=x.relu()
            x_dict=temp_1
        temp_2={}
        for node_type,x in x_dict.items():
            temp_2[node_type]=self.post(x)
        return temp_2


def objective(trial, data):

    hidden = trial.suggest_categorical("hidden_channels", [16, 32, 64, 128])
    head = trial.suggest_categorical("num_heads", [4, 8])
    layer = trial.suggest_int("num_layers", 1, 3)
    lr = trial.suggest_float("lr", 1e-4, 5e-3, log=True)
    latent = trial.suggest_categorical("latent_dim", [128, 256])
    epoch = trial.suggest_int("epochs", 20, 100)
    model = HGT(data.metadata(), hidden, latent, head, layer)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    return train_func(model,opt,epoch)

def train_func(model, opt,epoch):
    model.train()
    for _ in range(epoch):
        opt.zero_grad()
        out = model(data.x_dict, data.edge_index_dict)
        tot_loss = 0
        for edge_type, edge_index in data.edge_index_dict.items():
            src, rel, dst = edge_type
            z_u, z_v = out[src], out[dst]
            scores = (z_u[edge_index[0]] * z_v[edge_index[1]]).sum(dim=-1)
            if rel == "no_interaction":
                rloss = F.softplus(scores).mean()
            else:
                rloss = -F.logsigmoid(scores).mean()
            gloss = (z_u.mean(dim=0) * z_v.mean(dim=0)).sum()
            tot_loss += rloss + (0.01 * gloss)
        tot_loss.backward()
        opt.step()
    val_loss = tot_loss.item()
    return val_loss


study = optuna.create_study(direction="minimize")
study.optimize(lambda trial: objective(trial, data), n_trials=no_of_trials)
print(f"Best Params: {study.best_params}")
best_par = study.best_params
fmodel = HGT(data.metadata(), best_par['hidden_channels'], best_par['latent_dim'], best_par['num_heads'], best_par['num_layers'])
fmodel.eval()
with torch.no_grad():
    femb = fmodel(data.x_dict, data.edge_index_dict)
print(f"\nBest Hyperparameters: {best_par}")
print(f"Embedding Shape: {femb['drug'].shape}")


fs=open(r"trial_info_hgt.txt","w")
for trial in study.trials:
    if trial.state == optuna.trial.TrialState.PRUNED:
      continue
    print(f"Trial {trial.number}: Value={trial.value}, Params={trial.params}")
    s=""
    s+="trial="+str(trial.number+1)+"`"
    s+="value="+str(trial.value)+"`"
    s+="lr="+str(trial.params['lr'])+"`"
    s+="hidden_channels="+str(trial.params['hidden_channels'])+"`"
    s+="num_heads="+str(trial.params['num_heads'])+"`"
    s+="num_layers="+str(trial.params['num_layers'])+"`"
    s+="epochs="+str(trial.params['epochs'])+"`"
    s+="latent_dim="+str(trial.params['latent_dim'])+"`"
    s=s.strip()
    print(s)
    fs.write(s)
    fs.write("\n")
fs.close()

count=0
count_drug=0
symb="`"
fet_lst=[]
fs=open(r"merged_embeddings_hgt.txt","w")
for id,line in zip(ids_drug,femb['drug']):
    count+=1
    count_drug+=1
    s=""
    s+=str(id)+symb
    t=""
    for i in range(len(line)):
        s+=str(float(line[i]))+symb
        t+=str(float(line[i]))+symb
    s=s.strip()
    t=t.strip()
    fet_lst.append(t)
    fs.write(s)
    fs.write("\n")
    if count_drug<=10:
       print(s)
fet_lst=list(set(fet_lst))



count_go=0
for id,line in zip(ids_go,femb['go']):
    count+=1
    count_go+=1
    s=""
    s+=str(id)+symb
    t=""
    for i in range(len(line)):
        s+=str(float(line[i]))+symb
        t+=str(float(line[i]))+symb
    s=s.strip()
    t=t.strip()
    fet_lst.append(t)
    fs.write(s)
    fs.write("\n")
    if count_go<=10:
       print(s)
fet_lst=list(set(fet_lst))



count_path=0
for id,line in zip(ids_path,femb['pathway']):
    count+=1
    count_path+=1
    s=""
    s+=str(id)+symb
    t=""
    for i in range(len(line)):
        s+=str(float(line[i]))+symb
        t+=str(float(line[i]))+symb
    s=s.strip()
    t=t.strip()
    fet_lst.append(t)
    fs.write(s)
    fs.write("\n")
    if count_path<=10:
       print(s)
fet_lst=list(set(fet_lst))
print(count,len(fet_lst),len(fet_lst[0]))
fs.close()

import numpy as np

import pandas as pd
import pickle
from sklearn.preprocessing import Normalizer

df = pd.read_csv(r"merged_embeddings_hgt.txt", sep="`", header=None)
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
fp=open(r"reduced_merged_embeddings_hgt.txt","w")



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


