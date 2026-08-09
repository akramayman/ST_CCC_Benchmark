#!/usr/bin/env python
# coding: utf-8

# # spaCI : deciphering spatial cellular communication through adaptive graph model on S1R1 sample

# ### Import packages

# In[1]:


import os
import sys
import random
import numpy as np
import torch
import yaml

# Make spaCI code importable
sys.path.append('/home/akram/share/spaCI')

from main_yaml import train, predict 


# ### Create the configure file from my data

# In[2]:


cmd = 'python /home/akram/share/spaCI/configuration.py \
  --trainroot /home/akram/share/Result/spaCI_S1R1/triplet.csv \
  --testroot /home/akram/share/Result/spaCI_S1R1/test_pairs.csv \
  --predroot /home/akram/share/Result/spaCI_S1R1/test_lr_pairs.csv \
  --matrixroot /home/akram/share/Result/spaCI_S1R1/exp_data_LR.csv \
  --adjroot /home/akram/share/Result/spaCI_S1R1/spatial_graph.csv\
  --ymlname /home/akram/share/Result/spaCI_S1R1/result/configure_S1R1.yml \
  --threshold 0.5'
os.system(cmd)


# In[3]:


yaml_file = '/home/akram/share/Result/spaCI_S1R1/result/configure_S1R1.yml'
with open(yaml_file) as f:
    cfg = yaml.load(f, Loader=yaml.FullLoader)

cfg['use_cuda'] = 'cpu'
print("Using device:", cfg['use_cuda'])
# In[4]:


seed = 10
random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

########################################
# 4) Train spaCI using main_yaml.train
########################################

print("---- Training spaCI on S1R1 ----")
best_f1 = train(cfg)
print("Best F1 from training:", best_f1)
########################################
# 5) Predict & save embeddings using main_yaml.predict
########################################

print("---- Saving L-R embeddings and predictions ----")
predict(cfg, load_model='best_f1')
print("Done.")