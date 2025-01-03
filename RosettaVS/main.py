import os
from prepare import prepare
from GALigandDock import dock_vsx
from FFN import get_model, predict
import numpy as np
import pandas as pd
import shutil

druglike_centroid_path = 'druglike_centroid.csv'
screen_library_path = 'screen_library.csv'
train_size = 500000
test_size = 1000000
dG_start = 10
dG_end = 0.01
dG_intervals = 10
num_top_add = 250000
num_final = 500000

def get_dGs_cutoff(start=10, end=0.01, intervals=10):
    log_values = np.logspace(np.log10(start), np.log10(end), num=intervals)
    return log_values/100

def copy(task_name):
    task_dir = f'task/{task_name}'
    if not os.path.exists(task_dir): os.mkdir(task_dir)
    shutil.copy(f'task/dock_vsh.xml', task_dir)
    shutil.copy(f'task/dock_vsx.xml', task_dir)
    shutil.copy(f'task/receptor.pdb', task_dir)
    shutil.copy(f'task/ligand.mol2', task_dir)


df_cluster = pd.read_csv(druglike_centroid_path, sep='\t')
df_screen = pd.read_csv(screen_library_path, sep='\t')
dG_cutoffs = get_dGs_cutoff(start=dG_start, end=dG_end, intervals=dG_intervals)


for i,dG_cutoff in enumerate():
    if i == 0:
        sampled_df = df_cluster.sample(n=train_size+test_size, random_state=42)
        df_iter_train = sampled_df.iloc[:train_size, :]
        copy('train{i}')
        copy('test')
        df_iter_train.to_csv(f'task/train{i}/train{i}.csv', sep='\t')
        df_test = sampled_df.iloc[train_size:test_size, :]
        df_test.to_csv(f'task/test/test.csv', sep='\t')
        task_name = prepare(f'task/train{i}/train{i}.csv')
        df_train = dock_vsx(task_name)
        task_name = prepare(f'task/test/test.csv')
        df_test= dock_vsx(task_name)
    else:
        copy('train{i}')
        df_iter_train_add.to_csv(f'task/train{i}/train{i}.csv')
        task_name = prepare(f'task/train{i}/train{i}.csv')
        df_train = pd.concat([df_train, dock_vsx(task_name)], axis=0,)
    model = get_model(df_train, df_test, dG_cutoff)
    df_probs = predict(model, df_screen)
    if i < len(dG_cutoffs)-1 :
        df_tops = df_probs.nlargest(num_top_add, 'Probability')
        df_random = df_probs.sample(n=num_top_add)
        df_iter_train_add = pd.concat([df_tops,df_random], axis=0,)
    else:
        df_final = df_probs.nlargest(num_final, 'Probability')
        copy('final')
        df_final.to_csv(f'task/final/final.csv')
        task_name = prepare(f'task/final/final.csv')
        df_final= dock_vsx(task_name)
        df_final.to_csv(f'task/final/final_result.csv')

