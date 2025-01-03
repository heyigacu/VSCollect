
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from dimorphite_dl import DimorphiteDL
PARENT_DIR = os.path.abspath(os.path.dirname(__file__))
ROSETTA_DIR = '/home/hy/Softwares/rosetta/'
GALigandDock_app = ROSETTA_DIR+'source/bin/rosetta_scripts.linuxgccrelease'



def vsx_dock(task_name, output_dir, GALigandDock_app):
    cmd = f"cd {output_dir} ; {GALigandDock_app}\
    @ ./ligand_params.txt \
    -s ./complex.pdb \
    -extra_res_fa ./ligand.params \
    -gen_potential \
    -overwrite \
    -beta_cart \
    -parser:protocol ./dock_vsx.xml \
    -parser:script_vars liglist=./ligand_list.txt \
    -no_autogen_cart_improper \
    -multi_cool_annealer 10 \
    -missing_density_to_jump \
    -score::hb_don_strength hbdon_GENERIC_SC:1.45 \
    -score::hb_acc_strength hbacc_GENERIC_SP2SC:1.19 \
    -score::hb_acc_strength hbacc_GENERIC_SP3SC:1.19 \
    -score::hb_acc_strength hbacc_GENERIC_RINGSC:1.19 \
    -out:levels all:300 protocols.ligand_docking.GALigandDock:300 \
    -out:prefix {task_name}. \
    -out:file:silent ./output.out \
    -out:file:scorefile ./score.sc \
    -mute all \
    "
    os.system(cmd)


def vsh_dock(task_name, output_dir, GALigandDock_app):
    cmd = f"cd {output_dir} ; {GALigandDock_app}\
    @ ./ligand_params.txt \
    -s ./complex.pdb \
    -extra_res_fa ../ligand.params \
    -gen_potential \
    -overwrite \
    -beta_cart \
    -parser:protocol ./dock_vsh.xml \
    -parser:script_vars liglist=./ligand_list.txt \
    -no_autogen_cart_improper \
    -multi_cool_annealer 10 \
    -missing_density_to_jump \
    -score::hb_don_strength hbdon_GENERIC_SC:1.45 \
    -score::hb_acc_strength hbacc_GENERIC_SP2SC:1.19 \
    -score::hb_acc_strength hbacc_GENERIC_SP3SC:1.19 \
    -score::hb_acc_strength hbacc_GENERIC_RINGSC:1.19 \
    -out:levels all:300 protocols.ligand_docking.GALigandDock:300 \
    -out:prefix {task_name}. \
    -out:file:silent ./output.out \
    -out:file:scorefile ./score.sc \
    -mute all \
    "
    os.system(cmd)


def read_dock_result(output_dir):
    df_prepare = pd.read_csv(output_dir+'/prepare.csv', sep='\t')
    df_result = pd.read_csv(output_dir+'/score.sc', delim_whitespace=True, skiprows=1)
    df_result = df_result[['total_score', 'dG', 'ligandname']]
    df_result.columns = ['total_score', 'dG', 'Name']
    merged_df = pd.merge(df_prepare, df_result, on='Name', how='left')
    merged_df.to_csv(output_dir+'/dock.csv', sep='\t', index=False)
    return merged_df


def dock_vsx(task_name):
    task_dir = os.path.join(PARENT_DIR,'task',task_name)
    vsx_dock(task_name, task_dir, GALigandDock_app)
    return read_dock_result(output_dir=task_dir)

def dock_vsh(task_name):
    task_dir = os.path.join(PARENT_DIR,'task',task_name)
    vsh_dock(task_dir, GALigandDock_app)
    read_dock_result(output_dir=task_dir)

if __name__ == '__main__':
    dock_vsx('test1')