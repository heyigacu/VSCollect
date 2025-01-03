
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from dimorphite_dl import DimorphiteDL
import shutil

PARENT_DIR = os.path.abspath(os.path.dirname(__file__))
ROSETTA_DIR = '/home/hy/Softwares/rosetta/'
mol2genparams_app = ROSETTA_DIR+'source/scripts/python/public/generic_potential/mol2genparams.py'

# /home/hy/Softwares/rosetta/source/scripts/python/public/generic_potential/mol2genparams.py
def protonation(smileses, pH=7.4):
    rst_smileses = []
    dimorphite_dl = DimorphiteDL(
        min_ph=pH,
        max_ph=pH,
        max_variants=128,
        label_states=False,
        pka_precision=1.0
    )
    for smiles in smileses:
        rst_smiles = dimorphite_dl.protonate(smiles)
        if len(rst_smiles) == 0:
            rst_smileses.append(None)
        else:
            rst_smileses.append(rst_smiles[0])
    return rst_smileses


def gen_3d_mol2(smileses, names, output_dir):
    output_mol2_paths = [os.path.join(output_dir, name+'.mol2') for name in names]
    mol2_paths = []
    for index,smiles in enumerate(smileses):
        if smiles is not None:
            os.system(f'obabel -:"{smiles}" -O {output_mol2_paths[index]} --gen3d --partialcharge mmff94')
        if os.path.exists(output_mol2_paths[index]):
            mol2_paths.append(output_mol2_paths[index])
        else:
            mol2_paths.append(None)
    return mol2_paths


def gen_rosetta_parmas(mol2genparams_app, mol2_paths, names, output_dir):
    rosetta_parmas_ls = []
    for index,mol2_path in enumerate(mol2_paths):
        if mol2_path != None and os.path.exists(mol2_path):
            os.system(f'{mol2genparams_app} -s {mol2_path} --resname {names[index]} --outdir {output_dir} --no_pdb')
            rosetta_parmas_path = os.path.join(output_dir, os.path.splitext(os.path.basename(mol2_path))[0]+'.params')
            if os.path.exists(rosetta_parmas_path):
                rosetta_parmas_ls.append(rosetta_parmas_path)
            else:
                rosetta_parmas_ls.append(None)
        else:
            rosetta_parmas_ls.append(None)
    return rosetta_parmas_ls

def write_list(mol2_paths, rosetta_parmas_ls, names, output_dir):
    with open(output_dir+'/ligand_list.txt', 'w') as f1:
        with open(output_dir+'/ligand_params.txt', 'w') as f2:
            for index,mol2_path in enumerate(mol2_paths):
                if mol2_path != None and rosetta_parmas_ls[index] != None and os.path.exists(mol2_path) and os.path.exists(rosetta_parmas_ls[index]):
                    if os.path.splitext(os.path.basename(mol2_path))[0] == os.path.splitext(os.path.basename(rosetta_parmas_ls[index]))[0] == names[index]:
                        f1.write(names[index]+'\n')
                        f2.write('-extra_res_fa '+rosetta_parmas_ls[index]+'\n')



def prepare_complex(mol2genparams_app, ligand_mol2_path, receptor_pdb_path, output_dir):
    os.system(f'obabel {ligand_mol2_path} -O {ligand_mol2_path} --partialcharge mmff94')
    os.system(f'cd {output_dir}; {mol2genparams_app} -s {ligand_mol2_path} --resname ligand --outdir {output_dir}')
    ligand_pdb_path = output_dir+'/ligand_0001.pdb'
    ligand_rosetta_params_path = output_dir+'/ligand_0001.params' 
    with open(receptor_pdb_path, 'r') as f1:
        content1 = f1.read()
    with open(ligand_pdb_path, 'r') as f2:
        content2 = f2.read()
    combined_content = content1 + content2
    with open(output_dir+'/complex.pdb', 'w') as f_out:
        f_out.write(combined_content)


def prepare(df_path):
    # 0. parepare task dir
    df = pd.read_csv(df_path, sep='\t')
    smileses, names = list(df['Smiles']), list(df['Name'])
    task_name = os.path.splitext(os.path.basename(df_path))[0]
    
    task_dir = f'{PARENT_DIR}/task/{task_name}'
    if not os.path.exists(task_dir): os.mkdir(task_dir)
    for dir in [f'{task_dir}/mol2_dir', f'{task_dir}/parms_dir']:
        if not os.path.exists(dir): os.mkdir(dir)
    # 1. dimorphite_dl protonation
    rst_smileses = protonation(smileses)
    # 2. openbabel generate 3D mol2 with mmff94
    mol2_paths = gen_3d_mol2(rst_smileses, names, f'{task_dir}/mol2_dir')
    # 3. mol2genparams.py generate rosetta parmas
    rosetta_parmas_ls = gen_rosetta_parmas(mol2genparams_app, mol2_paths, names,  f'{task_dir}/parms_dir') 
    # 4. output prepared ligands csv
    write_list(mol2_paths, rosetta_parmas_ls, names, task_dir)
    df.insert(len(df.columns), 'mol2_path', mol2_paths)
    df.insert(len(df.columns), 'rosetta_parmas_path', rosetta_parmas_ls)
    df.to_csv(f'{task_dir}/prepare.csv', index=False, sep='\t')
    # 5. prepare reference complex
    prepare_complex(mol2genparams_app, ligand_mol2_path=task_dir+'/ligand.mol2', receptor_pdb_path=task_dir+'/receptor.pdb', output_dir=task_dir)

if __name__ == '__main__':
    prepare('task/test1/test1.csv')

