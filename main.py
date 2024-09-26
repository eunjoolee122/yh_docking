import os
import logging
import tempfile
import collections
from shutil import copyfile
from datetime import datetime
import mol_ngl_viewer as mol_viewer
import pandas as pd
import gradio as gr
from rdkit import Chem
from rdkit.Chem import PandasTools
from functions import generate_infile, ligprep_ligand, run_schrodinger, read_schrodinger, extract_pdbfile, parse_result

SAVE_PATH = './tmp'
os.makedirs(SAVE_PATH, exist_ok=True)

LOG_FILE = 'schdock.log'
logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


def run_docking(grid_in, lig_filein, lig_smilesin):
    
    _temp_dir = tempfile.TemporaryDirectory()
    temp_dir = _temp_dir.name
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    ## copy ligprep input file
    infile = os.path.join(temp_dir, 'ligprep.inp')
    copyfile('./data/ligprep.inp', infile)
    
    grid_file = os.path.join(temp_dir, os.path.basename(grid_in['path']))
    copyfile(grid_in['path'], grid_file)
    logging.info(f"grid file : {grid_file}")

    pdb_file = extract_pdbfile(grid_file, temp_dir)
    ## structure data load
    with open(pdb_file) as f:
        pdb_text = f.read()
    pdb_dict = {os.path.basename(grid_in['path']) : pdb_text}
    logging.info(f"grid file is converted to pdb file")
    if lig_filein:
        lig_file = os.path.join(temp_dir, os.path.basename(lig_filein['path']))
        copyfile(lig_filein['path'], lig_file)
        logging.info(f"ligand file : {lig_file}")

    elif lig_smilesin:
        m =Chem.MolFromSmiles(lig_smilesin)
        inchi = Chem.MolToInchiKey(m)
        lig_file = os.path.join(temp_dir, f'ligand_{inchi}.sdf')
        with Chem.SDWriter(lig_file) as w:
            w.write(m)
        logging.info(f"ligand smiles : {lig_smilesin}")

    # Validate file paths
    if not pdb_file or (not lig_file):
        logging.error("Receptor and ligand file paths are required")
        return "error : Receptor and ligand file paths are required"
    ligprep_path = os.path.join(temp_dir, 'ligand_ligprep.sdf')
    infile = os.path.join(temp_dir, 'glide.in')
    
    ### ligand preparation 
    ligprep_ligand(os.path.basename(lig_file), os.path.basename(ligprep_path), temp_dir)    
    logging.info(f"ligprep success")

    assert os.path.exists(ligprep_path)
    generate_infile(os.path.basename(grid_file), os.path.basename(ligprep_path), infile)
    logging.info(f"glide input file is generated")

    assert os.path.exists(infile)
    run_schrodinger(os.path.basename(infile), temp_dir)
    logging.info(f"glide docking ends")
    
    sdfgz_file = 'glide_lib.sdfgz'
    assert os.path.exists(os.path.join(temp_dir, sdfgz_file))
    result_file = os.path.join(temp_dir, 'glide_lib.sdf')
    logging.info(f"glide results...")
    ## copy result from temporary folder to static folder
    df = read_schrodinger(sdfgz_file, temp_dir)
    save_sdffile = os.path.join(SAVE_PATH, f"{current_time}_glide_result.sdf")
    copyfile(result_file, save_sdffile)
    df.sort_values(by='r_i_docking_score',ascending=False, inplace=True)

    view_selector_content = dict()

    for i, row in df.iterrows():
        docking_score = float(row['r_i_docking_score'])
        if docking_score is None:
            continue
        label = f"Docking Score {docking_score:.2f} {row['ID']}"
        sdf_text = row['molblock']
        output_viz = "Output visualisation unavailable"

        view_selector_content[label] = row['molblock'] ## output_viz
        
    labels = list(view_selector_content.keys())
    init_value = labels[0] if labels else None

    dropdown = gr.CheckboxGroup(interactive=True, label="Docking results",
                               choices=labels, value=init_value)
    return save_sdffile, view_selector_content, dropdown, pdb_dict
   
      
def update_view(view_selector_content, view_result_selector, pdb_check, pdb_dict):
    default_str="Output visualisation unavailable"
   
    sdf_dict = {}
    if view_selector_content and view_result_selector:
        sdf_texts = [view_selector_content.get(key) for key in view_result_selector]
        sdf_dict = {}

        for k in view_selector_content.keys():
            if k in view_result_selector:
                sdf_dict[k] = view_selector_content[k]

    if pdb_check:
        pdb_dict = pdb_dict
    else:
        pdb_dict = None
    if pdb_dict or sdf_dict:
        html_content = mol_viewer.gen_3dmol_vis(pdb_dict, sdf_dict)
        
        return html_content
    return default_str
    

def main_demo():
    
    with gr.Blocks(title="Schrodinger Glide Web") as demo:
        
        gr.Markdown("# Schrodinger Docking")

        with gr.Row():
            grid_in = gr.File(label="Protein Grid", height=256)
            lig_filein = gr.File(label="Ligand SDF")
            lig_smilesin = gr.Textbox(label="Ligand Smiles", placeholder="Provide smiles input")
            
        run_btn = gr.Button("Run")
        with gr.Blocks():
            gr.Markdown("# Output")
            with gr.Row():
                with gr.Column(scale=1, min_width=300):
                    output_file = gr.File(label="Output Files")
                    pdb_dict=gr.State()
                    pdbcheck = gr.Checkbox(interactive=True, label="PDB visualization")
                    init_value = "Schrodinger docking visualization"
                    view_result_selector = gr.CheckboxGroup(interactive=True, label="docking results")
                with gr.Column(scale=2, min_width=300):
                    viewer = gr.HTML(value=init_value, label="Protein Viewer", show_label=False)
        view_selector_content = gr.State()
        
        _inputs = [grid_in, lig_filein, lig_smilesin]
        _outputs = [output_file, view_selector_content, view_result_selector, pdb_dict]
        run_btn.click(fn=run_docking, inputs=_inputs, outputs=_outputs, preprocess=False)

        view_result_selector.change(fn=update_view,
                                    inputs=[view_selector_content, view_result_selector, pdbcheck, pdb_dict],
                                    outputs=viewer)
        pdbcheck.change(fn=update_view,
                 inputs=[view_selector_content, view_result_selector, pdbcheck, pdb_dict],
                 outputs=viewer)  
    demo.launch(
        server_name='0.0.0.0'
        )

if __name__ =='__main__':
    main_demo()