import os
import gzip
import glob
import rdkit 
from shutil import copyfile
from rdkit import Chem
from rdkit.Chem import PandasTools


# Define your template with placeholders
template = """
CALC_INPUT_RMS   True
FORCEFIELD   OPLS_2005
GRIDFILE   {receptor}
LIGANDFILE   {file}
POSE_OUTTYPE   ligandlib_sd
POSES_PER_LIG   1
POSTDOCK_NPOSE   1
PRECISION   SP
"""


def generate_infile(receptor_path, ligand_path, output_file):
    content = template.format(receptor=receptor_path, file=ligand_path)
    with open(output_file, 'w') as f:
        f.write(content)
    return output_file


def ligprep_ligand(ligand_path, ligand_out_path, data_dir):
    ligprep_command = f"{os.getenv('SCHRODINGER')}/ligprep -inp ligprep.inp "\
                      f"-isd {ligand_path} -osd {ligand_out_path} -NJOBS 1 -WAIT -LOCAL"
    os.system(f"cd {data_dir}; {ligprep_command}")
    #result = subprocess.run(ligprep_command, cwd = data_dir, capture_output=True, text=True)
    return #result.stdout, result.stderr


def run_schrodinger(infile, data_dir):
    schrodinger_command = f"{os.getenv('SCHRODINGER')}/glide {infile} -OVERWRITE -NJOBS 1 -TMPLAUNCHDIR -WAIT"
    os.system(f"cd {data_dir}; {schrodinger_command}")
    #result = subprocess.run(schrodinger_command, cwd = './data', capture_output=True, text=True)
    return #result.stdout, result.stderr


def read_schrodinger(sdfgz_file, data_dir):
    assert os.path.exists(os.path.join(data_dir, sdfgz_file))
    sdfgz_filename= sdfgz_file[:-5] + 'sdf.gz'
    copyfile(os.path.join(data_dir, sdfgz_file), os.path.join(data_dir, sdfgz_filename))
    os.system(f"cd {data_dir}; gzip -f -d {sdfgz_filename}")
    result_file = sdfgz_filename[:-3]
    df = parse_result(os.path.join(data_dir, result_file))
    return df
    
    
def parse_result(result_file):
    df = PandasTools.LoadSDF(result_file)
    df['molblock'] = df['ROMol'].apply(lambda mol: Chem.MolToMolBlock(mol) if mol else None)
    cols = ['i_m_source_file_index', 'r_i_docking_score', 'ID', 'molblock']
    return df[cols]


def extract_pdbfile(zipfile, data_dir):
    assert zipfile[-3:] == 'zip'

    os.system(f"cd {data_dir}; unzip {os.path.basename(zipfile)}")
    mae_file = glob.glob(data_dir+"/*.mae")[0]
    
    pdb_file = mae_file[:-3]+"pdb"
    os.system(f"{os.getenv('SCHRODINGER')}/utilities/structconvert {mae_file} {pdb_file}")
    return pdb_file
              
    
def view_in_ngl(molblock):
    # Placeholder: You will need to replace this with the NGL code that visualizes the molblock
    ngl_html = f"""
    <div id="viewport" style="width: 500px; height: 500px;"></div>
    <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.36/dist/ngl.js"></script>
    <script>
        var stage = new NGL.Stage("viewport");
        var stringBlob = new Blob([{molblock}], {{type: 'text/plain'}});
        stage.loadFile(stringBlob, {{ ext: 'sdf' }}).then(function (component) {{
            component.addRepresentation("ball+stick");
            stage.autoView();
        }});
    </script>
    """
    return ngl_html

def view_multiple_in_ngl(molblocks):
    # Generates the NGLView HTML for multiple molblocks
    ngl_html = """
    <div id="viewport" style="width: 500px; height: 500px;"></div>
    <script src="https://unpkg.com/ngl@2.3.1/dist/ngl.js"></script>
    <script>
        var stage = new NGL.Stage("viewport");
    """
    
    for molblock in molblocks:
        ngl_html += f"""
        var stringBlob = new Blob(['{molblock}'], {{type: 'text/plain'}});
        stage.loadFile(stringBlob, {{ ext: 'mol' }}).then(function (component) {{
            component.addRepresentation("ball+stick");
            stage.autoView();
        }});
        """
    
    ngl_html += "</script>"
    return ngl_html