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
              


def extract_structure_info(file_path):
    # Determine the file type based on the extension
    file_extension = file_path.split('.')[-1].lower()

    if file_extension == 'pdb':
        return extract_pdb_info(file_path)
    elif file_extension in ['sdf']:
        return extract_mol_sdf_info(file_path)
    else:
        return "Unsupported file format.Please insert SDF or PDB file"

def extract_pdb_info(file_path):
    with open(file_path) as f:
        pdb_text = f.read()
    return {os.path.basename(file_path): pdb_text}

def extract_mol_sdf_info(file_path):

    with open(file_path, 'r') as file:
        content = file.read()
    # Split the content into separate molecules based on SDF/MOL format
    #if file_path.endswith('sdf'):
        # Split by the "$$$$" delimiter used in SDF files
    molecules = content.split('$$$$\n')
    #else:
        # For MOL files, typically, there are no strict delimiters, but we will assume
        # each MOL file contains a single molecule and return as a list
    #    molecules = [content]

    result = {}
    for i, mol_text in enumerate(molecules):
        if mol_text.strip():  # Check if the text is not empty
            # Get the molecule name (usually the first line in MOL/SDF format)
            lines = mol_text.strip().splitlines()
            if lines:
                # For SDF, the name is often defined on the first line or the first line after a blank line
                name_line = lines[0].strip()
                molecule_name = name_line.split()[0]  # Extract the first word as the name
                if molecule_name is None:
                    molecule_name = f"os.path.basename(filepath)_{str(i)}"
                # Store the molecule name and content in the dictionary
                result[molecule_name] = mol_text.strip()  # Store the entire content

    return result


