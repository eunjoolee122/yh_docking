from rdkit import Chem

def ketcherdraw():       

    return f"""<iframe id="ifKetcher" src="/static/ketcher/index.html" width="100%" height="400px"></iframe>"""


# JavaScript code to fetch SMILES from the Ketcher iframe
js_code = """
async () => {
    function getMolfileFromKetcher() {
        var ketcherFrame = document.querySelector('#ifKetcher');
        var ketcher = null;

        if (ketcherFrame && 'contentDocument' in ketcherFrame) {
            ketcher = ketcherFrame.contentWindow.ketcher;
        } else if (ketcherFrame) {
            ketcher = document.frames[ketcherFrame].window.ketcher;
        }

        if (ketcher) {
            return ketcher.getMolfile();
            //const smiles = ketcher.getSmiles();
            //return {"mol" : molfile, "smiles" : smiles};  // Return molfile and smiles as an array
        } else {
            return "Ketcher instance not found";
        }
    }
    return getMolfileFromKetcher(); // Return the molfile text
    }
"""

def run(hidden_state):
    molfile =  f"{hidden_state}"
    print(molfile)
    m = Chem.MolFromMolBlock(molfile)
    smiles = Chem.MolToSmiles(m)
    return [molfile, smiles]