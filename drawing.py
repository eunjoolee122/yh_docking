

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

js_code_backup = """
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
            const molfile = ketcher.getMolfile();
            const smiles = ketcher.getSmiles();
            return [molfile, smiles];  // Fetch molfile from Ketcher
        } else {
            return "Ketcher instance not found";
        }
    }
    return getMolfileFromKetcher(); // Return the molfile text
    }
"""


get_js = """
async () => {
    function getSmilesFromKetcher() {
        var ketcherFrame = document.querySelector('#ifKetcher');
        // Return a promise to wait for the iframe to load
        return new Promise(function (resolve, reject) {
            if (!ketcherFrame) {
                reject("Ketcher iframe not found.");
                return;
            }

            // Wait until the iframe content is fully loaded
            ketcherFrame.onload = function () {
                var ketcher = null;

                if (ketcherFrame && 'contentWindow' in ketcherFrame) {
                    ketcher = ketcherFrame.contentWindow.ketcher;
                }

                if (ketcher) {
                    try {
                        // Fetch the SMILES string from Ketcher
                        resolve(ketcher.getSmiles());
                    } catch (err) {
                        reject("Error calling getSmiles: " + err);
                    }
                } else {
                    reject("Ketcher instance not found.");
                }
            };
        });
        }

    // Call the function and return the result to Gradio
    return getSmilesFromKetcher()
        .then(smiles => smiles)
        .catch(error => error)
        }
"""

def run(hidden_state):
    return f"{hidden_state}"