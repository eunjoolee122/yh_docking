

def gen_3dmol_vis(pdb_dict: dict, sdf_dict: dict):
    pdb_html = ""
    if pdb_dict:
        for name, pdb_text in pdb_dict.items():
            pdb_html = """var pdb = `"""+ pdb_text + """`
                    var stringBlob = new Blob([pdb], { type: "text/plain" });
                    stage.loadFile(stringBlob, { ext: "pdb", name :\""""+ name +"""\" }).then(function (component) {
                    component.addRepresentation("ball+stick", { radiusScale: 0.4, opacity: 0.4, colorScheme: "element" });
                    component.autoView();  // Adjust the view to focus on the molecule
                    component.addRepresentation("cartoon");
                }).catch(function (error) {
                    console.error("Error loading the file:", error);
                }); 
                """
    sdf_html = ""
    if sdf_dict:
        for name, sdf_text in sdf_dict.items():
            tmp_html = """var molblock = `"""+ sdf_text + """`
                var stringBlob = new Blob([molblock], { type: "text/plain" });
                stage.loadFile(stringBlob, { ext: "sdf", name :\""""+ name +"""\" }).then(function (component) {
                    component.addRepresentation("ball+stick", { radiusScale: 1.1, colorScheme: "element" });
                    component.autoView();  // Adjust the view to focus on the molecule
                }).catch(function (error) {
                    console.error("Error loading the file:", error);
                }); 
                """
            sdf_html += tmp_html

    x = (
        """<!DOCTYPE html>
            <html>
            <head>    
        <meta http-equiv="content-type" content="text/html; charset=UTF-8" />
        <style>
        body{
            font-family:sans-serif
        }
        .mol-container {
        width: 100%;
        height: 600px;
        position: relative;
        }
        .mol-container select{
            background-image:None;
        }
        </style>
         <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js" integrity="sha512-STof4xm1wgkfm7heWqFJVn58Hm3EtS31XFaagaa8VMReCXAkQnJZ+jEy8PCC/iT18dFy95WcExNHFTqLyp72eQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
        <script src="https://unpkg.com/ngl@2.3.1/dist/ngl.js"></script>
        </head>
        <body>  
        <div id="container" class="mol-container"></div>
            <script>
            document.addEventListener("DOMContentLoaded", function () {
                var stage = new NGL.Stage("container");
                """ + pdb_html + """
                """ + sdf_html + """
                })
        </script>
        </body></html>"""
        )
        

    return f"""<iframe style="width: 100%; height: 600px" name="result" sandbox="allow-modals allow-forms 
    allow-scripts allow-same-origin allow-popups 
    allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""

