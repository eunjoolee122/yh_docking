

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
            html, body {
                height: 100%;
                margin: 0;
                padding: 0;
                }
            body{
                font-family:sans-serif
            }
            .mol-container {
            width: 100%;
            height: 100%;
            position: relative;
            }
            #viewport {
                width: 100%;
                height: 100%;
            }
            .mol-container select{
                background-image:None;
            }
            </style>
         <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js" integrity="sha512-STof4xm1wgkfm7heWqFJVn58Hm3EtS31XFaagaa8VMReCXAkQnJZ+jEy8PCC/iT18dFy95WcExNHFTqLyp72eQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
        <script src="https://unpkg.com/ngl@2.3.1/dist/ngl.js"></script>
        </head>
        <body>  
        <div id="container" class="mol-container"><div id="viewport" class="viewport"></div></div>
            <script>
            // Function to resize the viewport
            function resizeViewport() {
                var viewport = document.getElementById("viewport");
                var container = viewport.parentElement;

                // Adjust the width and height of #viewport to match the container size
                viewport.style.width = container.clientWidth + "px";
                viewport.style.height = container.clientHeight + "px";
            }
    
            document.addEventListener("DOMContentLoaded", function () {
                // Call the function initially to set the size
                resizeViewport();
                
                var stage = new NGL.Stage("viewport");
                """ + pdb_html + """
                """ + sdf_html + """
                // Function to handle both resizing of the viewport and NGL stage
                
                function handleResize() {
                    resizeViewport();  // Resize the #viewport container
                    stage.handleResize();  // Let NGL know it needs to adjust the view
                }

                // Listen to the window resize event and adjust viewport and stage
                window.addEventListener("resize", handleResize);
                })
        </script>
        </body></html>"""
        )
        

    return f"""<iframe style="width: 100%; height: 600px" name="result" sandbox="allow-modals allow-forms 
    allow-scripts allow-same-origin allow-popups 
    allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""

