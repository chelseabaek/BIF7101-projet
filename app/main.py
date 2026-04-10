from flask import Flask, render_template, request, Response, url_for, render_template, send_from_directory
from Bio import Phylo, AlignIO
from io import StringIO
from werkzeug.utils import secure_filename
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

import subprocess, threading
import tempfile
import os
import matplotlib
import matplotlib.pyplot as plt
import uuid
import shutil
matplotlib.use("Agg")

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)

# TOOL: Other Tools -> Newick Tree Viewer
@app.route("/", methods=["GET", "POST"])
def tree_viewer():
    tree_image = None
    tree_error = None
    active_tab = request.form.get("active_tab") or request.args.get("tab") or "main-page"

    if request.method == "POST":
        newick_str = request.form.get("newick")
        active_tab = "newick_tree_viewer"  # keep tab active
        
        width = float(request.form.get("width") or 6.4)
        height = float(request.form.get("height") or 4.8)

        if not newick_str or newick_str.strip() == "":
            tree_error = "Please provide a Newick string."
        else:
            try:
                handle = StringIO(newick_str)
                tree = Phylo.read(handle, "newick")
                fig = plt.figure(figsize=(width, height))
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree, axes=ax, do_show=False)

                filename = f"tree_{uuid.uuid4().hex}.png"
                image_path = os.path.join(UPLOAD_FOLDER, filename)
                plt.savefig(image_path, bbox_inches="tight")
                plt.close()

                tree_image = url_for("uploaded_file", filename=filename)
                print(os.path.exists(image_path))

            except Exception as e:
                tree_error = f"Failed to generate tree: {str(e)}"

    return render_template("index.html", tree_image=tree_image, tree_error=tree_error, active_tab=active_tab)
    
# TOOL: Other Tools -> Conversion
@app.route("/convert", methods=["POST"])
def convert():
    file = request.files.get("fasta_file")
    molecule_type = request.form.get("molecule_type")

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.")

    if not molecule_type:
        return render_template("index.html", error="Please select a molecule type.")

    try:
        fasta_text = file.read().decode("utf-8")
        fasta_io = StringIO(fasta_text)
        alignment = AlignIO.read(fasta_io, "fasta")

        for record in alignment:
            record.annotations["molecule_type"] = molecule_type

        nexus_io = StringIO()
        AlignIO.write(alignment, nexus_io, "nexus")

        result = nexus_io.getvalue()
        
        fasta_filename = file.filename.rsplit('.', 1)[0]
        nexus_filename = f"{fasta_filename}.nex"

        return Response(
            result,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={nexus_filename}"
            }
        )
        
    except Exception as e:
        return render_template(
            "index.html",
            error=f"Conversion failed: {str(e)}"
        )

# https://ena01.uqam.ca/pluginfile.php/9775398/mod_resource/content/2/analyse_phylogenetique_alignement_distances_arbres.html?embed=1
def parse_iqtree_distance_matrix(file_path):
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    n = int(lines[0])
    names = []
    matrix = []
    
    for i in range(1, n + 1):
        parts = lines[i].split()
        names.append(parts[0])
        distances = []
        for j in range(1, i + 1):
            distances.append(float(parts[j]))
        
        matrix.append(distances)
    
    return DistanceMatrix(names, matrix)

# TOOL: Sequence Alignment -> MUSCLE 
@app.route("/align_muscle", methods=["POST"])
def align_muscle():
    file = request.files.get("fasta_file")
    mode = request.form.get("computation_mode", "align")
    active_tab = "muscle"

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    muscle_exe = "muscle"

    try:
        fasta_text = file.read()

        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_text)
            temp_in_path = temp_in.name

        temp_out_path = temp_in_path + "_aligned.fasta"

        # MUSCLE v5 Linux Syntax
        cmd = [muscle_exe, f"-{mode}", temp_in_path, "-output", temp_out_path]

        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if process.returncode != 0:
            raise Exception(f"MUSCLE Error: {process.stderr}")
        
        # Move the aligned file to the uploads folder
        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_name = f"{safe_original_name}_muscle_{uuid.uuid4().hex[:6]}.fasta"
        download_path = os.path.join(UPLOAD_FOLDER, download_name)

        shutil.move(temp_out_path, download_path)
        muscle_download_file = url_for("uploaded_file", filename=download_name)
        
        # Cleanup
        os.remove(temp_in_path)
        if os.path.exists(temp_out_path):
            os.remove(temp_out_path)

        return render_template("index.html", active_tab=active_tab, muscle_download_file=muscle_download_file)

    except Exception as e:
        if 'temp_in_path' in locals() and os.path.exists(temp_in_path): os.remove(temp_in_path)
        if 'temp_out_path' in locals() and os.path.exists(temp_out_path): os.remove(temp_out_path)
        return render_template("index.html", error=str(e), active_tab=active_tab)
    
# TOOL: Sequence Alignment -> MAFFT
@app.route("/align_mafft", methods=["POST"])
def align_mafft():
    file = request.files.get("fasta_file")
    
    # Capture the new strategy from the form
    strategy = request.form.get("mafft_strategy", "--auto")
    active_tab = "mafft"

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    # Path to the MAFFT binary in your /app folder
    # mafft_exe = os.path.join(os.getcwd(), "mafft")
    
    try:
        fasta_content = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_content)
            temp_in_path = temp_in.name
        
        temp_out_path = temp_in_path + "_mafft_aligned.fasta"

        if strategy == "--linsi":
            # L-INS-i: Local pair alignment (most accurate)
            cmd = ["mafft", "--localpair", "--maxiterate", "1000", temp_in_path]
        elif strategy == "--einsi":
            # E-INS-i: Generalized affine gap costs
            cmd = ["mafft", "--genafpair", "--ep", "0", "--maxiterate", "1000", temp_in_path]
        elif strategy == "--ginsi":
            # G-INS-i: Global pair alignment
            cmd = ["mafft", "--globalpair", "--maxiterate", "1000", temp_in_path]
        else:
            # Fallback to standard auto mode
            cmd = ["mafft", "--auto", temp_in_path]

        # Building command with the selected strategy
        # cmd = ["mafft", "--auto", temp_in_path]
        
        with open(temp_out_path, "w") as out_file:
            process = subprocess.run(cmd, 
                stdout=out_file, 
                stderr=subprocess.PIPE, 
                text=True)

        if process.returncode != 0:
            raise Exception(f"MAFFT Error: {process.stderr}")

        #with open(temp_out_path, "r") as f:
        #    aligned_result = f.read()

        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_name = f"{safe_original_name}_mafft_{uuid.uuid4().hex[:6]}.fasta"
        download_path = os.path.join(UPLOAD_FOLDER, download_name)

        shutil.move(temp_out_path, download_path)
        mafft_download_file = url_for("uploaded_file", filename=download_name)

        # Cleanup temporary files
        os.remove(temp_in_path)
        if os.path.exists(temp_out_path):
            os.remove(temp_out_path)           
        '''
        original_name = file.filename.rsplit('.', 1)[0]
        return Response(
            aligned_result,
            mimetype="text/plain",
            headers={"Content-Disposition": f"attachment; filename={original_name}_mafft.fasta"}
        )
        '''
        return render_template("index.html", active_tab=active_tab, mafft_download_file=mafft_download_file)

    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)

# TOOL: Sequence Alignment -> Clustal Omega 
@app.route("/align_clustalo", methods=["POST"])
def align_clustalo():
    file = request.files.get("fasta_file")
    active_tab = "clustalo"

    # Capture options from the form
    outfmt = request.form.get("outfmt", "fasta")
    iterations = request.form.get("iterations", "0")
    outorder = request.form.get("outorder", "input-order")

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    try:
        fasta_content = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_content)
            temp_in_path = temp_in.name
        
        temp_out_path = temp_in_path + "_clustalo_aligned.txt"

        # Build the Clustal Omega Command dynamically
        cmd = [
            "clustalo", 
            "-i", temp_in_path, 
            "-o", temp_out_path, 
            "--outfmt", outfmt,
            "--iter", iterations,
            "--output-order", outorder,
            "--force"
        ]
        
        process = subprocess.run(cmd, capture_output=True, text=True)

        if process.returncode != 0:
            raise Exception(f"Clustal Omega Error: {process.stderr}")

        # Match the downloaded file extension to the requested format
        ext_map = {"fasta": "fasta", "clustal": "aln", "phylip": "phy"}
        file_extension = ext_map.get(outfmt, "txt")
        
        # Move the aligned file to the uploads folder
        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_name = f"{safe_original_name}_clustalo_{uuid.uuid4().hex[:6]}.{file_extension}"
        download_path = os.path.join(UPLOAD_FOLDER, download_name)

        shutil.move(temp_out_path, download_path)
        clustalo_download_file = url_for("uploaded_file", filename=download_name)

        #with open(temp_out_path, "r") as f:
        #    aligned_result = f.read()

        # Cleanup
        os.remove(temp_in_path)
        if os.path.exists(temp_out_path):
            os.remove(temp_out_path)
        '''
        # Match the downloaded file extension to the requested format
        ext_map = {"fasta": "fasta", "clustal": "aln", "phylip": "phy"}
        file_extension = ext_map.get(outfmt, "txt")
        original_name = file.filename.rsplit('.', 1)[0]
        
        return Response(
            aligned_result,
            mimetype="text/plain",
            headers={"Content-Disposition": f"attachment; filename={original_name}_clustalo.{file_extension}"}
        )
        '''
        return render_template("index.html", active_tab=active_tab, clustalo_download_file=clustalo_download_file)
        
    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)
'''
# TOOL: Tree Inference -> Parsimony (MPBoot)
@app.route("/run_mpboot", methods=["POST"])
def run_mpboot():
    file = request.files.get("fasta_file")
    active_tab = "parsimony"

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    mpboot_exe = os.path.join(os.getcwd(), "mpboot")

    try:
        fasta_text = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_text)
            temp_in_path = temp_in.name
        
        output_prefix = temp_in_path + "_mpb"

        # -s: input, -bb: 1000 bootstrap replicates, -pre: output prefix
        cmd = [mpboot_exe, "-s", temp_in_path, "-bb", "1000", "-pre", output_prefix]
        process = subprocess.run(cmd, capture_output=True, text=True)

        if process.returncode != 0:
            raise Exception(f"MPBoot Error: {process.stderr}")

        treefile_path = output_prefix + ".treefile"
        with open(treefile_path, "r") as f:
            tree_data = f.read()

        # Cleanup intermediate files
        for ext in [".log", ".treefile", ".iqtree", ".ckp.gz"]:
            if os.path.exists(output_prefix + ext):
                os.remove(output_prefix + ext)
        os.remove(temp_in_path)

        return Response(
            tree_data,
            mimetype="text/plain",
            headers={"Content-Disposition": f"attachment; filename={file.filename}_mpboot.tree"}
        )

    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)
'''
# TOOL: Tree Inference -> Parsimony (MPBoot)
@app.route("/run_mpboot", methods=["POST"])
def run_mpboot():
    file = request.files.get("fasta_file")
    active_tab = "parsimony"

    # 1. Capture the advanced parameters from the UI
    bootstraps = request.form.get("bootstraps", "1000")
    iterations = request.form.get("iterations", "1000")
    nstop = request.form.get("nstop", "100")

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    mpboot_exe = os.path.join(os.getcwd(), "mpboot") if os.path.exists(os.path.join(os.getcwd(), "mpboot")) else "mpboot"

    try:
        fasta_text = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_text)
            temp_in_path = temp_in.name
        
        output_prefix = temp_in_path + "_mpb"

        # 2. Build command dynamically
        cmd = [mpboot_exe, "-s", temp_in_path, "-pre", output_prefix]
        
        if bootstraps != "0":
            cmd.extend(["-bb", bootstraps])
            
        cmd.extend(["-n", iterations, "-nstop", nstop])
            
        process = subprocess.run(cmd, capture_output=True, text=True)

        if process.returncode != 0:
            raise Exception(f"MPBoot Error: {process.stderr}")

        treefile_path = output_prefix + ".treefile"
        
        if not os.path.exists(treefile_path):
            raise Exception("Treefile was not generated. Check if the FASTA alignment is valid.")

        # 3. READ AND DRAW THE TREE
        with open(treefile_path, "r") as f:
            newick_str = f.read().strip()

        tree = Phylo.read(StringIO(newick_str), "newick")
        num_leaves = tree.count_terminals()
        set_height = max(5, num_leaves * 1.0)
        
        fig = plt.figure(figsize=(10, set_height))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        
        # Save the tree image to the uploads folder
        image_name = f"tree_mpboot_{uuid.uuid4().hex[:6]}.png"
        image_path = os.path.join(UPLOAD_FOLDER, image_name)
        plt.savefig(image_path, bbox_inches="tight")
        plt.close()

        tree_image_mpboot = url_for("uploaded_file", filename=image_name)

        # 4. PREPARE THE TREEFILE FOR DOWNLOAD
        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_tree_name = f"{safe_original_name}_mpboot_{uuid.uuid4().hex[:6]}.tree"
        download_tree_path = os.path.join(UPLOAD_FOLDER, download_tree_name)
        
        shutil.move(treefile_path, download_tree_path)
        mpboot_tree_file = url_for("uploaded_file", filename=download_tree_name)

        # Cleanup intermediate files
        for ext in [".log", ".iqtree", ".ckp.gz"]:
            if os.path.exists(output_prefix + ext):
                os.remove(output_prefix + ext)
        os.remove(temp_in_path)

        return render_template(
            "index.html", 
            active_tab=active_tab, 
            tree_image_mpboot=tree_image_mpboot,
            mpboot_tree_file=mpboot_tree_file
        )

    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)

# TOOL: Tree inference -> Maximum Likelihood (IQ-TREE)
@app.route("/run_iqtree", methods=["POST"])
def run_iqtree():
    file = request.files.get("fasta_file")
    active_tab = "ml"

    # 1. Capture parameters from the UI
    model = request.form.get("model", "MFP")
    bootstraps = request.form.get("bootstraps", "1000")
    alrt = request.form.get("alrt", "1000")

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    # Use the iqtree3 binary installed via Docker
    iqtree_exe = "iqtree3"

    try:
        fasta_text = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_text)
            temp_in_path = temp_in.name
        
        output_prefix = temp_in_path + "_iq"

        # 2. Build command dynamically
        cmd = [iqtree_exe, "-s", temp_in_path, "-m", model, "-nt", "AUTO", "-pre", output_prefix]
        
        if bootstraps != "0":
            cmd.extend(["-bb", bootstraps])
        if alrt != "0":
            cmd.extend(["-alrt", alrt])
            
        process = subprocess.run(cmd, capture_output=True, text=True)

        if process.returncode != 0:
            raise Exception(f"IQ-TREE Error: {process.stderr}")

        treefile_path = output_prefix + ".treefile"
        
        if not os.path.exists(treefile_path):
            raise Exception("Treefile was not generated. Check if the FASTA alignment is valid.")

        # 3. READ AND DRAW THE TREE
        with open(treefile_path, "r") as f:
            newick_str = f.read().strip()

        tree = Phylo.read(StringIO(newick_str), "newick")
        num_leaves = tree.count_terminals()
        set_height = max(5, num_leaves * 1.0)
        
        fig = plt.figure(figsize=(10, set_height))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        
        # Save image to uploads folder
        image_name = f"tree_iqtree_{uuid.uuid4().hex[:6]}.png"
        image_path = os.path.join(UPLOAD_FOLDER, image_name)
        plt.savefig(image_path, bbox_inches="tight")
        plt.close()

        tree_image_iqtree = url_for("uploaded_file", filename=image_name)

        # 4. PREPARE DOWNLOAD
        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_tree_name = f"{safe_original_name}_iqtree_{uuid.uuid4().hex[:6]}.tree"
        download_tree_path = os.path.join(UPLOAD_FOLDER, download_tree_name)
        
        shutil.move(treefile_path, download_tree_path)
        iqtree_download_file = url_for("uploaded_file", filename=download_tree_name)

        # Cleanup intermediate files (IQ-TREE makes a lot of them!)
        extensions = [".log", ".treefile", ".iqtree", ".ckp.gz", ".bionj", ".mldist", ".model.gz", ".splits.nex", ".contree"]
        for ext in extensions:
            if os.path.exists(output_prefix + ext):
                os.remove(output_prefix + ext)
        os.remove(temp_in_path)

        return render_template(
            "index.html", 
            active_tab=active_tab, 
            tree_image_iqtree=tree_image_iqtree,
            iqtree_download_file=iqtree_download_file
        )

    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)
'''
# TOOL: Tree inference -> Maximum Liklihoof (IQ-TREE)
@app.route("/run_iqtree", methods=["POST"])
def run_iqtree():
    file = request.files.get("alignment_file")
    active_tab = "ml" # Maximum Likelihood tab

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload an alignment file.", active_tab=active_tab)

    iqtree_exe = os.path.join(os.getcwd(), "iqtree")

    try:
        content = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(content)
            temp_in_path = temp_in.name
        
        output_prefix = temp_in_path + "_iq"

        # -s: input, -m MFP: ModelFinder Plus (finds best substitution model)
        # -bb: 1000 ultrafast bootstraps, -nt AUTO: use optimal threads
        cmd = [iqtree_exe, "-s", temp_in_path, "-m", "MFP", "-bb", "1000", "-nt", "AUTO", "-pre", output_prefix]
        
        process = subprocess.run(cmd, capture_output=True, text=True)

        if process.returncode != 0:
            raise Exception(f"IQ-TREE Error: {process.stderr}")

        # The ML tree with bootstrap supports is in the .treefile
        treefile_path = output_prefix + ".treefile"
        with open(treefile_path, "r") as f:
            tree_data = f.read()

        # Cleanup (IQ-TREE generates MANY files)
        extensions = [".log", ".treefile", ".iqtree", ".ckp.gz", ".bionj", ".mldist", ".model.gz", ".splits.nex"]
        for ext in extensions:
            if os.path.exists(output_prefix + ext):
                os.remove(output_prefix + ext)
        os.remove(temp_in_path)

        # Return as JSON so your new visualizer can print it on screen
        return jsonify({"newick": tree_data, "method": "Maximum Likelihood"})

    except Exception as e:
        return jsonify({"error": str(e)})
'''    
# TOOL: Tree Inference -> Distance Methods (Neighbour-Joining and UPGMA)
@app.route("/distance-methods", methods=["POST"])
def NJ():
    active_tab = "distance"
    file = request.files.get("fasta_file_distance")
    method = request.form.get("method")
    
    if not file or file.filename == "":
        return render_template(
            "index.html",
            error="Please upload a file.",
            active_tab="distance"
        )
    
    try:
        unique_folder = os.path.join(UPLOAD_FOLDER, str(uuid.uuid4()))
        os.makedirs(unique_folder, exist_ok=True)

        safe_filename = secure_filename(file.filename)
        filepath = os.path.join(unique_folder, safe_filename)

        file.save(filepath)
        print("Saved file to:", filepath)

        # run IQ-TREE to find best model
        cmd = [
            "iqtree3",
            "-s", filepath,
            "-m", "MFP",
            "-nt", "1",
            "-redo"
        ]     
        # print("Running command:", " ".join(cmd))

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        # print("STDOUT:", result.stdout)
        # print("STDERR:", result.stderr)
        
        # read .iqtree file to get the best model
        iqtree_file = filepath + ".iqtree"
        best_model = None

        if os.path.exists(iqtree_file):
            with open(iqtree_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if "Best-fit model according to BIC:" in line:
                        best_model = line.split(":")[1].strip()
                        break
        else :
            return render_template(
                    "index.html",
                    error="Could not find best-fit model",
                    active_tab="distance"
                )
        
        prefix = "best_model_tree"
        
        model_cmd = [
            "iqtree3",
            "-s", filepath,
            "-m", best_model,
            "-nt", "1",
            "-redo", 
            "-pre", os.path.join(os.path.dirname(filepath), prefix)
        ]
        
        result2 = subprocess.run(
            model_cmd,
            capture_output=True,
            text=True
        )
        
        # print("STDOUT:", result2.stdout)
        # print("STDERR:", result2.stderr)
        
        mldist_file = prefix + ".mldist"
        filepath_mldist_file = os.path.join(unique_folder, mldist_file)
        
        if not os.path.exists(filepath_mldist_file):
            return render_template(
                "index.html",
                error="Distance matrix (.mldist) not found.",
                active_tab="distance"
            )

        dm = parse_iqtree_distance_matrix(filepath_mldist_file)

        constructor = DistanceTreeConstructor()
        if method == "nj":
            tree_distance = constructor.nj(dm)
            image_name = "tree_nj.png"
        elif method == "upgma":
            tree_distance = constructor.upgma(dm)
            image_name = "tree_upgma.png"
        else:
            raise ValueError("Invalid method selected")

        nj_tree_file = os.path.join(unique_folder, "nj_tree.nwk")
        Phylo.write(tree_distance, nj_tree_file, "newick")
        
        with open(nj_tree_file, "r") as f:
            newick_str = f.read().strip()
            
        # TODO: Add lengths, rename root and nodes
        tree = Phylo.read(StringIO(newick_str), "newick")
        num_leaves = tree.count_terminals()
        set_height = max(5, num_leaves * 1.0)
        fig = plt.figure(figsize=(10, set_height))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        
        # TODO: eventually to be put into unique_folder
        image_path = os.path.join("static", image_name)
        plt.savefig(image_path, bbox_inches="tight")
        plt.close()

        tree_image_distance = url_for("static", filename=image_name)

    except Exception as e:
        return render_template("index.html", distance_error=str(e), active_tab="distance")
        
    return render_template("index.html", tree_image_distance=tree_image_distance, active_tab=active_tab)
        
# TOOL: Tree Inference -> Bayesian Inference (MrBayes)
def run_mrbayes(nexus_path, working_dir):
    try:
        print("MrBayes start")
        result = subprocess.run(
            ["mb", os.path.basename(nexus_path)],
            cwd=working_dir,
            capture_output=True,
            text=True
        )
        print("MrBayes end")

        log_path = os.path.join(working_dir, "mrbayes_log.txt")
        with open(log_path, "w") as f:
            f.write(result.stdout + "\n\n" + result.stderr)

    except Exception as e:
        print("MrBayes error:", str(e))
        
        
@app.route("/bayesian_inference", methods=["POST"])
def bayesian_inference():
    file = request.files.get("sequence_file")
    
    print("File object:", file)
    print("Filename:", file.filename if file else "No file")

    print("ngen:", request.form.get("ngen"))
    print("samplefreq:", request.form.get("samplefreq"))
    print("printfreq:", request.form.get("printfreq"))
    print("burnin:", request.form.get("burnin"))

    if not file or file.filename == "":
        return render_template(
            "index.html",
            error="Please upload a file.",
            active_tab="bayesian-inference"
        )

    try:
        filename = file.filename.lower()
        file_text = file.read().decode("utf-8")

        # Default MrBayes parameters
        ngen = request.form.get("ngen", 1000000)
        samplefreq = request.form.get("samplefreq", 100)
        printfreq = request.form.get("printfreq", 100000)
        burnin = request.form.get("burnin", 2500)

        if filename.endswith(".nex"):
            nexus_text = file_text
        else:
            return render_template(
                "index.html",
                error="Unsupported file type.",
                active_tab="bayesian-inference"
            )

        # Add MrBayes block
        nexus_text += "\nBEGIN MRBAYES;\n"
        nexus_text += "  set autoclose=yes nowarn=yes;\n"
        nexus_text += "  lset nst=6 rates=invgamma;\n"
        nexus_text += f"  mcmc ngen={ngen} samplefreq={samplefreq} printfreq={printfreq} nchains=4 nruns=2;\n"
        nexus_text += f"  sumt burnin={burnin};\n"
        nexus_text += f"  sump burnin={burnin};\n"
        nexus_text += "END;\n"

        # Save file
        # nexus_filename = "static/mrbayes_input.nex"
        
        unique_folder = os.path.join(UPLOAD_FOLDER, str(uuid.uuid4()))
        os.makedirs(unique_folder, exist_ok=True)

        nexus_filename = os.path.join(unique_folder, "mrbayes.nex")
        
        with open(nexus_filename, "w") as f:
            f.write(nexus_text)

        # Run MrBayes in background
        threading.Thread(target=run_mrbayes, args=(nexus_filename, unique_folder)).start()
        
        # run_mrbayes(nexus_filename, unique_folder)

        return render_template(
            "index.html",
            message="MrBayes is running... refresh in a few seconds.",
            active_tab="bayesian-inference"
        )

    except Exception as e:
        return render_template(
            "index.html",
            error=f"Error: {str(e)}",
            active_tab="bayesian-inference"
        )

# if __name__ == "__main__":
#     app.run(debug=True)
