from flask import Flask, render_template, request, Response, url_for, send_from_directory
from Bio import Phylo, AlignIO, SeqIO
from io import StringIO
from werkzeug.utils import secure_filename
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

import subprocess
import tempfile
import os
import glob
import re
import uuid
import shutil

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)

# TOOL: Other Tools -> Newick Tree Viewer
@app.route("/", methods=["GET", "POST"])
def tree_viewer():
    tree_error = None
    tree_newick = None
    newick_input = ""
    newick_branch_color = "#111827"
    newick_leaf_color = "#111827"
    newick_branch_length_color = "#111827"
    newick_internal_node_color = "#111827"
    active_tab = request.form.get("active_tab") or request.args.get("tab") or "main-page"

    if request.method == "POST":
        newick_str = (request.form.get("newick", "") or "").strip()
        newick_file = request.files.get("newick_file")
        active_tab = "newick_tree_viewer"  # keep tab active
        newick_branch_color = request.form.get("newick_branch_color", newick_branch_color)
        newick_leaf_color = request.form.get("newick_leaf_color", newick_leaf_color)
        newick_branch_length_color = request.form.get("newick_branch_length_color", newick_branch_length_color)
        newick_internal_node_color = request.form.get("newick_internal_node_color", newick_internal_node_color)

        # Prefer textarea input; fallback to uploaded text tree file.
        if not newick_str and newick_file and newick_file.filename:
            try:
                newick_str = newick_file.read().decode("utf-8", errors="ignore").strip()
            except Exception:
                newick_str = ""

        newick_input = newick_str

        if not newick_str or newick_str.strip() == "":
            tree_error = "Please provide a Newick string or upload a tree file (.nwk/.tree/.tre/.txt)."
        else:
            try:
                handle = StringIO(newick_str)
                tree = Phylo.read(handle, "newick")
                # Normalize parsed Newick so the frontend viewer gets a clean tree string.
                # normalized_newick = StringIO()
                # Phylo.write(tree, normalized_newick, "newick")
                # tree_newick = normalized_newick.getvalue().strip()
                
                # THE FIX: Check if tree is topology-only. 
                # If so, inject 1.0 lengths so WebGL doesn't divide by zero!
                total_len = sum(c.branch_length for c in tree.find_clades() if c.branch_length)
                if total_len == 0:
                    for clade in tree.find_clades():
                        clade.branch_length = 1.0
                        
                safe_io = StringIO()
                Phylo.write(tree, safe_io, "newick")
                tree_newick = safe_io.getvalue().strip()

            except Exception as e:
                tree_error = f"Failed to generate tree: {str(e)}"

    return render_template(
        "index.html",
        tree_error=tree_error,
        tree_newick=tree_newick,
        newick_input=newick_input,
        newick_branch_color=newick_branch_color,
        newick_leaf_color=newick_leaf_color,
        newick_branch_length_color=newick_branch_length_color,
        newick_internal_node_color=newick_internal_node_color,
        active_tab=active_tab,
    )
    
# TOOL: Other Tools -> Conversion
@app.route("/convert", methods=["POST"])
def convert():
    file = request.files.get("sequence_file")
    molecule_type = request.form.get("molecule_type")
    input_format = request.form.get("input_format", "fasta")
    output_format = request.form.get("output_format", "nexus")

    # Curated SeqIO formats with robust read/write support for this web workflow.
    conversion_formats = {
        "fasta": "FASTA",
        "fasta-2line": "FASTA (2-line)",
        "genbank": "GenBank",
        "embl": "EMBL",
        "nexus": "NEXUS",
        "phylip": "PHYLIP",
        "phylip-relaxed": "PHYLIP Relaxed",
        "stockholm": "Stockholm",
        "tab": "Tab-separated",
    }

    output_extensions = {
        "fasta": "fasta",
        "fasta-2line": "fasta",
        "genbank": "gb",
        "embl": "embl",
        "nexus": "nex",
        "phylip": "phy",
        "phylip-relaxed": "phy",
        "stockholm": "sto",
        "tab": "tab",
    }

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a sequence file.", active_tab="conversion")

    if input_format not in conversion_formats or output_format not in conversion_formats:
        return render_template("index.html", error="Unsupported conversion format selected.", active_tab="conversion")

    if output_format in {"nexus", "genbank", "embl"} and not molecule_type:
        return render_template("index.html", error="Please select a molecule type for this output format.", active_tab="conversion")

    try:
        sequence_text = file.read().decode("utf-8")
        input_io = StringIO(sequence_text)
        records = list(SeqIO.parse(input_io, input_format))

        if not records:
            return render_template("index.html", error="No records found. Check file content and selected input format.", active_tab="conversion")

        if molecule_type:
            for record in records:
                record.annotations["molecule_type"] = molecule_type

        output_io = StringIO()
        converted_count = SeqIO.write(records, output_io, output_format)

        if converted_count == 0:
            return render_template("index.html", error="Conversion produced no output records.", active_tab="conversion")

        result = output_io.getvalue()
        
        source_filename = file.filename.rsplit('.', 1)[0]
        out_ext = output_extensions.get(output_format, "txt")
        converted_filename = f"{source_filename}.{out_ext}"

        return Response(
            result,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={converted_filename}"
            }
        )
        
    except Exception as e:
        return render_template(
            "index.html",
            error=f"Conversion failed: {str(e)}",
            active_tab="conversion"
        )

# https://ena01.uqam.ca/pluginfile.php/9775398/mod_resource/content/2/analyse_phylogenetique_alignement_distances_arbres.html?embed=1
def parse_iqtree_distance_matrix(file_path):
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValueError("Distance matrix file is empty.")

    n = int(lines[0])
    names = []
    matrix = []

    if len(lines) < n + 1:
        raise ValueError("Distance matrix file is incomplete.")

    for i in range(1, n + 1):
        parts = lines[i].split()
        if len(parts) < i + 1:
            raise ValueError(f"Malformed distance matrix row for taxon: {parts[0] if parts else 'unknown'}")
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

        with open(temp_out_path, "w") as out_file:
            process = subprocess.run(cmd, 
                stdout=out_file, 
                stderr=subprocess.PIPE, 
                text=True)

        if process.returncode != 0:
            raise Exception(f"MAFFT Error: {process.stderr}")

        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_name = f"{safe_original_name}_mafft_{uuid.uuid4().hex[:6]}.fasta"
        download_path = os.path.join(UPLOAD_FOLDER, download_name)

        shutil.move(temp_out_path, download_path)
        mafft_download_file = url_for("uploaded_file", filename=download_name)

        # Cleanup temporary files
        os.remove(temp_in_path)
        if os.path.exists(temp_out_path):
            os.remove(temp_out_path)
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

        # Cleanup
        os.remove(temp_in_path)
        if os.path.exists(temp_out_path):
            os.remove(temp_out_path)
        return render_template("index.html", active_tab=active_tab, clustalo_download_file=clustalo_download_file)
        
    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)

# TOOL: Tree Inference -> Parsimony (MPBoot)
@app.route("/run_mpboot", methods=["POST"])
def run_mpboot():
    active_tab = "parsimony"
    file = request.files.get("fasta_file")
    
    # 1. Capture Advanced Parameters
    bootstraps = request.form.get("bootstraps", "1000")
    spr_rad = request.form.get("spr_rad", "6")
    ratchet_iter = request.form.get("ratchet_iter", "1")
    ratchet_percent = request.form.get("ratchet_percent", "50")
    nni_pars = request.form.get("nni_pars")
    mulhits = request.form.get("mulhits")

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.", active_tab=active_tab)

    mpboot_exe = "mpboot" # Assumes it is in your PATH or local bin as per Dockerfile [cite: 1, 2]

    try:
        fasta_text = file.read()
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_in:
            temp_in.write(fasta_text)
            temp_in_path = temp_in.name
        
        output_prefix = temp_in_path + "_mpb"

        # 2. Build the command dynamically
        cmd = [mpboot_exe, "-s", temp_in_path, "-pre", output_prefix]
        if bootstraps != "0":
            cmd.extend(["-bb", bootstraps])
        if spr_rad and spr_rad != "6":
            cmd.extend(["-spr_rad", spr_rad])
        if ratchet_iter and ratchet_iter != "1":
            cmd.extend(["-ratchet_iter", ratchet_iter])
        if ratchet_percent and ratchet_percent != "50":
            cmd.extend(["-ratchet_percent", ratchet_percent])
        if nni_pars:
            cmd.append("-nni_pars")
        if mulhits:
            cmd.append("-mulhits")

        # Execute MPBoot
        subprocess.run(cmd, capture_output=True, text=True, check=True)

        # 3. PREPARE THE VIEWER (Prioritizing the Consensus Tree)
        # The .contree file contains the bootstrap values seen in your successful tests
        newick_str = ""
        contree_path = output_prefix + ".contree"
        treefile_path = output_prefix + ".treefile"
        
        source_path = contree_path if os.path.exists(contree_path) else treefile_path
        
        if os.path.exists(source_path):
            tree = Phylo.read(source_path, "newick")
            
            # THE TOPOLOGY FIX: Inject 1.0 lengths if total length is 0
            # This prevents the 'invisible tree' WebGL error
            total_len = sum(c.branch_length for c in tree.find_clades() if c.branch_length)
            if total_len == 0:
                for clade in tree.find_clades():
                    clade.branch_length = 1.0
            
            safe_io = StringIO()
            Phylo.write(tree, safe_io, "newick")
            newick_str = safe_io.getvalue().strip()

        # 4. PREPARE DOWNLOADS
        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        run_id = uuid.uuid4().hex[:6]

        # MPBoot Report File (Summary)
        mpboot_report_file = None
        report_path = output_prefix + ".mpboot"
        if os.path.exists(report_path):
            report_name = f"{safe_original_name}_summary_{run_id}.txt"
            shutil.move(report_path, os.path.join(UPLOAD_FOLDER, report_name))
            mpboot_report_file = url_for("uploaded_file", filename=report_name)
        
        # Best Treefile (Maximum Parsimony)
        mpboot_tree_file = None
        if os.path.exists(treefile_path):
            tree_name = f"{safe_original_name}_best_{run_id}.tree"
            shutil.move(treefile_path, os.path.join(UPLOAD_FOLDER, tree_name))
            mpboot_tree_file = url_for("uploaded_file", filename=tree_name)

        # Consensus Tree (With Bootstraps)
        mpboot_contree_file = None
        if os.path.exists(contree_path):
            con_name = f"{safe_original_name}_consensus_{run_id}.contree"
            shutil.move(contree_path, os.path.join(UPLOAD_FOLDER, con_name))
            mpboot_contree_file = url_for("uploaded_file", filename=con_name)

        # Run Log
        mpboot_log_file = None
        logfile_path = output_prefix + ".log"
        if os.path.exists(logfile_path):
            log_name = f"{safe_original_name}_log_{run_id}.log"
            shutil.move(logfile_path, os.path.join(UPLOAD_FOLDER, log_name))
            mpboot_log_file = url_for("uploaded_file", filename=log_name)

        # Cleanup
        for ext in [".iqtree", ".ckp.gz", ".splits.nex", ".mpboot"]:
            path = output_prefix + ext
            if os.path.exists(path):
                os.remove(path)
        if os.path.exists(temp_in_path):
            os.remove(temp_in_path)

        return render_template(
            "index.html", 
            active_tab=active_tab, 
            mpboot_newick=newick_str,
            mpboot_tree_file=mpboot_tree_file,
            mpboot_contree_file=mpboot_contree_file,
            mpboot_log_file=mpboot_log_file
            mpboot_report_file=mpboot_report_file
        )

    except Exception as e:
        return render_template("index.html", error=f"MPBoot Error: {str(e)}", active_tab=active_tab)

# TOOL: Tree inference -> Maximum Likelihood (IQ-TREE)
@app.route("/run_iqtree", methods=["POST"])
def run_iqtree():
    file = request.files.get("fasta_file")
    active_tab = "maximum-likelihood"

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

        # 3. READ THE RAW NEWICK STRING
        with open(treefile_path, "r") as f:
            newick_str = f.read().strip()

        # 4. EXTRACT THE BEST-FIT MODEL
        best_model = model  # Default to whatever the user typed
        iqtree_report = output_prefix + ".iqtree"
        
        # If the report exists, search it for the winning model
        if os.path.exists(iqtree_report):
            with open(iqtree_report, 'r', encoding='utf-8') as f:
                for line in f:
                    if "Best-fit model according to BIC:" in line:
                        # Grabs the model name (e.g. "GTR+F+I+G4")
                        best_model = line.split(":")[1].strip()
                        break
                    elif "Best-fit model:" in line:
                        # Fallback for some versions of IQ-TREE
                        best_model = line.split("Best-fit model:")[1].split()[0]

        # 5. PREPARE DOWNLOAD
        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        download_tree_name = f"{safe_original_name}_iqtree_{uuid.uuid4().hex[:6]}.tree"
        download_tree_path = os.path.join(UPLOAD_FOLDER, download_tree_name)
        
        shutil.move(treefile_path, download_tree_path)
        iqtree_download_file = url_for("uploaded_file", filename=download_tree_name)

        # Cleanup intermediate files
        extensions = [".log", ".treefile", ".iqtree", ".ckp.gz", ".bionj", ".mldist", ".model.gz", ".splits.nex", ".contree"]
        for ext in extensions:
            if os.path.exists(output_prefix + ext):
                os.remove(output_prefix + ext)
        os.remove(temp_in_path)

        # PASS THE STRING AND THE MODEL TO THE TEMPLATE
        return render_template(
            "index.html", 
            active_tab=active_tab, 
            iqtree_newick=newick_str,
            iqtree_download_file=iqtree_download_file,
            iqtree_best_model=best_model # <--- Pass the new variable here!
        )

    except Exception as e:
        return render_template("index.html", error=str(e), active_tab=active_tab)
# TOOL: Tree Inference -> Distance Methods (Neighbour-Joining and UPGMA)
@app.route("/distance-methods", methods=["POST"])
def distance_methods():
    active_tab = "distance"
    file = request.files.get("fasta_file_distance")
    method = (request.form.get("method") or "").strip().lower()

    if not file or file.filename == "":
        return render_template(
            "index.html",
            error="Please upload a file.",
            active_tab="distance"
        )

    if method not in {"nj", "upgma"}:
        return render_template(
            "index.html",
            distance_error="Please choose a valid distance method (NJ or UPGMA).",
            active_tab=active_tab,
        )

    unique_folder = None
    try:
        unique_folder = os.path.join(UPLOAD_FOLDER, str(uuid.uuid4()))
        os.makedirs(unique_folder, exist_ok=True)

        safe_filename = secure_filename(file.filename)
        filepath = os.path.join(unique_folder, safe_filename)

        file.save(filepath)
        output_prefix = os.path.join(unique_folder, "distance_model")

        # One IQ-TREE run is enough: MFP selects the model and writes the ML distance matrix.
        cmd = ["iqtree3", "-s", filepath, "-m", "MFP", "-nt", "AUTO", "-pre", output_prefix, "-redo"]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "IQ-TREE failed while selecting a model.")

        iqtree_file = output_prefix + ".iqtree"
        best_model = None

        if os.path.exists(iqtree_file):
            with open(iqtree_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if "Best-fit model according to BIC:" in line:
                        best_model = line.split(":")[1].strip()
                        break
                    elif "Best-fit model:" in line:
                        best_model = line.split("Best-fit model:")[1].split()[0]
                        break

        if not best_model:
            return render_template(
                    "index.html",
                    error="Could not find best-fit model",
                    active_tab="distance"
                )

        filepath_mldist_file = output_prefix + ".mldist"

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
            method_label = "Neighbor-Joining"
        else:
            tree_distance = constructor.upgma(dm)
            method_label = "UPGMA"

        output_tree_file = os.path.join(unique_folder, f"{method}_tree.nwk")
        Phylo.write(tree_distance, output_tree_file, "newick")

        with open(output_tree_file, "r", encoding="utf-8") as f:
            newick_str = f.read().strip()

    except Exception as e:
        return render_template("index.html", distance_error=str(e), active_tab="distance")
    finally:
        if unique_folder and os.path.isdir(unique_folder):
            shutil.rmtree(unique_folder, ignore_errors=True)
        
    return render_template(
        "index.html",
        distance_newick=newick_str,
        distance_method_label=method_label,
        distance_substitution_model=best_model,
        active_tab=active_tab,
    )
        
# TOOL: Tree Inference -> Bayesian Inference (MrBayes)
def run_mrbayes(nexus_path, working_dir):
    try:
        # print("MrBayes start")
        result = subprocess.run(
            ["mb", os.path.basename(nexus_path)],
            cwd=working_dir,
            capture_output=True,
            text=True
        )
        # print("MrBayes end")

        log_path = os.path.join(working_dir, "mrbayes_log.txt")
        with open(log_path, "w") as f:
            f.write(result.stdout + "\n\n" + result.stderr)

        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "MrBayes exited with a non-zero status.")

        return log_path

    except Exception as e:
        print("MrBayes error:", str(e))
        raise


def find_mrbayes_consensus_file(working_dir):
    patterns = [
        "*.con.tre",
        "*.con",
        "*.t",
        "*.trprobs"
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.join(working_dir, pattern)))
        if matches:
            return matches[0]
    return None


def read_consensus_tree(consensus_path):
    for fmt in ("nexus", "newick"):
        try:
            trees = list(Phylo.parse(consensus_path, fmt))
            if trees:
                return trees[0]
        except Exception:
            continue
    raise ValueError("Could not parse MrBayes consensus tree file.")


def extract_clean_newick_from_mrbayes_text(raw_text):
    # Parse optional NEXUS translate block so numeric taxa IDs become real labels.
    translation = {}
    translate_match = re.search(r"translate\s+(.*?);", raw_text, flags=re.IGNORECASE | re.DOTALL)
    if translate_match:
        translate_block = translate_match.group(1)
        for key, value in re.findall(r"(\d+)\s+([^,;\n\r]+)", translate_block):
            translation[key] = value.strip().strip("'").strip('"')

    # Extract the tree assignment (everything after '=' up to ';').
    tree_match = re.search(r"tree\s+[^=]*=\s*(?:\[&U\]\s*)?(\(.*?\));", raw_text, flags=re.IGNORECASE | re.DOTALL)
    if not tree_match:
        return None

    newick = tree_match.group(1) + ";"

    # Remove all MrBayes annotation blocks, e.g. [&prob=...].
    newick = re.sub(r"\[&[^\]]*\]", "", newick)

    # Replace translated numeric taxon IDs with labels.
    if translation:
        pattern = re.compile(r"(?<=[(,])\s*(\d+)\s*(?=[:),])")

        def _replace_taxon(match):
            taxon_id = match.group(1)
            return translation.get(taxon_id, taxon_id)

        newick = pattern.sub(_replace_taxon, newick)

    # Normalize repeated whitespace.
    newick = re.sub(r"\s+", "", newick)
    return newick


def _guess_molecule_type_from_alignment(alignment):
    dna_chars = set("ACGTN-?")
    rna_chars = set("ACGUN-?")
    observed = set()

    for record in alignment:
        observed.update(str(record.seq).upper())

    if observed and observed.issubset(rna_chars):
        return "RNA"
    if observed and observed.issubset(dna_chars):
        return "DNA"
    return "protein"
        
        
@app.route("/bayesian_inference", methods=["POST"])
def bayesian_inference():
    file = request.files.get("sequence_file")
    
    # print("File object:", file)
    # print("Filename:", file.filename if file else "No file")

    # print("ngen:", request.form.get("ngen"))
    # print("samplefreq:", request.form.get("samplefreq"))
    # print("printfreq:", request.form.get("printfreq"))
    # print("burnin:", request.form.get("burnin"))

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
        elif filename.endswith((".fasta", ".fa", ".txt")):
            try:
                fasta_io = StringIO(file_text)
                alignment = AlignIO.read(fasta_io, "fasta")
                molecule_type = _guess_molecule_type_from_alignment(alignment)

                for record in alignment:
                    record.annotations["molecule_type"] = molecule_type

                nexus_io = StringIO()
                AlignIO.write(alignment, nexus_io, "nexus")
                nexus_text = nexus_io.getvalue()
            except Exception as convert_err:
                return render_template(
                    "index.html",
                    error=f"Could not convert FASTA to NEXUS: {str(convert_err)}",
                    active_tab="bayesian-inference"
                )
        else:
            return render_template(
                "index.html",
                error="Unsupported file type. Please upload a .nex, .fasta, .fa, or .txt FASTA file.",
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

        # Run MrBayes and wait until all outputs are generated.
        run_mrbayes(nexus_filename, unique_folder)

        safe_original_name = secure_filename(file.filename).rsplit('.', 1)[0]
        run_suffix = uuid.uuid4().hex[:6]

        # Package all MrBayes outputs from this run for automatic download.
        archive_basename = f"{safe_original_name}_mrbayes_{run_suffix}"
        archive_path_no_ext = os.path.join(UPLOAD_FOLDER, archive_basename)
        archive_path = shutil.make_archive(archive_path_no_ext, "zip", root_dir=unique_folder)
        bayes_results_zip = url_for("uploaded_file", filename=os.path.basename(archive_path))

        bayes_download_file = None
        bayes_newick = None
        consensus_path = find_mrbayes_consensus_file(unique_folder)
        if consensus_path:
            tree = read_consensus_tree(consensus_path)
            with open(consensus_path, "r", encoding="utf-8", errors="ignore") as f:
                consensus_text = f.read()
            bayes_newick = extract_clean_newick_from_mrbayes_text(consensus_text)
            if not bayes_newick:
                newick_io = StringIO()
                Phylo.write(tree, newick_io, "newick")
                bayes_newick = newick_io.getvalue().strip()

            treefile_name = f"{safe_original_name}_mrbayes_{run_suffix}.tree"
            treefile_path = os.path.join(UPLOAD_FOLDER, treefile_name)
            with open(treefile_path, "w", encoding="utf-8") as treefile_out:
                treefile_out.write(bayes_newick if bayes_newick.endswith(";") else bayes_newick + ";")
            bayes_download_file = url_for("uploaded_file", filename=treefile_name)

        shutil.rmtree(unique_folder, ignore_errors=True)

        return render_template(
            "index.html",
            active_tab="bayesian-inference",
            message="MrBayes completed.",
            bayes_newick=bayes_newick,
            bayes_download_file=bayes_download_file,
            auto_download_url=bayes_results_zip
        )

    except Exception as e:
        return render_template(
            "index.html",
            error=f"Error: {str(e)}",
            active_tab="bayesian-inference"
        )


@app.route("/view_mrbayes_tre", methods=["POST"])
def view_mrbayes_tre():
    file = request.files.get("mrbayes_tre_file")

    if not file or file.filename == "":
        return render_template(
            "index.html",
            error="Please upload a MrBayes .tre/.con.tre file.",
            active_tab="bayesian-inference"
        )

    try:
        suffix = os.path.splitext(file.filename)[1] or ".tre"
        raw_text = file.read().decode("utf-8", errors="ignore")
        with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as temp_file:
            temp_file.write(raw_text.encode("utf-8"))
            temp_path = temp_file.name

        tree = read_consensus_tree(temp_path)

        bayes_newick = extract_clean_newick_from_mrbayes_text(raw_text)
        if not bayes_newick:
            newick_io = StringIO()
            Phylo.write(tree, newick_io, "newick")
            bayes_newick = newick_io.getvalue().strip()

        uploaded_suffix = uuid.uuid4().hex[:6]
        treefile_name = f"mrbayes_uploaded_{uploaded_suffix}.tree"
        treefile_path = os.path.join(UPLOAD_FOLDER, treefile_name)
        with open(treefile_path, "w", encoding="utf-8") as treefile_out:
            treefile_out.write(bayes_newick if bayes_newick.endswith(";") else bayes_newick + ";")

        os.remove(temp_path)

        return render_template(
            "index.html",
            active_tab="bayesian-inference",
            message="Consensus tree loaded from uploaded file.",
            bayes_newick=bayes_newick,
            bayes_download_file=url_for("uploaded_file", filename=treefile_name)
        )

    except Exception as e:
        if 'temp_path' in locals() and os.path.exists(temp_path):
            os.remove(temp_path)
        return render_template(
            "index.html",
            error=f"Could not parse uploaded tree: {str(e)}",
            active_tab="bayesian-inference"
        )

# if __name__ == "__main__":
#     app.run(debug=True)
