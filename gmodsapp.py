#!/usr/bin/env python
"""
    Requirements: Python 3 or higher
                  Antechamber and related AmberTools
                  OpenBabel (strongly recommended for use with acpype)
				  acpype (latest version recommended with all its requirements)
				  Gromacs (Compulsory)
                  flask (Compulsory)
                  flaskwebgui (recommended)
                  pyfladesk (recommended)

    Daniyan, Michael Oluwatoyin, B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) Biochemistry
    Department of Pharmacology, Faculty of Pharmacy
    Obafemi Awolowo Universiy, Ile-Ife, Nigeria.
    >>http://www.oauife.edu.ng<<

    mdaniyan@oauife.edu.ng; toyinpharm@gmail.com
"""
import sys

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

import os
import subprocess
from pathlib import Path
import time
import shutil
import random
import string
import math

try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass

import gmodsScripts.gmodsRLmulti
import gmodsScripts.gmodsRLmany
import gmodsScripts.gmodsRLsingle
import gmodsScripts.gmodsSCmds
import gmodsScripts.gmodsPPmore
from gmodsScripts.gmodsHelpers import printWarning, printNote, cleanup

from flask import Flask, flash, request, redirect, url_for, render_template, g
from functools import wraps
from werkzeug.utils import secure_filename

appDIR = os.path.abspath(os.path.dirname(sys.argv[0]))
gmxDIR = Path.cwd()

if getattr(sys, 'frozen', False):
    template_folder = os.path.join(appDIR, 'templates')
    static_folder = os.path.join(appDIR, 'static')
    scripts_folder = os.path.join(appDIR, 'gmodsScripts')
    app = Flask(__name__, template_folder=template_folder, static_folder=static_folder, scripts_folder=scripts_folder)
else:
    app = Flask(__name__)

# Let's get templates to auto-reload
app.jinja_env.auto_reload = True
app.config["TEMPLATES_AUTO_RELOAD"] = True

try:
    from flaskwebgui import FlaskUI
    ui = FlaskUI(app, width=800, height=600, port=1380) # add app and parameters
except ImportError:
    printWarning("'flaskwebgui' is not installed or has errors, trying 'pyfladesk'.....")
    try:
        from pyfladesk import init_gui
        ui = init_gui
    except ImportError:
        printWarning("'pyfladesk' is not installed or has errors. GUI can not be launched")
        printNote("To access GUI interface, copy the generated web link to your default browser instead")
        time.sleep(10)
        ui = app

@app.errorhandler(Exception)
def all_exceptions(e):
    os.chdir(gmxDIR)
    return render_template("errors.html", e=e, N=N)

KEY = os.urandom(14)
app.config['SECRET_KEY'] = KEY

UPLOAD_FOLDER = os.path.join(appDIR, 'Uploads')
if not os.path.isdir(UPLOAD_FOLDER):
    os.mkdir(UPLOAD_FOLDER)
else:
    shutil.rmtree(UPLOAD_FOLDER)
    os.mkdir(UPLOAD_FOLDER)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ALLOWED_EXTENSIONS = {'pdb', 'mol2', 'gro', 'top', 'itp', 'mdp', 'h', 'cl', 'clh', 'in'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# Set some global varibales
global N, A, F, selected, ffselect
N = 0
A = 0
F = 0
selected = " "
ffselect = " "

def setup_route(f):
    """
    Decorate routes to require making choice of setup route to use.

    https://flask.palletsprojects.com/en/1.1.x/patterns/viewdecorators/
    """
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if N == 0:
            return redirect("/setroute")
        return f(*args, **kwargs)
    return decorated_function

@app.route('/')
def index():
    return render_template("index.html", N=N)

@app.route('/readme')
def readme():
    return render_template("readme.html", N=N)

@app.route("/setroute", methods=["GET", "POST"])
def setup_selection():
    """Determine the Setup Script to Use by selecting appropriate route"""
    global N
    global selected
    # User reached route via POST (as by submitting a form via POST)
    if request.method == "POST":
        # Remove all used subfolders present in upload directory
        listuploads = os.listdir(app.config['UPLOAD_FOLDER'])
        for item in listuploads:
            if not os.path.isdir(os.path.join(app.config['UPLOAD_FOLDER'], item)):
                os.remove(os.path.join(app.config['UPLOAD_FOLDER'], item))
            else:
                shutil.rmtree(os.path.join(app.config['UPLOAD_FOLDER'], item))

        # Ensure user choose from the dropdown menu
        if not request.form.get("damas"):
            flash('You must select a valid setup route to continue')
            return render_template("setroute.html", N=N, category='warning')

		# Store the submitted symbol as varibale
        selected = request.form.get("damas")

		# Check if selected route is valid
        if not (selected == "RLmulti" or selected == "RLmany" or selected == "RLsingle" or selected == "PPmore"):
            flash('You must select a valid setup route to continue')
            return render_template("setroute.html", N=N, category='warning')

		# If selection is valid, update the global value of N as follows
        if selected == "RLmulti":
            N = 1
        elif selected == "RLmany":
            N = 2
        elif selected == "RLsingle":
            N = 3
        elif selected == "PPmore":
            N = 4

        if N == 1 or N == 2 or N == 3:
            ligand_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands')
            if not os.path.isdir(ligand_dir):
                os.mkdir(ligand_dir)
            else:
                shutil.rmtree(ligand_dir)
                os.mkdir(ligand_dir)

            receptor_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors')
            if not os.path.isdir(receptor_dir):
                os.mkdir(receptor_dir)
            else:
                shutil.rmtree(receptor_dir)
                os.mkdir(receptor_dir)

        elif N == 4:
            peptide_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Peptides')
            if not os.path.isdir(peptide_dir):
                os.mkdir(peptide_dir)
            else:
                shutil.rmtree(peptide_dir)
                os.mkdir(peptide_dir)

            protein_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Proteins')
            if not os.path.isdir(protein_dir):
                os.mkdir(protein_dir)
            else:
                shutil.rmtree(protein_dir)
                os.mkdir(protein_dir)

        # Render the selectffgroup.html to select forcefields group
        return render_template("selectffgroup.html", N=N)

    # User reached route via GET (as by clicking a link or via redirect)
    else:
        return render_template("setroute.html", N=N)

@app.route("/selectffgroup", methods=["GET", "POST"])
def ff_selection():
    """Determine the group of forcefields preferred by the user"""
    """This will help to request additional files if needed for complex generation"""
    global N
    global F
    global selected
    global ffselect
	# User reached route via POST (as by submitting a form via POST)
    if request.method == "POST":
        # Cleanup the existing Ligsff subdirectory or file
        for item in os.listdir(app.config['UPLOAD_FOLDER']):
            if not (item == 'Ligands' or item == 'Receptors' or item == 'Peptides' or item == 'Proteins'):
                try:
                    shutil.rmtree(os.path.join(app.config['UPLOAD_FOLDER'], item))
                except:
                    os.remove(os.path.join(app.config['UPLOAD_FOLDER'], item))

        # Check if user choose a forcefield group from the dropdown menu or ignore
        if not request.form.get("gmodsff"):
            flash('You must select a forcefileds group')
            return redirect("/selectffgroup")

        if not N == 4:
            if not request.form.get("gmodsyn"):
                flash('You must select either YES or NO')
                return redirect("/selectffgroup")

        else:
            if not request.form.get("gmodsyn"):
                flash('You must select either Proteins or Peptides')
                return redirect("/selectffgroup")

		# Store the submitted symbol as varibale
        ffselect = request.form.get("gmodsff")
        ffupload = request.form.get("gmodsyn")

		# Check if selected route is valid
        if not (ffselect == "Amber" or ffselect == "Charmm" or ffselect == "Gromos" or ffselect == "Opls"):
            flash('You must select a forcefileds group')
            return redirect("/selectffgroup")

        if not (ffupload == "YES" or ffupload == "NO") and not N == 4:
            flash('You must select either YES or NO')
            return redirect("/selectffgroup")

        elif not (ffupload == "Proteins" or ffupload == "Peptides") and N == 4:
            flash('You must select either Proteins or Peptides for upload')
            return redirect("/selectffgroup")

		# If selection is valid, update the global value of F as follows
        if N == 1 or N == 2:
            if ffupload == "YES":
                F = 1
                if (ffselect == "Amber" or ffselect == "Charmm" or ffselect == "Gromos" or ffselect == "Opls"):
                    ligsff_dir = os.path.join(app.config['UPLOAD_FOLDER'], ffselect)
                    if not (os.path.isdir(ligsff_dir) or os.path.isfile(ligsff_dir)):
                        os.mkdir(ligsff_dir)
                    else:
                        try:
                            shutil.rmtree(ligsff_dir)
                        except:
                            os.remove(ligsff_dir)
                        os.mkdir(ligsff_dir)
            else:
                F = 0
                if (ffselect == "Amber" or ffselect == "Charmm" or ffselect == "Gromos" or ffselect == "Opls"):
                    ligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], ffselect)
                    if not (os.path.isdir(ligsff_file) or os.path.isfile(ligsff_file)):
                        fileff = open(ligsff_file, "w")
                        fileff.close()
                    else:
                        try:
                            shutil.rmtree(ligsff_file)
                        except:
                            os.remove(ligsff_file)
                        fileff = open(ligsff_file, "w")
                        fileff.close()

        elif N == 3:
            if ffupload == "YES":
                F = 1
            else:
                F = 0

            ligsff_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff')
            if not (os.path.isdir(ligsff_dir) or os.path.isfile(ligsff_dir)):
                os.mkdir(ligsff_dir)
            else:
                try:
                    shutil.rmtree(ligsff_dir)
                except:
                    os.remove(ligsff_dir)
                os.mkdir(ligsff_dir)

        elif N == 4:
            if (ffselect == "Amber" or ffselect == "Charmm" or ffselect == "Gromos" or ffselect == "Opls"):
                ligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], ffselect)
                if not (os.path.isdir(ligsff_file) or os.path.isfile(ligsff_file)):
                    fileff = open(ligsff_file, "w")
                    fileff.close()
                else:
                    try:
                        shutil.rmtree(ligsff_file)
                    except:
                        os.remove(ligsff_file)
                    fileff = open(ligsff_file, "w")
                    fileff.close()

            listUdir = os.listdir(app.config['UPLOAD_FOLDER'])
            if ffupload == "Peptides":
                F = 2
                if "Proteins" in listUdir:
                    shutil.rmtree(os.path.join(app.config['UPLOAD_FOLDER'], "Proteins"))

            elif ffupload == "Proteins":
                F = 3
                if "Peptides" in listUdir:
                    shutil.rmtree(os.path.join(app.config['UPLOAD_FOLDER'], "Peptides"))

        # Render the quoted values in quoted.html
        return render_template("routetopselected.html", selected=selected, ffselect=ffselect, N=N, F=F)

    # User reached route via GET (as by clicking a link or via redirect)
    else:
        return render_template("selectffgroup.html", N=N)

@app.route('/runscripts')
@setup_route
def runscripts():
    return render_template("runscripts.html", N=N)

@app.route('/rlfilesupload', methods=['GET', 'POST'])
@setup_route
def rlupload():
    global N
    global F
    global ffselect
    if request.method == 'POST':
        if (N == 1 and F == 0) or (N == 2 and F == 0):
            listRdir = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors'))
            listLdir = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands'))
            if len(listLdir) > 0 or len(listRdir) > 0:
                flash('Additional files uploads not allowed. Setroute and upload all files once')
                return redirect("/")

        # check if the post request has the file part
            if 'Lfiles[]' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Rfile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            # if valid file(s) submitted, safe ligands
            Lfiles = request.files.getlist('Lfiles[]')
            for fileL in Lfiles:
                if fileL.filename == '':
                    flash('No selected files. Click upload menu to upload valid files')
                    return redirect("/")

                if fileL and allowed_file(fileL.filename):
                    Lfilename = secure_filename(fileL.filename)
                    fileL.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands', Lfilename))

            # if valid file(s) submitted, safe receptor
            Rfile = request.files['Rfile']
            if Rfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Rfile and allowed_file(Rfile.filename):
                Rfilename = secure_filename(Rfile.filename)
                Rfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfilename))

            flash('Successful Receptor and Ligands files Upload')
            return render_template("upload.html", Lfilelist=Lfiles, Rfilelist=Rfile, N=N, F=F, category='success')

        elif (N == 1 and F == 1) or (N == 2 and F == 1):
            listLdir = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands'))
            listTdir = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], ffselect))
            listRdir = os.listdir(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors'))
            if len(listLdir) > 0 or len(listTdir) > 0 or len(listRdir) > 0:
                flash('Additional files uploads not allowed. Setroute and upload all files once')
                return redirect("/")

            # check if the post request has the file part
            if 'Lfiles[]' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Ffiles[]' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Rfile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            # if valid file(s) submitted, safe ligands
            Lfiles = request.files.getlist('Lfiles[]')
            for fileL in Lfiles:
                if fileL.filename == '':
                    flash('No selected file. Click upload menu to upload valid files')
                    return redirect("/")

                if fileL and allowed_file(fileL.filename):
                    Lfilename = secure_filename(fileL.filename)
                    fileL.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands', Lfilename))

            # if valid file(s) submitted, safe ligands topologies
            Ffiles = request.files.getlist('Ffiles[]')
            for fileF in Ffiles:
                if fileF.filename == '':
                    flash('No selected file. Click upload menu to upload valid files')
                    return redirect("/")

                if fileF and allowed_file(fileF.filename):
                    Ffilename = secure_filename(fileF.filename)
                    fileF.save(os.path.join(app.config['UPLOAD_FOLDER'], ffselect, Ffilename))

            # if valid file(s) submitted, safe receptor
            Rfile = request.files['Rfile']
            if Rfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Rfile and allowed_file(Rfile.filename):
                Rfilename = secure_filename(Rfile.filename)
                Rfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfilename))

            flash('Successful Receptor and Ligands files Upload')
            return render_template("upload.html", Lfilelist=Lfiles, Ffilelist=Ffiles, Rfilelist=Rfile, N=N, F=F, category='success')

        elif (N == 3 and F == 0):
            # check if the post request has the file part
            if 'Lfile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Rfile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if not request.form.get("gmodsff"):
                flash('You must select a forcefileds group')
                return redirect("/")

            # Create subfolder for each ligand and receptor uploads
            x = 1
            while x > 0:
                Lfolder = "ligand" + str(x)
                Ldir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands', Lfolder)
                if not os.path.isdir(Ldir):
                    os.mkdir(Ldir)
                    break
                else:
                    x += 1

            z = 1
            while z > 0:
                Ffolder = "ligtop" + str(z)
                Fdir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder)
                if not os.path.isdir(Fdir):
                    os.mkdir(Fdir)
                    break
                else:
                    z += 1

            y = 1
            while y > 0:
                Rfolder = "receptor" + str(y)
                Rdir = os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfolder)
                if not os.path.isdir(Rdir):
                    os.mkdir(Rdir)
                    break
                else:
                    y += 1

			# if valid file(s) submitted, safe ligands
            Lfile = request.files['Lfile']
            if Lfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Lfile and allowed_file(Lfile.filename):
                Lfilename = secure_filename(Lfile.filename)
                Lfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands', Lfolder, Lfilename))

			# if valid file(s) submitted, safe receptor
            Rfile = request.files['Rfile']
            if Rfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Rfile and allowed_file(Rfile.filename):
                Rfilename = secure_filename(Rfile.filename)
                Rfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfolder, Rfilename))

		    # Store the submitted symbol as varibale and check if valid
            ffselect = request.form.get("gmodsff")
            if not (ffselect == "Amber" or ffselect == "Charmm" or ffselect == "Gromos" or ffselect == "Opls"):
                flash('You must select a forcefileds group')
                return redirect("/")

            ligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffselect)
            if not (os.path.isdir(ligsff_file) or os.path.isfile(ligsff_file)):
                fileff = open(ligsff_file, "w")
                fileff.close()

            flash('Successful Receptor and Ligands files Upload. Click Uploads Menu to upload additional pair of receptor and ligand')
            return render_template("upload.html", Lfilelist=Lfile, Rfilelist=Rfile, N=N, F=F, category='success')

        elif (N == 3 and F == 1):
            # check if the post request has the file part
            if 'Lfile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Ffile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Rfile' not in request.files:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if not request.form.get("gmodsff"):
                flash('You must select a forcefileds group')
                return redirect("/")

            # Create subfolder for each ligand and receptor uploads
            x = 1
            while x > 0:
                Lfolder = "ligand" + str(x)
                Ldir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands', Lfolder)
                if not os.path.isdir(Ldir):
                    os.mkdir(Ldir)
                    break
                else:
                    x += 1

            z = 1
            while z > 0:
                Ffolder = "ligtop" + str(z)
                Fdir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder)
                if not os.path.isdir(Fdir):
                    os.mkdir(Fdir)
                    break
                else:
                    z += 1

            y = 1
            while y > 0:
                Rfolder = "receptor" + str(y)
                Rdir = os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfolder)
                if not os.path.isdir(Rdir):
                    os.mkdir(Rdir)
                    break
                else:
                    y += 1

		    # Store the submitted symbol as varibale and check if valid
            ffselect = request.form.get("gmodsff")
            if not (ffselect == "Amber" or ffselect == "Charmm" or ffselect == "Gromos" or ffselect == "Opls"):
                flash('You must select a forcefileds group')
                return redirect("/")

            ligsff_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffselect)
            if not os.path.isdir(ligsff_dir):
                os.mkdir(ligsff_dir)

			# if valid file(s) submitted, safe ligands
            Lfile = request.files['Lfile']
            if Lfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Lfile and allowed_file(Lfile.filename):
                Lfilename = secure_filename(Lfile.filename)
                Lfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligands', Lfolder, Lfilename))

			# if valid file(s) submitted, safe ligand topology
            Ffile = request.files['Ffile']
            if Ffile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Ffile and allowed_file(Ffile.filename):
                Ffilename = secure_filename(Ffile.filename)
                Ffile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffselect, Ffilename))

			# if valid file(s) submitted, safe receptor
            Rfile = request.files['Rfile']
            if Rfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Rfile and allowed_file(Rfile.filename):
                Rfilename = secure_filename(Rfile.filename)
                Rfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfolder, Rfilename))

            flash('Successful Receptor, Ligands and Topology files Upload. Click Uploads Menu to add new set')
            return render_template("upload.html", Lfilelist=Lfile, Ffilelist=Ffile, Rfilelist=Rfile, N=N, F=F, category='success')

        elif (N == 4 and F > 1):
            # check if the post request has the file part
            if 'Ptfiles[]' not in request.files and not F == 3:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            if 'Prfiles[]' not in request.files and not F == 2:
                flash('No file part. Click upload menu to upload valid files')
                return redirect("/")

            # if valid file(s) submitted, safe peptides
            if F == 2:
                Ptfiles = request.files.getlist('Ptfiles[]')
                for filePt in Ptfiles:
                    if filePt.filename == '':
                        flash('No selected files. Click upload menu to upload valid files')
                        return redirect("/")

                    if filePt and allowed_file(filePt.filename):
                        Ptfilename = secure_filename(filePt.filename)
                        filePt.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Peptides', Ptfilename))

                flash('Successful Peptide Sequence files Upload. Click Uploads Menu to add more Peptides')
                return render_template("upload.html", Ptfilelist=Ptfiles, N=N, F=F, category='success')

            # if valid file(s) submitted, safe proteins
            if F == 3:
                Prfiles = request.files.getlist('Prfiles[]')
                for filePr in Prfiles:
                    if filePr.filename == '':
                        flash('No selected files. Click upload menu to upload valid files')
                        return redirect("/")

                    if filePr and allowed_file(filePr.filename):
                        Prfilename = secure_filename(filePr.filename)
                        filePr.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Proteins', Prfilename))

                flash('Successful Protein Structure files Upload. Click Uploads Menu to add more Proteins')
                return render_template("upload.html", Prfilelist=Prfiles, N=N, F=F, category='success')

    else:
        return render_template("rlfilesupload.html", N=N, F=F)

@app.route('/rlmulti')
def rlmulti_script():
    gmodsScripts.gmodsRLmulti.RLmulti()
    flash('Protein - Ligand Complex Generation was Successful')
    return render_template("runscripts.html", N=N, category='success')

@app.route('/rlmany')
def rlmany_script():
    gmodsScripts.gmodsRLmany.RLmany()
    flash('Protein - Ligand Complex Generation was Successful')
    return render_template("runscripts.html", N=N, category='success')

@app.route('/rlsingle')
def rlsingle_script():
    gmodsScripts.gmodsRLsingle.RLsingle()
    flash('Protein - Ligand Complex Generation was Successful', 'success')
    return render_template("runscripts.html", N=N)

@app.route('/ppmore')
def ppmore_script():
    gmodsScripts.gmodsPPmore.PPmore()
    flash('Protein(s) and/or Peptide(s) MDS Setup was Successful')
    return render_template("runscripts.html", N=N, category='success')

@app.route('/scmds', methods=['GET', 'POST'])
def scmds_setup():
    global N
    if request.method == 'POST':
        # check if the post request has the file part
        if 'mdsfiles[]' not in request.files:
            flash('No file part. If loading files manually, scroll down for instruction')
            return redirect("/scmds")

        mdpids = ['minzsd', 'minzcg', 'equnvt', 'equnpt', 'equpmd', 'pmds']
        for id in mdpids:
            if id not in request.files:
                flash('No file part. If loading files manually, scroll down for instruction')
                return redirect("/scmds")

        # Create MDP and MDSHOME folders. Subfolders will then be created inside MDSHOME for each simulation can be created and populated
        mdp_dir = os.path.join(appDIR, 'MDP')
        if not os.path.isdir(mdp_dir):
            os.mkdir(mdp_dir)
        else:
            shutil.rmtree(mdp_dir)
            os.mkdir(mdp_dir)

        mdsdir = os.path.join(gmxDIR, 'MDSHOME')
        if not os.path.isdir(mdsdir):
            os.mkdir(mdsdir)

        subdir = "mdsgmx"
        mdsubdir = os.path.join(mdsdir, subdir)
        while True:
            if not os.path.isdir(mdsubdir):
                os.mkdir(mdsubdir)
                break
            elif os.path.isdir(mdsubdir):
                if len(os.listdir(mdsubdir)) == 0:
                    break
                else:
                    x = 1
                    while x > 0:
                        backup = subdir + str(x) + ".bk"
                        backupf = os.path.join(mdsdir, backup)
                        if os.path.isdir(backupf):
                            x += 1
                            continue
                        else:
                            os.rename(mdsubdir, backupf)
                            break
                    continue

        # if valid file(s) submitted, safe needed gmx files in the subdirectory
        gmxfiles = request.files.getlist('mdsfiles[]')
        for gfile in gmxfiles:
            if gfile.filename == '':
                flash('No selected file. If loading files manually, scroll down for instruction')
                return redirect("/scmds")

            if gfile and allowed_file(gfile.filename):
                gfilename = secure_filename(gfile.filename)
                gfile.save(os.path.join(mdsubdir, gfilename))

        # if valid file(s) submitted, safe mdp files in the MDP directory
        mdpfiles = ['mzsdfile', 'mzcgfile', 'envtfile', 'enptfile', 'epmdfile', 'pmdsfile']
        mdpidrename = ['minzsd.mdp', 'minzcg.mdp', 'equnvt.mdp', 'equnpt.mdp', 'equpmd.mdp', 'pmds.mdp']
        mdpids = ['minzsd', 'minzcg', 'equnvt', 'equnpt', 'equpmd', 'pmds']

        mdpI = 0
        for file in mdpfiles:
            file = request.files[mdpids[mdpI]]
            if file.filename == '':
                flash('No selected file. If loading files manually, scroll down for instruction')
                return redirect("/scmds")

            if file and allowed_file(file.filename):
                mdpfilename = secure_filename(file.filename)
                file.save(os.path.join(mdp_dir, mdpfilename))
                savedfile = os.path.join(mdp_dir, mdpfilename)
                os.rename(savedfile, os.path.join(mdp_dir, mdpidrename[mdpI]))
                mdpI += 1

        # Render template or redirect to run gmodsSCmds.py script
        flash('Setting up files for MDS was Successful')
        return render_template("scmdsrun.html", N=N, category='success')
    else:
        return render_template("scmds.html", N=N)

@app.route('/openclheaders', methods=['GET', 'POST'])
def oclupload():
    global N
    global A

    if request.method == 'POST':
        #Check to be sure copying OpenCL Header was activated
        if not A > 0:
            return render_template("openclheaders.html", N=N)

        # check if the post request has the file part
        if 'oclheaders[]' not in request.files:
            flash('No file part detected')
            return redirect("/scmdsrun")

        #Set variables and generate directories id not yet exist
        mdsmain = os.path.join(gmxDIR, 'MDSHOME')
        if not os.path.isdir(mdsmain):
            os.mkdir(mdsmain)

        oclfolder = os.path.join(mdsmain, 'opencl')
        if not os.path.isdir(oclfolder):
            os.mkdir(oclfolder)
        else:
            shutil.rmtree(oclfolder)
            os.mkdir(oclfolder)

        # if valid file(s) submitted, safe ligands
        oclfiles = request.files.getlist('oclheaders[]')
        for header in oclfiles:
            if header.filename == '':
                flash('No selected file. Please upload OpenCl Headers')
                return redirect("/scmdsrun")

            if header and allowed_file(header.filename):
                oclfilename = secure_filename(header.filename)
                header.save(os.path.join(oclfolder, oclfilename))

        flash('Successful OpenCL Headers files Upload')
        return render_template("scmdsrun.html", N=N, category='success')

    else:
        A += 1
        return render_template("openclheaders.html", N=N)

@app.route('/runscmds')
def runmds_script():
    gmodsScripts.gmodsSCmds.SCmds()
    flash('MDS of solvated complex was Successful')
    return render_template("index.html", N=N, category='success')

@app.route('/scmdsrun')
def mdsrun_page():
    return render_template("scmdsrun.html", N=N)

if __name__ == "__main__":
    # app.run() for debug
    try:
        ui.run()
    except:
        ui(app, width=800, height=600, port=1380, window_title="pyGROMODS")
