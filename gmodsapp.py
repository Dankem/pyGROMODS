#!/usr/bin/env python
"""
    pyGROMODS-v2024.02 Release

                <<<  NO WARRANTY AT ALL!!!  >>>
    
    Requirements: Python 3.5 or higher
                  Antechamber and related AmberTools
                  OpenBabel (strongly recommended for use with acpype)
				  acpype (latest version recommended with all its requirements)
				  GROMACS (Compulsory)

    Daniyan, Michael Oluwatoyin, B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) Biochemistry
    Department of Pharmacology, Faculty of Pharmacy
    Obafemi Awolowo Universiy, Ile-Ife, Nigeria.
    >>http://www.oauife.edu.ng<<

    mdaniyan@oauife.edu.ng; toyinpharm@gmail.com
"""
import sys

if sys.version_info < (3, 5):
    raise Exception("Python 3.5 or a more recent version is required.")

from pathlib import Path
import os
import subprocess
import shutil
import platform
import time
import random

try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass

import gmodsScripts.gmodsRLmulti
import gmodsScripts.gmodsRLmany
import gmodsScripts.gmodsRLsingle
import gmodsScripts.gmodsSCmds
import gmodsScripts.gmodsCONmds
import gmodsScripts.gmodsPPmore
from gmodsScripts.gmodsCRPackages import reqPackages_Check
from gmodsScripts.gmodsHelpers import printWarning, printNote, select_folder, gmxtop, gmxmdsFChecks

from flask import Flask, flash, request, redirect, render_template, g
from functools import wraps
from werkzeug.utils import secure_filename
from tkinter import Tk, messagebox

gmxDIR = Path.cwd()
appDIR = os.path.abspath(os.path.dirname(sys.argv[0]))

if getattr(sys, 'frozen', False):
    template_folder = os.path.join(appDIR, 'templates')
    static_folder = os.path.join(appDIR, 'static')
    scripts_folder = os.path.join(appDIR, 'gmodsScripts')
    app = Flask(__name__, template_folder=template_folder, static_folder=static_folder)
else:
    app = Flask(__name__)

# Let's get templates to auto-reload
app.jinja_env.auto_reload = True
app.config["TEMPLATES_AUTO_RELOAD"] = True

# Let's check your system to be sure all required packages are installed
printNote("System checks for required packages in progress ...")
check_packages = reqPackages_Check()
printNote("Please check the above results and Respond to Popup Question")
time.sleep(5)

crp = Tk()
crp.title("Platform Check Results")
crp.withdraw()
crp.attributes('-topmost', True)

checkmessage = messagebox.askquestion("Platform Check Results", f"{check_packages} \n\n DO YOU REALLY WANT TO CONTINUE?")
crp.destroy()

if not checkmessage.lower() == "yes":
    raise Exception("Operation Interrupted by User. Make necessary corrections and restart")
else:
    print(check_packages)
    printNote("You have choosen to continue with the process")
time.sleep(5)

# # Let's setup our GUI interface and add app and other parameters
# First, Let's generate unique port and check it availability
gmodsport = ""
while True:
    tryport = ""
    for rdport in random.sample(range(3, 10), 4):
        tryport += str(rdport)

    try:
        import socket as sk
    except ImportError:
        gmodsport = int(tryport)
        break

    s = sk.socket(sk.AF_INET, sk.SOCK_STREAM)
    result = s.connect_ex(("127.0.0.1", int(tryport)))
    if result == 0:
        continue
    else:
        gmodsport = int(tryport)
        break

# Tested and trusted versions of flaskwebgui are bundled with the package to avoid breaking
plat4m = platform.system()
if plat4m.lower() == "linux" or plat4m.lower() == "posix":
    try:
        from gmodsScripts.flaskwebgui.flaskwebgui037 import FlaskUI
        ui = FlaskUI(app, width=900, height=650, port=gmodsport) # for v037
    except Exception as e1:
        print(f"'flaskwebgui' failed with error: {e1}")
        print("Trying again ...")
        try:
            from gmodsScripts.flaskwebgui.flaskwebgui112 import FlaskUI
            ui = FlaskUI(app=app, width=900, height=650, port=gmodsport, server="flask") # for v112
        except Exception as e2:
            printWarning(f"'flaskwebgui' failed with error: {e2}")
            print("Please check and correct the errors")
            printNote("To access GUI interface, copy the generated web link below to your default browser instead")
            ui = app
else:
    try:
        from gmodsScripts.flaskwebgui.flaskwebgui112 import FlaskUI
        ui = FlaskUI(app=app, width=900, height=650, port=gmodsport, server="flask") # for v112
    except Exception as e1:
        print(f"'flaskwebgui' failed with error: {e1}")
        print("Trying again ...")
        try:
            from gmodsScripts.flaskwebgui.flaskwebgui037 import FlaskUI
            ui = FlaskUI(app, width=900, height=650, port=gmodsport) # for v037
        except Exception as e2:
            printWarning(f"'flaskwebgui' failed with error: {e2}")
            print("Check and correct the errors")
            printNote("To access GUI interface, copy the generated web link below to your default browser instead")
            ui = app

# Setup how errors will be handled
@app.errorhandler(Exception)
def all_exceptions(e):
    os.chdir(gmxDIR)
    return render_template("errors.html", e=e, N=N)

# Auto generate the required secret key
KEY = os.urandom(14)
app.config['SECRET_KEY'] = KEY

# Performing some folder setup
UPLOAD_FOLDER = os.path.join(appDIR, 'Uploads')
if not os.path.isdir(UPLOAD_FOLDER):
    os.mkdir(UPLOAD_FOLDER)
else:
    shutil.rmtree(UPLOAD_FOLDER)
    os.mkdir(UPLOAD_FOLDER)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

ALLOWED_EXTENSIONS = {'pdb', 'mol2', 'gro', 'top', 'itp', 'mdp', 'h', 'cl', 'clh', 'in', 'ndx'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# Choose the preferred Working directory and change to it
title = "Select Working Directory"
wFolder = select_folder(title)
gmxDIR = os.path.join(os.getcwd(), wFolder)
os.chdir(gmxDIR)
print(f"Selected Working Directory is {gmxDIR}")

# Set some global varibales
global N, A, F, selected, ffselect, gmxtoplist, gmxtopdir, dictff
N = 0
A = 0
F = 0
selected = " "
ffselect = " "
gmxtoplist = " "
gmxtopdir = " "
dictff = {}

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

@app.route('/quit')
def quitme():
    qws = Tk()
    qws.title('Quiting pyGROMODS')
    qws.withdraw()
    qws.attributes('-topmost', True)

    quitmessage = messagebox.askquestion('Quiting pyGROMODS', 'DO YOU REALLY WANT TO QUIT pyGROMODS?')
    qws.destroy()

    if quitmessage.lower() == "yes":
        printNote("Closing the application. To restart, run again")
        print("Manually close the web GUI")
        os.kill(os.getpid(), 9)
    else:
        return render_template("index.html", N=N)

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

        # Setup directories to store ligands and protein receptor
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
    global gmxtoplist
    global gmxtopdir
    global dictff

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
            flash('You must select a GROMACS compatible forcefileds group')
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

        # Get the list of forcefields compatible with ffselect
        gmxtoplist, gmxtopdir = gmxtop()
        seltoplist = []
        for gtl in gmxtoplist:
            if ffselect == "Amber" and gtl[0:5].lower() == ffselect.lower():
                seltoplist.append(gtl)
            elif ffselect == "Charmm" and gtl[0:6].lower() == ffselect.lower():
                seltoplist.append(gtl)
            elif ffselect == "Gromos" and gtl[0:6].lower() == ffselect.lower():
                seltoplist.append(gtl)
            elif ffselect == "Opls" and gtl[0:4].lower() == ffselect.lower():
                seltoplist.append(gtl)

        # Setup a Dictionary with ff as key and water models as values
        dictff = {}
        for ftl in gmxtoplist:
            ftlfolder = ftl + ".ff"
            ftldir = os.path.join(gmxtopdir, ftlfolder)
            listftldir = os.listdir(ftldir)
            wmodels = []
            if not "watermodels.dat" in listftldir:
                wmodels.append("none")
            else:
                wmodels.append("none")
                wmfile = os.path.join(ftldir, "watermodels.dat")
                wmfileopen = open(wmfile, "r")
                wmfileread = wmfileopen.readlines()
                for wline in wmfileread:
                    wmlist = wline.split()
                    wmodels.append(wmlist[0])
                wmfileopen.close()
            dictff[ftl] = wmodels

        # Render the collected values in routetopselected.html
        return render_template("routetopselected.html", selected=selected, ffselect=ffselect, seltoplist=seltoplist, N=N, F=F)

    # User reached route via GET (as by clicking a link or via redirect)
    else:
        return render_template("selectffgroup.html", N=N)

@app.route('/rlfilesupload', methods=['GET', 'POST'])
@setup_route
def rlupload():
    global N
    global F
    global ffselect
    global gmxtoplist
    global gmxtopdir
    global dictff

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
                flash('You must select a forcefiled or opt for interactive mode')
                return redirect("/")

            if not request.form.get("watertype"):
                flash('You must select a water model or opt for interactive mode')
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
            ffwater = request.form.get("watertype")
            ffselect = request.form.get("gmodsff")
            if not (ffselect[0:5].capitalize() == "Amber" or ffselect[0:6].capitalize() == "Charmm" or ffselect[0:6].capitalize() == "Gromos" or ffselect[0:4].capitalize() == "Opls" or ffselect.capitalize() == "Select"):
                flash('You must select a forcefiled. Try upload again')
                return redirect("/")

            ligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffselect)
            if not (os.path.isdir(ligsff_file) or os.path.isfile(ligsff_file)) and not ffselect == "select":
                fileff = open(ligsff_file, "w")
                fileff.close()

            elif ffselect == "select":
                newffselect = ffselect + "_ff"
                nligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, newffselect)
                if not (os.path.isdir(nligsff_file) or os.path.isfile(nligsff_file)):
                    fileff = open(nligsff_file, "w")
                    fileff.close()

            ligsww_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffwater)
            if not (os.path.isdir(ligsww_file) or os.path.isfile(ligsww_file)) and not ffwater == "select":
                fileww = open(ligsww_file, "w")
                fileww.close()

            elif ffwater == "select":
                newffwater = ffwater + "_ww"
                nligsww_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, newffwater)
                if not (os.path.isdir(nligsww_file) or os.path.isfile(nligsww_file)):
                    fileww = open(nligsww_file, "w")
                    fileww.close()

            flash('Successful Receptor and Ligands files Upload. Click Uploads Menu to upload additional pair of receptor and ligand')
            return render_template("upload.html", Lfilelist=Lfile, Rfilelist=Rfile, ffselect=ffselect, ffwater=ffwater, N=N, F=F, category='success')

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
                flash('You must select a forcefiled or opt for interactive mode')
                return redirect("/")

            if not request.form.get("watertype"):
                flash('You must select a water model or opt for interactive mode')
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
            ffwater = request.form.get("watertype")
            ffselect = request.form.get("gmodsff")

            if not (ffselect[0:5].capitalize() == "Amber" or ffselect[0:6].capitalize() == "Charmm" or ffselect[0:6].capitalize() == "Gromos" or ffselect[0:4].capitalize() == "Opls" or ffselect[0:6].capitalize() == "Select"):
                flash('You must select a forcefiled. Try upload again')
                return redirect("/")

            if ffselect[0:5].capitalize() == "Amber":
                ff_dir = "Amber"
            elif ffselect[0:6].capitalize() == "Charmm":
                ff_dir = "Charmm"
            elif ffselect[0:6].capitalize() == "Gromos":
                ff_dir = "Gromos"
            elif ffselect[0:4].capitalize() == "Opls":
                ff_dir = "Opls"
            elif ffselect[0:6].capitalize() == "Select":
                ff_dir = "Select"

            ligsff_dir = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ff_dir)
            if not os.path.isdir(ligsff_dir):
                os.mkdir(ligsff_dir)

            ligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffselect)
            if not (os.path.isdir(ligsff_file) or os.path.isfile(ligsff_file)):
                fileff = open(ligsff_file, "w")
                fileff.close()

            elif os.path.isdir(ligsff_file) and ffselect == "select":
                #rename ffselect to create a file that can be tracked
                newffselect = ffselect + "_ff"
                nligsff_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, newffselect)
                if not os.path.isfile(nligsff_file):
                    fileff = open(nligsff_file, "w")
                    fileff.close()

            ligsww_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ffwater)
            if not (os.path.isdir(ligsww_file) or os.path.isfile(ligsww_file)):
                fileww = open(ligsww_file, "w")
                fileww.close()

            elif os.path.isdir(ligsww_file) and ffwater == "select":
                # rename ffwater to create a file that can be tracked
                newffwater = ffwater + "_ww"
                nligsww_file = os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, newffwater)
                if not os.path.isfile(nligsww_file):
                    fileff = open(nligsww_file, "w")
                    fileff.close()

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
                Ffile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Ligsff', Ffolder, ff_dir, Ffilename))

			# if valid file(s) submitted, safe receptor
            Rfile = request.files['Rfile']
            if Rfile.filename == '':
                flash('No selected file. Click upload menu to upload valid files')
                return redirect("/")

            if Rfile and allowed_file(Rfile.filename):
                Rfilename = secure_filename(Rfile.filename)
                Rfile.save(os.path.join(app.config['UPLOAD_FOLDER'], 'Receptors', Rfolder, Rfilename))

            flash('Successful Receptor, Ligands and Topology files Upload. Click Uploads Menu to add new set')
            return render_template("upload.html", Lfilelist=Lfile, Ffilelist=Ffile, Rfilelist=Rfile, ffselect=ffselect, ffwater=ffwater, N=N, F=F, category='success')

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
        return render_template("rlfilesupload.html", N=N, F=F, gmxtoplist=gmxtoplist, dictff=dictff)

@app.route('/runscripts', methods=['GET', 'POST'])
@setup_route
def runscripts():
    global selected
    global gmxtoplist
    global gmxtopdir
    global dictff

    if request.method == 'POST':
        # Check if user choose the correct default values
        if not request.form.get("timeout"):
            flash('You must supply a valid value for timeout')
            return redirect("/runscripts")

        if not request.form.get("changed"):
            flash('You must supply a valid value for -d box size')
            return redirect("/runscripts")

        if not request.form.get("changebt"):
            flash('You must select a valid value for -bt box type')
            return redirect("/runscripts")

        if not request.form.get("approach"):
            flash('You must select a valid approach to use')
            return redirect("/runscripts")

        if not request.form.get("forcefield"):
            flash('You must select a valid forcefield to use')
            return redirect("/runscripts")

        if not request.form.get("watertype"):
            flash('You must select a valid water model to use')
            return redirect("/runscripts")

        # If the submitted data are valid, store them in varibales
        ftimeout = int(request.form.get("timeout"))
        fchanged = float(request.form.get("changed"))
        fchangebt = request.form.get("changebt")
        fapproach = request.form.get("approach")
        fforcefield = request.form.get("forcefield")
        fwatertype = request.form.get("watertype")

        if not (fapproach == "Interactive" and fforcefield == "select" and fwatertype == "select"):
            approach = "B"
        else:
            approach = "A"

		# Set the variable up as list of default values
        fdefaults = [fforcefield, fwatertype, fchangebt, fchanged, ftimeout, approach]

		# Pass the defaults to the script
        if selected == "RLmulti":
            gmodsScripts.gmodsRLmulti.RLmulti(appDIR, gmxDIR, fdefaults)
            print("Performing post setup analysis...")
            rlstatus = ""
            if "gmxmds" in os.listdir():
                if "fsolvated.gro" in os.listdir(os.path.join(os.getcwd(), "gmxmds")) or "ufsolvate.gro" in os.listdir(os.path.join(os.getcwd(), "gmxmds")):
                    rlstatus = "pass"
                else:
                    rlstatus = "fail"
            else:
                print(f"No gmxmds subfolder was detected in {os.getcwd()}")
                print("Assuming process completed successfully")
                rlstatus = "pass"
            print("Analysis completed")

            os.chdir(gmxDIR)

            if rlstatus == "pass":
                flash('Protein - Ligand Complex Generation was Successful')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='success')
            else:
                flash('Protein - Ligand Complex Generation was Unsuccessful')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='warning')

        elif selected == "RLmany":
            gmodsScripts.gmodsRLmany.RLmany(appDIR, gmxDIR, fdefaults)
            print("Performing post job analysis...")
            nmpass = 0
            nmfail = 0
            for lmany in os.listdir():
                if "gmxmds" in os.listdir(os.path.join(os.getcwd(), lmany)):
                    if "fsolvated.gro" in os.listdir(os.path.join(os.getcwd(), lmany, "gmxmds")) or "ufsolvate.gro" in os.listdir(os.path.join(os.getcwd(), lmany, "gmxmds")):
                        nmpass += 1
                    else:
                        nmfail += 1
            print("Analysis completed")

            os.chdir(gmxDIR)

            if nmpass > 0 and nmfail == 0:
                flash(f'Protein - Ligand Complex Generation Completed with: {nmpass} PASSED / {nmfail} FAILED')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='success')
            else:
                flash(f'Protein - Ligand Complex Generation Completed with: {nmpass} PASSED / {nmfail} FAILED')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='warning')

        elif selected == "RLsingle":
            gmodsScripts.gmodsRLsingle.RLsingle(appDIR, gmxDIR, fdefaults, dictff)
            print("Performing post job analysis...")
            nspass = 0
            nsfail = 0
            for lsingle in os.listdir():
                if "gmxmds" in os.listdir(os.path.join(os.getcwd(), lsingle)):
                    if "fsolvated.gro" in os.listdir(os.path.join(os.getcwd(), lsingle, "gmxmds")) or "ufsolvate.gro" in os.listdir(os.path.join(os.getcwd(), lsingle, "gmxmds")):
                        nspass += 1
                    else:
                        nsfail += 1
            print("Analysis completed")

            os.chdir(gmxDIR)

            if nspass > 0 and nsfail == 0:
                flash(f'Protein - Ligand Complex Generation Completed with: {nspass} PASSED / {nsfail} FAILED')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='success')
            else:
                flash(f'Protein - Ligand Complex Generation Completed with: {nspass} PASSED / {nsfail} FAILED')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='warning')

        elif selected == "PPmore":
            gmodsScripts.gmodsPPmore.PPmore(appDIR, gmxDIR, fdefaults)
            print("Performing post job analysis...")
            nppass = 0
            npfail = 0
            for pmore in os.listdir():
                if "gmxmds" in os.listdir(os.path.join(os.getcwd(), pmore)):
                    if "fsolvated.gro" in os.listdir(os.path.join(os.getcwd(), pmore, "gmxmds")) or "ufsolvate.gro" in os.listdir(os.path.join(os.getcwd(), pmore, "gmxmds")):
                        nppass += 1
                    else:
                        npfail += 1
            print("Analysis completed")

            os.chdir(gmxDIR)

            if nppass > 0 and npfail == 0:
                flash(f'Protein - Ligand Complex Generation Completed with: {nppass} PASSED / {npfail} FAILED')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='success')
            else:
                flash(f'Protein - Ligand Complex Generation Completed with: {nppass} PASSED / {npfail} FAILED')
                return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='warning')

        else:
            os.chdir(gmxDIR)
            flash('No recognized route for generating MDS input files was detected')
            return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff, category='warning')

    else:
        return render_template("runscripts.html", N=N, selected=selected, gmxtoplist=gmxtoplist, dictff=dictff)

@app.route('/scmds', methods=['GET', 'POST'])
def scmds_setup():
    global N
    global A

    # Reset the value of A
    A = 0

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
        return render_template("scmdsrun.html", N=N, A=A, category='success')
    else:
        return render_template("scmds.html", N=N)

@app.route('/openclactivate', methods=['GET', 'POST'])
def oclactivate():
    global N
    global A

    if request.method == 'POST':
        # Check if user choose the correct default values
        if not request.form.get("activate"):
            flash('You must supply YES or NO for OpenCL activation')
            return redirect("/scmdsrun.html")

        # If valid value was supplied
        if request.form.get("activate") == "NO":
            A += 1
            return render_template("scmdsrun.html", N=N, A=A)

        if request.form.get("activate") == "YES":
            A += 1
            return render_template("openclheaders.html", N=N)

    else:
        return render_template("scmdsrun.html", N=N)

@app.route('/openclheaders', methods=['GET', 'POST'])
def oclupload():
    global N
    global A

    if request.method == 'POST':
        #Check to be sure copying OpenCL Header was activated
        if not A > 0:
            return render_template("scmdsrun.html", N=N)

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
        return render_template("scmdsrun.html", N=N, A=A, category='success')

    else:
        A += 1
        return render_template("openclheaders.html", N=N)

@app.route('/runscmds')
def runmds_script():
    gmodsScripts.gmodsSCmds.SCmds(appDIR, gmxDIR)
    os.chdir(gmxDIR)
    flash('MDS of solvated complex was Successful')
    return render_template("index.html", N=N, category='success')

@app.route('/scmdsrun')
def mdsrun_page():
    global A
    return render_template("scmdsrun.html", N=N, A=A)

@app.route("/selectmdstype", methods=["GET", "POST"])
def mdstype_selection():
    """Determine the approach to use for the MDS run"""
    global N
    global F
    global selected
    global ffselect

    # User reached route via POST (as by submitting a form via POST)
    if request.method == "POST":
        # Check if user choose a MDS type from the dropdown menu or ignore
        if not request.form.get("mdstype"):
            flash('You must select an approach to use for MDS run')
            return redirect("/selectmdstype")

        # Store the submitted symbol as varibale
        mdtype = request.form.get("mdstype")

        # Check if selected route is valid
        if not (mdtype == "Fresh" or mdtype == "Continuation"):
            flash('You must select an approach to use for MDS run')
            return redirect("/selectmdstype")

        # If selection is valid, perform appropraite operation
        if mdtype == "Fresh":
            return render_template("scmds.html", N=N)
        elif mdtype == "Continuation":
            return render_template("scmdcontinuation.html", N=N)

    # User reached route via GET (as by clicking a link or via redirect)
    else:
        return render_template("selectmdstype.html", N=N)

@app.route('/runconmds')
def conmds_script():
    # Select the mdsgmx_* directory meant for continuation and change to it
    title = "Select Continuation Working Directory"
    cFolder = select_folder(title)
    conDIR = os.path.join(os.getcwd(), cFolder)
    print("Selected Continuation Working Directory is ", conDIR)

    os.chdir(conDIR)
    cwdirall = os.listdir()

    # Checking currently available files and folders in the directory
    printNote("Checking the previous files ....")

    if not ('fsolvated.gro' in cwdirall and 'topol.top' in cwdirall):
        flash('The required compulsory files are missing from the selected directory')
        return render_template("scmdcontinuation.html", N=N)

    if not ('posre.itp' in cwdirall and 'posre_udp.itp' in cwdirall):
        printNote("No restraint file was found. If defined -DP0SRE OR -DPOSRE_UDP in .mdp file, you must include posre.itp or posre_udp.itp respectively in the working folder")

    elif 'posre_udp.itp' in cwdirall and 'posre.itp' in cwdirall:
        printNote("Two restraint files are found. Only the one defined .mdp file will be used")

    elif 'posre_udp.itp' in cwdirall or 'posre.itp' in cwdirall:
        for sre in cwdirall:
            if sre == 'posre_udp.itp':
                print("Only user defined restriant file was found. To use it, you must define -DPOSRE_UDP in .mdp file")
                print("However, if already defined -DPOSRE, rename the file to posre.itp")
                response = input("Type YES to rename. Otherwise, press ENTER to continue: ")
                if response.lower() == "yes" or response.lower() == "y":
                    os.rename('posre_udp.itp', 'posre.itp')
            elif sre == 'posre.itp':
                print("Only posre.itp restraint file was found. To use it, please check that -DPOSRE is defined in .mdp file")

    if 'indexfile.ndx' in cwdirall:
        printNote("Found index file and will be supplied to grompp")
        
    if not ('LIGS_at.itp' in cwdirall or 'LIGS_mt.itp' in cwdirall):
        printNote("One or more Ligand(s) related files are missing. They may not be needed if all needed parameters are contained in topol.top file. Otherwise, the process may fail")

    if not 'mdps' in cwdirall:
        flash('The mdps folder is missing from the selected directory')
        return render_template("scmdcontinuation.html", N=N)
    else:
        mdpfilename = ['minzsd.mdp', 'minzcg.mdp', 'equnvt.mdp', 'equnpt.mdp', 'equpmd.mdp', 'pmds.mdp']
        MDPfiles = os.listdir('mdps')
        if len(MDPfiles) == 0:
            flash('The mdps folder, meant to contain required .mdp files, is empty')
            return render_template("scmdcontinuation.html", N=N)

        mdpf = 0
        for file in MDPfiles:
            if file in mdpfilename:
                mdpf += 1

        if not mdpf == int(len(mdpfilename)):
            flash('One or more required mdps file(s) missing or renamed')
            return render_template("scmdcontinuation.html", N=N)

        print("found needed .mdp files")

    # Checking the accuracy of generated folders and files to determine point of continuation
    printNote("Checking to determine the point of continuation....")

    if 'EM' in cwdirall:
        emdir = os.listdir('EM')
        if not 'minzcg.gro' in emdir:
            flash("This appears to be a Fresh run. Please consider Fresh run instead")
            return render_template("selectmdstype.html", N=N)
        else:
            printNote("Noted: Minimization was successfully completed. Nothing to do")
    else:
        printNote("This appears to be a Fresh run. Please consider Fresh run instead")
        flash("This appears to be a Fresh run. Please consider Fresh run instead")
        return render_template("selectmdstype.html", N=N)

    if 'NVT' in cwdirall:
        nvtdir = os.listdir('NVT')
        if not ('equnvt.gro' in nvtdir and 'equnvt.cpt' in nvtdir):
            printNote("MDS continuation will begin with NVT Equilibration")
            response = input("To continue, type YES/y. Otherwise press ENTER to abort: ")
            if not (response.lower() == "yes" or response.lower() == "y"):
                os.chdir(gmxDIR)
                flash("MDS Continuation aborted")
                return render_template("selectmdstype.html", N=N)
            else:
                gmodsScripts.gmodsCONmds.nvtSCmds()
                os.chdir(gmxDIR)
                flash('MDS Continuation was Successfully completed')
                return render_template("index.html", N=N, category='success')
        else:
            printNote("Noted: NVT Equilibration was successfully completed. Nothing to do")
    else:
        printNote("MDS continuation will begin with NVT Equilibration")
        response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            os.chdir(gmxDIR)
            flash("MDS Continuation aborted")
            return render_template("selectmdstype.html", N=N)
        else:
            gmodsScripts.gmodsCONmds.nvtSCmds()
            os.chdir(gmxDIR)
            flash('MDS Continuation was Successful completed')
            return render_template("index.html", N=N, category='success')

    if 'NPT' in cwdirall:
        nptdir = os.listdir('NPT')
        if not ('equnpt.gro' in nptdir and 'equnpt.cpt' in nptdir):
            printNote("MDS continuation will begin with NPT Equilibration")
            response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
            if not (response.lower() == "yes" or response.lower() == "y"):
                os.chdir(gmxDIR)
                flash("MDS Continuation aborted")
                return render_template("selectmdstype.html", N=N)
            else:
                gmodsScripts.gmodsCONmds.nptSCmds()
                os.chdir(gmxDIR)
                flash('MDS Continuation was Successful completed')
                return render_template("index.html", N=N, category='success')
        else:
            printNote("Noted: NPT Equilibration was successfully completed. Nothing to do")
    else:
        printNote("MDS continuation will begin with NPT Equilibration")
        response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            os.chdir(gmxDIR)
            flash("MDS Continuation aborted")
            return render_template("selectmdstype.html", N=N)
        else:
            gmodsScripts.gmodsCONmds.nptSCmds()
            os.chdir(gmxDIR)
            flash('MDS Continuation was Successfully completed')
            return render_template("index.html", N=N, category='success')

    fcer = open("fconerror.txt", "a")
    if 'EMD' in cwdirall:
        emddir = os.listdir('EMD')
        if not 'equpmd.gro' in emddir:
            if not ('equpmd.tpr' in emddir and 'equpmd.cpt' in emddir):
                printNote("MDS continuation will begin with Equilibration Production run")
                response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
                if not (response.lower() == "yes" or response.lower() == "y"):
                    os.chdir(gmxDIR)
                    flash("MDS Continuation aborted")
                    fcer.close()
                    return render_template("selectmdstype.html", N=N, category='warning')
                else:
                    gmodsScripts.gmodsCONmds.emdSCmds()
                    os.chdir(gmxDIR)
                    flash('MDS Continuation was Successfully completed')
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
            else:
                printNote("Its appears Equilibration Production run was interrupted")
                response = input("To continue from where its ended, type YES/y. Otherwise press ENTER to continue: ")
                if not (response.lower() == "yes" or response.lower() == "y"):
                    os.chdir(gmxDIR)
                    flash("MDS Continuation aborted")
                    fcer.close()
                    return render_template("selectmdstype.html", N=N, category='warning')
                else:
                    print("Continuation assume interruption in previous run")
                    os.chdir('EMD')
                    print("Executing: Continuation Equilibration Production MD Simulation...")
                    try:
                        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equpmd', '-cpi', 'equpmd.cpt'], check=True, text=True, stderr=subprocess.STDOUT, stdout=fcer)
                    except subprocess.CalledProcessError as e:
                        print(e)
                        fcer.close()
                        fcer.close()
                        raise Exception("Something was not right with EMD continuation. Please check")

                    printNote("Interrupted EMD was successfuly completed")
                    os.chdir('../')

                    gmodsScripts.gmodsCONmds.pmdSCmds()
                    os.chdir(gmxDIR)
                    flash('MDS Continuation was Successfully completed')
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
        else:
            print("Checking your EMD folder for required files ....")
            if not 'equpmd.cpt' in emddir:
                printWarning("equpmd.cpt file is missing in EMD directory. Please check and correct")
                os.chdir(gmxDIR)
                flash("MDS Continuation aborted. Missing equpmd.cpt file")
                fcer.close()
                return render_template("selectmdstype.html", N=N, category='warning')

            if "extend.tpr" in emddir:
                print("Please backup any existing PMD folder if you want to extend equilibration step")
                response = input("To extend Equilibration MDS, type YES/y. Otherwise press ENTER to continue: ")
                if not (response.lower() == "yes" or response.lower() == "y"):
                    printNote("Noted: Equilibration production was previously completed. Nothing to do")
                else:
                    print("Continuation assume extension of previous run")
                    os.chdir('EMD')
                    print("Executing: Extension of Equilibration Production MD Simulation...")
                    try:
                        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equpmd', '-s', 'extend.tpr', '-cpi', 'equpmd.cpt'], check=True, text=True, stderr=subprocess.STDOUT, stdout=fcer)
                    except subprocess.CalledProcessError as e:
                        print(e)
                        fcer.close()
                        fcer.close()
                        raise Exception("Something was not right with EMD continuation. Please check")

                    printNote("EMD Extension was successfuly completed")
                    os.chdir('../')

                    gmodsScripts.gmodsCONmds.pmdSCmds()
                    os.chdir(gmxDIR)
                    flash('MDS Continuation was Successfully completed')
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
            else:
                printNote("Noted: Equilibration production was previously completed. Nothing to do")
    else:
        printNote("MDS continuation will begin with Equilibration Production run")
        response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            os.chdir(gmxDIR)
            flash("MDS Continuation aborted")
            fcer.close()
            return render_template("selectmdstype.html", N=N, category='warning')
        else:
            gmodsScripts.gmodsCONmds.emdSCmds()
            os.chdir(gmxDIR)
            flash('MDS Continuation was Successfully completed')
            fcer.close()
            return render_template("index.html", N=N, category='success')

    if 'PMD' in cwdirall:
        pmddir = os.listdir('PMD')
        if not 'pmds.gro' in pmddir:
            if not ('pmds.tpr' in pmddir and 'pmds.cpt' in pmddir):
                printNote("MDS continuation will begin with Production MD run")
                response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
                if not (response.lower() == "yes" or response.lower() == "y"):
                    os.chdir(gmxDIR)
                    flash("MDS Continuation aborted")
                    fcer.close()
                    return render_template("selectmdstype.html", N=N, category='warning')
                else:
                    gmodsScripts.gmodsCONmds.pmdSCmds()
                    os.chdir(gmxDIR)
                    flash('PMD Continuation was Successfully completed')
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
            else:
                printNote("Its appears Production MDS run was interrupted")
                response = input("To continue from where its ended, type YES/y. Otherwise press ENTER to abort: ")
                if not (response.lower() == "yes" or response.lower() == "y"):
                    os.chdir(gmxDIR)
                    flash("MDS Continuation aborted")
                    fcer.close()
                    return render_template("selectmdstype.html", N=N, category='warning')
                else:
                    print("Continuation assume interruption in previous run")
                    os.chdir('PMD')
                    print("Executing: Continuation Production MD Simulation...")
                    try:
                        subprocess.run(['gmx', 'mdrun', '-deffnm', 'pmds', '-cpi', 'pmds.cpt'], check=True, text=True, stderr=subprocess.STDOUT, stdout=fcer)
                    except subprocess.CalledProcessError as e:
                        print(e)
                        fcer.close()
                        fcer.close()
                        raise Exception("Something was not right with Production MDS continuation. Please check")

                    printNote("Interrupted Production MDS was successfuly completed")
                    os.chdir(gmxDIR)
                    flash('MDS Continuation was Successfully completed')
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
        else:
            print("Checking your PMD folder for required files ....")
            if not 'pmds.cpt' in pmddir:
                printWarning("pmds.cpt file is missing in the PMD directory. Please check and correct")
                os.chdir(gmxDIR)
                flash('Production MDS Continuation aborted. Missing pmds.cpt file')
                fcer.close()
                return render_template("selectmdstype.html", N=N, category='warning')

            if "extend.tpr" in pmddir:
                printNote("Continuation assume extension of previous run")
                response = input("To continue from where its ended, type YES/y. Otherwise press ENTER to end: ")
                if not (response.lower() == "yes" or response.lower() == "y"):
                    os.chdir(gmxDIR)
                    flash("MDS process already successfuly completed")
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
                else:
                    os.chdir('PMD')
                    print("Executing: Equilibration Production MD Simulation...")
                    try:
                        subprocess.run(['gmx', 'mdrun', '-deffnm', 'pmds', '-s', 'extend.tpr', '-cpi', 'pmds.cpt'], check=True, text=True, stderr=subprocess.STDOUT, stdout=fcer)
                    except subprocess.CalledProcessError as e:
                        print(e)
                        fcer.close()
                        raise Exception("Something was not right with PMD continuation. Please check")

                    os.chdir(gmxDIR)
                    flash("Production MDS Extension was successfuly completed")
                    fcer.close()
                    return render_template("index.html", N=N, category='success')
            else:
                os.chdir(gmxDIR)
                flash("MDS process already successfuly completed")
                fcer.close()
                return render_template("index.html", N=N, category='success')
    else:
        printNote("MDS continuation will proceed with Production MD run")
        response = input("To continue, type YES/y. Otherwise press ENTER to end: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            os.chdir(gmxDIR)
            flash("MDS Continuation aborted")
            fcer.close()
            return render_template("selectmdstype.html", N=N, category='warning')
        else:
            gmodsScripts.gmodsCONmds.pmdSCmds()
            os.chdir(gmxDIR)
            flash('MDS Continuation was Successfully completed')
            fcer.close()
            return render_template("index.html", N=N, category='success')

if __name__ == "__main__":
    # app.run() for debug
    from waitress import serve
    gmods = ui.run()
    serve(gmods)
