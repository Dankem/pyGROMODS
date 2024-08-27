#!/usr/bin/env python

"""
    pyGROMODS-v2024.02 Release

          <<<  NO WARRANTY AT ALL!!!  >>>

	Daniyan, Oluwatoyin Michael, B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) Biochemistry
    Department of Pharmacology, Faculty of Pharmacy
    Obafemi Awolowo Universiy, Ile-Ife, Nigeria.
    >>http://www.oauife.edu.ng<<

    mdaniyan@oauife.edu.ng; toyinpharm@gmail.com
"""
import sys

if sys.version_info < (3, 5):
	raise Exception("Python 3.5 or a more recent version is required.")

import os
import subprocess
from pathlib import Path
import time
import shutil
import random
import string
from datetime import datetime
from tkinter import Tk, messagebox

from gmodsScripts.gmodsHelpers import ligtopol, receptopol, topolsplit, indexoflines, complexgen, pdbcatogro, solvation, insertdetails, printWarning, printNote, tinput, defaults1, defaults2, udrestraint, gmxmdsFChecks, gmxmdsFClean, gmxmdsFEChecks

from gmodsScripts.gmodsTLptopol import TLtopol, tlpfinal
from gmodsScripts.gmodsTScheck import Checkligtop
from gmodsScripts.gmodsOPLStop import OPLStop, OPLSacpype

def RLmany(appDIR, gmxDIR, fdefaults):
	# Send out popup message to capture user's attention to a need for final setup
	messages = "WE SHALL CHECK YOUR UPLOADE FILES FOR CORRECTNESS. AT THE END, CONFIRMATION IS REQUIRED TO PROCEED"

	crp = Tk()
	crp.title("MDS INPUT FILES GENERATION!!!")
	crp.withdraw()
	crp.attributes("-topmost", True)
	crp.geometry("250x750")

	checkmessage = messagebox.askokcancel("IMPORTANT:", f"{messages} \n\n CLICK OK PLEASE")
	crp.destroy()
	if not checkmessage == True:
		raise Exception("Operation Interrupted by User. Make necessary corrections and restart")
	else:
		printNote("Please respond to questions as appropriate")
	time.sleep(5)
	print('\n')

	# Set some environment variables
	scriptDIR = appDIR
	GMX_MDS = gmxDIR

	print(f"User working directory set to: {GMX_MDS}")
	print('\n')

	# Set global variable to access user supplied topology file(s)
	TFF = 0
	Ligsff = " "
	tff = " "

	# Set global variable to automatically run pdb2gmx, editconf and others
	# forcefields & water (select), -bt (triclinic), -d (0.1), and timeout (60)
	defaults = fdefaults 
	
	printNote("Let us check again your selected default values ..... ")
	
	print(f"Default forcefield is {defaults[0]}")
	print(f"Default water model is {defaults[1]}")
	print(f"Default editconf -bt option is {defaults[2]}")
	print(f"Default editconf -d option is {defaults[3]}")
	print(f"Default timeout for input() request is {defaults[4]}")

	if defaults[5] == "A":
		print("Default mode for generating input file is Interractive")
		response = tinput("To revert to Noninteractive mode type YES/y: ", 30, "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[5] = "B"
			printNote("You have activated non-interactive mode")
			print("Your preferred forcefield and water model will be autodetected following your first interactive selection")
	else:
		print("Default mode for generating input file is Noninterractive")
		response = tinput("To revert back to Interactive mode type YES/y: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[0] = "select"
			defaults[1] = "select"
			defaults[5] = "A"
		else:
			defaults[5] = "C"
	print('\n')

	response = tinput("To adjust further the selected default values of -d, -bt and timeout type YES/y: ", defaults[4], "n")
	if response.lower() == "yes" or response.lower() == "y":
		defaults[2], defaults[3], defaults[4] = defaults1() 

	# Set some global parameter variables
	PLD = os.path.join(scriptDIR, 'gmodsTSF')
	PLDfiles = os.listdir(PLD)

	fSamples = os.path.join(scriptDIR, 'sampleFiles')
	fPDB = os.path.join(scriptDIR, 'Uploads')

	LIGfiles = os.listdir(os.path.join(fPDB, 'Ligands'))
	RECfiles = os.listdir(os.path.join(fPDB, 'Receptors'))
	TOPfiles = []

	# Check preferred forcefield selection, forcefield folder and its content
	fPDBfolders = os.listdir(fPDB)
	for itemf in fPDBfolders:
		if itemf == 'Gromos' or itemf == 'Opls' or itemf == 'Amber' or itemf == 'Charmm':
			Ligsff = itemf
			
	if not Ligsff == " ":
		printNote("PLEASE NOTE")
		print(f"Your preselected forcefield group is {Ligsff}")

	ligsff_dir = os.path.join(fPDB, Ligsff)
	if os.path.isdir(ligsff_dir):
		TOPfiles = os.listdir(ligsff_dir)
		if len(TOPfiles) > 0:
			TFF = 1
			printNote("Uploaded ligand(s) topology(ies) detected. Checks in progress ...")
		else:
			TFF = 0		
			printNote("No uploaded ligand(s) topology(ies) detected. Platform default will be used")
	else:
		TFF = 0
		printNote("No uploaded ligand(s) topology(ies) detected. Platform default will be used")

	# Check that number and format of uploaded ligand topologies are correct
	if TFF > 0:
		if not len(LIGfiles) == len(os.listdir(ligsff_dir)):
			printWarning("Number of supplied topology files does not match the number of ligands")
			printNote("By default, user uploaded ligand topologies will be ignored")
			response = input("To ignore uploaded ligand topology and continue, type YES/y, Otherwise press ENTER to abort: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Number of topology files does not match the number of ligands. Check and restart")
			else:
				TFF = 0

		else:
			ntopfile = 0
			ntopformat = 0
			for file in os.listdir(ligsff_dir):
				if not (Path(file).suffix == ".itp" or Path(file).suffix == ".top"):
					ntopformat += 1
				else:
					checkindex = indexoflines(os.path.join(fPDB, Ligsff, file))
					filecheck = open(os.path.join(fPDB, Ligsff, file), "r")
					readcheck = filecheck.readlines()
					checklist = ['atomtypes', 'moleculetype', 'system']
					nckl = 0
					for ckl in checklist:
						try:
							if not ckl in readcheck[checkindex[ckl]].split():
								nckl += 1
						except Exception as e:
							nckl += 1
					filecheck.close()
					if nckl > 0:
						ntopfile += 1

			if not ntopformat == 0:
				print(f"{ntopformat} Topology file(s) lacked required .itp or .top format")
				printNote("Default is to continue without the topology files")
				response = tinput("To abort and make correction, type YES/y: ", defaults[4], "n")
				if response.lower() == "yes" or response.lower() == "y":
					raise Exception("Sorry, unacceptable topology file format detected")
				else:
					TFF = 0

			elif not ntopfile == 0:
				print(f"{ntopfile} Topology file(s) lacked required ['atomtypes'], ['moleculetype'] and/or [ system ] subheading")
				print("Check and include the subheadings, with or without expected accompanied values")
				printNote("Default is to continue without the topology files")
				response = tinput("To abort and make correction, type YES/y: ", defaults[4], "n")
				if response.lower() == "yes" or response.lower() == "y":
					raise Exception("Checking ligand topology failed. Check, correc and restart")
				else:
					TFF = 0

	# Check the suitability of the uploaded files
	for file in LIGfiles:
		if not (Path(file).suffix == ".pdb" or Path(file).suffix == ".mol2"):
			raise Exception("Sorry, only pdb or mol2 file formats is acceptable. Setroute again and upload required files")

	for file in RECfiles:
		if not Path(file).suffix == ".pdb":
			raise Exception("Sorry, Receptor file must be pdb files. Setroute again and upload required files")

	if not (len(RECfiles) > 0 and len(LIGfiles) > 0):
		raise Exception("Needed ligands and/or receptor files is/are missing. Make sure required files are uploaded")

	if len(RECfiles) > 1:
		raise Exception("Sorry, RLmany route allows only one Receptor per ligand or set of ligands. You have " + str(len(RECfiles)) + " receptors. Try RLsingle to upload pairs of receptor and ligands")

	if len(LIGfiles) > 1:
		print("You have more than one ligand. Each ligand will form separate complex with receptor. Are you sure you want to continue?")
		response = tinput("Type YES/y to continue: ", defaults[4], "y")
		if not (response.lower() == "yes" or response.lower() == "y"):
			raise Exception("You have choosen to Abort!!!. Check and restart the process")

	# Check that appropraite tleap source file is available for use
	while True:
		PLDfiles = os.listdir(PLD)
		numberL = 1
		tlpresent = 0
		for nL in PLDfiles:
			tlsflist = nL.split("-")
			if str(numberL) in tlsflist:
				tlpresent += 1

		if tlpresent > 1:
			raise Exception("Multiple tleap source file with similar identifiers found. Check samples in", PLD)

		elif tlpresent < 1:
			printWarning("WARNING: You have no tleap source files with identifier that match number of uploaded ligands")
			print(f"Please refer to README.md file for instruction on how to prepare a suitable tleap source file, appropriate for your actual number of uploaded ligands, and saved appropriately in {PLD}")
			print("Alternatively, you may also reduce the number of uploaded ligands, to make use of relevant pre-installed tleap source file")
			printNote("No alternative GROMACS compatible tolopogy file will be generated")
			break

		elif tlpresent == 1:
			printNote("Found the appropriate tleap source file for your work")
			break
	print('\n')

	# Set the starting date and time
	printNote("IF YOU GET HERE, FILES CHECK WAS MOST PROBABLY SUCCESSFUL!!!")
	printNote(
		"CAREFULLY CHECK THE CHECK RESULTS ABOVE AND TYPE 'YES/y' TO CONFIRM OR 'NO/n' TO ABORT")
	confirmation = input("Response is: ")
	if not (confirmation.lower() == "yes" or confirmation.lower() == "y"):
		raise Exception("Process aborted by the User. Please rerun from SetRoute")
	else:
		rlmanytime = datetime.now()
		Tstart = rlmanytime.strftime("%B %d, %Y %H:%M:%S")
		print(f'Generation of MDS input files begins: {Tstart}')
	print('\n')

	# Get user imput for project name
	while True:
		Xname = tinput("Suppy a name for the current project: ", defaults[4], "RLmany")
		if Xname == " ":
			print("You must supply a name for the project")
			continue
		elif Xname.isalpha() == False:
			print("Supplied name must be all alphabets")
			continue
		else:
			break

	# Generate unique id number for the project
	while True:
		idnumber = random.randint(0, 9)
		idalpha1 = random.choice(string.ascii_uppercase)
		idalpha2 = random.choice(string.ascii_uppercase)
		ID = idalpha1 + str(idnumber) + idalpha2
		foldername = Xname + "_" + str(ID)
		if not os.path.isdir(foldername):
			print(f"Your Unique ID is: {ID}")
			print(f"Your current work will be stored in {foldername}")
			break
		else:
			continue

	# Create working folder using the generated name
	os.mkdir(foldername)
	os.chdir(foldername)

	work_dir = os.path.join(GMX_MDS, foldername)
	print(f"pyGROMODS current workspace directory set to {work_dir}")

	# Copy sample mdp and peptide files to the working directory
	os.mkdir('samples')
	os.chdir('samples')
	try:
		for sample in os.listdir(fSamples):
			sfolder = os.path.join(scriptDIR, 'sampleFiles', sample)
			if os.path.isdir(sfolder):
				os.system('cp -r ' + sfolder + ' ./')
			else:
				shutil.copy(sfolder, './')
		printNote("Samples of .mdp, peptide and Tleap source files have been copied to your workspace directory")
	except:
		pass

	os.chdir('../')
	print('\n')
	time.sleep(5)

	RLname = "RLmany_"
	npairs = 0
	mdsfgerrors = 0
	while npairs < len(LIGfiles):
		# Create host directory for each pair of receptor and ligand
		RLcount = npairs + 1
		lig = LIGfiles[npairs]
		M = TFF
		rls_dir = RLname + str(RLcount)

		if os.path.isdir(os.path.join(work_dir, rls_dir)):
			shutil.rmtree(os.path.join(work_dir, rls_dir))
			time.sleep(5)

		os.mkdir(rls_dir)
		os.chdir(rls_dir)

		workhost_dir = os.path.join(work_dir, rls_dir)
		print(f"Current project host directory set to {workhost_dir}")
		print('\n')

		# Adjust defaults values of -bt, -b and timeout if ff and water are 'select'
		if defaults[0] == "select" and defaults[1] == "select" and mdsfgerrors == 0:
			printNote("It appears you're working Interactively")
			print(f"Current default values for -bt is: {defaults[2]}")
			print(f"Current default values for -d is: {defaults[3]}")
			print(f"Current default values for timeout is: {defaults[4]}")
			response = tinput("To adjust these values for current protein - ligand complex, type YES/y: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				defaults[2], defaults[3], defaults[4] = defaults1() 

		elif defaults[0] == "select" and defaults[1] == "select" and mdsfgerrors > 0:
			printNote("This appears to be a rerun of failed MDS input files generation")
			time.sleep(5)

		###############################################
		# GENERATING PROTEIN TOPOLOGIES AND STRUCTURE #
		###############################################
		# Create Receptor directory and populate it with required files
		os.mkdir('Receptor')
		os.chdir('Receptor')

		receptor_dir = os.path.join(workhost_dir, 'Receptor')
		print(f"Protein data directory set to {receptor_dir}")

		printNote("Generating Protein topology and parameter files...")

		# Generate Protein topologies and associated parameters
		rname = "receptor"
		rep = RECfiles[0]
		selff = defaults[0]
		selwater = defaults[1]

		shutil.copy(os.path.join(fPDB, 'Receptors', rep), './')

		notproceed = 0
		while True:
			RFtop, RFpdb, RFposre = receptopol(rep, rname, selff, selwater)

			# Let us find out the forcefield choosen by the user for protein topology
			t = open(RFtop, "r")
			tread = t.readlines()
			tf = 1
			for tl in tread:
				if not tl.split() == []:
					if 'Include' in tl.split() and 'forcefield' in tl.split() and 'parameters' in tl.split():
						tln1 = tread[tf].split()
						tln2 = tln1[1].split('"')[1].split('/')
						for tfd in tln2:
							if Path(tfd).suffix == ".ff":
								tff = tfd
					else:
						tf += 1
				else:
					tf += 1
			t.close()

			if not (tff[0:4].capitalize() == Ligsff or tff[0:5].capitalize() == Ligsff or tff[0:6].capitalize() == Ligsff):
				print(f"{tff} forcefield in topol.top file does not match the preselected forcefield group: {Ligsff}")
				print(f"By default forcefield group will be changed to {tff} and uploaded ligand topology(ies) will be ignored")

				printNote("To rerun, Type YES/y. To continue with current selection press ENTER")
				response = tinput("Response: ", defaults[4], "n")
				if response.lower() == "yes" or response.lower() == "y":
					selff = "select"
					selwater = "select"
					continue
				else:
					if tff[0:4].capitalize() == "Opls":
						Ligsff = "Opls"
					elif tff[0:5].capitalize() == "Amber":
						Ligsff = "Amber"
					elif tff[0:6].capitalize() == "Gromos":
						Ligsff = "Gromos"
					elif tff[0:6].capitalize() == "Charmm":
						Ligsff = "Charmm"
					else:
						print(f"{tff} does not match any known forcefield group")
						printNote("If using a self created or modified forcefield, adopt standard naming convention")
						response = input("Type YES/y to select new forcefield interactively. Otherwise press ENTER: ")
						if not (response.lower() == "yes" or response.lower() == "y"):
							printNote("Generation of MDS input files for the current pairs will not proceed")
							notproceed += 1
							break
						else:
							selff = "select"
							selwater = "select"
							continue

					Ufolders = os.listdir(fPDB)
					if not ('Amber' in Ufolders or 'Charmm' in Ufolders or 'Gromos' in Ufolders or 'Opls' in Ufolders):
						notproceed += 1
						break
					else:
						for sff in Ufolders:
							if sff == 'Amber' or sff == 'Charmm' or sff == 'Gromos' or sff == 'Opls':
								os.rename(os.path.join(fPDB, sff), os.path.join(fPDB, Ligsff))

					print(f"Your preselected forcefield group has been changed to {Ligsff}")
					ligsff_dir = os.path.join(fPDB, Ligsff)
					if TFF > 0:
						print(f"Subdirectory for uploaded ligand topology is now: {ligsff_dir}")
						printNote("To use with uploaded ligand topology, type YES/y. Otherwise press ENTER to ignore")
						response = tinput("Response: ", defaults[4], "y")
						if not (response.lower() == "yes" or response.lower() == "y"):
							TFF = 0
							break
						else:
							TFF = 1
							break
					else:
						break
			else:
				print(f"Your forcefiled as contained in topol.top file is {tff}")
				break

		if not notproceed == 0:
			os.chdir('../../')
			npairs += 1
			TFF = M
			continue

		try:
			shutil.copy(rep, 'xreceptor.pdb')
			shutil.copy(RFtop, 'nReceptor.top')
			shutil.copy(RFpdb, 'nReceptor.pdb')
		except Exception as e:
			printWarning(e)
			pass

		# Check if regenerating and decide either to continue interactively or not
		if mdsfgerrors > 0:
			print('\n')
			printNote("This is a rerun of failed MDS input files generation with Interactive mode")
			print("Type YES/y to continue interactviely, Or press ENTER to change to Noninteractive mode")
			mode = input("Response is: ")
			if not (mode.lower() == "yes" or mode.lower() == "y"):
				print("We are changing to Noninteractive mode")
				defaults[5] = "B"
			else:
				defaults[5] = "A"
				defaults[0] = "select"
				defaults[1] = "select"
				mdsfgerrors = 0

		# If needed, we now change the default parameters and lock the new default mode
		if not (defaults[0] == "select" and defaults[1] == "select" and mdsfgerrors == 0):
			if defaults[5] == "B":
				if Path(tff).suffix == ".ff":
					defaults[0] = Path(tff).stem
				else:
					defaults[0] = tff
				
				defaults[1] = defaults2(RFtop)
				if defaults[1] == "none":
					print("No water model was detected for your system")
				else:
					print(f"Your water model as contained in topol.top file is {defaults[1]}")
				
				printNote("Your selected default values are as follows: ")
				print(f"--------Default forcefield: {defaults[0]}")
				print(f"--------Default water model: {defaults[1]}")
				print(f"--------Default editconf -bt: {defaults[2]}")
				print(f"--------Default editconf -d: {defaults[3]}")
				print(f"--------Default input timeout: {defaults[4]}")
				print("---------Default mode: non-interactive")
				defaults[5] = "C"
			else:
				defaults[1] = defaults2(RFtop)
				if defaults[1] == "none":
					print("No water model was detected for your system")
				else:
					print(f"Your water model as contained in topol.top file is {defaults[1]}")
			mdsfgerrors = 0

		# Checking for posre.itp file or create one to include relevant posre files
		print("Updating topology file...")
		if not "posre.itp" in os.listdir():
			# Create new posre file to accomodate all generated chains posres
			posrefileopen = open("posre.itp", "+a")
			
			# Add headers to show it was not created by gromacs			
			posreheaders = ["; This file was generated to captured all pdb2gmx generated posre files", "; This file is needed when multiple protein chains are involved, each having its own posre file", "; The file should be checked for correctness and used with caution, and may be edited as approapriate if needed", "; Also, details here can be transfered to topol.top directly under '#ifdef POSRE' and ignore this file", "; Note that pdb2gmx prepared posre for each chain and included '#ifdef POSRE' in chain topology that has now been removed", "; A backup of the original chain topology files can be found in Receptor and gmxmds/backup folders", "; To use the backup chains topology files instead, comment in '; ' or delete the lines below"]
			for pheader in posreheaders:
				posrefileopen.write(pheader)
				posrefileopen.write('\n')
			posrefileopen.write('\n')

			# Check for and include all generated posres file here, make backup copy
			for posrefile in os.listdir():
				if posrefile[0:5].lower() == "posre" and not posrefile == "posre.itp":
					includeposre = '#include' + ' ' + '"' + posrefile + '"'
					posrefileopen.write(includeposre)
					posrefileopen.write('\n')
			posrefileopen.write('\n')
			posrefileopen.close()

			# Create backup copy of the chain topology files and remove 'Include Position restraint file' statement
			for chain in os.listdir():
				if Path(chain).suffix == ".itp" and not chain[0:5] == "posre":
					bkposre = "bk" + chain
					shutil.copy(chain, bkposre)

					ciindex1 = 0
					readL = 0
					
					with open(chain, "r") as pfrc:
						chainreadlines = pfrc.readlines()
						readL = len(chainreadlines) - 1

						for ciline in chainreadlines:
							if "include" in ciline.lower().split() and "position" in ciline.lower().split() and "restraint" in ciline.lower().split():
								break
							else:
								ciindex1 += 1

					ciindex2 = ciindex1 + 4
					if ciindex1 < readL and not ciindex2 > readL:
						with open(chain, "r") as pfr:
							chainlines = pfr.readlines()

						with open(chain, "w") as pfw:
							for pline in chainlines:
								if not pline in chainlines[ciindex1:ciindex2]:
									pfw.write(pline)

		# Checking for inclusion of restraint files
		nrecopen = open("nReceptor.top", "r")
		nrecreadlines = nrecopen.readlines()
		includelines = []
		for nline in nrecreadlines:
			xnline = []
			for xnl in nline.split():
				xnline.append(xnl.lower())

			if "#ifdef" in xnline and "posres" in xnline:
				includelines.append("#ifdef")
				includelines.append("posres")
		nrecopen.close()

		if "#ifdef" in includelines and "posres" in includelines:
			shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')
			insertUDP = "; Include water topology"
			insertdetails('nReceptor.top', 'mt-topol2.itp', insertUDP)
			os.remove('mt-topol2.itp')
		else:
			shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol3.itp'), './')
			insertP = "; Include water topology"
			insertdetails('nReceptor.top', 'mt-topol3.itp', insertP)
			os.remove('mt-topol3.itp')

		os.chdir('../')
		time.sleep(10)

		##########################################################
		# DETERMINE ROUTE FOR OPLS LIGANDS TOPOLOGIES GENERATION #
		##########################################################
		# Determine and choose preferred route for platform generated opls ligand topology
		opls_route = 0
		if Ligsff == 'Opls':
			print('\n')
			printNote("The following options are available to generate OPLS compatible topology:")
			print(">>>>> 1. Using platform opls compatible ligand topology generation - RECOMMENDED")
			print(">>>>> 2. Using acpype opls compatible ligand topology generation")
			print(">>>>> 3. Using platform default - ignore opls compatibility - NOT RECOMMENDED")
			while True:
				response = tinput("Choose your preferred option: ", defaults[4], "1")
				if response == '1' or response == '2' or response == '3':
					print(f"Option {response} Selected!")
					printNote("Type YES/y to confirm or press ENTER to choose a different option")
					confirm = tinput("Response YES/y or press ENTER: ", defaults[4], "y")
					if not (confirm.lower() == "yes" or confirm.lower() == "y"):
						continue
					else:
						if response == '1':
							printNote("Confirmed: The platform opls compatible ligand topology generation") 
						elif response == '2':
							printNote("Confirmed: The acpype opls compatible ligand topology generation") 
						elif response == '3':
							printNote("Confirmed: The platform default - ignore opls compatibility") 
						opls_route = int(response)
						break
				else:
					print("Wrong selection detected. Please type 1, 2 or 3")
					continue

			printNote("Other options will be attempted if your preferred option fails")
		print('\n')

		#################################
		# GENERATING LIGANDS TOPOLOGIES #
		#################################

		# Create Ligand directory and populate it with required files
		printNote("Generating topology and parameter files for ligands...")
		time.sleep(2)

		os.mkdir('Ligands')
		os.chdir('Ligands')

		lig_dir = os.path.join(workhost_dir, 'Ligands')
		print(f"Ligands data directory set to {lig_dir}")

		Lname = "LIG"
		name = Lname + str(RLcount)

		shutil.copy(os.path.join(PLD, 'ligtoptleap.in'), './')
		shutil.copy(os.path.join(fPDB, 'Ligands', lig), './')

		LFgro, LFtop = ligtopol(lig, name)

		# Identify the pdb file to use subsequently, depending on uploaded ligand format
		print("Generating optimized ligand structure...")
		ligE = open('ligcleanerror.txt', 'a')

		nligname = name + "up.pdb"
		if Path(lig).suffix == ".pdb":
			try:
				subprocess.run(['gmx', 'editconf', '-f', lig, '-o', nligname], check=True, stderr=subprocess.STDOUT, stdout=ligE, text=True)
			except subprocess.SubprocessError as e:
				printNote("Trying again...")
				shutil.copy(lig, nligname)

		elif Path(lig).suffix == ".mol2":
			ligtopol_pdb = name + "new" + ".pdb"
			try:
				subprocess.run(['gmx', 'editconf', '-f', ligtopol_pdb, '-o', nligname], check=True, stderr=subprocess.STDOUT, stdout=ligE, text=True)
			except subprocess.SubprocessError as e:
				printNote("Trying again...")
				shutil.copy(ligtopol_pdb, nligname)
		ligE.close()

		# Determine the molecule identifier
		ligsmnfile = name + "_mn.itp"
		file_mn = open(ligsmnfile, "+a")
		LFtopfile = open(LFtop, "r")
		readLFtopfile = LFtopfile.readlines()
		lastline = len(readLFtopfile) - 1
		while lastline < len(readLFtopfile):
			if readLFtopfile[lastline].split() == "":
				lastline += -1
				continue
			else:
				file_mn.write(readLFtopfile[lastline])
				break
		file_mn.close()
		LFtopfile.close()

		# Prepare atomtypes and moleculetype files for user supplied topology files
		if TFF > 0:
			os.mkdir('usertops')
			os.chdir('usertops')

			u = RLcount - 1
			utop = TOPfiles[u]
			shutil.copy(os.path.join(fPDB, Ligsff, utop), './')
			ulig = nligname
			shutil.copy(os.path.join(lig_dir, ulig), './')

			utindex = indexoflines(utop)
			uml = utindex['moleculetype'] + 2

			utopopen = open(utop, "r")
			utopread = utopopen.readlines()
			umoltype = utopread[uml].split()[0]

			# Check if the supplied topology file match the ligand pdb
			utopopen.seek(0)
			ml = 0
			for mline in utopread:
				try:
					if umoltype in mline.split():
						ml += 1
				except:
					pass

			uligopen = open(ulig, "r")
			uligread = uligopen.readlines()
			ul = 0
			for uline in uligread:
				try:
					if umoltype in uline.split():
						ul += 1
				except:
					pass
			
			utopopen.close()
			uligopen.close()

			time.sleep(5)

			if ul > ml:
				dul_ml = ul - ml
			else:
				dul_ml = ml - ul

			if not (dul_ml > ul or dul_ml > ml):
				print(f"{utop} the user supplied topology match the ligand named {ulig}")
				time.sleep(5)

				# Check the file for duplicate atom types with forcefield standard atomtypes
				uligcheck = Checkligtop(utop, tff)
				uctindex = indexoflines(uligcheck)

				# Generate the uLIGS_at.itp (atomtypes) and uLIGS_mt.itp (moleculetypes) file
				uTname = "u" + name
				uligand_at, uligand_mt = topolsplit(uligcheck, uTname, uctindex)
				utopname = uTname + ".top"
				os.rename(uligcheck, utopname)
				shutil.copy(utopname, '../')
				shutil.copy(uligand_at, '../')
				shutil.copy(uligand_mt, '../')

				os.chdir('../')

			else:
				print(f"{utop} user supplied topology does not match the ligand named {ulig}")
				print("By default, this topology file will be ignored for this pair")
				print("If you're sure the topology is correct, type YES/y to continue with it")
				response = tinput("Response YES/y or press ENTER: ", defaults[4], "n")
				if not (response.lower() == "yes" or response.lower() == "y"):
					TFF = 0
					os.chdir('../')
				else:
					# Check the file for duplicate atom types with forcefield standard atomtypes
					uligcheck = Checkligtop(utop, tff)
					uctindex = indexoflines(uligcheck)

					# Generate the uLIGS_at.itp (atomtypes) and uLIGS_mt.itp (moleculetypes) file
					uTname = "u" + name
					uligand_at, uligand_mt = topolsplit(uligcheck, uTname, uctindex)
					utopname = uTname + ".top"
					os.rename(uligcheck, utopname)
					shutil.copy(utopname, '../')
					shutil.copy(uligand_at, '../')
					shutil.copy(uligand_mt, '../')

					os.chdir('../')

		# If need be, generate opls compatible ligand topology using the selected approach above
		oplstop = " "

		if opls_route == 1:
			mergename = "opls" + name + ".top"
			oplstop = OPLStop(LFtop, tff, name)

			if not oplstop == mergename:
				printWarning("Platform OPLS compatible ligand topology generation failed")
				printNote("Trying alternative option...")
				oplstop = OPLSacpype(LFtop, lig, name)
				if not oplstop == mergename:
					printWarning("Alternative generation of OPLS compatible ligand topology failed")
					printNote("Platform default, ignoring OPLS compatibility, will now be used")
					oplstop = LFtop		

		elif opls_route == 2:
			mergename = "opls" + name + ".top"
			oplstop = OPLSacpype(LFtop, lig, name)

			if not oplstop == mergename:
				printWarning("Alternative OPLS compatible ligand topology generation failed")
				printNote("Trying default Platform OPLS compatible option ...")
				oplstop = OPLStop(LFtop, tff, name)
				if not oplstop == mergename:
					printWarning("Platform OPLS compatible ligand topology generation failed")
					printNote("Platform default, ignoring OPLS compatibility, will now be used")
					oplstop = LFtop		

		else:
			oplstop = LFtop

		# Check the generated ligand topology file for duplicate atom types with forcefield standard atomtypes
		ligcheck = Checkligtop(oplstop, tff)

		nametop = " "
		if not ligcheck == oplstop:
			nametop = "ck" + name + ".top"
			os.rename(ligcheck, nametop)
		else:
			nametop = oplstop

		findex = indexoflines(nametop)
		ligand_at, ligand_mt = topolsplit(nametop, name, findex)
	
		# Check ligand directory for posre file and back it up
		for posrefile in os.listdir():
			if posrefile[0:5] == "posre":
				posrebk = posrefile + ".bk"
				os.rename(posrefile, posrebk)

		os.chdir('../')
		print('\n')

		###########################################################
		# GENERATING PROTEIN-LIG COMPLEX TOPOLOGIES AND STRUCTURE #
		###########################################################
		# Create Complex directory and populate it with required files
		os.mkdir('Complex')
		os.chdir('Complex')

		complex_dir = os.path.join(workhost_dir, 'Complex')
		print(f"Protein-Ligand Complex directory set to {complex_dir}")
		time.sleep(2)

		cname = Lname + str(RLcount)
		lp_gro = cname + ".gro"
		lp_mol2 = cname + ".mol2"
		lp_prm = cname + ".prmtop"
		lp_inp = cname + ".inpcrd"
		lp_frc = cname + ".frcmod"
		lp_top = cname + ".top"
		LIG_at = cname + "_at.itp"
		LIG_mt = cname + "_mt.itp"
		LIG_mn = cname + "_mn.itp"

		if TFF > 0:
			ulp_top = "u" + cname + ".top"
			uLIG_at = "u" + cname + "_at.itp"
			uLIG_mt = "u" + cname + "_mt.itp"

		# Gather needed files to generate complex topologies
		listLdir = os.listdir(os.path.join(workhost_dir, 'Ligands'))
		for L in listLdir:
			if (Path(L).suffix == ".itp" or Path(L).suffix == ".gro" or Path(L).suffix == ".mol2" or Path(L).suffix == ".frcmod"):
				shutil.copy(os.path.join(workhost_dir, 'Ligands', L), './')

		listRdir = os.listdir(os.path.join(workhost_dir, 'Receptor'))
		for R in listRdir:
			if R == 'xreceptor.pdb' or R == 'nReceptor.pdb' or R == 'nReceptor.top':
				shutil.copy(os.path.join(workhost_dir, 'Receptor', R), './')

			elif Path(R).suffix == ".itp":
				shutil.copy(os.path.join(workhost_dir, 'Receptor', R), './')

		shutil.copy(os.path.join(PLD, 'tleap-1-ligs.in'), './')
		shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'at-topol.itp'), './')
		shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol.itp'), './')

		try:
			os.rename(lp_mol2, "xLIG1.mol2")
			os.rename(lp_frc, "yLIG1.frcmod")
			os.rename(LIG_at, "LIGS_at.itp")
			os.rename(LIG_mt, "LIGS_mt.itp")
			os.rename(LIG_mn, "LIGS_mn.itp")
		except Exception as e:
			printWarning("Something is not right. Check the error below:")
			printWarning(e)
			time.sleep(5)

		if Ligsff == "Opls":
			shutil.copy("LIGS_at.itp", "oplsLIGS_at.itp")

		if TFF > 0:
			try:
				os.rename(uLIG_at, "uLIGS_at.itp")
				os.rename(uLIG_mt, "uLIGS_mt.itp")
			except Exception as e:
				printWarning("Something is not right. Check the error below:")
				printWarning(e)
				time.sleep(5)

		# Now we shall insert details into receptor topol file, if correctly generated
		listCpldir = os.listdir()
		if 'nReceptor.top' in listCpldir:
			print(f"Generating complex topologies and parameters for Receptor and {cname} ....")
			shutil.copy('nReceptor.top', 'topol.top')
			insertAT = "[ moleculetype ]"
			insertdetails('topol.top', 'at-topol.itp', insertAT)

			insertMT = "; Include water topology"
			insertdetails('topol.top', 'mt-topol.itp', insertMT)

			topolfile = open("topol.top", "a")
			ligmnfile = open("LIGS_mn.itp", "r")
			readtops = ligmnfile.read()
			topolfile.write(readtops)
			topolfile.close()
			ligmnfile.close()

			if TFF > 0:
				shutil.copy('nReceptor.top', 'utopol.top')
				uinsertAT = "[ moleculetype ]"
				insertdetails('utopol.top', 'at-topol.itp', uinsertAT)

				uinsertMT = "; Include water topology"
				insertdetails('utopol.top', 'mt-topol.itp', uinsertMT)

				utopolfile = open("utopol.top", "a")
				uligmnfile = open("LIGS_mn.itp", "r")
				ureadtops = uligmnfile.read()
				utopolfile.write(ureadtops)
				utopolfile.close()
				uligmnfile.close()

		# Make directory and copy lig pdb into it for pdbcatogro() to access for alternative complex generation
		os.mkdir('ligpdb')
		os.chdir('ligpdb')

		tlplig = cname + "up.pdb"
		try:
			shutil.copy(os.path.join(workhost_dir, 'Ligands', tlplig), './')
		except Exception as e:
			printWarning(e)
			printWarning("The required ligand file was not found")
		os.chdir('../')

		# Check if tleap source file is present and attempt to generate complex via tleap or use alternative
		complexdir = os.listdir()
		tleapfile = "tleap-1-ligs.in"
		if not tleapfile in complexdir:
			printNote("No matching tLeap input file is found for your number of ligands")
			printNote("As such, alternative GROMACS compatible topology file can't be generated")

			print('\n')
			print("Trying alterenative complex generation...")

			catcomplx_gro, catcomplx_pdb = pdbcatogro()

			if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
				print("Alternative approach failed to generate complex")
				printWarning("Process aborted for this pair. Make necessary corrections and rerun")
				os.chdir('../../')
				npairs += 1
				TFF = M
				continue
	
			else:
				os.rename(catcomplx_gro, 'catComplex.gro')
				os.rename(catcomplx_pdb, 'catComplex.pdb')
				print("Complex generation with alternative approach was successful")

		else:
			complx_gro, complx_top = complexgen(tleapfile)

			if not (complx_gro == 'Complex.gro' and complx_top == 'Complex.top'):
				printWarning("Generation of complex structure and topologies failed")
				print("Please check error file for details, and make necessary corrections")

				print("Trying alternative complex generation...")

				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach also failed to generate complex")
					printWarning("Process aborted for this pair. Make necessary corrections and rerun")
					os.chdir('../../')
					npairs += 1
					TFF = M
					continue

				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')
					print("Complex generation with alternative approach was successful")

			else:
				os.rename(complx_gro, 'tlpComplex.gro')
				os.rename(complx_top, 'tlpComplex.top')
				print("Complex generation with tleap was successful")

				print("Generating a backup complex....")
				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach failed to generate complex to serve as fallback")
	
				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')
					print("Backup Complex generation was successful")

		try:
			os.remove("LIGS_mn.itp")
			os.remove("uLIGS_mn.itp")
		except:
			pass

		os.chdir('../')
		print('\n')

	    ####################################
	    # SOLVATION OF PROTEIN-LIG COMPLEX #
	    ####################################
		# Create solvation directory and populate it with required files
		print(f"Generating complex Solvation for Receptor and {cname} ....")
		time.sleep(2)

		solname = "Solvation"
		os.mkdir(solname)
		os.chdir(solname)

		listCdir = os.listdir(os.path.join(workhost_dir, 'Complex'))
		for C in listCdir:
			if Path(C).suffix == ".itp" or Path(C).suffix == ".top" or Path(C).suffix == ".gro" or Path(C).suffix == ".pdb":
				shutil.copy(os.path.join(workhost_dir, 'Complex', C), './')
	
		try:
			os.remove('mt-topol.itp')
			os.remove('Receptor_new.pdb')
		except:
			pass

		if not "at-topol.itp" in os.listdir():
			shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'at-topol.itp'), './')

		# Determine if the forcefield belong to the selected group - amber, charmm, gromos and opls
		if TFF > 0:
			if tff[0:4].lower() == Ligsff.lower() or tff[0:5].lower() == Ligsff.lower() or tff[0:6].lower() == Ligsff.lower():
				printNote("Platform generated and user uploaded ligand topologies are available")
				printNote("PLEASE NOTE THAT THE DEFAULT ORDER ARE:")
				print("------1. User uploaded ligand topologies will be attempted first")
				print("------2. Platform generated ligand topologies will serve as fallback")
				print('\n')
				printNote("To reverse the order, type YES/y, Otherwise press ENTER to continue")
				response = tinput("Response: ", defaults[4], "n")
				if not (response.lower() == "yes" or response.lower() == "y"):
					os.rename('LIGS_at.itp', 'bkLIGS_at.itp')
					os.rename('LIGS_mt.itp', 'bkLIGS_mt.itp')
					os.rename('topol.top', 'bktopol.top')
					os.rename('uLIGS_at.itp', 'LIGS_at.itp')
					os.rename('uLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('utopol.top', 'topol.top')

				else:
					print("The order has been reversed")
					os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
					os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
					os.rename('utopol.top', 'bktopol.top')

			else:
				print(f"{tff} in topol.top does not match the selected forcefield group, {Ligsff}")
				printNote("Solvation will most likely failed. Please check")
				os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('utopol.top', 'bktopol.top')

		else:
			if tff[0:5].lower() == 'amber' or tff[0:6].lower() == 'charmm' or tff[0:6].lower() == 'gromos' or tff[0:4].lower() == 'opls':
				print(f"{tff} is a GROMACS compatible forcefield")
				printNote("Default auto-generated ligand topology will be used") 
				printNote("For Gromos and OPLS, uploading compatible ligand topology is recommended") 
			else:
				print(f"{tff} in topol.top does not match any known forcefield group in GROMACS")

		# Time to prepared solvated complex
		if not defaults[1] == "select":
			selwater = defaults[1]
		else:
			selwater = defaults2("topol.top")
		selbt = defaults[2]
		seld = defaults[3]

		chkSoldir = os.listdir()
		if ('tlpComplex.gro' in chkSoldir and 'tlpComplex.top' in chkSoldir):

			# Prepare a version of amber tleap generated topology for suitable use with Gromacs
			print("Preparing tleap topology file, for possible use with Gromacs...")
			tlpgmxtopol = TLtopol('nReceptor.top', 'tlpComplex.top', tff)
			if not tlpgmxtopol == "tlptopol.top":
				printNote("Returning tleap topology as originally generated")
			else:
				shutil.copy(tlpgmxtopol, 'xtlptopol.top')
				if "posre_complex.itp" in os.listdir():
					os.rename("posre_complex.itp", "posre_tlptopol.itp")

			# Now it's time to solvate and add ions
			printNote("Solvation of complex in progress.....")
			grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)
			
			if not grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir:
				printNote("Solvation of complex unsuccessful")

				printNote("Trying with backup complex.....")
				grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro' and TFF > 0:
					printNote("Solvation failed with backup complex")

					print("Reverting to the backup topology files...")

					os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
					os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
					os.rename('topol.top', '#topol.top.bk#')
					os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
					os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('bktopol.top', 'topol.top')
					if "tlptopol.top" in os.listdir() and "xtlptopol.top" in os.listdir():
						os.rename('tlptopol.top', '#tlptopol.top.bk#')
						shutil.copy('xtlptopol.top', 'tlptopol.top')
	
					printNote("Repeating Complex Solvation ...")
					grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

					if not grosolvated == 'fsolvated.gro':
						printNote("Solvation of complex failed again")

						printNote("Trying Solvation with backup complex ...")
						grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

						if not grosolvated == 'fsolvated.gro':
							printNote("Solvation using backup complex unsuccessful")
							print("Solvation can not continue. You may wish to do solvation manually")

						else:
							printNote("Solvation of with backup complex and topologies was successful")

					else:
						printNote("Solvation of complex with backup topologies was successful")

				elif not grosolvated == 'fsolvated.gro' and TFF == 0:
					printNote("Solvation with backup complex unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")

				else:
					printNote("Solvation with backup complex was successful")

			elif not (grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir):
				if TFF > 0:
					os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
					os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
					os.rename('topol.top', '#topol.top.bk#')
					os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
					os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('bktopol.top', 'topol.top')
					if "tlptopol.top" in os.listdir() and "xtlptopol.top" in os.listdir():
						os.rename('tlptopol.top', '#tlptopol.top.bk#')
						shutil.copy('xtlptopol.top', 'tlptopol.top')
	
					printNote("Solvation of complex in progress.....")
					grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

					if not grosolvated == 'fsolvated.gro':
						printNote("Solvation of complex failed. Please check error file")
						print("Solvation can not continue. You may wish to do solvation manually")

					else:
						printNote("Solvation of complex was successful")

				else:
					printNote("Solvation of complex unsuccessful. Check error file")
					print("Solvation can not continue. You may wish to do solvation manually")
				
			else:
				printNote("Solvation of complex was successful")

		elif 'catComplex.gro' in chkSoldir:
			printNote("Solvation with backup complex in progress.....")
			grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

			if not grosolvated == 'fsolvated.gro' and TFF > 0:
				printNote("Solvation with backup complex unsuccessful")

				print("Trying Solvation eith backup topologies ...")
				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bktopol.top', 'topol.top')
				if "tlptopol.top" in os.listdir() and "xtlptopol.top" in os.listdir():
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				printNote("Repeating Solvation with backup in progress.....")
				grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro':
					printNote("Solvation with backup topologies failed again")
					print("Solvation can not continue. You may wish to do solvation manually")
				else:
					printNote("Solvation of backup complex with backup topologies was successful")

			elif not grosolvated == 'fsolvated.gro' and TFF == 0:
				printNote("Solvation with backup complex topologies unsuccessful")
				print("Solvation can not continue. You may wish to do solvation manually")

			else:
				printNote("Solvation of with backup complex and topologies was successful")

		else:
			printNote("Required gro file for solvation not found / generated. Solvation can not continue")
			printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

		os.chdir('../')
		print('\n')

		######################################
		# GATHERING FILES FOR MD SIMULATIONS #
		######################################
		# Create gmxmds directory and populate it with needed data
		print("Gathering and sorting files needed for MDS in progress ...")
		time.sleep(2)

		gmxname = "gmxmds"
		os.mkdir(gmxname)
		os.chdir(gmxname)

		listsoldir = os.listdir(os.path.join(workhost_dir, 'Solvation'))
		for file in listsoldir:
			if Path(file).suffix == ".itp" or Path(file).suffix == ".top" or Path(file).suffix == ".gro":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')

		# Removing files not needed
		gmxmdsrequired = ['fsolvated.gro', 'ufsolvate.gro', 'tlpSolvated.gro', 'tlpSolvated.top', 'tlptopol.top', 'topol.top', 'utopol.top', 'bktopol.top']
		for itpf in os.listdir():
			if Path(itpf).suffix == ".itp":
				gmxmdsrequired.append(itpf)
		 
		for rmf in os.listdir():
			if not rmf in gmxmdsrequired: 
				os.remove(rmf)

		if 'fsolvated.gro' in os.listdir():
			try:
				os.remove('tlpSolvated.gro')
				os.remove('tlpSolvated.top')
			except:
				pass

		# Getting Include files ready and up-to-date
		neededIncludeFiles = gmxmdsFChecks(os.listdir())
		notNeededIncludeFiles = []

		for includefile in os.listdir():
			if Path(includefile).suffix == ".itp":
				if not includefile[0:5] == "posre":
					if not includefile in neededIncludeFiles:
						if not includefile in notNeededIncludeFiles:
							notNeededIncludeFiles.append(includefile)
		
		os.mkdir("not4mds")
		os.mkdir("backup")

		for notfile in notNeededIncludeFiles:
			if notfile[0:2].lower() == "bk":
				shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'backup'))
			else:
				shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'not4mds'))

		for bkfile in os.listdir():
			if bkfile[0:2].lower() == "bk":
				shutil.move(bkfile, os.path.join(workhost_dir, 'gmxmds', 'backup'))

		# Performing final processing for tlptopol.top file, if it has not been used instead of topol.top file
		if "topol.top" in os.listdir() and "tlptopol.top" in os.listdir() and not defaults[1] == "none":
			print("Final processing of tlptopol.top file in progress .....")
			tlpfinal("tlptopol.top", "topol.top")

		# If required, new restraint file can now be generated
		if "fsolvated.gro" in os.listdir():
			grosolvated = "fsolvated.gro"

			print('\n')
			# Generating new restraint file if needed
			printNote("To generate a new restraint file interactively, type YES/y, otherwise press ENTER")
			response = tinput("Response YES/y or press ENTER: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				restrfile = udrestraint(grosolvated)
				if not (restrfile == "posre_udp.itp" or "posre_udp.itp" in os.listdir()):
					printNote("No user defined restriant file has been generated")
				else:
					printNote("User defined restraint file, posre_udp.itp, was genreated successfully")
					printNote("To use it in MDS, define -DPOSRES_UDP in .mdp files OR if already defined -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")

			else:
				printNote("No user defined restraint file has been generated")

		elif "tlpSolvated.gro" in os.listdir() and "tlpSolvated.top" in os.listdir():
			print("TLeap generated solvated files have been moved to gmxmds folder. Please read README.md file for further guidance on how to use it")
    
		else:
			print("To generate solvated files using amber tleap, follow instruction in README.md file")

		# Final Processing of generated files for better formating
		print("Performing Final Checking of Processed MDS input files...")
		if "fsolvated.gro" in os.listdir() or "ufsolvate.gro" in os.listdir():
			for mdsfile in os.listdir():
				if not os.path.isdir(mdsfile):
					gmxmdsFClean(mdsfile)
			os.chdir('../../')
			print('\n')
			npairs += 1
			TFF = M

		else:
			mdsfgerrors += 1
			decision, newdefaults = gmxmdsFEChecks(defaults)
			if decision == "abort":
				raise Exception("Process aborted to allow for user to checks and make corrections")
			elif decision == "continue":
				printWarning("You have chosen to continue with the process despite errors. Please check")
				os.chdir('../../')
				print('\n')
				npairs += 1
				TFF = M
			elif decision == "regenerate":
				printNote("Getting files and defauts values ready for Regeneration...")
				defaults = newdefaults
				os.chdir('../../')
				print('\n')
				TFF = M

	printNote("PLEASE NOTE:")
	print(f"{rls_dir} gmxmds subfolder has been populated and ready for molecular dynamic simulation")
	print("The files in subfolder 'not4mds' of 'gmxmds' are considered not necessary for MDS. Please check")
	print("The files in subfolder 'backup' of 'gmxmds' can be used to replace relevant ones in gmxmds folder for MDS")

	print('\n')
	printNote("Setup with RLmany route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS") 

	# Set the starting date and time
	print('\n')
	rlmanytime2 = datetime.now()
	Tend = rlmanytime2.strftime("%B %d, %Y %H:%M:%S")
	print(f'Generation of MDS input files ends: {Tend}')
