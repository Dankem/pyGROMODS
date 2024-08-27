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

from gmodsScripts.gmodsHelpers import ligtopol, receptopol, topolsplit, indexoflines, complexgen, pdbcatogro, solvation, insertdetails, printWarning, printNote, tinput, defaults1, defaults2, udrestraint, gmxmdsFChecks, gmxmdsFClean

from gmodsScripts.gmodsTLptopol import TLtopol, tlpfinal
from gmodsScripts.gmodsTScheck import checktopsimilarities, Checkligtop
from gmodsScripts.gmodsOPLStop import OPLStop, OPLSacpype

def RLmulti(appDIR, gmxDIR, fdefaults):
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

	# Set global variable to access user supplied topology file(s)
	TFF = 0
	Ligsff = " "
	tff = " "

	# Set global variable to automatically run pdb2gmx, editconf and others
	# forcefields & water (select), -bt (triclinic), -d (0.1), and timeout (60)
	defaults = fdefaults 
	
	print('\n')
	printNote("Let us check again your selected default values ..... ")

	print(f"Default forcefield is {defaults[0]}")
	print(f"Default water model is {defaults[1]}")
	print(f"Default editconf -bt option is {defaults[2]}")
	print(f"Default editconf -d option is {defaults[3]}")
	print(f"Default timeout for input() request is {defaults[4]}")

	print('\n')
	if defaults[5] == "A":
		printNote("Your selected default mode for generating input file is Interractive")
		response = tinput("To revert to Noninteractive mode type YES/y: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[5] = "B"
			printNote("Noninteractive mode: Effective following first interactive selection of forcefield and water")
	else:
		printNote("Your selected default mode for generating input file is Noninterractive")
		response = tinput("To revert back to Interactive mode type YES/y: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[0] = "select"
			defaults[1] = "select"
			defaults[5] = "A"
			printNote("Interactive mode: Effective immediately")
		else:
			defaults[5] = "C"

	response = tinput("To adjust further the selected default values of -d, -bt and timeout type YES/y: ", defaults[4], "n")
	if response.lower() == "yes" or response.lower() == "y":
		defaults[2], defaults[3], defaults[4] = defaults1() 
	print('\n')

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
			printNote("If the topology is important, please abort, crosscheck and rerun")
			printNote("Otherwise, by default user uploaded ligand topology files will be ignored")
			response = input("To continue and ignore uploaded topology files, type YES/y, or press ENTER to abort: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Incomplete topology files!!!. Check and restart the process")
			else:
				printNote("Uploaded topology files will be ignored")
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
				printNote("Default setting is to continue without the topology files")
				response = tinput("To abort and make correction, type YES/y: ", defaults[4], "n")
				if response.lower() == "yes" or response.lower() == "y":
					raise Exception("Sorry, unacceptable topology file format detected")
				else:
					TFF = 0

			elif not ntopfile == 0:
				print(f"{ntopfile} Topology file(s) lacked required ['atomtypes'], ['moleculetype'] and/or [ system ] subheading")
				print("Check and include the subheadings, with or without expected accompanied values")
				printNote("Default setting is to continue without the topology files")
				response = tinput("To abort and make correction, type YES/y: ", defaults[4], "n")
				if response.lower() == "yes" or response.lower() == "y":
					raise Exception("Checking ligand topology failed. Check, correc and restart")
				else:
					TFF = 0

	# Check that all files have been correctly uploaded
	for file in LIGfiles:
		if not (Path(file).suffix == ".pdb" or Path(file).suffix == ".mol2"):
			raise Exception("Sorry, ligand files must be pdb or mol2 files. Setroute again and upload required files")

	for file in RECfiles:
		if not Path(file).suffix == ".pdb":
			raise Exception("Sorry, Receptor file must be pdb file. Setroute again and upload required files")

	if not (len(RECfiles) and len(LIGfiles)) > 0:
		raise Exception("Needed ligands and/or receptor files is/are missing. Please upload all required files")

	if len(RECfiles) > 1:
		raise Exception("Sorry, only one Receptor is allowed per set of ligands. You have " + str(len(RECfiles)) + " receptors")

	if len(LIGfiles) > 1:
		printNote("You have more than one ligand. All the ligands will be treated together to form just one multi-ligands complex with the receptor")
		response = tinput("Type YES/y to continue or press ENTER to abort: ", defaults[4], "y")
		if not (response.lower() == "yes" or response.lower() == "y"):
			raise Exception("You have choosen to Abort!!!. Check and restart the process")

	# Check that appropraite tleap source file is available for use
	while True:
		PLDfiles = os.listdir(PLD)
		LIGfiles = os.listdir(os.path.join(fPDB, 'Ligands'))
		numberL = len(LIGfiles)
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
	printNote("CAREFULLY CHECK THE ABOVE RESULTS AND TYPE 'YES/y' TO CONFIRM OR 'NO/n' TO ABORT")
	confirmation = input("Response is: ")
	if not (confirmation.lower() == "yes" or confirmation.lower() == "y"):
		raise Exception("Process aborted by the User. Please rerun from SetRoute")
	else:
		rlmultitime = datetime.now()
		Tstart = rlmultitime.strftime("%B %d, %Y %H:%M:%S")
		print(f'Generation of MDS input files begins: {Tstart}')
	print('\n')
	time.sleep(5)

	# Get user imput for project name
	while True:
		name = tinput("Suppy a name for the current project: ", defaults[4], "RLmulti")
		if name == " ":
			print("You must supply a name for the project")
			continue
		elif name.isalpha() == False:
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
		foldername = name + "_" + str(ID)
		if not os.path.isdir(foldername):
			print("Your Unique ID is: ", ID)
			print("Your current work will be stored in ", foldername)
			break
		else:
			continue

	# Create working folder using the generated project name
	os.mkdir(foldername)
	os.chdir(foldername)

	work_dir = os.path.join(GMX_MDS, foldername)
	print(f"pyGROMODS Current workspace directory set to {work_dir}")

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

	###############################################
	# GENERATING PROTEIN TOPOLOGIES AND STRUCTURE #
	###############################################

	os.mkdir('Receptor')
	os.chdir('Receptor')

	receptor_dir = os.path.join(work_dir, 'Receptor')
	print(f"Protein data directory set to {receptor_dir}")
	time.sleep(2)

	printNote("Generating Protein topology and parameter files...")

	# Generate Protein topologies and associated parameters
	rname = "receptor"
	rep = RECfiles[0]
	selff = defaults[0]
	selwater = defaults[1]

	shutil.copy(os.path.join(fPDB, 'Receptors', rep), './')

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
			print(f"{tff} forcefield does not match the preselected forcefield group: {Ligsff}")
			printNote("PLEASE NOTE - If you choose to continue: ")
			print(f"A). Your forcefield group will be changed to {tff}")
			print("B). By default, any uploaded ligand topology(ies) will be ignored")

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
					print(f"{tff} does not match any known forcefield group. Please rerun")
					printNote("If using a self created or modified forcefield, adopt standard naming convention")
					response = input("Type YES/y to select new forcefield interactively. Otherwise press ENTER: ")
					if not (response.lower() == "yes" or response.lower() == "y"):
						printNote("Generation of MDS input files will not proceed")
						raise Exception("Process aborted. Make necessary corrections and Rerun setup")
					else:
						selff = "select"
						selwater = "select"
						continue

				Ufolders = os.listdir(fPDB)
				if not ('Amber' in Ufolders or 'Charmm' in Ufolders or 'Gromos' in Ufolders or 'Opls' in Ufolders):
					raise Exception("Expected forcefield folder or file was not found. Please set route correctly")
				else:
					for sff in Ufolders:
						if sff == 'Amber' or sff == 'Charmm' or sff == 'Gromos' or sff == 'Opls':
							os.rename(os.path.join(fPDB, sff), os.path.join(fPDB, Ligsff))

				print(f"Your preselected forcefield group has been changed to {Ligsff}")
				ligsff_dir = os.path.join(fPDB, Ligsff)
				if TFF > 0:
					print(f"Subdirectory for uploaded ligand topology is now: {ligsff_dir}")
					printNote("To use with uploaded ligand topology, type YES/y. Otherwise press ENTER to ignore")
					response = tinput("Response: ", 30, "n")
					if not (response.lower() == "yes" or response.lower() == "y"):
						TFF = 0
						break
					else:
						TFF = 1
						break
				else:
					break
		else:
			print("Your forcefiled as contained in topol.top file is", tff)
			break

	try:
		shutil.copy(rep, 'xreceptor.pdb')
		shutil.copy(RFtop, 'nReceptor.top')
		shutil.copy(RFpdb, 'nReceptor.pdb')
	except Exception as e:
		printWarning(e)
		pass

	# If needed, we now change the default parameters and lock the new default mode
	if not (defaults[0] == "select" or defaults[1] == "select"):
		if defaults[5] == "B":
			if Path(tff).suffix == ".ff":
				defaults[0] = Path(tff).stem
			else:
				defaults[0] = tff
			
			defaults[1] = defaults2(RFtop)
			if defaults[1] == "none":
				print("No water model was detected for your system")
			else:
				print("Your water model as contained in topol.top file is ", defaults[1])
			
			printNote("Your selected default values are as follows: ")
			print(f"		Default forcefield: {defaults[0]}")
			print(f"		Default water model: {defaults[1]}")
			print(f"		Default editconf -bt: {defaults[2]}")
			print(f"		Default editconf -d: {defaults[3]}")
			print(f"		Default input timeout: {defaults[4]}")
			print("		Default mode: non-interactive")
			defaults[5] = "C"
		else:
			defaults[1] = defaults2(RFtop)
			if defaults[1] == "none":
				print("No water model was detected for your system")
			else:
				print(f"Your water model as contained in topol.top file is {defaults[1]}")

	# Checking for posre.itp file or create one to include relevant posre files
	print("Updating topology file .....")
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
		printNote("The following OPLS options are available to generate ligand topology:")
		print("	1. Using platform OPLS compatible ligand topology generation - RECOMMENDED")
		print("	2. Using acpype OPLS compatible ligand topology generation")
		print("	3. Using platform default - ignore opls compatibility - NOT RECOMMENDED")
		while True:
			response = tinput("Choose your preferred option: ", defaults[4], "1")
			if response == '1' or response == '2' or response == '3':
				print(f"Option {response} was Selected!")
				confirm = tinput("Type YES/y to confirm or press ENTER to choose a different option: ", defaults[4], "y")
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
		time.sleep(5)

	#################################
	# GENERATING LIGANDS TOPOLOGIES #
	#################################
	print('\n')
	printNote("Generating topology and parameter files for ligands...")

	os.mkdir('Ligands')
	os.chdir('Ligands')

	lig_dir = os.path.join(work_dir, 'Ligands')
	print("Ligands data directory set to ", lig_dir)

	os.mkdir('at_ligs')
	os.mkdir('mt_ligs')
	os.mkdir('top_ligs')
	os.mkdir('st_ligs')
	os.mkdir('para_ligs')

	if TFF > 0:
		os.mkdir('user_top')
		os.mkdir('user_at')
		os.mkdir('user_mt')

	Lcount = 0
	for lig in LIGfiles:
		Lcount += 1
		Lname = "LIG" + str(Lcount)
		os.mkdir(Lname)
		os.chdir(Lname)
		shutil.copy(os.path.join(PLD, 'ligtoptleap.in'), './')
		shutil.copy(os.path.join(fPDB, 'Ligands', lig), './')

		LFgro, LFtop = ligtopol(lig, Lname)

		# Identify the pdb file to use subsequently, depending on uploaded ligand format
		print(f"Generating optimized {lig} structure...")
		ligE = open('ligcleanerror.txt', 'a')
		nligname = Lname + "up.pdb"

		if Path(lig).suffix == ".pdb":
			ligtopol_pdb = lig
		elif Path(lig).suffix == ".mol2":
			ligtopol_pdb = Lname + "new" + ".pdb"

		try:
			subprocess.run(['gmx', 'editconf', '-f', ligtopol_pdb, '-o', nligname], check=True, stderr=subprocess.STDOUT, stdout=ligE, text=True)
		except subprocess.SubprocessError as e:
			shutil.copy(ligtopol_pdb, nligname)
		ligE.close()

		# Prepare atomtypes and moleculetype files for user supplied topology files
		if TFF > 0:
			os.mkdir('usertops')
			os.chdir('usertops')

			u = Lcount - 1
			utop = TOPfiles[u]
			shutil.copy(os.path.join(fPDB, Ligsff, utop), './')
			ulig = nligname
			shutil.copy(os.path.join(lig_dir, Lname, ulig), './')

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

			ignore = 0
			if not (dul_ml > ul or dul_ml > ml):
				print(f"{utop}, the user supplied topology match the ligand named {ulig}")
				time.sleep(3)
			else:
				printWarning("WARNINGS TO NOTE:")
				print(f"{utop} user supplied topology does not match the ligand named {ulig}")
				print("This may happen if the topology(ies) is/are not named similar to the corresponding ligand(s)")
				print("Please compare the two files and if you're sure it's correct, continue, otherwise abort")
				printNote("PLEASE TYPE YES/y TO USE THE TOPOLOGY FILES AS UPLOADED")
				printNote("PLEASE TYPE NO/n TO IGNORE THE UPLOADED TOPOLOGY FILES")
				printNote("PLEASE PRESS ENTER TO ABORT THE PROCESS")
				response = input("Response is: ")
				if response.lower() == "no" or response.lower() == "n":
					TFF = 0
					ignore += 1
				if not (response.lower() == "yes" or response.lower() == "y" or response.lower() == "no" or response.lower() == "n"):
					raise Exception("Process aborted. Restart and rerun")

			# Check the file for duplicate atom types with forcefield standard atomtypes
			if ignore == 0:
				uligcheck = Checkligtop(utop, tff)
				uctindex = indexoflines(uligcheck)

				# Generate the atomtypes and moleculetypes file
				uligand_at, uligand_mt = topolsplit(uligcheck, Lname, uctindex)

				shutil.copy(uligand_at, '../../user_at')
				utopname = Lname + ".top"
				os.rename(uligcheck, utopname)
				shutil.copy(utopname, '../../user_top')
				os.chdir('../')
			else:
				os.chdir('../')

		# Create variables for easy copy other needed files into created folders
		st_mol2 = Lname + ".mol2"
		lp_prm = Lname + ".prmtop"
		lp_inp = Lname + ".inpcrd"
		lp_frc = Lname + ".frcmod"
		st_pdb = Lname + ".pdb"

		# If need be, generate opls compatible ligand topology using acpype
		oplstop = " "

		if opls_route == 1:
			mergename = "opls" + Lname + ".top"
			oplstop = OPLStop(LFtop, tff, Lname)

			if not oplstop == mergename:
				printWarning("Platform opls compatible ligand topology generation failed")
				printNote("Trying acpype...")
				oplstop = OPLSacpype(LFtop, lig, Lname)
				if not oplstop == mergename:
					printWarning("acpype opls compatible ligand topology generation also failed")
					print("Check that you have correctly installed the latest acpype")
					printNote("Platform default will now be used")
					oplstop = LFtop		

		elif opls_route == 2:
			mergename = "opls" + Lname + ".top"

			oplstop = OPLSacpype(LFtop, lig, Lname)
			if not oplstop == mergename:
				printWarning("acpype opls compatible ligand topology generation failed")
				printNote("Trying Platform opls compatible ligand topology generation ...")
				oplstop = OPLStop(LFtop, tff, Lname)
				if not oplstop == mergename:
					printWarning("Platform opls compatible ligand topology generation also failed")
					printNote("Platform default will now be used")
					oplstop = LFtop		

		else:
			oplstop = LFtop

		# Check the generated ligand topology file for duplicate atom types with forcefield standard atomtypes
		ligcheck = Checkligtop(oplstop, tff)

		nametop = " "
		if not ligcheck == oplstop:
			nametop = "ck" + Lname + ".top"
			os.rename(ligcheck, nametop)
		else:
			nametop = oplstop

		topindex = indexoflines(nametop)
		ligand_at, ligand_mt = topolsplit(nametop, Lname, topindex)

		shutil.copy(nametop, '../top_ligs')
		shutil.copy(LFgro, '../st_ligs')
		shutil.copy(st_mol2, '../st_ligs')
		shutil.copy(lp_prm, '../para_ligs')
		shutil.copy(lp_inp, '../para_ligs')
		shutil.copy(lp_frc, '../para_ligs')
		shutil.copy(ligand_at, '../at_ligs')

		shutil.copy(nligname, st_pdb)
		shutil.move(st_pdb, '../st_ligs')

		# Check ligand directory for posre file and back it up
		for posrefile in os.listdir():
			if posrefile[0:5] == "posre":
				posrebk = posrefile + ".bk"
				os.rename(posrefile, posrebk)

		os.chdir('../')

	# Check topology files for similarities
	os.chdir('top_ligs')
	mtdir = os.path.join(work_dir, 'Ligands', 'mt_ligs')
	ligand_mn = checktopsimilarities(mtdir)
	shutil.copy(ligand_mn, '../')
	os.chdir('../')

	if TFF > 0:
		os.chdir('user_top')
		umtdir = os.path.join(work_dir, 'Ligands', 'user_mt')
		uligand_mn = checktopsimilarities(umtdir)
		os.rename(uligand_mn, 'uLIGS_mn.itp')
		shutil.copy('uLIGS_mn.itp', '../')
		os.chdir('../')

	# Combine the atomtypes files into one file called LIGS_at.itp
	file_at = open("LIGS_at.itp", "a+")
	at_list = os.listdir('at_ligs')
	for file in at_list:
		newfile = os.path.join(lig_dir, 'at_ligs', file)
		atfile = open(newfile, "r")
		readatfile = atfile.readlines()

		if not "atomtypes" in (readatfile[0].split() or readatfile[1].split() or readatfile[2].split() or readatfile[3].split() or readatfile[4].split()):
			print(file, "lack the [atomtypes] header signature. This file will be skipped")
			continue

		for line in readatfile:
			file_at.seek(0)
			listligs_at = []
			readfile_at = file_at.readlines()
			for atline in readfile_at:
				if not (atline.split() == [] or atline.split()[0] in listligs_at):
					listligs_at.append(atline.split()[0])

			if not (line.split() == [] or line.split()[0] in listligs_at):
				file_at.write(line)
			else:
				continue
		atfile.close()
	file_at.close()
	if Ligsff == "Opls":
		shutil.copy("LIGS_at.itp", "oplsLIGS_at.itp")

	# Combine the user uploaded topology(ies) atomtypes files into one file called uLIGS_at.itp
	if TFF > 0:
		ufile_at = open("uLIGS_at.itp", "a+")
		userat_list = os.listdir('user_at')
		for ufile in userat_list:
			unewfile = os.path.join(lig_dir, 'user_at', ufile)
			uatfile = open(unewfile, "r")
			ureadatfile = uatfile.readlines()

			if not "atomtypes" in (ureadatfile[0].split() or ureadatfile[1].split() or ureadatfile[2].split() or ureadatfile[3].split() or ureadatfile[4].split()):
				print(ufile, "lack the [atomtypes] header signature. This file will be skipped")
				continue

			for uline in ureadatfile:
				ufile_at.seek(0)
				ulistligs_at = []
				ureadfile_at = ufile_at.readlines()
				for uatline in ureadfile_at:
					if not (uatline.split() == [] or uatline.split()[0] in ulistligs_at):
						ulistligs_at.append(uatline.split()[0])

				if not (uline.split() == [] or uline.split()[0] in ulistligs_at):
					ufile_at.write(uline)
				else:
					continue
			uatfile.close()
		ufile_at.close()

	# Combine the moleculetype files into one file called LIGS_mt.itp
	file_mt = open("LIGS_mt.itp", "a+")
	mt_list = os.listdir('mt_ligs')
	for file in mt_list:
		nextfile = os.path.join(lig_dir, 'mt_ligs', file)
		mtfile = open(nextfile, "r")
		readmtfile = mtfile.readlines()

		if not "moleculetype" in (readmtfile[0].split() or readmtfile[1].split() or readmtfile[2].split() or readmtfile[3].split() or readmtfile[4].split()):
			print(file, "lack the [moleculetype] header signature. This file will be skipped")
			continue

		for line in readmtfile:
			file_mt.write(line)
		file_mt.write('\n')
		mtfile.close()
	file_mt.close()

	# Combine the moleculetype from user uploaded files into one file called LIGSopls_mt.itp
	if TFF > 0:
		ufile_mt = open("uLIGS_mt.itp", "a+")
		usermt_list = os.listdir('user_mt')
		for ufmt in usermt_list:
			unextfile = os.path.join(lig_dir, 'user_mt', ufmt)
			umtfile = open(unextfile, "r")
			ureadmtfile = umtfile.readlines()

			if not "moleculetype" in (ureadmtfile[0].split() or ureadmtfile[1].split() or ureadmtfile[2].split() or ureadmtfile[3].split() or ureadmtfile[4].split()):
				print(ufmt, "lack the [moleculetype] header signature. This file will be skipped")
				continue

			for lineu in ureadmtfile:
				ufile_mt.write(lineu)
			ufile_mt.write('\n')
			umtfile.close()
		ufile_mt.close()

	os.chdir('../')

	###########################################################
	# GENERATING PROTEIN-LIG COMPLEX TOPOLOGIES AND STRUCTURE #
	###########################################################
	print('\n')
	os.mkdir('Complex')
	os.chdir('Complex')

	complex_dir = os.path.join(work_dir, 'Complex')
	print(f"Protein-Ligand Complex data directory set to {complex_dir}")
	time.sleep(2)

	# Generate topologies and associated ligand-complex for each receptor
	printNote("Generating Protein-Ligands complex topologies and parameter files...")

	# Gather other necessary files into the folder to generate complex
	listLdir = os.listdir(os.path.join(work_dir, 'Ligands'))
	for L in listLdir:
		if Path(L).suffix == ".itp":
			shutil.copy(os.path.join(work_dir, 'Ligands', L), './')

	stligs = os.listdir(os.path.join(work_dir, 'Ligands', 'st_ligs'))
	for lig in stligs:
		if Path(lig).suffix == ".gro":
			shutil.copy(os.path.join(work_dir, 'Ligands', 'st_ligs', lig), './')

	listRdir = os.listdir(os.path.join(work_dir, 'Receptor'))
	for R in listRdir:
		if R == 'xreceptor.pdb' or R == 'nReceptor.top' or R == 'nReceptor.pdb':
			shutil.copy(os.path.join(work_dir, 'Receptor', R), './')
		elif Path(R).suffix == ".itp":
			shutil.copy(os.path.join(work_dir, 'Receptor', R), './')

	shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'at-topol.itp'), './')
	shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol.itp'), './')

	# Now we shall insert details into pdb2gmx topol file, if successfully generated
	listCpldir = os.listdir()
	if 'nReceptor.top' in listCpldir:
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
			insertUAT = "[ moleculetype ]"
			insertdetails('utopol.top', 'at-topol.itp', insertUAT)

			insertUMT = "; Include water topology"
			insertdetails('utopol.top', 'mt-topol.itp', insertUMT)

			utopolfile = open("utopol.top", "a")
			uligmnfile = open("uLIGS_mn.itp", "r")
			ureadtops = uligmnfile.read()
			utopolfile.write(ureadtops)
			utopolfile.close()
			uligmnfile.close()

	print("Generating Protein-Ligands complex structures...")

	# Generate receptor - ligand complex with amber
	os.mkdir('ligmol')
	os.mkdir('ligfrcmod')
	os.mkdir('ligpdb')

	os.chdir('ligmol')
	try:
		subprocess.run('cp ../../Ligands/st_ligs/*.mol2 ./', shell=True, stderr=subprocess.STDOUT, check=True, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning("Something is not right. Check the above error message")
		time.sleep(5)			

	os.chdir('../ligfrcmod')
	try:
		subprocess.run('cp ../../Ligands/para_ligs/*.frcmod ./', shell=True, stderr=subprocess.STDOUT, check=True, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning("Something is not right. Check the above error message")
		time.sleep(5)			
		
	os.chdir('../ligpdb')
	try:
		subprocess.run('cp ../../Ligands/st_ligs/*.pdb ./', shell=True, stderr=subprocess.STDOUT, check=True, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning("Something is not right. Check the above error message")			
		time.sleep(5)

	os.chdir('../')

	listmol = os.listdir('ligmol')
	listfrcmod = os.listdir('ligfrcmod')
	if not len(listmol) == len(listfrcmod):
		printWarning("Something went wrong, number of .mol2 files didn't match number of .frcmod files")
		printNote("Defualt complex generation will not continue. Please crosscheck")
	else:
		fx = 1
		for filex in listmol:
			xname = "xLIG" + str(fx) + ".mol2"
			shutil.copy(os.path.join(complex_dir, 'ligmol', filex), './')
			os.rename(filex, xname)
			fx += 1

		fy = 1
		for filey in listfrcmod:
			yname = "yLIG" + str(fy) + ".frcmod"
			shutil.copy(os.path.join(complex_dir, 'ligfrcmod', filey), './')
			os.rename(filey, yname)
			fy += 1

		numL = len(listmol)
		tleapfile = "tlpname.in"
		for filez in PLDfiles:
			tlplist = filez.split("-")
			if str(numL) in tlplist:
				shutil.copy(os.path.join(PLD, filez), './')
				os.rename(filez, tleapfile)
				break

	# Check if tleap source file is present and attempt to generate complex via tleap or use alternative
	complexdir = os.listdir()
	if not tleapfile in complexdir:
		print('\n')
		printNote("Missing tleap input files! Can't proceed with default complex generation. Check README.md")
		print("Trying alterenative complex generation...")

		catcomplx_gro, catcomplx_pdb = pdbcatogro()

		if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
			print("Alternative approach failed to generate complex")
			raise Exception("Process aborted. Make necessary corrections and rerun")
	
		else:
			os.rename(catcomplx_gro, 'catComplex.gro')
			os.rename(catcomplx_pdb, 'catComplex.pdb')

	else:
		complx_gro, complx_top = complexgen(tleapfile)

		if not (complx_gro == 'Complex.gro' and complx_top == 'Complex.top'):
			printWarning("Default generation of complex structure and topologies failed")
			printNote("Please check the outputed errors message for detail and README.md for help")

			print("Trying alterenative complex generation...")
			time.sleep(5)

			catcomplx_gro, catcomplx_pdb = pdbcatogro()

			if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
				printWarning("Alternative approach also failed to generate complex")
				raise Exception("Process aborted. Make necessary corrections and rerun")
	
			else:
				os.rename(catcomplx_gro, 'catComplex.gro')
				os.rename(catcomplx_pdb, 'catComplex.pdb')

		else:
			printNote("Default complex generation was successful")
			os.rename(complx_gro, 'tlpComplex.gro')
			os.rename(complx_top, 'tlpComplex.top')

			print("Attempting alternative complex generation as fall back ...")
			catcomplx_gro, catcomplx_pdb = pdbcatogro()

			if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
				printNote("Alternative approach failed, No complex to serve as fallback")
	
			else:
				printNote("Alternative complex generation successful")
				os.rename(catcomplx_gro, 'catComplex.gro')
				os.rename(catcomplx_pdb, 'catComplex.pdb')

	try:
		os.remove("LIGS_mn.itp")
		os.remove("uLIGS_mn.itp")
	except:
		pass

	os.chdir('../')

	#####################################################
	# GENERATING SOLVATED PROTEIN-LIG COMPLEX STRUCTURE #
	#####################################################
	print('\n')
	os.mkdir('Solvation')
	os.chdir('Solvation')

	solvation_dir = os.path.join(work_dir, 'Solvation')
	print(f"Solvated data directory set to {solvation_dir}")
	time.sleep(2)

	grosolvated = " "

	printNote("Generating Solvated Protein-Ligands complex structures...")

	listLdir = os.listdir(os.path.join(work_dir, 'Ligands'))
	for L in listLdir:
		if Path(L).suffix == ".itp":
			shutil.copy(os.path.join(work_dir, 'Ligands', L), './')

	listRdir = os.listdir(os.path.join(work_dir, 'Receptor'))
	for R in listRdir:
		if R == 'xreceptor.pdb' or R == 'nReceptor.top' or R == 'nReceptor.pdb':
			shutil.copy(os.path.join(work_dir, 'Receptor', R), './')
		elif Path(R).suffix == ".itp":
			shutil.copy(os.path.join(work_dir, 'Receptor', R), './')

	listTdir = os.listdir(os.path.join(work_dir, 'Complex'))
	for T in listTdir:
		if T == 'tlpComplex.gro' or T == 'tlpComplex.top' or T == 'catComplex.gro' or T == 'catComplex.pdb' or T == 'tlpSolvated.gro' or T == 'tlpSolvated.top' or T == 'topol.top' or T == 'utopol.top':
			shutil.copy(os.path.join(work_dir, 'Complex', T), './')
		elif T[0:5] == 'posre' and not T in os.listdir():
			shutil.copy(os.path.join(work_dir, 'Complex', T), './')

	if not "at-topol.itp" in os.listdir():
		shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'at-topol.itp'), './')

	# Determine if the forcefield belong to the selected group - amber, charmm, gromos and opls
	if TFF > 0:
		print('\n')
		if tff[0:4].lower() == Ligsff.lower() or tff[0:5].lower() == Ligsff.lower() or tff[0:6].lower() == Ligsff.lower():
			printNote("PLEASE NOTE THAT THE DEFAULT ORDER ARE:")
			print(">>>>>1. User uploaded ligand topologies will be attempted first")
			print(">>>>>2. Platform generated ligand topologies will serve as fallback")
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
				printNote("DEFAULT ORDER IS MAINTAINED")

			else:
				os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('utopol.top', 'bktopol.top')
				printNote("THE DEFAULT ORDER HAS BEEN REVERSED")

		else:
			print('\n')
			print(f"{tff} in topol.top doesn't match the selected forcefield group, {Ligsff}")
			printNote("Solvation will most likely fail. It is advisable to check and rerun the process")
			# Backing off some relevant files
			try:
				os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('utopol.top', 'bktopol.top')
			except Exception as e:
				print(e)
				pass

	else:
		print('\n')
		if tff[0:5].lower() == 'amber':
			print(f"{tff} is an amber forcefield. Auto-generated ligand topology will be used") 

		elif tff[0:6].lower() == 'charmm':
			print(f"{tff} is a charmm forcefield. Auto-generated ligand topology will be used") 

		elif tff[0:6].lower() == 'gromos':
			print(f"{tff} is a gromos forcefield. Auto-generated ligand topology will be used") 
			printNote("If this fails, please try uploading gromos compatible ligand topology")

		elif tff[0:4].lower() == 'opls':
			print(f"{tff} is a opls forcefield. Auto-generated ligand topology will be used") 
			printNote("If this fails, please try uploading opls compatible ligand topology")
	
	# Time to prepared solvated complex
	print('\n')
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
		printNote("Solvation in progress ...")
		grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)
			
		if not grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir:
			printNote("Solvation failed with default complex")

			printNote("Trying Solvation with the backup complex ...")
			grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

			if not grosolvated == 'fsolvated.gro' and TFF > 0:
				printNote("Solvation with backup complex unsuccessful")
				print("Re-trying default complex with backup topology files ...")

				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bktopol.top', 'topol.top')
				if "tlptopol.top" in os.listdir() and "xtlptopol.top" in os.listdir():
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro':
					printWarning("Solvation with backup topologies failed")

					print("Re-trying Solvation with the backup complex and topologies ...")
					grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

					if not grosolvated == 'fsolvated.gro':
						printWarning("Solvation using backup complex and topologies unsuccessful")

					else:
						printNote("Solvation of backup complex with backup topologies was successful")

				else:
					printNote("Solvation of default complex with backup topologies was successful")

			elif not grosolvated == 'fsolvated.gro' and TFF == 0:
				printWarning("Solvation with backup complex unsuccessful")

			else:
				printNote("Solvation with backup complex was successful")

		elif not (grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir):
			if TFF > 0:
				print("Trying solvation of default complex with backup topology files ...")

				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bktopol.top', 'topol.top')
				if "tlptopol.top" in os.listdir() and "xtlptopol.top" in os.listdir():
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro':
					printWarning("Solvation of default complex with backup topologies failed")

				else:
					printNote("Solvation of default complex with backup topologies was successful")

			else:
				printWarning("Solvation with default complex unsuccessful")
				
		else:
			printNote("Solvation with default complex was successful")

	elif 'catComplex.gro' in chkSoldir:
		printNote("Solvation with alternative complex in progress.....")
		grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

		if not grosolvated == 'fsolvated.gro' and TFF > 0:
			printWarning("Solvation with backup complex unsuccessful")

			print("Re-trying with backup topologies ...")

			os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
			os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
			os.rename('topol.top', '#topol.top.bk#')
			os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
			os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
			os.rename('bktopol.top', 'topol.top')
			if "tlptopol.top" in os.listdir() and "xtlptopol.top" in os.listdir():
				os.rename('tlptopol.top', '#tlptopol.top.bk#')
				shutil.copy('xtlptopol.top', 'tlptopol.top')
	
			grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

			if not grosolvated == 'fsolvated.gro':
				printWarning("Solvation with backup topologies failed")

			else:
				printNote("Solvation with backup topologies was successful")

		elif not grosolvated == 'fsolvated.gro' and TFF == 0:
			printWarning("Solvation with backup complex unsuccessful")

		else:
			printNote("Solvation with backup complex was successful")

	else:
		printNote("Required gro file for solvation not found / generated. Solvation can not continue")
		printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

	os.chdir('../')

	######################################
	# GATHERING FILES FOR MD SIMULATIONS #
	######################################
	# Make gmxmds directory and populate it with needed data
	print('\n')
	os.mkdir('gmxmds')
	os.chdir('gmxmds')

	gmxmds_dir = os.path.join(work_dir, 'gmxmds')
	print(f"GMX MDS data directory set to {gmxmds_dir}")
	time.sleep(2)

	printNote("Gathering require files for MDS run into gmxmds directory...")

	listsoldir = os.listdir(os.path.join(work_dir, 'Solvation'))
	for file in listsoldir:
		if Path(file).suffix == ".itp" or Path(file).suffix == ".gro" or Path(file).suffix == ".top":
			shutil.copy(os.path.join(work_dir, 'Solvation', file), './')

	# Removing files that are not needed for MDS with Gromacs
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
		else:
			pass
	
	os.mkdir("not4mds")
	os.mkdir("backup")

	for notfile in notNeededIncludeFiles:
		if notfile[0:2].lower() == "bk":
			shutil.move(notfile, os.path.join(work_dir, 'gmxmds', 'backup'))
		else:
			shutil.move(notfile, os.path.join(work_dir, 'gmxmds', 'not4mds'))

	for bkfile in os.listdir():
		if bkfile[0:2].lower() == "bk":
			shutil.move(bkfile, os.path.join(work_dir, 'gmxmds', 'backup'))

	# Performing final processing for tlptopol.top file, if it has not been used instead of topol.top file
	if "topol.top" in os.listdir() and "tlptopol.top" in os.listdir() and not defaults[1] == "none":
		print("Final processing of tlptopol.top file in progress .....")
		time.sleep(3)
		tlpfinal("tlptopol.top", "topol.top")
		print('\n')

	# Now if need be, let's generate a new restraint file
	if "fsolvated.gro" in os.listdir():
		grosolvated = "fsolvated.gro"
		restrfile = ""

		printNote("To generate a new restraint file, type YES/y, otherwise press ENTER")
		response = tinput("Response: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			restrfile = udrestraint(grosolvated)
			if restrfile == "posre_udp.itp" or "posre_udp.itp" in os.listdir():
				printNote("To use posre_udp.itp in MDS, define -DPOSRES_UDP in .mdp files")
				printNote("To use with defined -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")

	elif "tlpSolvated.gro" in os.listdir() and "tlpSolvated.top" in os.listdir():
		print("TLeap generated solvated files were found and have been moved to gmxmds folder. Please read README.md file for further gudiance on how to use it")

	time.sleep(5)

	# Final Processing of generated files for better formating
	print("Performing Final Checking of Processed MDS input files...")
	if "fsolvated.gro" in os.listdir() or "ufsolvate.gro" in os.listdir():
		for mdsfile in os.listdir():
			if not os.path.isdir(mdsfile):
				gmxmdsFClean(mdsfile)
	else:
		print('\n')
		printWarning("It's appears the process failed to produced all the needed MDS input files")
		print("For more information, Please check the relevant error files as listed below: ")
		print(">>>>>>>> ligerror.txt in Ligand Directory or Subdirectories")
		print(">>>>>>>> recerror.txt in Receptor Directory")
		print(">>>>>>>> cplxerror.txt or cplxBKerror.txt in Complext Directory")
		print(">>>>>>>> checklog.txt in Solvation Directory")
		print('\n')
		printNote("If the error has to do the default values, do the following: ")
		print(">>>>>1). Retrun back to Generate Menu")
		print(">>>>>2). Re-select or re-input a new set of default values")
		print(">>>>>3). Rerun the process again")
		time.sleep(5)

	os.chdir('../')
	print('\n')
	printNote("PLEASE NOTE:")
	printNote("RLmulti gmxmds subfolder has been populated and ready for use for simulation")
	print("The files in folder 'not4mds' of 'gmxmds' are considered not necessary for MDS. Please check")
	print("The files in subfolder 'backup' of 'gmxmds' can be used to replace relevant ones in gmxmds folder for MDS")

	print('\n')
	printNote("Setup with RLmulti route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS") 

	# Set the starting date and time
	print('\n')
	rlmultitime2 = datetime.now()
	Tend = rlmultitime2.strftime("%B %d, %Y %H:%M:%S")
	print(f'Generation of MDS input files ends: {Tend}')
