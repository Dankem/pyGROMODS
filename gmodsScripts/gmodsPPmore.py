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

from gmodsScripts.gmodsHelpers import receptopol, complexgen, solvation, insertdetails, printWarning, printNote, tinput, defaults1, defaults2, udrestraint, gmxmdsFChecks, gmxmdsFClean, gmxmdsFEChecks

from gmodsScripts.gmodsTLptopol import TLtopol, tlpfinal

def PPmore(appDIR, gmxDIR, fdefaults):
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
	global Ligsff
	global tff
	global fPDB
	global PLD

	Ligsff = " "
	tff = " "

	# Set global variable to automatically run pdb2gmx, editconf and others
	# forcefields & water (select), -bt (triclinic), -d (0.1), and timeout (60)
	defaults = fdefaults 
	
	printNote("Here are your selected default values ..... ")
	
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
			print("You have changed to pdb2gmx non-interactive mode")
			print("Your preferred forcefield and water model will be autodetected following your first interactive selection")
	else:
		printNote("Your selected default mode for generating input file is Noninterractive")
		response = tinput("To revert back to Interactive mode type YES/y: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[0] = "select"
			defaults[1] = "select"
			defaults[5] = "A"
		else:
			defaults[5] = "C"

	response = tinput("To adjust further the selected default values of -d, -bt and timeout type YES/y: ", defaults[4], "n")
	if response.lower() == "yes" or response.lower() == "y":
		defaults[2], defaults[3], defaults[4] = defaults1() 
	print('\n')

	# Set some global parameter variables
	PLD = os.path.join(scriptDIR, 'gmodsTSF')

	fSamples = os.path.join(scriptDIR, 'sampleFiles')
	fPDB = os.path.join(scriptDIR, 'Uploads')

	PEPfiles = []
	PROfiles = []

	# Perform initial checks of the ligands and receptor files
	listfPDB = os.listdir(fPDB)
	if 'Peptides' in listfPDB:
		PEPfiles = os.listdir(os.path.join(fPDB, 'Peptides'))
		if not len(PEPfiles) > 0:
			raise Exception("Needed peptide(s) file(s) is/are missing. Please upload required files")
		for file in PEPfiles:
			if not Path(file).suffix == ".in":
				raise Exception("Sorry, Peptide sequnece must be saved as .in files. Check sample for acceptable format, then setroute again and upload files")

	elif 'Proteins' in listfPDB:
		PROfiles = os.listdir(os.path.join(fPDB, 'Proteins'))
		if not len(PROfiles) > 0:
			raise Exception("Needed protein(s) file(s) is/are missing. Please upload required files")
		for file in PROfiles:
			if not Path(file).suffix == ".pdb":
				raise Exception("Sorry, Protein file must be pdb files. Check sample for acceptable format, then setroute again and upload required files")

	# Check preferred forcefield selection, forcefield folder and its content
	for itemf in listfPDB:
		if itemf == 'Gromos' or itemf == 'Opls' or itemf == 'Amber' or itemf == 'Charmm':
			Ligsff = itemf
			
	if not (Ligsff == " " or Ligsff == "Select"):
		printNote("PLEASE NOTE")
		print(f"Your preselected forcefield group is {Ligsff}")
	elif Ligsff == "Select":
		print("Your forcefield will be selected interactively")
	print('\n')

	# Set the starting date and time
	printNote("IF YOU GET HERE, FILES CHECK WAS MOST PROBABLY SUCCESSFUL!!!")
	printNote(
		"CAREFULLY CHECK THE CHECK RESULTS ABOVE AND TYPE 'YES/y' TO CONFIRM OR 'NO/n' TO ABORT")
	confirmation = input("Response is: ")
	if not (confirmation.lower() == "yes" or confirmation.lower() == "y"):
		raise Exception("Process aborted by the User. Please rerun from SetRoute")
	else:
		ppmoretime = datetime.now()
		Tstart = ppmoretime.strftime("%B %d, %Y %H:%M:%S")
		if 'Peptides' in os.listdir(os.path.join(scriptDIR, 'Uploads')):
			print(f'Generation of MDS input files for Peptides begins: {Tstart}')
		else:
			print(f'Generation of MDS input files for Proteins begins: {Tstart}')

	# Get user imput for project name
	while True:
		name = tinput("Suppy a name for the current project: ", defaults[4], "PPmore")
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
			print(f"Your Unique ID is: {ID}")
			print(f"Your current work will be stored in {foldername}")
			break
		else:
			continue

	def pptopsolgen(pname, rname, tleapfile, rundefaults, mdsfgerrors):
		# Get needed variables as global
		global Ligsff
		global tff
		global PLD
		global fPDB

		f = open("pperror.txt", "a")
		udefaults = rundefaults
		errors = mdsfgerrors

		# Set workhost directory similar to the main function
		workhost_dir = Path.cwd()

		#######################################################
		# GENERATING GROMACS AND ALTERNATIVE AMBER TOPOLOGIES #
		#######################################################
		print('\n')
		os.mkdir("ppTOPS")
		os.chdir("ppTOPS")

		pptops_dir = os.path.join(workhost_dir, "ppTOPS")
		print(f"Protein topology data directory set to {pptops_dir}")

		printNote("Generating topology and parameter files...")

		selff = udefaults[0]
		selwater = udefaults[1]

		shutil.copy(os.path.join(workhost_dir, pname), './')

		# Generating Protein topology and parameter files

		while True:
			RFtop, RFpdb, RFposre = receptopol(pname, rname, selff, selwater)

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
				print('\n')
				print(f"{tff} forcefield in your topol file does not match the preselected group: {Ligsff}")

				printNote("To rerun, Type YES/y. Otherwise press ENTER to continue with current selection")
				response = tinput("Response: ", 30, "n")
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
						print('\n')
						print(f"{tff} does not match any known forcefield group. Please check")
						printNote("If using a self created or modified forcefield, ensure to use standard naming convention, similar to GROMACS naming")
						response = tinput("To abort, Type YES/y. To continue anyway, press ENTER: ", udefaults[4], "y")
						if response.lower() == "yes" or response.lower() == "y":
							raise Exception("Process aborted. Make necessary corrections and Rerun setup")
						else:
							break

					Ufolders = os.listdir(fPDB)
					if not ('Amber' in Ufolders or 'Charmm' in Ufolders or 'Gromos' in Ufolders or 'Opls' in Ufolders):
						raise Exception("No folder or file matching forcefield that can be renamed was found")
					else:
						for sff in Ufolders:
							if sff == 'Amber' or sff == 'Charmm' or sff == 'Gromos' or sff == 'Opls':
								os.rename(os.path.join(fPDB, sff), os.path.join(fPDB, Ligsff))

					print(f"Your preferred forcefield group has been changed to {Ligsff}")
					break
			else:
				print(f"Your forcefiled as contained in topol.top file is {tff}")
				break

		rftopfile = "n" + RFtop
		rfpdbfile = "n" + RFpdb
		try:
			os.rename(RFtop, rftopfile)
			os.rename(RFpdb, rfpdbfile)
		except Exception as e:
			printWarning(e)

		# Check if regenerating and decide either to continue interactively or not
		if errors > 0:
			print('\n')
			printNote("This is a rerun of failed MDS input files generation with Interactive mode")
			print("Type YEs/y to continue interactviely, Or press ENTER to change to Noninteractive mode")
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
		if not (udefaults[0] == "select" and udefaults[1] == "select" and errors == 0):
			if udefaults[5] == "B":
				print('\n')
				if Path(tff).suffix == ".ff":
					udefaults[0] = Path(tff).stem
				else:
					udefaults[0] = tff
				
				udefaults[1] = defaults2(rftopfile)
				if udefaults[1] == "none":
					print("No water model was detected for your system")
				else:
					print(f"Your water model as contained in topol.top file is {udefaults[1]}")
				
				printNote("Your selected default values are as follows: ")
				print(f"		Default forcefield: {udefaults[0]}")
				print(f"		Default water model: {udefaults[1]}")
				print(f"		Default editconf -bt: {udefaults[2]}")
				print(f"		Default editconf -d: {udefaults[3]}")
				print(f"		Default input timeout: {udefaults[4]}")
				print("		Default mode: non-interactive")
				udefaults[5] = "C"
			else:
				udefaults[1] = defaults2(rftopfile)
				if udefaults[1] == "none":
					print("No water model was detected for your system")
				else:
					print(f"Your water model as contained in topol.top file is {udefaults[1]}")

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
		nrecopen = open(rftopfile, "r")
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
			insertdetails(rftopfile, 'mt-topol2.itp', insertUDP)
			os.remove('mt-topol2.itp')
		else:
			shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol3.itp'), './')
			insertP = "; Include water topology"
			insertdetails(rftopfile, 'mt-topol3.itp', insertP)
			os.remove('mt-topol3.itp')

		# Now, let us genreated alternative amber topology file for use with gromacs
		time.sleep(5)
		print('\n')
		print("Generating alternative topology file...")

		# Gather needed files to generate amber topologies
		shutil.copy(os.path.join(workhost_dir, tleapfile), './')

		if 'Peptides' in os.listdir(fPDB):		
			shutil.copy(os.path.join(workhost_dir, 'ppmore.in'), './')
		else:
			shutil.copy(os.path.join(workhost_dir, pname), './')
			os.rename(pname, 'ppmore.pdb')

		tlpgro = "tlp" + rname + ".gro"
		tlptop = "tlp" + rname + ".top"
		tlpsolvgro = "tlpSolv" + rname + ".gro"
		tlpsolvtop = "tlpSolv" + rname + ".top"

		# Check if tleap source file is present and attempt to generate complex via tleap or use alternative
		ppatopsdir = os.listdir()
		if not tleapfile in ppatopsdir:
			printNote("An alternative amber generated GROMACS compatible topology file can not be generated")

		else:
			print("Generating an alternative GROMACS compatible topology file ...")
			ppmore_gro, ppmore_top = complexgen(tleapfile)

			if not (ppmore_gro == 'Complex.gro' and ppmore_top == 'Complex.top'):
				printNote("Generation of alternative topology failed")

			else:
				printNote("Generation of alternative topology successful")
				try:
					os.rename(ppmore_gro, tlpgro)
					os.rename(ppmore_top, tlptop)
				except Exception as e:
					pass

		for file in ppatopsdir:
			if file == 'tlpSolvated.gro':
				os.rename('tlpSolvated.gro', tlpsolvgro)
			elif file == 'tlpSolvated.top':
				os.rename('tlpSolvated.top', tlpsolvtop)
			elif file == 'posre_complex.itp':
				os.rename('posre_complex.itp', 'posre_tlptopol.itp')

		os.chdir('../')
		time.sleep(5)

		####################################
		# SOLVATION OF PROTEIN STRUCTURE   #
		####################################
		# Creating solvation directory and populating it with required files
		print('\n')
		print("Generating Solvated Structure ...")

		solname = "Solvation"
		os.mkdir(solname)
		os.chdir(solname)

		listTpdir = os.listdir(pptops_dir)
		for Tp in listTpdir:
			if Path(Tp).suffix == ".itp" or Path(Tp).suffix == ".top" or Path(Tp).suffix == ".gro" or Path(Tp).suffix == ".pdb":
				shutil.copy(os.path.join(workhost_dir, "ppTOPS", Tp), './')

		# Determine if the forcefield belong to the selected group - amber, charmm, gromos and opls
		if not (tff[0:4].lower() == Ligsff.lower() or tff[0:5].lower() == Ligsff.lower() or tff[0:6].lower() == Ligsff.lower()):
			print(f"The {tff} forecfield in topol.top did not match the preferred forcefield group, {Ligsff}")
			
		if tff[0:5].lower() == 'amber':
			print(f"{tff} is an amber forcefield")

		elif tff[0:6].lower() == 'charmm':
			print(f"{tff} is a charmm forcefield")

		elif tff[0:6].lower() == 'gromos':
			print(f"{tff} is a gromos forcefield")
			
		elif tff[0:4].lower() == 'opls':
			print(f"{tff} is an opls forcefield")

		# Prepare a version of amber tleap generated topology for suitable use with Gromacs
		chkSoldir = os.listdir()
		if (tlpgro in chkSoldir and tlptop in chkSoldir):
			tlpgmxtopol = TLtopol(rftopfile, tlptop, tff)
			if not tlpgmxtopol == "tlptopol.top":
				print("Reverting to the generated tleap topology file")
			elif "tlptopol.top" in chkSoldir and "posre_complex.itp" in chkSoldir:
				shutil.copy(tlpgmxtopol, 'xtlptopol.top')
				os.rename('posre_complex.itp', 'posre_tlptopol.itp')
			else:
				shutil.copy(tlpgmxtopol, 'xtlptopol.top')

		# Time to prepared solvated complex
		if not udefaults[1] == "select":
			selwater = udefaults[1]
		else:
			selwater = defaults2("topol.top")
		selbt = udefaults[2]
		seld = udefaults[3]

		chkSoldir = os.listdir()
		if (rftopfile in chkSoldir and rfpdbfile in chkSoldir):
			print(f"{rname} Solvation in progress.....")

			shutil.copy(rftopfile, 'topol.top')
			grosolvated = solvation(rfpdbfile, 'topol.top', selwater, selbt, seld)
			
			if not grosolvated == 'fsolvated.gro':
				print(f"{rname} Solvation unsuccessful")

			else:
				print(f"{rname} Solvation was successful")

		else:
			print(f"Required file(s) for {rname} solvation not found / generated. Solvation can not continue")
			printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

		os.chdir('../')
		time.sleep(5)

		######################################
		# GATHERING FILES FOR MD SIMULATIONS #
		######################################
		# Making gmxmds directory and populate it with needed data
		print('\n')
		print("Gathering and sorting files needed for MDS run ...")

		gmxname = "gmxmds"
		os.mkdir(gmxname)
		os.chdir(gmxname)

		soldirlist = os.listdir(os.path.join(workhost_dir, 'Solvation'))
		for file in soldirlist:
			if Path(file).suffix == ".itp":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == "topol.top" or file == "tlptopol.top" or file == "fsolvated.gro" or file == "ufsolvate.gro":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif not ("fsolvated.gro" in soldirlist and defaults[1] == "none"):
				if file == 'tlpSolvated.gro' or file == 'tlpSolvated.top':
					shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')

		# Removing files not needed for MDS
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

		if "posre_complex.itp" in os.listdir():
			os.rename('posre_complex.itp', 'posre_tlptopol.itp')
			
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
				shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'backup'))
			else:
				shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'not4mds'))

		# Performing final processing for tlptopol.top file, if it has not been used instead of topol.top file
		if "topol.top" in os.listdir() and "tlptopol.top" in os.listdir() and not defaults[1] == "none":
			print("Final processing of alternative topology file in progress .....")
			tlpfinal("tlptopol.top", "topol.top")

		for bkfile in os.listdir():
			if bkfile[0:2].lower() == "bk":
				shutil.move(bkfile, os.path.join(work_dir, 'gmxmds', 'backup'))

		# If required, new restraint file can now be generated
		listgmxmds = os.listdir()
		grosolvated = "none"

		if "fsolvated.gro" in listgmxmds:
			print('\n')
			grosolvated = "fsolvated.gro"

			# Generating new restraint file if needed
			response = tinput("To generate a new restraint interactively, type YES/y, or press ENTER: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				restrfile = udrestraint(grosolvated)
				if not (restrfile == "posre_udp.itp" or "posre_udp.itp" in os.listdir()):
					printNote("No user defined restraint file has been generated")
				else:
					printNote("Your generated restraint file is 'posre_udp.itp'")
					printNote("To use it in MDS, define -DPOSRES_UDP in .mdp files OR if already define -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")

			else:
				printNote("No user defined restraint file has been generated. If need be, generate it manually")

		elif "ufsolvate.gro" in listgmxmds:
			grosolvated = "ufsolvate.gro"
			print('\n')
			print("Manual solvation and/or generation of restraint is required if necessary")

		elif "tlpSolvated.gro" in listgmxmds:
			grosolvated = "tlpSolvated.gro"
			print('\n')
			print("Manual generation of restraint for tleap generated solvated complex is required, if necessary")

		else:
			grosolvated = "none"
			print('\n')
			print("No solvated or unsolvated complex was detected in gmxmds subfolder. Please check")

		# Final Processing of generated files for better formating
		print("Performing Final Checking of Processed MDS files...")
		for mdsfile in os.listdir():
			if not os.path.isdir(mdsfile):
				gmxmdsFClean(mdsfile)
		
		print("gmxmds subfolder has been populated")
		f.close()
		return grosolvated, udefaults

	time.sleep(5)

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

	if 'Peptides' in listfPDB:		
		print('\n')
		print(f"pyGROMODS Peptides project directory set to {work_dir}")
	
		mdsfgerrors = 0
		pepname = "Peptide_"
		ppcount = 0
		while ppcount < len(PEPfiles):
			# Create host directory for each pair of receptor and ligand
			pep = PEPfiles[ppcount]
			count = ppcount + 1
			pep_dir = pepname + str(count)

			if os.path.isdir(os.path.join(work_dir, pep_dir)):
				shutil.rmtree(os.path.join(work_dir, pep_dir))
				time.sleep(5)

			os.mkdir(pep_dir)
			os.chdir(pep_dir)
		
			workhost_dir = os.path.join(work_dir, pep_dir)
			print(f"Current project host directory set to {workhost_dir}")

			####################################################
			# GENERATE PROTEIN STRUCTURE FROM PEPTIDE SEQUENCE #
			####################################################

			rname = 'peptide' + str(count)
			pname = 'new' + rname + '.pdb'

			os.mkdir('pepPDB')
			os.chdir('pepPDB')

			# Copy needed files into the directory
			shutil.copy(os.path.join(fPDB, 'Peptides', pep), './')
			shutil.copy(os.path.join(PLD, 'peptopdbtleap.in'), './')
			shutil.copy(pep, 'pep.in')

			peppdb_dir = os.path.join(workhost_dir, 'pepPDB')
			print(f"Peptide structure data directory set to {peppdb_dir}")

			printNote("Convertion of Peptide sequence into pdb structure file...")

			# Run amber tleap to convert to pdb
			fpep = open("ppdberror.txt", "a")
			try:
				subprocess.run(['tleap', '-s', '-f', 'peptopdbtleap.in'], stderr=subprocess.STDOUT, stdout=fpep, check=True, text=True)
			except subprocess.SubprocessError as e:
				printWarning("Peptide preparation failed with above error. The process will be aborted")
				printNote("Check the above error and make corrections and rerun. Read README.md for help")
				fpep.close()
				raise Exception(f"Process Aborted with error: {e}. Make necessary corrections and restart")

			shutil.copy('sequence.pdb', '../')
			fpep.close()
			os.chdir('../')
			os.rename('sequence.pdb', pname)

			shutil.copy(os.path.join(fPDB, 'Peptides', pep), './')
			os.rename(pep, 'ppmore.in')
			shutil.copy(os.path.join(PLD, 'peptoptleap.in'), './')
			tleapfile = "peptoptleap.in"

			pepsolvated, udefaults = pptopsolgen(pname, rname, tleapfile, defaults, mdsfgerrors)
			defaults = udefaults

			# Final Processing of generated files for better formating
			print("Performing Final Checking of Processed MDS input files...")
			if pepsolvated == "fsolvated.gro" or pepsolvated == "ufsolvate.gro":
				os.chdir('../')
				for fitem in os.listdir():
					if not os.path.isdir(fitem):
						os.remove(fitem)
				ppcount += 1
				os.chdir('../')

			else:
				mdsfgerrors += 1
				decision, newdefaults = gmxmdsFEChecks(defaults)
				if decision == "abort":
					raise Exception("Process aborted to allow for user to checks and make corrections")
				elif decision == "continue":
					printWarning("You have chosen to continue with the process despite errors. Please check")
					os.chdir('../')
					for fitem in os.listdir():
						if not os.path.isdir(fitem):
							os.remove(fitem)
					ppcount += 1
					os.chdir('../')
				elif decision == "regenerate":
					printNote("Getting files and defauts values ready for Regeneration...")
					defaults = newdefaults
					os.chdir('../../')

	if 'Proteins' in listfPDB:
		print('\n')
		print(f"pyGROMODS Proteins project directory set to {work_dir}")

		mdsfgerrors = 0
		proname = "Protein_"
		prcount = 0
		while prcount < len(PROfiles):
			# Create host directory for each protein to be setup for MDS
			pro = PROfiles[prcount]
			count = prcount + 1
			pro_dir = proname + str(count)

			if os.path.isdir(os.path.join(work_dir, pro_dir)):
				shutil.rmtree(os.path.join(work_dir, pro_dir))
				time.sleep(5)

			os.mkdir(pro_dir)
			os.chdir(pro_dir)

			workhost_dir = os.path.join(work_dir, pro_dir)
			print(f"Current project host directory set to {workhost_dir}")
			print("Copying needed files into the working directory...")

			rrname = 'protein' + str(count)
			prname = 'new' + rrname + '.pdb'
			shutil.copy(os.path.join(fPDB, 'Proteins', pro), './')
			shutil.copy(pro, prname)
			shutil.copy(os.path.join(PLD, 'protoptleap.in'), './')
			tleapfile = "protoptleap.in"

			protsolvated, udefaults = pptopsolgen(prname, rrname, tleapfile, defaults, mdsfgerrors)
			defaults = udefaults

			# Final Processing of generated files for better formating
			print("Performing Final Checking of Processed MDS input files...")
			if protsolvated == "fsolvated.gro" or protsolvated == "ufsolvate.gro":
				os.chdir('../')
				for fitem in os.listdir():
					if not os.path.isdir(fitem):
						os.remove(fitem)
				prcount += 1
				os.chdir('../')

			else:
				mdsfgerrors += 1
				decision, newdefaults = gmxmdsFEChecks(defaults)
				if decision == "abort":
					raise Exception("Process aborted to allow for user to checks and make corrections")
				elif decision == "continue":
					printWarning("You have chosen to continue with the process despite errors. Please check")
					os.chdir('../')
					for fitem in os.listdir():
						if not os.path.isdir(fitem):
							os.remove(fitem)
					prcount += 1
					os.chdir('../')
				elif decision == "regenerate":
					printNote("Getting files and defauts values ready for Regeneration...")
					defaults = newdefaults
					os.chdir('../../')

	print('\n')
	printNote("PLEASE NOTE:")
	print("The files in folder 'not4mds' of 'gmxmds' are considered not necessary for MDS. Please check")
	print("The files in subfolder 'backup' of 'gmxmds' can be used to replace relevant ones in gmxmds folder for MDS")

	print('\n')
	printNote("Setup with PPmore route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS")

	# Set the ending date and time
	print('\n')
	ppmoretime2 = datetime.now()
	Tend = ppmoretime2.strftime("%B %d, %Y %H:%M:%S")
	if 'Peptides' in os.listdir(os.path.join(scriptDIR, 'Uploads')):
		print(f'Generation of MDS input files for Peptides ends: {Tend}')
	else:
		print(f'Generation of MDS input files for Proteins ends: {Tend}')
