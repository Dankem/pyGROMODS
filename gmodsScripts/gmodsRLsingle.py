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

def RLsingle(appDIR, gmxDIR, fdefaults, dictff):
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

	# Set global variable to automatically run pdb2gmx, editconf and others
	# forcefields & water (select), -bt (triclinic), -d (0.1), and timeout (60)
	defaults = fdefaults 
	
	printNote("Let us check again your selected default values ..... ")	
	print(f"Default forcefield is {defaults[0]}")
	print(f"Default water model is {defaults[1]}")
	print(f"Default editconf -bt option is {defaults[2]}")
	print(f"Default editconf -d option is {defaults[3]}")
	print(f"Default timeout for input() request is {defaults[4]}")
	print('\n')

	if defaults[5] == "A":
		printNote("Your selected default mode: Interractive")
		response = tinput("To revert to Noninteractive type YES/y: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[5] = "B"
			printNote("You have changed to Noninteractive mode")
			print("Your preferred forcefield and water model will be autodetected following your first interactive selection")
	else:
		printNote("Your selected default mode: Noninterractive")
		response = tinput("To revert back to Interactive type YES/y: ", defaults[4], "n")
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

	LIGfolders = os.listdir(os.path.join(fPDB, 'Ligands'))
	RECfolders = os.listdir(os.path.join(fPDB, 'Receptors'))
	TOPfolders = os.listdir(os.path.join(fPDB, 'Ligsff'))

	# Check receptor, ligand and, if available, topology files for correctness
	if not (len(RECfolders) > 0 and len(LIGfolders) > 0 and len(TOPfolders) > 0):
		raise Exception("Needed ligand, receptor and topology subfolder(s) is/are missing. Rerun the setup")

	for fileL in LIGfolders:
		Lsubdir = os.path.join(fPDB, 'Ligands', fileL)
		if os.path.isdir(Lsubdir):
			Lsdirlist = os.listdir(Lsubdir)
			if not len(Lsdirlist) == 1:
				print(f"Expected an instance of one ligand per subdirectory, but detected {len(Lsdirlist)}") 
				raise Exception("Sorry, unacceptable number of ligand detected. Set route and rerun")
			for Lfile in Lsdirlist:
				if not (Path(Lfile).suffix == ".pdb" or Path(Lfile).suffix == ".mol2"):
					raise Exception("Sorry, ligand files must be pdb or mol2 files. Upload required files")
		else:
			raise Exception("Sorry, file detected where subdirectory should be, you might need to select different route")

	for fileR in RECfolders:
		Rsubdir = os.path.join(fPDB, 'Receptors', fileR)
		if os.path.isdir(Rsubdir):
			Rsdirlist = os.listdir(Rsubdir)
			if not len(Rsdirlist) == 1:
				print("Expected an instance of one receptor per subdirectory, but detected", len(Rsdirlist)) 
				raise Exception("Sorry, unacceptable numbers of receptors detected. Set route and rerun")
			for Rfile in Rsdirlist:
				if not (Path(Rfile).suffix == ".pdb"):
					raise Exception("Sorry, Receptor files must be pdb file. Upload required files")
		else:
			raise Exception("Sorry, file detected where subdirectory should be, you might need to select different route")

	numfc = 0
	numtp = 0
	for fileT in TOPfolders:
		Tsubdir = os.path.join(fPDB, 'Ligsff', fileT)
		if os.path.isdir(Tsubdir):
			Tsdirlist = os.listdir(Tsubdir)
			wmodel = []
			numff = 0
			numfw = 0
			for Tfile in Tsdirlist:
				Tdir = os.path.join(fPDB, 'Ligsff', fileT, Tfile)
				if os.path.isdir(Tdir):
					if Tfile.lower() == "amber" or Tfile.lower() == "charmm" or Tfile.lower() == "gromos" or Tfile.lower() == "opls":
						printNote("PLEASE NOTE!!!")
						print(f">>>>> {Tfile} is your preferred forcefield group for this pair of receptor - ligand")
						print(">>>>> Make sure the uploaded ligand topology was generated using similar forcefield")

					elif Tfile.lower() == "select":
						printNote("PLEASE NOTE!!!")
						print(f"Interactively {Tfile} forcefield for this pair of receptor - ligand")
						print("You must chose forcefield compatible with your uploaded ligand topology")

					else:
						printWarning("NOTE WARNING!!!")
						print(f">>>>> {Tfile} could not be identify with any forcefields in gromacs")
						print(">>>>> The uploaded ligand topology will be ignored for this pair")
						shutil.rmtree(Tdir)
						continue

					listTdir = os.listdir(Tdir)
					if not len(listTdir) == 1:
						printWarning("NOTE WARNING!!!")
						print(">>>>> Expecting one topology file per subdirectory, but found", len(listTdir))
						print(">>>>> The uploaded ligand topology will be ignored for this pair")
						shutil.rmtree(Tdir)
						continue
 
					for Topf in listTdir:
						if not (Path(Topf).suffix == ".itp" or Path(Topf).suffix == ".top"):
							printWarning("NOTE WARNING!!!")
							print(">>>>> Sorry, unacceptable topology file format detected. Must be .itp or .top")
							print(">>>>> The uploaded ligand topology will be ignored for this pair")
							shutil.rmtree(Tdir)
							break

						else:
							checkindex = indexoflines(os.path.join(fPDB, 'Ligsff', fileT, Tfile, Topf))
							filecheck = open(os.path.join(fPDB, 'Ligsff', fileT, Tfile, Topf), "r")
							readcheck = filecheck.readlines()
							checklist = ['atomtypes', 'moleculetype', 'system']
							nckl = 0
							eckl = 0
							for ckl in checklist:
								try:
									if not ckl in readcheck[checkindex[ckl]].split():
										nckl += 1
								except Exception as e:
									eckl += 1
							filecheck.close()

							if nckl == 0 and eckl == 0:
								print(f"{Topf}, an acceptable uploaded ligand topology file detected")
								numtp += 1
								break
							else:
								printWarning("NOTE WARNING!!!")
								print(">>>>> Checking a user supplied ligand topology file failed with errors")
								print(">>>>> This file may lack ['atomtypes'], ['moleculetype'] and/or [ system ] subheader")
								print(">>>>> The uploaded ligand topology will be ignored for this pair")
								shutil.rmtree(Tdir)
								break
								
				elif os.path.isfile(Tdir):
					# Check for forcefield selection and generate water model list
					if Tfile[0:5].capitalize() == "Amber" or Tfile[0:6].capitalize() == "Charmm" or Tfile[0:6].capitalize() == "Gromos" or Tfile[0:4].capitalize() == "Opls":
						numff += 1
						wmodel = dictff[Tfile]
						print(f"{Tfile} forcefiled detected for this pair of receptor - ligand")

					elif Tfile.lower() == "select_ff":
						numff += 1
						wmodel.append("select_ww")
						print(f"Interactively '{Tfile[0:6]}' forcefield detected for this pair of receptor - ligand")

					else:
						pass	

			# Check if water model is present and correspond to selected ff
			for Wfile in Tsdirlist:
				if Wfile in wmodel and not Wfile.lower() == "select_ww":
					numfw += 1
					print(f"{Wfile} water model detected for this pair of receptor - ligand")

				elif Wfile in wmodel and Wfile.lower() == "select_ww":
					numfw += 1
					print(f"Interactively '{Wfile[0:6]}' water model detected for this pair of receptor - ligand")

			revert_mode = 0
			if numff > 1 or numfw > 1:
				printWarning("NOTE WARNING!!!")
				print(">>>>> More than one forcefield or water model detected for this pair")
				print(">>>>> We shall attempt to revert to interactive mode for this pair")
				revert_mode += 1

			elif numff == 0 or numfw == 0:
				printWarning("NOTE WARNING!!!")
				print(">>>>> Needed forcefield and/or its corresponding water model is missing for this pair")
				print(">>>>> We shall attempt to revert to interactive mode for this pair")
				revert_mode += 1

			if revert_mode > 0:			
				for Tfile in Tsdirlist:
					Tdir = os.path.join(Tsubdir, Tfile)
					if os.path.isdir(Tdir):
						os.rename(Tdir, os.path.join(Tsubdir, "Select"))
					else:
						os.remove(Tdir)

				fileff = open(os.path.join(Tsubdir, 'select_ff'), "w")
				numff = 0
				fileww = open(os.path.join(Tsubdir, 'select_ww'), "w")
				numfw = 0
				fileff.close()
				fileww.close()
				numfc += 1

			else:
				numfc += 1
		else:
			raise Exception("Sorry, a file detected where subdirectory should be. Please setroute again and rerun")
		print('\n')

	if not len(LIGfolders) == len(TOPfolders):
		raise Exception("Matching pairs of ligand and forcefield/ligand topology not detected. Please cross check")

	if numfc == len(LIGfolders) and numtp == 0:
		TFF = 0
		printNote("No user supplied topology file(s) was detected for all pairs of Receptor:Ligand")

	elif numtp == len(LIGfolders) and numfc == len(LIGfolders):
		TFF = 1
		printNote("User supplied topology file(s) was detected for all pairs of Receptor:Ligand")

	else:
		TFF = 3
		printNote("One or more user supplied topology file(s) will be ignored")
		response = tinput("Do you want to proceed anyway? YES/y or NO/n: ", defaults[4], "y")
		if not (response.lower() == "yes" or response.lower() == "y"):
			raise Exception("Process Abort!!!. Check your files and rerun")

	if not len(RECfolders) == len(LIGfolders):
		raise Exception("Matching pairs of receptor and ligand could not be detected. Please cross check")
	else:
		printNote("Detected matching pairs of receptor and ligand. Solvated complex will be generated for each pair")
	print('\n')

	# Set the starting date and time
	printNote("IF YOU GET HERE, FILES CHECK WAS MOST PROBABLY SUCCESSFUL!!!")
	printNote("CAREFULLY CHECK THE CHECK RESULTS ABOVE AND TYPE 'YES/y' TO CONFIRM OR 'NO/n' TO ABORT")
	confirmation = input("Response is: ")
	if not (confirmation.lower() == "yes" or confirmation.lower() == "y"):
		raise Exception("Process aborted by the User. Please rerun from SetRoute")
	else:
		rlsingletime = datetime.now()
		Tstart = rlsingletime.strftime("%B %d, %Y %H:%M:%S")
		print(f'Generation of MDS input files begins: {Tstart}')
	print('\n')
	time.sleep(5)

	# Get user imput for project name
	while True:
		name = tinput("Suppy a name for the current project: ", defaults[4], "RLsingle")
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
	print('\n')

	# Create working directory using genrated project name
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
		printNote("Samples of .mdp, peptide, and Tleap source files have been copied to your workspace directory")
	except:
		pass
	print('\n')

	os.chdir('../')
	time.sleep(2)

	RLname = "RLsingle_"
	npairs = 0
	mdsfgerrors = 0
	while npairs < len(RECfolders):
		# Set global variable to access user supplied topology file(s)
		RLcount = npairs + 1
		RL = RECfolders[npairs]
		Ligsff = " "
		tff = " "
		ligsff_dir = " "
		M = TFF

		# Create host directory for each pair of receptor and ligand
		rls_dir = RLname + str(RLcount)
		if os.path.isdir(os.path.join(work_dir, rls_dir)):
			shutil.rmtree(os.path.join(work_dir, rls_dir))
			time.sleep(5)

		workhost_dir = os.path.join(work_dir, rls_dir)
		print(f"Current project host directory set to {workhost_dir}")

		os.mkdir(rls_dir)
		os.chdir(rls_dir)

		# Check preferred forcefield selection for the pair under consideration
		ligtop = "ligtop" + str(RLcount)
		ligtopdir = os.path.join(fPDB, 'Ligsff', ligtop)
		listligtopdir = os.listdir(ligtopdir)
		for itemf in listligtopdir:
			if itemf == 'Gromos' or itemf == 'Opls' or itemf == 'Amber' or itemf == 'Charmm' or itemf == 'Select':
				Ligsff = itemf
				break
			elif itemf[0:4].capitalize() == "Opls":
				Ligsff = "Opls"
				break
			elif itemf[0:5].capitalize() == "Amber":
				Ligsff = "Amber"
				break
			elif itemf[0:6].capitalize() == "Charmm":
				Ligsff = "Charmm"
				break
			elif itemf[0:6].capitalize() == "Gromos":
				Ligsff = "Gromos"
				break
			elif itemf[0:6].capitalize() == "Select":
				Ligsff = "Select"
				break

		ligsff_dir = os.path.join(fPDB, 'Ligsff', ligtop, Ligsff)
		if TFF == 3:
			if os.path.isdir(ligsff_dir) and len(os.listdir(ligsff_dir)) > 0:
				TFF = 1
				printNote("Detected uploaded ligand topology file for this pair")
			else:
				TFF = 0		
				printNote("No uploaded ligand(s) topology(ies) was detected for this pair")

		if not (Ligsff == " " or Ligsff == "Select"):
			printNote("PLEASE NOTE")
			print(f"Your preselected forcefield group for this pair is {Ligsff}")

		elif Ligsff == "Select":
			printNote("PLEASE NOTE")
			print("Forcefield will be selected interactively for this pair")

		################################################
		# GENERATING RECEPTOR TOPOLOGIES AND STRUCTURE #
		################################################
		print('\n')
		os.mkdir('Receptors')
		os.chdir('Receptors')

		receptor_dir = os.path.join(workhost_dir, 'Receptors')
		print(f"Receptor data directory set to {receptor_dir}")

		printNote("Generating Protein topology and parameter files...")
		time.sleep(2)

		# Set variables for generating Protein topologies and parameter files
		Rname = "receptor"
		drname = Rname + str(RLcount)
		selff = defaults[0]
		selwater = defaults[1]
		newff = ""
		newWater = ""

		# Get the list of all the avaialble water models
		fwmodels = []
		ffields = []
		for key in dictff:
			if not (key in ffields and key == "select"):
				ffields.append(key)
			for fw in dictff[key]:
				if not (fw in fwmodels and fw == "select"):
					fwmodels.append(fw)

		# Get the preselected ff and water model for current pair
		preff = ""
		prewater = ""
		for itemfw in listligtopdir:
			ffww_dir = os.path.join(ligtopdir, itemfw)
			if not os.path.isdir(ffww_dir):
				if itemfw in ffields:
					preff = itemfw
				elif itemfw in fwmodels:
					prewater = itemfw
				elif itemfw == "select_ff":
					preff = "select"
				elif itemfw == "select_ww":
					prewater = "select"
			else:
				pass

		# Overide the default generate menu selections if ff and water are not 'select'
		if selff == "select" and selwater == "select" and mdsfgerrors == 0:
			printNote("This appears to be an interactive mode")
			print(f"Current default values for -bt is: {defaults[2]}")
			print(f"Current default values for -d is: {defaults[3]}")
			print(f"Current default values for timeout is: {defaults[4]}")
			response = tinput("To adjust these values for current protein - ligand complex, type YES/y: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				defaults[2], defaults[3], defaults[4] = defaults1() 
		elif selff == "select" and selwater == "select" and mdsfgerrors > 0:
			printNote("This appears to be a rerun of failed MDS input files generation")
			printNote("The following are the current default values")
			print(f"Current default values for forcefield is: {defaults[0]}")
			print(f"Current default values for water model is: {defaults[1]}")
			print(f"Current default values for -bt is: {defaults[2]}")
			print(f"Current default values for -d is: {defaults[3]}")
			print(f"Current default values for timeout is: {defaults[4]}")
			if defaults[5] == "A":
				print("Current default mode of generation is: Interactive")
			else:
				print("Current default mode of generation is: Noninteractive")
		else:
			selff = preff
			selwater = prewater
    
		# Populate the Receptor folder and generate topologies
		Rdir = os.path.join(fPDB, 'Receptors', RL)
		Rdirlist = os.listdir(Rdir)
		for rep in Rdirlist:
			shutil.copy(os.path.join(fPDB, 'Receptors', RL, rep), './')

		while True:
			RFtop, RFpdb, RFposre = receptopol(rep, drname, selff, selwater)
				
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

			sff_folders = os.listdir(os.path.join(fPDB, 'Ligsff', ligtop))
			newff = Path(tff).stem
			newWater = defaults2(RFtop)

			if not newff in sff_folders:
				if not (tff[0:4].capitalize() == Ligsff or tff[0:5].capitalize() == Ligsff or tff[0:6].capitalize() == Ligsff):
					print(f"{tff} forcefield does not match the preselected forcefield group: {Ligsff}")
					print(f"By default, your forcefield group will be changed to match {tff} any uploaded ligand topology will be ignored")

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
							print(f"{tff} does not match any known forcefield group. A rerun is recommended")
							printNote("If using a self created or modified forcefield, use standard naming convention for forcefields")
							print("To abort, Type YES/y. To continue anyway, press ENTER")
							response = tinput("Response: ", 30, "n")
							if response.lower() == "yes" or response.lower() == "y":
								raise Exception("Process aborted. Make necessary corrections and Rerun setup")
							else:
								break

						for sff in sff_folders:
							if sff == 'Amber' or sff == 'Charmm' or sff == 'Gromos' or sff == 'Opls' or sff == 'Select':
								os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, Ligsff))
								print(f"Your selected forcefield group has been changed to {Ligsff}")

							elif (sff in dictff and sff == preff) or sff == "select_ff":
								os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, newff))
								print(f"Your selected forcefield has been changed to {tff}")

							elif (sff == prewater and newWater in dictff[Path(tff).stem]) or sff == "select_ww":
								os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, newWater))
								print(f"Your selected water model has been changed to {newWater}")

							else:
								pass

						ligsff_dir = os.path.join(fPDB, 'Ligsff', ligtop, Ligsff)
						if TFF == 1:
							print(f"Subdirectory for uploaded ligand topology is now: {ligsff_dir}")
							printNote("Type YES/y to use with uploaded ligand topology")
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
					print(f"{tff} forcefield match the preselected forcefield group: {Ligsff}")
					print(f"However, {tff} is different from the preselected forcefield for this pair")

					printNote("To rerun, Type YES/y. Otherwise press ENTER to continue with current selection")
					response = tinput("Response: ", defaults[4], "n")
					if response.lower() == "yes" or response.lower() == "y":
						selff = "select"
						selwater = "select"
						continue
					else:
						for sff in sff_folders:
							if sff in dictff and sff == preff:
								os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, newff))
								print(f"Your selected forcefield has been changed to {tff}")

							elif sff == prewater and newWater in dictff[Path(tff).stem]:
								os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, newWater))
								print(f"Your selected water model has been changed to {newWater}")

							else:
								pass

						ligsff_dir = os.path.join(fPDB, 'Ligsff', ligtop, Ligsff)
						if TFF == 1:
							print(f"Subdirectory for uploaded ligand topology is: {ligsff_dir}")
							printNote("Type YES/y to use with uploaded ligand topology")
							response = tinput("Response: ", defaults[4], "y")
							if not (response.lower() == "yes" or response.lower() == "y"):
								TFF = 0
								break
							else:
								TFF = 1
								break
						else:
							break

			elif not newWater in sff_folders:
				print(f"The detected {newWater} does not match the preselected water model: {selwater}")

				printNote("To rerun, Type YES/y. Otherwise press ENTER to continue with current selection")
				response = tinput("Response: ", defaults[4], "n")
				if response.lower() == "yes" or response.lower() == "y":
					selff = "select"
					selwater = "select"
					continue
				else:
					for sff in sff_folders:
						if sff == prewater and newWater in dictff[Path(tff).stem]:
							os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, newWater))
							print(f"Your selected water model has been changed to {newWater}")
							break
						else:
							pass
					break

			else:
				print(f"Your forcefiled as contained in topol.top file is {tff}")
				break

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
			print("Type YES/y to continue interactviely, Or press ENTER for Noninteractive mode")
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
			printNote("The following options are available to generate OPLS compatible ligand topology:")
			print("	1. Using platform OPLS compatible ligand topology generation - RECOMMENDED")
			print("	2. Using acpype OPLS compatible ligand topology generation")
			print("	3. Using platform default - ignore OPLS compatibility - NOT RECOMMENDED")

			while True:
				response = tinput("Choose your preferred option: ", defaults[4], "1")
				if response == '1' or response == '2' or response == '3':
					print(f"Option {response} Selected!")
					confirm = tinput("Type YES/y to confirm or press ENTER to choose a different option: ", 5, "y")
					if not (confirm.lower() == "yes" or confirm.lower() == "y"):
						continue
					else:
						if response == '1':
							printNote("Confirmed: The platform opls compatible ligand topology generation") 
						elif response == '2':
							printNote("Confirmed: The acpype opls compatible ligand topology generation") 
						elif response == '3':
							printNote("Confirmed: The platform default - not opls") 
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

		lig_dir = os.path.join(workhost_dir, 'Ligands')
		print(f"Ligands data directory set to {lig_dir}")
		time.sleep(2)

		Lname = "LIG"
		Lfolder = LIGfolders[npairs]
		dlname = Lname + str(RLcount)

		shutil.copy(os.path.join(PLD, 'ligtoptleap.in'), './')

		Ldir = os.path.join(fPDB, 'Ligands', Lfolder)
		Ldirlist = os.listdir(Ldir)
		lig = Ldirlist[0]
		shutil.copy(os.path.join(fPDB, 'Ligands', Lfolder, lig), './')

		LFgro, LFtop = ligtopol(lig, dlname)

		# Identify the pdb file to use subsequently, depending on uploaded ligand format
		print("Generating optimized ligand structure...")
		ligE = open('ligcleanerror.txt', 'a')

		nligname = dlname + "up.pdb"
		if Path(lig).suffix == ".pdb":
			try:
				subprocess.run(['gmx', 'editconf', '-f', lig, '-o', nligname], check=True, stderr=subprocess.STDOUT, stdout=ligE, text=True)
			except subprocess.SubprocessError as e:
				printNote("Trying again...")
				shutil.copy(lig, nligname)

		elif Path(lig).suffix == ".mol2":
			ligtopol_pdb = dlname + "new" + ".pdb"
			try:
				subprocess.run(['gmx', 'editconf', '-f', ligtopol_pdb, '-o', nligname], check=True, stderr=subprocess.STDOUT, stdout=ligE, text=True)
			except subprocess.SubprocessError as e:
				printNote("Trying again...")
				shutil.copy(ligtopol_pdb, nligname)
		ligE.close()

		# Prepare molecule number files
		ligsmnfile = dlname + "_mn.itp"
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

			uTdirlist = os.listdir(ligsff_dir)
			for utop in uTdirlist:
				shutil.copy(os.path.join(ligsff_dir, utop), './')

			uLdirlist = os.listdir(Ldir)
			for ulig in uLdirlist:
				shutil.copy(os.path.join(Ldir, ulig), './')

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
				print(utop,", the user supplied topology match the ligand named,", ulig)
				time.sleep(5)

				# Check the file for duplicate atom types with forcefield standard atomtypes
				uligcheck = Checkligtop(utop, tff)
				uctindex = indexoflines(uligcheck)

				# Generate the uLIGS_at.itp (atomtypes) and uLIGS_mt.itp (moleculetypes) file
				uTname = "u" + dlname
				uligand_at, uligand_mt = topolsplit(uligcheck, uTname, uctindex)
				utopname = uTname + ".top"
				os.rename(uligcheck, utopname)
				shutil.copy(utopname, '../')
				shutil.copy(uligand_at, '../')
				shutil.copy(uligand_mt, '../')

				os.chdir('../')

			else:
				print(f"{utop} user supplied topology does not match the ligand named {ulig}")
				printNote("By default, uploaded ligand topology for this pair will be ignored")
				printNote("Type YES/y to continue with the topology file. Otherwise press ENTER")
				response = tinput("Response: ", defaults[4], "n")
				if not (response.lower() == "yes" or response.lower() == "y"):
					TFF = 0
					printNote("Defaut used: The uploaded topology for the current ligand will be ignored")
					os.chdir('../')

				else:
					# Check the file for duplicate atom types with forcefield standard atomtypes
					uligcheck = Checkligtop(utop, tff)
					uctindex = indexoflines(uligcheck)

					# Generate the uLIGS_at.itp (atomtypes) and uLIGS_mt.itp (moleculetypes) file
					uTname = "u" + dlname
					uligand_at, uligand_mt = topolsplit(uligcheck, uTname, uctindex)
					utopname = uTname + ".top"
					os.rename(uligcheck, utopname)
					shutil.copy(utopname, '../')
					shutil.copy(uligand_at, '../')
					shutil.copy(uligand_mt, '../')

					os.chdir('../')

		# If need be, generate opls compatible ligand topology using acpype
		oplstop = " "

		if opls_route == 1:
			mergename = "opls" + dlname + ".top"
			oplstop = OPLStop(LFtop, tff, dlname)

			if not oplstop == mergename:
				printWarning("Selelcted option 1 for OPLS compatible ligand topology generation failed")
				printNote("Trying alternative option 2...")
				oplstop = OPLSacpype(LFtop, lig, dlname)
				if not oplstop == mergename:
					printWarning("Alternative option 2 for ligand topology generation also failed")
					printNote("Platform default will now be used")
					oplstop = LFtop		

		elif opls_route == 2:
			mergename = "opls" + dlname + ".top"
			oplstop = OPLSacpype(LFtop, lig, dlname)

			if not oplstop == mergename:
				printWarning("Selected option 2 for OPLS ligand topology generation failed")
				print("Trying alternative option...")
				oplstop = OPLStop(LFtop, tff, dlname)
				if not oplstop == mergename:
					printWarning("Alternative option 1 for ligand topology generation also failed")
					printNote("Platform default will now be used")
					oplstop = LFtop		

		else:
			oplstop = LFtop

		# Check the generated ligand topology file for duplicate atom types with forcefield standard atomtypes
		ligcheck = Checkligtop(oplstop, tff)

		nametop = " "
		if not ligcheck == oplstop:
			nametop = "ck" + dlname + ".top"
			os.rename(ligcheck, nametop)
		else:
			nametop = oplstop

		findex = indexoflines(nametop)
		ligand_at, ligand_mt = topolsplit(nametop, dlname, findex)
	
		# Check ligand directory for posre file and back it up
		for posrefile in os.listdir():
			if posrefile[0:5] == "posre":
				posrebk = posrefile + ".bk"
				os.rename(posrefile, posrebk)

		os.chdir('../')

		###########################################################
		# GENERATING PROTEIN-LIG COMPLEX TOPOLOGIES AND STRUCTURE #
		###########################################################
		# Creating needed directories and file identifiers
		print('\n')
		os.mkdir('Complex')
		os.chdir('Complex')

		complex_dir = os.path.join(workhost_dir, 'Complex')
		print(f"Protein-Ligand Complex data directory set to {complex_dir}")

		printNote("Generating Protein-Ligands complex topologies and parameters...")
		time.sleep(2)

		clname = Lname + str(RLcount)
		lp_gro = clname + ".gro"
		lp_mol2 = clname + ".mol2"
		lp_prm = clname + ".prmtop"
		lp_inp = clname + ".inpcrd"
		lp_frc = clname + ".frcmod"
		lp_top = clname + ".top"
		LIG_at = clname + "_at.itp"
		LIG_mt = clname + "_mt.itp"
		LIG_mn = clname + "_mn.itp"

		if TFF > 0:
			ulp_top = "u" + clname + ".top"
			uLIG_at = "u" + clname + "_at.itp"
			uLIG_mt = "u" + clname + "_mt.itp"

		# Gather needed files to generate complex topologies
		listcln = os.listdir(lig_dir)
		for cln in listcln:
			if (Path(cln).suffix == ".itp" or Path(cln).suffix == ".gro" or Path(cln).suffix == ".mol2" or Path(cln).suffix == ".frcmod"):
				shutil.copy(os.path.join(lig_dir, cln), './')

		listRdir = os.listdir(receptor_dir)
		for R in listRdir:
			if R == 'xreceptor.pdb' or R == 'nReceptor.pdb' or R == 'nReceptor.top':
				shutil.copy(os.path.join(receptor_dir, R), './')
			elif Path(R).suffix == ".itp":
				shutil.copy(os.path.join(receptor_dir, R), './')

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

		# Now we shall insert details into receptor topol file
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
		
		# Make dir and copy lig pdb into it for pdbcatogro() to access for alternative complex generation
		os.mkdir('ligpdb')
		os.chdir('ligpdb')

		tlplig = clname + "up.pdb"
		shutil.copy(os.path.join(workhost_dir, 'Ligands', tlplig), './')
		os.chdir('../')

		# Check if tleap source file is present and attempt to generate complex via tleap or use alternative
		complexdir = os.listdir()
		tleapfile = "tleap-1-ligs.in"
		if not tleapfile in complexdir:
			print('\n')
			printWarning("GROMACS compatible tleap generated complex structure and topology - Not available")
			print("To use tleap, check README.md for guide on how to create a matching input file")

			printNote("Trying alternative approach to complex generation...")
			catcomplx_gro, catcomplx_pdb = pdbcatogro()

			if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
				print("Alternative approach also failed to generate complex")
				raise Exception("Process aborted. Make necessary corrections and rerun")
	
			else:
				os.rename(catcomplx_gro, 'catComplex.gro')
				os.rename(catcomplx_pdb, 'catComplex.pdb')

		else:
			complx_gro, complx_top = complexgen(tleapfile)

			if not (complx_gro == 'Complex.gro' and complx_top == 'Complex.top'):
				print('\n')
				printWarning("Generation of GROMACS compatible tleap generated complex structure and topology FAILED")

				printNote("Trying alterenative complex generation...")

				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach also failed to generate complex")
					raise Exception("Process aborted. Make necessary corrections and rerun")
	
				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')

			else:
				os.rename(complx_gro, 'tlpComplex.gro')
				os.rename(complx_top, 'tlpComplex.top')

				print("Generating additional alternative complex structure....")
				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach failed to generate complex to serve as fallback")
	
				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')

		try:
			os.remove("LIGS_mn.itp")
			os.remove("uLIGS_mn.itp")
		except:
			pass

		os.chdir('../')

	    #################################
	    # SOLVATING PROTEIN-LIG COMPLEX #
	    #################################
		print('\n')
		printNote("Generating Solvated Protein-Ligands complex structures...")

		os.mkdir("Solvation")
		os.chdir("Solvation")
		solvation_dir = os.path.join(workhost_dir, 'Solvation')
		print("Solvation directory set to: ", solvation_dir)

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
				print('\n')
				printNote("Your uploaded ligand topologies have been prepared for use")
				printNote("PLEASE NOTE THE DEFAULT ORDER OF USE:")
				print(">>>>> 1. User uploaded ligand topologies will be attempted first")
				print(">>>>> 2. Platform generated ligand topologies will serve as fallback")
				response = tinput("To reverse the order, type YES/y, Otherwise press ENTER: ", defaults[4], "n")
				if not (response.lower() == "yes" or response.lower() == "y"):
					os.rename('LIGS_at.itp', 'bkLIGS_at.itp')
					os.rename('LIGS_mt.itp', 'bkLIGS_mt.itp')
					os.rename('topol.top', 'bktopol.top')
					os.rename('uLIGS_at.itp', 'LIGS_at.itp')
					os.rename('uLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('utopol.top', 'topol.top')
					printNote("THE DEFAULT ORDER WILL BE MAINTAINED")
					print('\n')

				else:
					printNote("THE DEFAULT ORDER HAS BEEN REVERSED")
					os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
					os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
					os.rename('utopol.top', 'bktopol.top')
					print('\n')

			else:
				print('\n')
				print(f"{tff} in topol.top does not match the selected forcefield group, {Ligsff}")
				printNote("Forcefield in topol.top will be used")

		else:
			if tff[0:5].lower() == 'amber' or tff[0:6].lower() == 'charmm' or tff[0:6].lower() == 'gromos' or tff[0:4].lower() == 'opls':
				print(f"{tff} is an acceptable GROMACS compatible forcefield")
				printNote("Default auto-generated ligand topology will be used") 
				if tff[0:6].lower() == 'gromos' or tff[0:4].lower() == 'opls':
					print(f"for {tff}, uploading prepared ligand topology file is recommended")
	
		# Time to prepared solvated complex
		if not defaults[1] == "select":
			selwater = newWater
		else:
			selwater = defaults2("topol.top")
		selbt = defaults[2]
		seld = defaults[3]

		chkSoldir = os.listdir()
		if ('tlpComplex.gro' in chkSoldir and 'tlpComplex.top' in chkSoldir):

			# Prepare a version of amber tleap generated topology for suitable use with Gromacs
			print("Preparing tleap topology file for possible use with Gromacs...")
			tlpgmxtopol = TLtopol('nReceptor.top', 'tlpComplex.top', tff)
			if not tlpgmxtopol == "tlptopol.top":
				printNote("Returning tleap topology as originally generated")
			else:
				shutil.copy(tlpgmxtopol, 'xtlptopol.top')
				if "posre_complex.itp" in os.listdir():
					os.rename("posre_complex.itp", "posre_tlptopol.itp")

			# Now it's time to solvate and add ions
			printNote("Solvation in progress.....")
			grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)
			
			if not grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir:
				printWarning("Solvation with default complex unsuccessful")

				printNote("Trying Solvation with backup complex...")
				grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro' and TFF > 0:
					printWarning("Solvation with backup complex unsuccessful")

					print("Trying solvation with backup topologies...")
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
						printWarning("Solvation failed again with default complex")

						printNote("Trying Solvation with the backup complex ...")
						grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

						if not grosolvated == 'fsolvated.gro':
							printWarning("Solvation failed again with backup complex")
							printNote("Solvation can not continue. You may wish to do solvation manually")

						else:
							printNote("Solvation with backup complex was successful")

					else:
						printNote("Solvation with default complex was successful")

				elif not grosolvated == 'fsolvated.gro' and TFF == 0:
					printNote("Solvation with backup complex unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")

				else:
					printNote("Solvation with backup complex was successful")

			elif not (grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir):
				if TFF > 0:
					print("Trying the backup topology files...")
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
						printNote("Solvation with backup topologies and default complex failed")
						print("Solvation can not continue. You may wish to do solvation manually")

					else:
						printNote("Solvation with default complex and backup topologies was successful")

				else:
					printNote("Solvation with default complex unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")
				
			else:
				printNote("Solvation with default complex was successful")

		elif 'catComplex.gro' in chkSoldir:
			printNote("Solvation with backup complex in progress.....")
			grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

			if not grosolvated == 'fsolvated.gro' and TFF > 0:
				printNote("Solvation with backup complex unsuccessful")

				print("Trying the backup topologies ...")
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
					printNote("Solvation with backup complex and topologies failed")
					print("Solvation can not continue. You may wish to do solvation manually")

				else:
					printNote("Solvation with backup complex and topologies was successful")

			elif not grosolvated == 'fsolvated.gro' and TFF == 0:
				printNote("Solvation with backup complex and topologies was unsuccessful")
				print("Solvation can not continue. You may wish to do solvation manually")

			else:
				printNote("Solvation with backup complex and topologies was successful")

		else:
			printNote("Required gro file for solvation not found / generated. Solvation can not continue")
			printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

		os.chdir('../')

	    ##############################################
	    # GATHERING FILES INTO gmxmds FOLDER FOR MDS #
	    ##############################################
		# Make gmxmds directory and populate it with needed data
		print('\n')
		print("Gathering files needed for MDS run ...")

		os.mkdir("gmxmds")
		os.chdir("gmxmds")

		gmxmds_dir = os.path.join(workhost_dir, 'gmxmds')
		print(f"gmxmds directory set to: {gmxmds_dir}")
		time.sleep(2)

		listsol = os.listdir(solvation_dir)
		for file in listsol:
			if Path(file).suffix == ".itp" or Path(file).suffix == ".top" or Path(file).suffix == ".gro":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')

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
				shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'backup'))
			else:
				shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'not4mds'))

		for bkfile in os.listdir():
			if bkfile[0:2].lower() == "bk":
				shutil.move(bkfile, os.path.join(workhost_dir, 'gmxmds', 'backup'))

		# Performing final processing for tlptopol.top file, if it has not been used instead of topol.top file
		if "topol.top" in os.listdir() and "tlptopol.top" in os.listdir() and not defaults[1] == "none":
			print("Final processing of tlptopol.top file in progress .....")
			time.sleep(5)
			tlpfinal("tlptopol.top", "topol.top")

		# If required, new restraint file can now be generated
		listgmxmds = os.listdir()

		if "fsolvated.gro" in listgmxmds:
			grosolvated = "fsolvated.gro"

			printNote("To generate a new restraint interactively, type YES/y, otherwise press ENTER")
			response = tinput("Response: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				restrfile = udrestraint(grosolvated)
				if restrfile == "posre_udp.itp" or "posre_udp.itp" in os.listdir():
					printNote("You have generated posre_udp.itp. To use it in MDS, define -DPOSRES_UDP in .mdp files")
					printNote("OR if already define -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")
			TFF = M

		elif "tlpSolvated.gro" in os.listdir() and "tlpSolvated.top" in os.listdir():
			print("TLeap generated solvated files were found and have been moved to gmxmds folder")
			print("Please read README.md file for further guides on how to use it")
			TFF = M

		else:
			print("To generate alternative tleap solvated files, follow instruction in README.md file to edit relevant tleap file and rerun the process")
			TFF = M

		# Final Processing of generated files for better formating
		print("Performing Final Checking of Processed MDS input files...")
		if "fsolvated.gro" in os.listdir() or "ufsolvate.gro" in os.listdir():
			for mdsfile in os.listdir():
				if not os.path.isdir(mdsfile):
					gmxmdsFClean(mdsfile)
			os.chdir('../../')
			print('\n')
			npairs += 1

		else:
			mdsfgerrors += 1
			decision, newdefaults = gmxmdsFEChecks(defaults)
			if decision == "abort":
				raise Exception("Process aborted. Please check and make corrections")
			elif decision == "continue":
				printWarning("You have chosen to continue with the process despite errors. Please check")
				os.chdir('../../')
				print('\n')
				npairs += 1
			elif decision == "regenerate":
				printNote("Getting files and defauts values ready for Regeneration...")
				defaults = newdefaults
				os.chdir('../../')
				print('\n')

	print(f"{rls_dir} gmxmds subfolder has been populated and ready for use for simulation")
	print('\n')

	printNote("PLEASE NOTE:")
	print("The files in subfolder 'not4mds' of 'gmxmds' are considered not necessary for MDS. Please check")
	print("The files in subfolder 'backup' of 'gmxmds' can be used to replace relevant ones in gmxmds folder for MDS")

	print('\n')
	printNote("Setup with RLsingle route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS")

	# Set the starting date and time
	print('\n')
	rlsingletime2 = datetime.now()
	Tend = rlsingletime2.strftime("%B %d, %Y %H:%M:%S")
	print(f'Generation of MDS input files ends: {Tend}')
