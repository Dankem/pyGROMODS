#!/usr/bin/env python

"""
    pyGROMODS-v2023.05.1 Release

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

from gmodsScripts.gmodsHelpers import ligtopol, receptopol, topolsplit, indexoflines, complexgen, pdbcatogro, solvation, insertdetails, printWarning, printNote, tinput, defaults1, defaults2

from gmodsScripts.gmodsTLptopol import TLtopol
from gmodsScripts.gmodsTScheck import Checkligtop
from gmodsScripts.gmodsOPLStop import OPLStop, OPLSacpype

def RLmany(appDIR, gmxDIR, fdefaults):
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
	print('\n')

	if defaults[5] == "A":
		printNote("Your selected default mode for generating input file is Interractive")
		response = tinput("To revert to Noninteractive mode type YES/y: ", 30, "n")
		if response.lower() == "yes" or response.lower() == "y":
			defaults[5] = "B"
			printNote("You have changed to pdb2gmx non-interactive mode")
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
			print("Number of supplied topology files does not match the number of ligands")
			print("It is highly recommended to abort, crosscheck the files and rerun")
			print("If you choose to continue, user uploaded ligand topologies will be ignored")
			response = tinput("To continue, type YES/y, Otherwise press ENTER to abort: ", defaults[4], "n")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Number of topology files does not match the number of ligands. Check and restart")
			else:
				TFF = 0

		else:
			for file in os.listdir(ligsff_dir):
				if not (Path(file).suffix == ".itp" or Path(file).suffix == ".top"):
					print("Topology file format must be .itp or .top")
					raise Exception("Sorry, unacceptable topology file format detected")
				else:
					checkindex = indexoflines(os.path.join(fPDB, Ligsff, file))
					filecheck = open(os.path.join(fPDB, Ligsff, file), "r")
					readcheck = filecheck.readlines()
					try:
						if not 'atomtypes' in readcheck[checkindex['atomtypes']].split():
							raise Exception("Sorry, one or more uploaded topology file(s) lack [ atomtypes ] subheading")
						elif not 'moleculetype' in readcheck[checkindex['moleculetype']].split():
							raise Exception("Sorry, one or more uploaded topology file(s) lack [ moleculetype ] subheading")
						elif not 'system' in readcheck[checkindex['system']].split():
							raise Exception("Sorry, one or more uploaded topology file(s) lack [ system ] subheading")
						else:
							print(f"{file} is an acceptable uploaded ligand topology file")
					except Exception as e:
						print(f"Checking a user supplied ligand topology file failed with error {e}")
						printNote("This file may lack ['atomtypes'], ['moleculetype'] and/or [ system ] subheaders")
						print("Check and include the subheadings, with or without expected accompanied values")
						print("To abort and restart, type YES/y. Otherwise the process will ignore uploaded ligand topology")
						response = tinput("Response: ", defaults[4], "y")						
						if not (response.lower() == "yes" or response.lower() == "y"):
							TFF = 0
						else:
							raise Exception("Checking supplied ligand topology failed. Check the file and restart")
					filecheck.close()

	# Check the suitability of the uploaded files
	for file in LIGfiles:
		if not (Path(file).suffix == ".pdb" or Path(file).suffix == ".mol2"):
			raise Exception("Sorry, only pdb or mol2 file formats is acceptable. Setroute again and upload required files")

	for file in RECfiles:
		if not Path(file).suffix == ".pdb":
			raise Exception("Sorry, Receptor file must be pdb files. Setroute again and upload required files")

	if not (len(RECfiles) and len(LIGfiles)) > 0:
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
			printWarning("WARNING: You have more than one tleap source files with similar identifier. This may be a case of duplications. Please check and make necessary corrections, and try again")
			raise Exception("Process aborted. Crosscheck the tleap source files found in", PLD)

		elif tlpresent < 1:
			printWarning("WARNING: The needed tleap source files is missing or may have been renamed with missing identifier in the file name")
			raise Exception("Process aborted. Crosscheck the tleap source file")

		elif tlpresent == 1:
			printNote("Found the needed tleap source file for your work")
			break
	print('\n')

	# Get user imput for project name
	while True:
		name = tinput("Suppy a name for the current project: ", defaults[4], "RLmany")
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
	time.sleep(10)

	RLname = "RLmany_"
	RLcount = 0
	for lig in LIGfiles:
		# Create host directory for each pair of receptor and ligand
		M = TFF
		RLcount += 1
		rls_dir = RLname + str(RLcount)
		os.mkdir(rls_dir)
		os.chdir(rls_dir)

		workhost_dir = os.path.join(work_dir, rls_dir)
		print(f"Current project host directory set to {workhost_dir}")
		print('\n')

		# Adjust defaults values of -bt, -b and timeout if ff and water are 'select'
		if defaults[0] == "select" or defaults[1] == "select":
			print(f"Current default values for -bt is: {defaults[2]}")
			print(f"Current default values for -d is: {defaults[3]}")
			print(f"Current default values for timeout is: {defaults[4]}")
			response = tinput("To adjust these values for current protein - ligand complex, type YES/y: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				defaults[2], defaults[3], defaults[4] = defaults1() 

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
				print(f"Your forcefiled as contained in receptor topology file is {tff}")
				print(f"{tff} forcefield does not match the preselected forcefield group: {Ligsff}")
				print("It is advisable to rerun and choose forcefield that match preselected")
				printNote("PLEASE NOTE - If you choose to continue: ")
				print(f"A). Your forcefield group will be changed to match {tff}")
				print("B). The default is that any uploaded ligand topology(ies) will be ignored")
				print("C). However, you may choose to keep your uploaded ligand topology(ies) if compatible with your current forcefield selection")

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
						print("This may happen if you used a self created or modified forcefield. As such, standard naming convention for forcefield should be used. E.g. Amber group of forcefields starts with amber, Gromos with gromos, etc. OR it may happen if generation of topol.top fails.")
						printNote("It is strongly recommended to check the uploaded file for correctness, and try again. Check README.md file for some troubleshooting tips")
						printNote("To abort, Type YES/y. To continue anyway, press ENTER")
						response = tinput("Response: ", defaults[4], "n")
						if response.lower() == "yes" or response.lower() == "y":
							raise Exception("Process aborted. Make necessary corrections and Rerun setup")
						else:
							break

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

		try:
			shutil.copy(rep, 'xreceptor.pdb')
			shutil.copy(RFtop, 'nReceptor.top')
			shutil.copy(RFpdb, 'nReceptor.pdb')
		except Exception as e:
			printWarning(e)
			pass

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
			print(f"		Default forcefield: {defaults[0]}")
			print(f"		Default water model: {defaults[1]}")
			print(f"		Default editconf -bt: {defaults[2]}")
			print(f"		Default editconf -d: {defaults[3]}")
			print(f"		Default input timeout: {defaults[4]}")
			print("		Default mode: non-interactive")

			# We will now lock these defaults by changing mode to C
			defaults[5] = "C"
		else:
			defaults[1] = defaults2(RFtop)
			if defaults[1] == "none":
				print("No water model was detected for your system")
			else:
				print(f"Your water model as contained in topol.top file is {defaults[1]}")

		os.chdir('../')
		time.sleep(5)

		# Determine and choose preferred route for platform generated opls ligand topology
		opls_route = 0
		if Ligsff == 'Opls':
			printNote("You have selected opls as your preferred forcefield")
			printNote("Attempt will be made to generate opls compatible topology for each ligand")
			printNote("The following options are available to generate ligand topology:")
			print("	1. Using platform opls compatible ligand topology generation - RECOMMENDED")
			print("	2. Using acpype opls compatible ligand topology generation")
			print("	3. Using platform default - ignore opls compatibility - NOT RECOMMENDED")
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
				print("Please compare the two files and if correct, continue, otherwise, it advisable to ignore the affected uploaded ligand topology file")
				print("If you need the topology, rerun the process and upload correct file")
				printNote("Type YES/y to continue, Otherwise topology will be ignored")
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
				printWarning("Platform opls compatible ligand topology generation failed")
				printNote("Trying acpype...")
				oplstop = OPLSacpype(LFtop, lig, name)
				if not oplstop == mergename:
					printWarning("acpype opls compatible ligand topology generation also failed")
					print("Check that you have correctly installed the latest acpype")
					printNote("Platform default will now be used")
					oplstop = LFtop		

		elif opls_route == 2:
			mergename = "opls" + name + ".top"
			oplstop = OPLSacpype(LFtop, lig, name)

			if not oplstop == mergename:
				printWarning("acpype opls compatible ligand topology generation also failed")
				print("Check that you have correctly installed the latest acpype")
				printNote("Trying Platform opls compatible ligand topology generation ...")
				oplstop = OPLStop(LFtop, tff, name)
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
			nametop = "ck" + name + ".top"
			os.rename(ligcheck, nametop)
		else:
			nametop = oplstop

		findex = indexoflines(nametop)
		ligand_at, ligand_mt = topolsplit(nametop, name, findex)
	
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
			if R == 'xreceptor.pdb' or R == 'nReceptor.pdb' or R == 'nReceptor.top' or R == 'posre.itp':
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
			printNote("#####################################################################")
			print("No matching tLeap input file is found for your number of ligands")
			print("As such, complex can not be generated using tleap")
			print("To use tleap, check README.md for guide on how to create a matching input file")
			print("Note that it's only through tleap that added useful alternative files can be generated")
			print("It is therefore advisable to abort and attempt to fix the errors")
			printNote("#####################################################################")

			print("However, you may wish to try an alternative approach to complex generation")
			printNote("Type YES/y to try an alternative approach. Otherwise press ENTER to abort")
			response = tinput("Response YES/y: ", defaults[4], "y")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Process aborted. Make necessary corrections and rerun")

			print("Trying alterenative complex generation...")

			catcomplx_gro, catcomplx_pdb = pdbcatogro()

			if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
				print("Alternative approach also failed to generate complex")
				raise Exception("Process aborted. Make necessary corrections and rerun")
	
			else:
				os.rename(catcomplx_gro, 'catComplex.gro')
				os.rename(catcomplx_pdb, 'catComplex.pdb')
				print("Complex generation with alternative approach was successful")

		else:
			complx_gro, complx_top = complexgen(tleapfile)

			if not (complx_gro == 'Complex.gro' and complx_top == 'Complex.top'):
				printNote("#####################################################################")
				printWarning("tleap could not successfully generate complex structure and topologies")
				print("This is understandable if the choosen forcefield did not match amber forcefields")
				print("If you use amber or charmm forcefields, it is advisable to abort and attempt to fix the error")
				printNote("#####################################################################")

				print("However, you may wish to try alternative approach to complex generation")
				printNote("Type YES/y to try an alternative approach. Otherwise press ENTER to abort")
				response = tinput("Response YES/y or press ENTER: ", defaults[4], "y")
				if not (response.lower() == "yes" or response.lower() == "y"):
					raise Exception("Process aborted. Make necessary corrections and rerun")

				print("Trying alterenative complex generation...")

				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach also failed to generate complex")
					raise Exception("Process aborted. Make necessary corrections and rerun")
	
				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')
					print("Complex generation with alternative approach was successful")

			else:
				os.rename(complx_gro, 'tlpComplex.gro')
				os.rename(complx_top, 'tlpComplex.top')
				print("Complex generation with tleap was successful")

				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach failed to generate complex to serve as fallback")
	
				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')
					print("Complex generation with alternative approach was successful")

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
			os.remove('at-topol.itp')
			os.remove('mt-topol.itp')
			os.remove('Receptor_new.pdb')
		except:
			pass

		# Determine if the forcefield belong to the selected group - amber, charmm, gromos and opls
		if TFF > 0:
			if tff[0:4].lower() == Ligsff.lower() or tff[0:5].lower() == Ligsff.lower() or tff[0:6].lower() == Ligsff.lower():
				printNote("The forecfield in topol.top match the preferred forcefield group selected at setup")
				printNote("Your uploaded ligand topologies have been prepared for use")
				print('\n')
				printNote("PLEASE NOTE THAT THE DEFAULT ORDER ARE:")
				print("1. User uploaded ligand topologies will be attempted first")
				print("2. Platform generated ligand topologies will serve as fallback")
				print("It is strongly recommended to maintain this order for Gromos and Opls forcefields")
				print('\n')
				printNote("To reverse the order, type YES/y, Otherwise press ENTER to continue")
				response = tinput("Response YES/y or press ENTER: ", defaults[4], "n")
				if not (response.lower() == "yes" or response.lower() == "y"):
					print("Backing off and renaming relevant files....")
					os.rename('LIGS_at.itp', 'bkLIGS_at.itp')
					os.rename('LIGS_mt.itp', 'bkLIGS_mt.itp')
					os.rename('topol.top', 'bktopol.top')
					os.rename('uLIGS_at.itp', 'LIGS_at.itp')
					os.rename('uLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('utopol.top', 'topol.top')

				else:
					print("The order has been reversed as follow:")
					print("1. Platform generated ligand topology(ies) will be used first")
					print("2. User uploaded ligand topologies will serve as fallback")
					print("Backing off some relevant files....")
					os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
					os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
					os.rename('utopol.top', 'bktopol.top')

			else:
				print(f"{tff} in topol.top does not match the selected forcefield group, {Ligsff}")
				printNote("Solvation will most likely failed. It is advisable to check and rerun the process")
				print("Backing off some relevant files....")
				os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('utopol.top', 'bktopol.top')

		else:
			if tff[0:5].lower() == 'amber':
				print(f"{tff} is an amber forcefield")
				printNote("Default auto-generated ligand topology will be used") 

			elif tff[0:6].lower() == 'charmm':
				print(f"{tff} is a charmm forcefield")
				printNote("Default auto-generated ligand topology will be used") 

			elif tff[0:6].lower() == 'gromos':
				print(f"{tff} is a gromos forcefield and you have not uploaded ligand topology files.")
				printNote("Default auto-generated ligand topology will be used, but may fail") 

			elif tff[0:4].lower() == 'opls':
				print(f"{tff} is an opls forcefield and you have not uploaded ligand topology files.")
				printNote("Auto-generated OPLS compatible ligand topology will be used, but may fail") 
	
		# Time to prepared solvated complex
		selwater = defaults[1]
		selbt = defaults[2]
		seld = defaults[3]

		chkSoldir = os.listdir()
		if ('tlpComplex.gro' in chkSoldir and 'tlpComplex.top' in chkSoldir):

			# Prepare a version of amber tleap generated topology for suitable use with Gromacs
			print("Preparing tleap topology file, for possible use with Gromacs...")
			tlpgmxtopol = TLtopol('nReceptor.top', 'tlpComplex.top', tff)
			shutil.copy(tlpgmxtopol, 'xtlptopol.top')

			# Now it's time to solvate and add ions

			printNote("Solvation with tlpComplex.gro in progress.....")
			grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)
			
			if not grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir:
				printNote("Solvation with tlpComplex.gro unsuccessful")

				printNote("Trying Solvation with the alternative complex in progress.....")
				grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro' and TFF > 0:
					printNote("Solvation with alternative complex unsuccessful")

					print("Trying the backup topology files with tlpComplex.gro")
					print("Backing up updated files...")

					os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
					os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
					os.rename('topol.top', '#topol.top.bk#')
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
					os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('bktopol.top', 'topol.top')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
					printNote("Repeating Solvation with tlpComplex.gro in progress.....")
					grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

					if not grosolvated == 'fsolvated.gro':
						printNote("Solvation with backup topologies failed with tlpComplex.gro")

						printNote("Trying Solvation with the alternative complex with backup topologies ...")
						grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

						if not grosolvated == 'fsolvated.gro':
							printNote("Solvation using backup topologies with alternative complex unsuccessful")
							print("Solvation can not continue. You may wish to do solvation manually")

						else:
							printNote("Solvation of catComplex.gro with backup topologies was successful")

					else:
						printNote("Solvation of tlpComplex.gro with backup topologies was successful")

				elif not grosolvated == 'fsolvated.gro' and TFF == 0:
					printNote("Solvation with alternative complex unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")

				else:
					printNote("Solvation with alternative complex was successful")

			elif not (grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir):
				if TFF > 0:
					print("Backing up relevant files...")

					os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
					os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
					os.rename('topol.top', '#topol.top.bk#')
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
					os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('bktopol.top', 'topol.top')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
					printNote("Repeating Solvation with tlpComplex.gro in progress.....")
					grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

					if not grosolvated == 'fsolvated.gro':
						printNote("Solvation with backup topologies failed with tlpComplex.gro")
						print("Solvation can not continue. You may wish to do solvation manually")

					else:
						printNote("Solvation of tlpComplex.gro with backup topologies was successful")

				else:
					printNote("Solvation with tlpComplex.gro unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")
				
			else:
				printNote("Solvation with tlpComplex.gro was successful")

		elif 'catComplex.gro' in chkSoldir:
			printNote("Solvation with catComplex.gro in progress.....")
			grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

			if not grosolvated == 'fsolvated.gro' and TFF > 0:
				printNote("Solvation with catComplex.gro unsuccessful")

				print("Trying the backup topologies with catComplex.gro")
				print("Backing up updated files..")

				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('tlptopol.top', '#tlptopol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bktopol.top', 'topol.top')
				shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				printNote("Repeating Solvation with catComplex.gro in progress.....")
				grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro':
					printNote("Solvation with backup topologies failed with catComplex.gro")
					print("Solvation can not continue. You may wish to do solvation manually")

				else:
					printNote("Solvation of catComplex.gro with backup topologies was successful")

			elif not grosolvated == 'fsolvated.gro' and TFF == 0:
				printNote("Solvation of catComplex.gro with backup topologies unsuccessful")
				print("Solvation can not continue. You may wish to do solvation manually")

			else:
				printNote("Solvation of catComplex.gro with backup topologies was successful")

		else:
			printNote("Required gro file for solvation not found / generated. Solvation can not continue")
			printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

		os.chdir('../')
		print('\n')

		######################################
		# GATHERING FILES FOR MD SIMULATIONS #
		######################################
		# Create gmxmds directory and populate it with needed data
		print("Gathering files needed for MDS run for Receptor and ", cname, "Solvated complex...")
		time.sleep(2)

		gmxname = "gmxmds"
		os.mkdir(gmxname)
		os.chdir(gmxname)

		listsoldir = os.listdir(os.path.join(workhost_dir, 'Solvation'))
		for file in listsoldir:
			if Path(file).suffix == ".itp":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == "topol.top" or file == "tlptopol.top" or file == "fsolvated.gro" or file == "ufsolvate.gro":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif not ("fsolvated.gro" in listsoldir and defaults[1] == "none"):
				if file == "tlpSolvated.gro" or file == "tlpSolvated.top":
					shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')

		# Removing files not needed
		print("Removing files not needed for MDS from gmxmds directory")
		print("These files can still be found in Solvation directory")

		gmxmdsrequired = ['fsolvated.gro', 'ufsolvate.gro', 'tlpSolvated.gro', 'tlpSolvated.top', 'tlptopol.top', 'topol.top', 'LIGS_at.itp', 'LIGS_mt.itp', 'posre.itp']
		for rmf in os.listdir():
			if not rmf in gmxmdsrequired: 
				os.remove(rmf)

		# If required, new restraint file can now be generated
		if "fsolvated.gro" in os.listdir():
			grosolvated = "fsolvated.gro"

			printNote("##############################################################################")
			print("# The current posre.itp restrain all heavy atoms which include Backbone atoms")
			print("# You can generate your desired restrain file if this does not meet your need")
			print("# This should be named posre_udp. To use posre.itp, define -DPOSRE in .mdp files")
			print("# To use posre_udp.itp instead, define -DPOSRE_UDP in .mdp files")
			printNote("##############################################################################")

			time.sleep(5)

			# Generating new restraint file if needed
			printNote("To generate a new restraint interactively, type YES/y, otherwise press ENTER")
			response = tinput("Response YES/y or press ENTER: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				success = 0
				try:
					subprocess.run('gmx genrestr -f ' + grosolvated + ' -o posre_udp.itp', shell=True)
				except subprocess.CalledProcessError as e:
					print(e)
					printWarning("Something went wrong with the above error message. Restrint was not successful")			
					success = 1
				
				if success == 0:
					printNote("You have generated posre_udp.itp. To use it in MDS, define -DPOSRES_UDP in .mdp files")
					printNote("OR if already define -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")
				else:
					printNote("No posre_udp.itp has been generated. If need be, generate it manually")

			else:
				printNote("No posre_udp.itp has been generated. If need be, generate it manually")

			shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')
			insertUDP = "; Include water topology"
			insertdetails('topol.top', 'mt-topol2.itp', insertUDP)
			os.remove('mt-topol2.itp')

			print(f"{rls_dir} gmxmds subfolder has been populated and ready for molecular dynamic simulation")
			time.sleep(5)
			os.chdir('../../')

		elif "tlpSolvated.gro" in os.listdir() and "tlpSolvated.top" in os.listdir():
			print("TLeap generated solvated files have been moved to gmxmds folder. Please read README.md file for further gudiance on how to use it")
			print(f"{rls_dir} gmxmds subfolder has been populated and ready for manual solvation")
			time.sleep(5)
			os.chdir('../../')
    
		else:
			print(f"{rls_dir} gmxmds subfolder has been populated and ready for manual solvation")
			print("To generate alternative tleap solvated files, follow instruction in README.md file to edit relevant tleap file and rerun the process")
			time.sleep(5)
			os.chdir('../../')

		print('\n')
		TFF = M

	printNote("Setup with RLmany route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS") 
