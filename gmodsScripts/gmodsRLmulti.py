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
from gmodsScripts.gmodsTScheck import checktopsimilarities, Checkligtop
from gmodsScripts.gmodsOPLStop import OPLStop, OPLSacpype

def RLmulti(appDIR, gmxDIR, fdefaults):
	# Set some environment variables
	print('\n')
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
	
	printNote("Let us check again your selected default values ..... ")

	print(f"Default forcefield is {defaults[0]}")
	print(f"Default water model is {defaults[1]}")
	print(f"Default editconf -bt option is {defaults[2]}")
	print(f"Default editconf -d option is {defaults[3]}")
	print(f"Default timeout for input() request is {defaults[4]}")

	if defaults[5] == "A":
		printNote("Your selected default mode for generating input file is Interractive")
		response = tinput("To revert to Noninteractive mode type YES/y: ", defaults[4], "n")
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
			printNote("Number of supplied topology files does not match the number of ligands")
			printNote("If the topology is important, please abort, crosscheck and rerun")
			printWarning("If you choose to continue, user uploaded topology will be ignored")
			response = tinput("To continue, type YES/y, or press ENTER: ", defaults[4], "n")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Incomplete topology files!!!. Check and restart the process")
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
						printNote("This file may lack ['atomtypes'], ['moleculetype'] and/or [ system ] subheading")
						print("Check and include the subheadings, with or without expected accompanied values")
						print("To restart, type YES/y. Otherwise the process will ignore all uploaded ligand topologies")
						response = tinput("Response: ", defaults[4], "y")						
						if not (response.lower() == "yes" or response.lower() == "y"):
							TFF = 0
						else:
							raise Exception("Checking ligand topology failed. Check, correc and restart")
					filecheck.close()
					
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
		response = tinput("Do you want to continue, type YES/y or press ENTER to abort: ", defaults[4], "y")
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
			printWarning("WARNING: You have more than one tleap source files with similar identifier")
			raise Exception("Multiple tleap source file with similar identifiers found. Check samples in", PLD)

		elif tlpresent < 1:
			printWarning("WARNING: You have no tleap source files with identifier that match number of uploaded ligands")
			print(f"Please refer to README.md file for instruction on how to prepare a suitable tleap source file, appropriate for your actual number of uploaded ligands, and saved appropriately in {PLD}")
			print("Alternatively, you may also reduce the number of uploaded ligands, to make use of relevant pre-installed tleap source file")
			raise Exception("No suitable tleap source file found")

		elif tlpresent == 1:
			printNote("Found the appropriate tleap source file for your work")
			break

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
			print(f"Your forcefiled as contained in receptor topology file is {tff}")
			print(f"{tff} forcefield does not match the preselected forcefield group: {Ligsff}")
			print("It is advisable to rerun and choose forcefield that match preselected")
			printNote("PLEASE NOTE - If you choose to continue: ")
			print(f"A). Your forcefield group will be changed to {tff}")
			print("B). By default, any uploaded ligand topology(ies) will be ignored")
			print("C). However, you may choose to keep your uploaded ligand topology(ies) if compatible with your current forcefield selection")

			printNote("To rerun, Type YES/y. To continue with current selection press ENTER")
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
					print(f"{tff} does not match any known forcefield group. Please rerun")
					printNote("This may happen if you used a self created or modified forcefield. As such standard naming convention for forcefield should be used. E.g. Amber group of forcefields starts with amber, Gromos with gromos, etc. OR it may happen if generation of topol.top fails.")
					printWarning("It is strongly recommended to check the uploaded file for correctness, and try again. Check README.md file for some troubleshooting tips")
					raise Exception("Process aborted. Make necessary corrections and Rerun setup")

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
					response = tinput("Response: ", defaults[4], "n")
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
		print("Generating optimized ligand structure...")
		ligE = open('ligcleanerror.txt', 'a')
		nligname = Lname + "up.pdb"
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

			if not (dul_ml > ul or dul_ml > ml):
				print(f"{utop} the user supplied topology match the ligand named {ulig}")
				time.sleep(10)
			else:
				print(f"{utop} user supplied topology does not match the ligand named {ulig}")
				print("This may happen if the topology(ies) is/are not named similar to the corresponding ligand(s)")
				print("Please compare the two files and if you're sure it's correct, continue, otherwise abort")
				response = tinput("Type YES/y to continue, Otherwise press ENTER to abort: ", defaults[4], "n")
				if not (response.lower() == "yes" or response.lower() == "y"):
					raise Exception("Process aborted. Restart and rerun")

			# Check the file for duplicate atom types with forcefield standard atomtypes
			uligcheck = Checkligtop(utop, tff)
			uctindex = indexoflines(uligcheck)

			# Generate the atomtypes and moleculetypes file
			uligand_at, uligand_mt = topolsplit(uligcheck, Lname, uctindex)

			shutil.copy(uligand_at, '../../user_at')
			utopname = Lname + ".top"
			os.rename(uligcheck, utopname)
			shutil.copy(utopname, '../../user_top')
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
				printWarning("acpype opls compatible ligand topology generation also failed")
				print("Check that you have correctly installed the latest acpype")
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
		if R == 'xreceptor.pdb' or R == 'nReceptor.top' or R == 'nReceptor.pdb' or R == 'posre.itp':
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
		printWarning("Defualt complex generation will not continue. Please crosscheck")
		printNote("Meanwhile, attempt will be made to use the alternative approach")
		time.sleep(7)
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
		printNote("#####################################################################")
		print("No matching tLeap input file is found for your number of ligands")
		print("As such, complex can not be generated using tleap")
		print("To use tleap, check README.md for guide on how to create a matching input file")
		print("Note that it's only through tleap that added useful alternative files can be generated")
		print("It is therefore advisable to attempt to fix the errors and rerun the process")
		printNote("#####################################################################")

		print("However, attempt will be made to use alternative approach to complex generation")
		time.sleep(5)
		print("Trying alterenative complex generation...")

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
			printNote("############################################################################")
			printWarning("tleap could not successfully generate complex structure and topologies")
			print("Check the outputed errors and README.md for troubleshooting tips to fix the errors")
			print("You can rerun the process afterwards")
			printNote("############################################################################")

			print("However, attempt will be made to use alternative approach to complex generation")
			time.sleep(5)
			print("Trying alterenative complex generation...")

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

			catcomplx_gro, catcomplx_pdb = pdbcatogro()

			if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
				print("Alternative approach failed to generate complex to serve as fallback")
				print("This may not affect your work")
	
			else:
				os.rename(catcomplx_gro, 'catComplex.gro')
				os.rename(catcomplx_pdb, 'catComplex.pdb')

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
	print("Gathering require files...")

	listLdir = os.listdir(os.path.join(work_dir, 'Ligands'))
	for L in listLdir:
		if Path(L).suffix == ".itp":
			shutil.copy(os.path.join(work_dir, 'Ligands', L), './')

	listRdir = os.listdir(os.path.join(work_dir, 'Receptor'))
	for R in listRdir:
		if R == 'xreceptor.pdb' or R == 'nReceptor.top' or R == 'nReceptor.pdb':
			shutil.copy(os.path.join(work_dir, 'Receptor', R), './')

	listTdir = os.listdir(os.path.join(work_dir, 'Complex'))
	for T in listTdir:
		if T == 'tlpComplex.gro' or T == 'tlpComplex.top' or T == 'catComplex.gro' or T == 'catComplex.pdb' or T == 'tlpSolvated.gro' or T == 'tlpSolvated.top' or T == 'posre.itp' or T == 'topol.top' or T == 'utopol.top':
			shutil.copy(os.path.join(work_dir, 'Complex', T), './')

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
			response = tinput("Response: ", defaults[4], "n")
			if not (response.lower() == "yes" or response.lower() == "y"):
				print("Backing off and renaming relevant files....")
				os.rename('LIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('LIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('LIGS_mn.itp', 'bkLIGS_mn.itp')
				os.rename('topol.top', 'bktopol.top')
				os.rename('uLIGS_at.itp', 'LIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('uLIGS_mn.itp', 'LIGS_mn.itp')
				os.rename('utopol.top', 'topol.top')

			else:
				print("The order has been reversed as follow:")
				print("1. Platform generated ligand topology(ies) will be used first")
				print("2. User uploaded ligand topologies will serve as fallback")
				print("Backing off some relevant files....")
				os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('uLIGS_mn.itp', 'bkLIGS_mn.itp')
				os.rename('utopol.top', 'bktopol.top')

		else:
			print(f"{tff} in topol.top does not match the earlier selected forcefield group, {Ligsff}")
			printNote("Solvation will most likely failed. It is advisable to check and rerun the process")
			print("Backing off some relevant files....")
			try:
				os.rename('uLIGS_at.itp', 'bkLIGS_at.itp')
				os.rename('uLIGS_mt.itp', 'bkLIGS_mt.itp')
				os.rename('uLIGS_mn.itp', 'bkLIGS_mn.itp')
				os.rename('utopol.top', 'bktopol.top')
			except Exception as e:
				print(e)
				pass

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
				print("Doing backup of updated files...")

				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('LIGS_mn.itp', '#LIGS_mn.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('tlptopol.top', '#tlptopol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bkLIGS_mn.itp', 'LIGS_mn.itp')
				os.rename('bktopol.top', 'topol.top')
				shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				printNote("Repeating Solvation with tlpComplex.gro in progress.....")
				grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro':
					printWarning("Solvation with backup topologies failed with tlpComplex.gro")

					print("Trying Solvation with the alternative complex with backup topologies ...")
					grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

					if not grosolvated == 'fsolvated.gro':
						printWarning("Solvation using backup topologies with alternative complex unsuccessful")

					else:
						printNote("Solvation of catComplex.gro with backup topologies was successful")

				else:
					printNote("Solvation of tlpComplex.gro with backup topologies was successful")

			elif not grosolvated == 'fsolvated.gro' and TFF == 0:
				printWarning("Solvation with alternative complex unsuccessful")

			else:
				printNote("Solvation with alternative complex was successful")

		elif not (grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir):
			if TFF > 0:
				print("Trying the backup topology files with tlpComplex.gro")
				print("Doing backup of updated files...")

				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('LIGS_mn.itp', '#LIGS_mn.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('tlptopol.top', '#tlptopol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bkLIGS_mn.itp', 'LIGS_mn.itp')
				os.rename('bktopol.top', 'topol.top')
				shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				printNote("Repeating Solvation with tlpComplex.gro in progress.....")
				grosolvated = solvation('tlpComplex.gro', 'topol.top', selwater, selbt, seld)

				if not grosolvated == 'fsolvated.gro':
					printWarning("Solvation with backup topologies failed with tlpComplex.gro")

				else:
					printNote("Solvation of tlpComplex.gro with backup topologies was successful")

			else:
				printWarning("Solvation with tlpComplex.gro unsuccessful")
				
		else:
			printNote("Solvation with tlpComplex.gro was successful")

	elif 'catComplex.gro' in chkSoldir:
		printNote("Solvation with catComplex.gro in progress.....")
		grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

		if not grosolvated == 'fsolvated.gro' and TFF > 0:
			printWarning("Solvation with catComplex.gro unsuccessful")

			print("Trying the backup topologies with catComplex.gro")
			print("Doing backup of updated files...")

			os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
			os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
			os.rename('LIGS_mn.itp', '#LIGS_mn.itp.bk#')
			os.rename('topol.top', '#topol.top.bk#')
			os.rename('tlptopol.top', '#tlptopol.top.bk#')
			os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
			os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
			os.rename('bkLIGS_mn.itp', 'LIGS_mn.itp')
			os.rename('bktopol.top', 'topol.top')
			shutil.copy('xtlptopol.top', 'tlptopol.top')
	
			printNote("Repeating Solvation with catComplex.gro in progress.....")
			grosolvated = solvation('catComplex.gro', 'topol.top', selwater, selbt, seld)

			if not grosolvated == 'fsolvated.gro':
				printWarning("Solvation with backup topologies failed with catComplex.gro")

			else:
				printNote("Solvation of catComplex.gro with backup topologies was successful")

		elif not grosolvated == 'fsolvated.gro' and TFF == 0:
			printWarning("Solvation with catComplex.gro unsuccessful")

		else:
			printNote("Solvation with catComplex.gro was successful")

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
		if (Path(file).suffix == ".itp" or Path(file).suffix == ".gro" or Path(file).suffix == ".top"):
			shutil.copy(os.path.join(work_dir, 'Solvation', file), './')

	# Removing files that are not needed for MDS with Gromacs
	print("Removing files not needed for MDS from gmxmds directory")
	print("They can be found in Solvation directory")

	gmxmdsrequired = ['fsolvated.gro', 'ufsolvate.gro', 'tlpSolvated.gro', 'tlpSolvated.top', 'tlptopol.top', 'topol.top', 'LIGS_at.itp', 'LIGS_mt.itp', 'posre.itp']
	for rmf in os.listdir():
		if not rmf in gmxmdsrequired: 
			os.remove(rmf)

	if 'fsolvated.gro' in os.listdir():
		try:
			os.remove('tlpSolvated.gro')
			os.remove('tlpSolvated.top')
		except:
			pass

	# Now if need be, let's generate a new restraint file
	if "fsolvated.gro" in os.listdir():
		grosolvated = "fsolvated.gro"

		printNote("##############################################################################")
		print("# The current posre.itp restrain all heavy atoms which include Backbone atoms")
		print("# You can generate your desired restrain file if this does not meet your need")
		print("# This should be named posre_udp. To use posre.itp, define -DPOSRE in .mdp files")
		print("# To use posre_udp.itp instead, define -DPOSRE_UDP in .mdp files")
		printNote("##############################################################################")

		printNote("To generate a new restraint interactively, type YES/y, otherwise press ENTER")
		response = tinput("Response: ", defaults[4], "n")
		if response.lower() == "yes" or response.lower() == "y":
			success = 0
			try:
				subprocess.run('gmx genrestr -f ' + grosolvated + ' -o posre_udp.itp', shell=True, stderr=subprocess.STDOUT, check=True, text=True)
			except subprocess.SubprocessError as e:
				print(e)
				printWarning("Something went wrong with the above error message. Please check")			
				time.sleep(5)
				success = 1
			
			if success == 0:
				printNote("You have generated posre_udp.itp. To use it in MDS, define -DPOSRES_UDP in .mdp files")
				printNote("OR if already define -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")
			else:
				printNote("No posre_udp.itp has been generated. If need be, generate it manually")
			time.sleep(5)

		else:
			printNote("No posre_udp.itp has been generated. If need be, generate it manually")
			time.sleep(5)

		shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')
		insertUDP = "; Include water topology"
		insertdetails('topol.top', 'mt-topol2.itp', insertUDP)
		os.remove('mt-topol2.itp')

		printNote("RLmulti gmxmds subfolder has been populated and ready for use for simulation")
		os.chdir('../')

	elif "tlpSolvated.gro" in os.listdir() and "tlpSolvated.top" in os.listdir():
		print("TLeap generated solvated files were found and have been moved to gmxmds folder. Please read README.md file for further gudiance on how to use it")
		printNote("RLmulti gmxmds subfolder has been populated and ready for use for simulation")
		time.sleep(5)
		os.chdir('../')

	else:
		printNote("RLmulti gmxmds subfolder has been populated and ready for manual solvation")
		print("To generate alternative tleap solvated files, follow instruction in README.md file to edit relevant tleap file and rerun the process")
		time.sleep(5)
		os.chdir('../')

	print('\n')
	printNote("Setup with RLmulti route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS") 
