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

    This code is released under GNU General Public License V3.

          <<<  NO WARRANTY AT ALL!!!  >>>

    It was inspired by:

    - CS50 online training for which this code serves as part of the final project

	- PLEASE Read the README.md file and also follow instructions on the GUI and/or Terminal

	Daniyan, Oluwatoyin Michael, B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) Biochemistry
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
import glob
from colored import fore, back, style

from gmodsScripts.gmodsHelpers import ligtopol, receptopol, topolsplit, indexoflines, complexgen, pdbcatogro, solvation, insertdetails, printWarning, printNote

from gmodsScripts.gmodsTLptopol import TLtopol
from gmodsScripts.gmodsTScheck import Checkligtop
from gmodsScripts.gmodsOPLStop import OPLStop, OPLSacpype

def RLsingle():
	# Set global variable to access user supplied topology file(s)
	TFF = 0

	# Set some environment variables
	scriptDIR = os.path.abspath(os.path.dirname(sys.argv[0]))

	GMX_MDS=Path.cwd()
	print("pyGROMODS working directory set to ", GMX_MDS)

	time.sleep(5)

	PLD = os.path.join(scriptDIR, 'gmodsTSF')
	PLDfiles = os.listdir(PLD)

	time.sleep(5)

	fSamples = os.path.join(scriptDIR, 'sampleFiles')
	fPDB = os.path.join(scriptDIR, 'Uploads')

	# Perform initial checks of the ligand, receptor and/or topology folders
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
				print("Expected an instance of one ligand per subdirectory, but detected", len(Lsdirlist)) 
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

	numff = 0
	numtp = 0
	for fileT in TOPfolders:
		Tsubdir = os.path.join(fPDB, 'Ligsff', fileT)
		if os.path.isdir(Tsubdir):
			Tsdirlist = os.listdir(Tsubdir)
			if not len(Tsdirlist) == 1:
				raise Exception("Sorry, expecting one file or folder but found", len(Tsdirlist))

			for Tfile in Tsdirlist:
				if not (Tfile == "Amber" or Tfile == "Charmm" or Tfile == "Gromos" or Tfile == "Opls"):
					raise Exception(Tfile, "could not be identify with any forcefields in gromacs")

				Tdir = os.path.join(fPDB, 'Ligsff', fileT, Tfile)
				if os.path.isdir(Tdir):
					listTdir = os.listdir(Tdir)
					if not len(listTdir) == 1:
						raise Exception("Sorry, expecting one topology file per subdirectory, but found", len(listTdir))

					for Topf in listTdir:
						if not (Path(Topf).suffix == ".itp" or Path(Topf).suffix == ".top"):
							raise Exception("Sorry, unacceptable topology file format detected. Must be .itp or .top")

						else:
							checkindex = indexoflines(os.path.join(fPDB, 'Ligsff', fileT, Tfile, Topf))
							filecheck = open(os.path.join(fPDB, 'Ligsff', fileT, Tfile, Topf), "r")
							readcheck = filecheck.readlines()
							try:
								if not 'atomtypes' in readcheck[checkindex['atomtypes']].split():
									raise Exception("Sorry, an uploaded topology file lack [ atomtypes ] subheading")
								elif not 'moleculetype' in readcheck[checkindex['moleculetype']].split():
									raise Exception("Sorry, an uploaded topology file lack [ moleculetype ] subheading")
								elif not 'system' in readcheck[checkindex['system']].split():
									raise Exception("Sorry, an uploaded topology file lack [ system ] subheading")
								else:
									print(Topf, "an acceptable uploaded ligand topology file detected")
							except Exception as e:
								print("Checking a user supplied ligand topology file failed with error", e)
								printNote("This file may lack ['atomtypes'], ['moleculetype'] and/or [ system ] subheader")
								print("Check and include the subheading with or without expected accompanied values")
								print("To abort, type YES/y. Otherwise the process will ignore uploaded ligand topology")
								response = input("Response: ")						
								if not (response.lower() == "yes" or response.lower() == "y"):
									print(Topf, "topology file will be ignored")
									shutil.rmtree(Tdir)
									changedF = open(Tdir, "w")
									changedF.close()
									numff += 1
									continue
								else:
									raise Exception("Make necessary corrections and restart")
							filecheck.close()

					numtp += 1
				elif os.path.isfile(Tdir):
					numff += 1
				else:
					raise Exception("Empty forcefields subfolder detected")
		else:
			raise Exception("Sorry, a file detected where subdirectory should be. Please setroute again and rerun")

	if not len(RECfolders) == len(LIGfolders):
		raise Exception("Matching pairs of receptor and ligand could not be detected. Please cross check")

	printNote("You have matching pairs of receptor and ligand. A solvated complex will now be generated for each pair. Are you sure you want to continue?")

	response = input("Type YES/y to continue or press ENTER to abort: ")
	if not(response.lower() == "yes" or response.lower() == "y"):
		raise Exception("Process Abort!!!. Setroute again and upload required files")

	time.sleep(5)

	if not len(LIGfolders) == len(TOPfolders):
		raise Exception("Matching pairs of ligand and forcefield/ligand topology not detected. Please cross check")

	if numff == len(LIGfolders) and numtp == 0:
		TFF = 0
		printNote("No user supplied topology file(s) was detected for all pairs of Receptor:Ligand")

	elif numtp == len(LIGfolders) and numff == 0:
		TFF = 1
		printNote("User supplied topology file(s) was detected for all pairs of Receptor:Ligand")

	else:
		TFF = 3
		printNote("One or more user supplied topology file(s) will be ignored")
		response = input("Do you want to proceed anyway? YES/y or NO/n: ")
		if not(response.lower() == "yes" or response.lower() == "y"):
			raise Exception("Process Abort!!!. Check your files and rerun")
	
	time.sleep(5)

	fMDP = os.path.join(scriptDIR, 'MDP')

	# Get user imput for project name
	while True:
		name = input("Suppy a name for the current project: ")
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
		idalpha1 = random.choice(string.ascii_letters)
		idalpha2 = random.choice(string.ascii_letters)
		ID = idalpha1 + str(idnumber) + idalpha2
		foldername = name + "_" + str(ID)
		if not os.path.isdir(foldername):
			print("Your Unique ID is: ", ID)
			print("Your current work will be stored in ", foldername)
			break
		else:
			continue

	# Create working directory using genrated project name
	os.mkdir(foldername)
	os.chdir(foldername)

	work_dir = os.path.join(GMX_MDS, foldername)
	print("pyGROMODS Current workspace directory set to ", work_dir)

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

	os.chdir('../')

	time.sleep(10)

	RLname = "RLsingle_"
	RLcount = 0
	for RL in RECfolders:
		# Set global variable to access user supplied topology file(s)
		Ligsff = " "
		tff = " "
		ligsff_dir = " "
		M = TFF

		# Create host directory for each pair of receptor and ligand
		RLcount += 1
		rls_dir = RLname + str(RLcount)
		workhost_dir = os.path.join(work_dir, rls_dir)
		print("Current project host directory set to ", workhost_dir)

		os.mkdir(rls_dir)
		os.chdir(rls_dir)

		# Check preferred forcefield selection for the pair under consideration
		ligtop = "ligtop" + str(RLcount)
		ligtopdir = os.path.join(fPDB, 'Ligsff', ligtop)
		listligtopdir = os.listdir(ligtopdir)
		for itemf in listligtopdir:
			if itemf == 'Gromos' or itemf == 'Opls' or itemf == 'Amber' or itemf == 'Charmm':
				Ligsff = itemf
			
		ligsff_dir = os.path.join(fPDB, 'Ligsff', ligtop, Ligsff)
		if TFF == 3:
			if os.path.isdir(ligsff_dir) and len(os.listdir(ligsff_dir)) > 0:
				TFF = 1
				printNote("Detected uploaded ligand topology file for this pair")
			else:
				TFF = 0		
				printNote("No uploaded ligand(s) topology(ies) was detected for this pair")

		if not Ligsff == " ":
			printNote("PLEASE NOTE")
			print("Your preselected forcefield group for this pair is", Ligsff)

		################################################
		# GENERATING RECEPTOR TOPOLOGIES AND STRUCTURE #
		################################################

		os.mkdir('Receptors')
		os.chdir('Receptors')

		receptor_dir = os.path.join(workhost_dir, 'Receptors')
		print("Receptor data directory set to ", receptor_dir)

		printNote("Generating Protein topology and parameter files...")

		# Generate Protein topologies and parameter files
		Rname = "receptor"
		drname = Rname + str(RLcount)

		Rdir = os.path.join(fPDB, 'Receptors', RL)
		Rdirlist = os.listdir(Rdir)
		for rep in Rdirlist:
			shutil.copy(os.path.join(fPDB, 'Receptors', RL, rep), './')

		while True:
			RFtop, RFpdb, RFposre = receptopol(rep, drname)

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
				print("Your forcefiled as contained in receptor topology file is", tff)
				print(tff, "forcefield does not match the preselected forcefield group:", Ligsff)
				print("It is advisable to rerun and choose forcefield that match preselected")
				printNote("PLEASE NOTE - If you choose to continue: ")
				print("A). Your forcefield group will be changed to match", tff)
				print("B). By default, any uploaded ligand topology will be ignored")
				print("C). However, if your uploaded ligand topology(ies) is/are compatible with your current forcefield selection, you may choose to keep it, when provied the options below")

				printNote("To rerun, Type YES/y. Otherwise press ENTER to continue with current selection")
				response = input("Response: ")
				if response.lower() == "yes" or response.lower() == "y":
					continue
				else:
					printNote("You have choosen to continue with current forcefield selection")					
					if tff[0:4].capitalize() == "Opls":
						Ligsff = "Opls"
					elif tff[0:5].capitalize() == "Amber":
						Ligsff = "Amber"
					elif tff[0:6].capitalize() == "Gromos":
						Ligsff = "Gromos"
					elif tff[0:6].capitalize() == "Charmm":
						Ligsff = "Charmm"
					else:
						print(tff, "does not match any known forcefield group. Please rerun")
						printNote("This may happen if you used a self created or modified forcefield. As such standard naming convention for forcefield should be used. E.g. Amber group of forcefields starts with amber, Gromos with gromos, etc. OR it may happen if generation of topol.top fails.")
						printWarning("It is strongly recommended to abort the process, check uploaded file for correctness, and try again. Check README.md file for some troubleshooting tips")
						printNote("To abort, Type YES/y. To continue anyway, press ENTER")
						response = input("Response: ")
						if response.lower() == "yes" or response.lower() == "y":
							raise Exception("Process aborted. Make necessary corrections and Rerun setup")
						else:
							break

					sff_folders = os.listdir(os.path.join(fPDB, 'Ligsff', ligtop))
					for sff in sff_folders:
						if not (sff == 'Amber' or sff == 'Charmm' or sff == 'Gromos' or sff == 'Opls'):
							raise Exception("No folder or file matching forcefield that can be renamed was found")
						else:
							os.rename(os.path.join(fPDB, 'Ligsff', ligtop, sff), os.path.join(fPDB, 'Ligsff', ligtop, Ligsff))

					print("Your selected forcefield group has been changed to ", Ligsff)
					ligsff_dir = os.path.join(fPDB, 'Ligsff', ligtop, Ligsff)

					if TFF == 1:
						print("Subdirectory for uploaded ligand topology is now: ", ligsff_dir)
						printNote("To use with uploaded ligand topology, type YES/y")
						printNote("To use without uploaded ligand topology, type NO/n")
						response = input("Response: ")
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
				response = input("Choose your preferred option: ")
				if response == '1' or response == '2' or response == '3':
					print("Option", response, "Selected!")
					confirm = input("Type YES/y to confirm or press ENTER to choose a different option: ")
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
		printNote("Generating topology and parameter files for ligands...")

		os.mkdir('Ligands')
		os.chdir('Ligands')

		lig_dir = os.path.join(workhost_dir, 'Ligands')
		print("Ligands data directory set to ", lig_dir)

		Lname = "LIG"
		Lfolder = LIGfolders[RLcount - 1]
		dlname = Lname + str(RLcount)

		shutil.copy(os.path.join(PLD, 'ligtoptleap.in'), './')

		Ldir = os.path.join(fPDB, 'Ligands', Lfolder)
		Ldirlist = os.listdir(Ldir)
		for lig in Ldirlist:
			shutil.copy(os.path.join(fPDB, 'Ligands', Lfolder, lig), './')

		nligname = dlname + "up.pdb"
		shutil.copy(lig, nligname)

		LFgro, LFtop = ligtopol(lig, dlname)

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
				print(utop, "user supplied topology does not match the ligand named", ulig)
				print("Please compare the two files and if correct, continue, otherwise, the uploaded topology for the current ligand will be ignored")
				print("If you need the topology, rerun the process and upload correct file")
				response = input("Type YES/y to continue, Otherwise topology will be ignored: ")
				if not (response.lower() == "yes" or response.lower() == "y"):
					TFF = 0
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
				printWarning("Platform opls compatible ligand topology generation failed")
				printNote("Trying acpype...")
				oplstop = OPLSacpype(LFtop, lig, dlname)
				if not oplstop == mergename:
					printWarning("acpype opls compatible ligand topology generation also failed")
					print("Check that you have correctly installed the latest acpype")
					printNote("Platform default will now be used")
					oplstop = LFtop		

		elif opls_route == 2:
			mergename = "opls" + dlname + ".top"
			oplstop = OPLSacpype(LFtop, lig, dlname)

			if not oplstop == mergename:
				printWarning("acpype opls compatible ligand topology generation also failed")
				print("Check that you have correctly installed the latest acpype")
				printNote("Trying Platform opls compatible ligand topology generation ...")
				oplstop = OPLStop(LFtop, tff, dlname)
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
			nametop = "ck" +dlname + ".top"
			os.rename(ligcheck, nametop)
		else:
			nametop = oplstop

		findex = indexoflines(nametop)
		ligand_at, ligand_mt = topolsplit(nametop, dlname, findex)
	
		os.chdir('../')

		time.sleep(10)

		###########################################################
		# GENERATING PROTEIN-LIG COMPLEX TOPOLOGIES AND STRUCTURE #
		###########################################################
		# Creating needed directories and file identifiers
		os.mkdir('Complex')
		os.chdir('Complex')

		complex_dir = os.path.join(workhost_dir, 'Complex')
		print("Protein-Ligand Complex data directory set to ", complex_dir)

		printNote("Generating Protein-Ligands complex topologies and parameters for each pair of receptor and ligand...")

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

		print("Generating Protein-Ligands complex structures...")
		print("Gathering needed files...")

		# Gather needed files to generate complex topologies
		listcln = os.listdir(lig_dir)
		for cln in listcln:
			if (Path(cln).suffix == ".itp" or Path(cln).suffix == ".gro" or Path(cln).suffix == ".mol2" or Path(cln).suffix == ".frcmod"):
				shutil.copy(os.path.join(lig_dir, cln), './')

		listRdir = os.listdir(receptor_dir)
		for R in listRdir:
			if R == 'xreceptor.pdb' or R == 'nReceptor.pdb' or R == 'nReceptor.top' or R == 'posre.itp':
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
			print("Generating complex topologies and parameters for Receptor - Ligand pair ...")
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
			printNote("#####################################################################")
			print("No matching tleap input file is found for your number of ligands")
			print("As such, complex can not be generated using tleap")
			print("To use tleap, check README.md for guide on how to create a matching input file")
			print("Note that it's only through tleap that added useful alternative files can be generated")
			print("It is therefore advisable to abort and attempt to fix the errors")
			printNote("#####################################################################")

			print("However, you may wish to try an alternative approach to complex generation")
			response = input("Type YES/y to try an alternative approach. Otherwise press ENTER to abort: ")
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

		else:
			complx_gro, complx_top = complexgen(tleapfile)

			if not (complx_gro == 'Complex.gro' and complx_top == 'Complex.top'):
				printNote("#####################################################################")
				printWarning("tleap could not successfully generate complex structure and topologies")
				printNote("This is understandable if the choosen forcefield did not match amber. It is advisable to abort and attempt to fix the error")
				printNote("#####################################################################")

				print("However, you may wish to try alternative approach to complex generation")
				response = input("Type YES/y to try an alternative approach. Otherwise press ENTER to abort: ")
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

			else:
				os.rename(complx_gro, 'tlpComplex.gro')
				os.rename(complx_top, 'tlpComplex.top')

				catcomplx_gro, catcomplx_pdb = pdbcatogro()

				if not (catcomplx_gro == 'conComplex.gro' and catcomplx_pdb == 'conComplex.pdb'):
					print("Alternative approach failed to generate complex to serve as fallback")
	
				else:
					os.rename(catcomplx_gro, 'catComplex.gro')
					os.rename(catcomplx_pdb, 'catComplex.pdb')

		os.chdir('../')

		time.sleep(10)

	    #################################
	    # SOLVATING PROTEIN-LIG COMPLEX #
	    #################################
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
				response = input("Response: ")
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
				print(tff, "in topol.top does not match the selected forcefield group,", Ligsff)
				raise Exception("Platform Aborted!!!. Check and restart the process")

		else:
			if tff[0:5].lower() == 'amber':
				print(tff, "is an amber forcefield")
				printNote("Default auto-generated ligand topology will be used") 

			elif tff[0:6].lower() == 'charmm':
				print(tff, "is a charmm forcefield")
				printNote("Default auto-generated ligand topology will be used") 

			elif tff[0:6].lower() == 'gromos':
				print(tff, "is a gromos forcefield and you have not uploaded ligand topology files.")
				printNote("Default auto-generated ligand topology will be used, but may fail") 

			elif tff[0:4].lower() == 'opls':
				print(tff, "is an opls forcefield and you have not uploaded ligand topology files.")
				printNote("Auto-generated OPLS compatible ligand topology will be used, but may fail") 
	
		time.sleep(5)

		# Time to prepared solvated complex
		chkSoldir = os.listdir()
		if ('tlpComplex.gro' in chkSoldir and 'tlpComplex.top' in chkSoldir):

			# Prepare a version of amber tleap generated topology for suitable use with Gromacs
			print("Preparing tleap topology file for possible use with Gromacs...")
			tlpgmxtopol = TLtopol('nReceptor.top', 'tlpComplex.top', tff)
			shutil.copy(tlpgmxtopol, 'xtlptopol.top')

			# Now it's time to solvate and add ions
			print('\n')

			printNote("Solvation with tlpComplex.gro in progress.....")
			grosolvated = solvation('tlpComplex.gro', 'topol.top')
			
			if not grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir:
				printNote("Solvation with tlpComplex.gro unsuccessful")

				printNote("Trying Solvation with the alternative complex in progress.....")
				grosolvated = solvation('catComplex.gro', 'topol.top')

				if not grosolvated == 'fsolvated.gro' and TFF > 0:
					printNote("Solvation with alternative complex unsuccessful")
					print('\n')

					print("Trying the backup topology files with tlpComplex.gro")
					print("Doing backup to avoid using already updated files..")

					os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
					os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
					os.rename('topol.top', '#topol.top.bk#')
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
					os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('bktopol.top', 'topol.top')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
					printNote("Repeating Solvation with tlpComplex.gro in progress.....")
					grosolvated = solvation('tlpComplex.gro', 'topol.top')

					if not grosolvated == 'fsolvated.gro':
						printNote("Solvation with backup topologies failed with tlpComplex.gro")
						print('\n')

						printNote("Trying Solvation with the alternative complex with backup topologies ...")
						grosolvated = solvation('catComplex.gro', 'topol.top')

						if not grosolvated == 'fsolvated.gro':
							printNote("Solvation using backup topologies with alternative complex unsuccessful")
							print("Solvation can not continue. You may wish to do solvation manually")
							print('\n')

						else:
							printNote("Solvation of catComplex.gro with backup topologies was successful")
							print('\n')

					else:
						printNote("Solvation of tlpComplex.gro with backup topologies was successful")

				elif not grosolvated == 'fsolvated.gro' and TFF == 0:
					printNote("Solvation with alternative complex unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")
					print('\n')

				else:
					printNote("Solvation with alternative complex was successful")

			elif not (grosolvated == 'fsolvated.gro' and 'catComplex.gro' in chkSoldir):
				if TFF > 0:
					print("Trying the backup topology files with tlpComplex.gro")
					print("Doing backup of updated files..")

					os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
					os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
					os.rename('topol.top', '#topol.top.bk#')
					os.rename('tlptopol.top', '#tlptopol.top.bk#')
					os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
					os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
					os.rename('bktopol.top', 'topol.top')
					shutil.copy('xtlptopol.top', 'tlptopol.top')
	
					printNote("Repeating Solvation with tlpComplex.gro in progress.....")
					grosolvated = solvation('tlpComplex.gro', 'topol.top')

					if not grosolvated == 'fsolvated.gro':
						printNote("Solvation with backup topologies failed with tlpComplex.gro")
						print("Solvation can not continue. You may wish to do solvation manually")
						print('\n')

					else:
						printNote("Solvation of tlpComplex.gro with backup topologies was successful")
						print('\n')

				else:
					printNote("Solvation with tlpComplex.gro unsuccessful")
					print("Solvation can not continue. You may wish to do solvation manually")
					print('\n')
				
			else:
				printNote("Solvation with tlpComplex.gro was successful")

		elif 'catComplex.gro' in chkSoldir:
			printNote("Solvation with catComplex.gro in progress.....")
			grosolvated = solvation('catComplex.gro', 'topol.top')

			if not grosolvated == 'fsolvated.gro' and TFF > 0:
				printNote("Solvation with catComplex.gro unsuccessful")

				print("Trying the backup topologies with catComplex.gro")
				print("Doing backup of updated files..")

				os.rename('LIGS_at.itp', '#LIGS_at.itp.bk#')
				os.rename('LIGS_mt.itp', '#LIGS_mt.itp.bk#')
				os.rename('topol.top', '#topol.top.bk#')
				os.rename('tlptopol.top', '#tlptopol.top.bk#')
				os.rename('bkLIGS_at.itp', 'LIGS_at.itp')
				os.rename('bkLIGS_mt.itp', 'LIGS_mt.itp')
				os.rename('bktopol.top', 'topol.top')
				shutil.copy('xtlptopol.top', 'tlptopol.top')
	
				printNote("Repeating Solvation with catComplex.gro in progress.....")
				grosolvated = solvation('catComplex.gro', 'topol.top')

				if not grosolvated == 'fsolvated.gro':
					printNote("Solvation with backup topologies failed with catComplex.gro")
					print("Solvation can not continue. You may wish to do solvation manually")
					print('\n')

				else:
					printNote("Solvation of catComplex.gro with backup topologies was successful")
					print('\n')

			elif not grosolvated == 'fsolvated.gro' and TFF == 0:
				printNote("Solvation of catComplex.gro with backup topologies unsuccessful")
				print("Solvation can not continue. You may wish to do solvation manually")
				print('\n')

			else:
				printNote("Solvation of catComplex.gro with backup topologies was successful")
				print('\n')

		else:
			printNote("Required gro file for solvation not found / generated. Solvation can not continue")
			printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

		if 'notipw' in os.listdir():
			try:
				os.rename('fsolvated.gro', 'ufsolvate.gro')
			except Exception as e:
				print("Something went wrong with error", e)

		os.chdir('../')

		time.sleep(10)

	    ##############################################
	    # GATHERING FILES INTO gmxmds FOLDER FOR MDS #
	    ##############################################
		# Make gmxmds directory and populate it with needed data
		print("Gathering files needed for MDS run ...")

		os.mkdir("gmxmds")
		os.chdir("gmxmds")
		gmxmds_dir = os.path.join(workhost_dir, 'gmxmds')
		print("gmxmds directory set to: ", gmxmds_dir)

		listsol = os.listdir(solvation_dir)
		for file in listsol:
			if (Path(file).suffix == ".itp" or Path(file).suffix == ".top"):
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == "fsolvated.gro" or file == "ufsolvate.gro":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif not ('fsolvated.gro' in listsol or 'ufsolvate.gro' in listsol):
				if file == 'tlpSolvated.gro' or file == 'tlpSolvated.top':
					shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
				elif not 'tlpSolvated.gro' in listsol and 'tlpSolvated.top' in listsol:
					if file == 'tlpComplex.gro' or file == 'tlpComplex.top' or file == 'catComplex.gro':
						shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')

		# Removing files that are not needed for MDS with Gromacs
		print("Removing files not needed for MDS from gmxmds directory")
		print("They can be found in Solvation directory")

		removal = ['nReceptor.top', 'tlpComplex.top', 'tlptopol_mn.itp', 'ffnonbonded.itp', 'xtlptopol.top', 'LIGS_mn.itp']
		for rmf in removal:
			try:
				os.remove(rmf)
			except:
				pass

		time.sleep(5)

		# If required, new restraint file can now be generated
		if not grosolvated == "fsolvated.gro":
			print(rls_dir, "gmxmds subfolder has been populated and ready for manual solvation or simulation")
			os.chdir('../../')
		else:
			printNote("##############################################################################")
			print("# The current posre.itp restrain all heavy atoms which include Backbone atoms")
			print("# You can generate your desired restrain file if this does not meet your need")
			print("# This should be named posre_udp. To use posre.itp, define -DPOSRE in .mdp files")
			print("# To use posre_udp.itp instead, define -DPOSRE_UDP in .mdp files")
			printNote("##############################################################################")

			response = input("To generate a new restraint interactively, type YES/y, otherwise the default will be used: ")
			if response.lower() == "yes" or response.lower() == "y":
				success = 0
				try:
					subprocess.run('gmx genrestr -f ' + grosolvated + ' -o posre_udp.itp', shell=True)
				except subprocess.CalledProcessError as e:
					print(e)
					printWarning("Something went wrong with the above error message, Trying again ....")			
					time.sleep(5)
					Err1 = os.system('gmx genrestr -f ' + grosolvated + ' -o posre_udp.itp')
					if not Err1 == 0:
						printWarning("No solvated complex found. Restraint can't be generated successfully")
						success = 1
				
				shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')
				insertUDP = "; Include water topology"
				insertdetails('topol.top', 'mt-topol2.itp', insertUDP)
				os.remove('mt-topol2.itp')

				if success == 0:
					printNote("You have generated posre_udp.itp. To use it in MDS, define -DPOSRES_UDP in .mdp files")
					printNote("OR if already define -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")
				else:
					printNote("No posre_udp.itp has been generated. If need be, generate it manually")
				time.sleep(5)

			else:
				shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')

				insertUDP = "; Include water topology"
				insertdetails('topol.top', 'mt-topol2.itp', insertUDP)
				os.remove('mt-topol2.itp')

				printNote("No posre_udp.itp has been generated. If need be, generate it manually")
				time.sleep(5)

			print(rls_dir, "gmxmds subfolder has been populated and ready for use for simulation")
			TFF = M
			time.sleep(5)
			os.chdir('../../')

	print('\n')
	printNote("Setup with RLmulti route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS") 

	os.chdir('../')
