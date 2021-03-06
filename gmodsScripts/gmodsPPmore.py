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

from gmodsScripts.gmodsHelpers import receptopol, topolsplit, indexoflines, complexgen, pdbcatogro, solvation, insertdetails, printWarning, printNote

from gmodsScripts.gmodsTLptopol import TLtopol

def PPmore():
	# Set global variable to access user supplied topology file(s)
	global Ligsff
	global tff
	global fPDB
	global PLD

	Ligsff = " "
	tff = " "

	# Set some environment variables
	scriptDIR = os.path.abspath(os.path.dirname(sys.argv[0]))

	GMX_MDS=Path.cwd()
	print("User working directory set to ", GMX_MDS)

	time.sleep(5)

	PLD = os.path.join(scriptDIR, 'gmodsTSF')
	PLDfiles = os.listdir(PLD)

	time.sleep(10)

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
				raise Exception("Sorry, Peptide sequnece must be saved as .in files. Setroute again and upload files")

	elif 'Proteins' in listfPDB:
		PROfiles = os.listdir(os.path.join(fPDB, 'Proteins'))
		if not len(PROfiles) > 0:
			raise Exception("Needed protein(s) file(s) is/are missing. Please upload required files")
		for file in PROfiles:
			if not Path(file).suffix == ".pdb":
				raise Exception("Sorry, Protein file must be pdb files. Setroute again and upload required files")

	# Check preferred forcefield selection, forcefield folder and its content
	for itemf in listfPDB:
		if itemf == 'Gromos' or itemf == 'Opls' or itemf == 'Amber' or itemf == 'Charmm':
			Ligsff = itemf
			
	if not Ligsff == " ":
		printNote("PLEASE NOTE")
		print("Your preselected forcefield group is", Ligsff)

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

	def pptopsolgen(pname, rname, tleapfile):
		# Get needed variables as global
		global Ligsff
		global tff
		global PLD
		global fPDB

		# Set workhost directory similar to the main function
		workhost_dir = Path.cwd()

		#######################################################
		# GENERATING GROMACS AND ALTERNATIVE AMBER TOPOLOGIES #
		#######################################################

		os.mkdir("ppTOPS")
		os.chdir("ppTOPS")

		pptops_dir = os.path.join(workhost_dir, "ppTOPS")
		print("Protein topology data directory set to ", pptops_dir)

		printNote("Generating topology and parameter files...")

		shutil.copy(os.path.join(workhost_dir, pname), './')

		# Generating Protein topology and parameter files

		while True:
			RFtop, RFpdb, RFposre = receptopol(pname, rname)

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
				print("B). By default, any actions related to former forcefield will be ignored")

				printNote("To rerun, Type YES/y. Otherwise press ENTER to continue with current selection")
				response = input("Response: ")
				if response.lower() == "yes" or response.lower() == "y":
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
						print(tff, "does not match any known forcefield group. Please rerun")
						printNote("This may happen if you used a self created or modified forcefield. As such standard naming convention for forcefield should be used. E.g. Amber group of forcefields starts with amber, Gromos with gromos, etc. OR it may happen if generation of topol.top fails.")
						printWarning("It is strongly recommended to abort the process, check uploaded file for correctness, and try again. Check README.md file for some troubleshooting tips")
						printNote("To abort, Type YES/y. To continue anyway, press ENTER")
						response = input("Response: ")
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

					print("Your preferred forcefield group has been changed to", Ligsff)
					ligsff_dir = os.path.join(fPDB, Ligsff)
					break
			else:
				print("Your forcefiled as contained in topol.top file is", tff)
				break

		rftopfile = "n" + RFtop
		rfpdbfile = "n" + RFpdb
		try:
			os.rename(RFtop, rftopfile)
			os.rename(RFpdb, rfpdbfile)
		except Exception as e:
			printWarning(e)

		# Now, let us genreated alternative amber topology file for use with gromacs

		print("Generating alternative topology file...")

		time.sleep(5)

		# Gather needed files to generate complex topologies

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
			printNote("#####################################################################")
			print("No matching tLeap input file is found for your peptide/protein setup can be found")
			print("As such, alternative amber genrated topology file can not be generated")
			print("It is therefore advisable to abort and attempt to fix the errors")
			printNote("#####################################################################")

			response = input("Type YES/y to continue. Otherwise press ENTER to abort: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Process aborted. Make necessary corrections and rerun")

		else:
			ppmore_gro, ppmore_top = complexgen(tleapfile)

			if not (ppmore_gro == 'Complex.gro' and ppmore_top == 'Complex.top'):
				printWarning("tleap could not successfully generate alternative topologies")
				print("Please check, and if need be abort and attempt to fix the error")

				response = input("Type YES/y to continue. Otherwise press ENTER to abort: ")
				if not (response.lower() == "yes" or response.lower() == "y"):
					raise Exception("Process aborted. Make necessary corrections and rerun")

			else:
				try:
					os.rename(ppmore_gro, tlpgro)
					os.rename(ppmore_top, tlptop)
				except Exception as e:
					print("One or more files could not be found. Renaming failed with error", e)

		for file in ppatopsdir:
			if file == 'tlpSolvated.gro':
				os.rename('tlpSolvated.gro', tlpsolvgro)
			elif file == 'tlpSolvated.top':
				os.rename('tlpSolvated.top', tlpsolvtop)

		os.chdir('../')

		time.sleep(10)

		####################################
		# SOLVATION OF PROTEIN STRUCTURE   #
		####################################
		# Creating solvation directory and populating it with required files
		print("Generating Solvated Structure ...")

		solname = "Solvation"
		os.mkdir(solname)
		os.chdir(solname)

		listTpdir = os.listdir(pptops_dir)
		for Tp in listTpdir:
			if Path(Tp).suffix == ".itp" or Path(Tp).suffix == ".top" or Path(Tp).suffix == ".gro" or Path(Tp).suffix == ".pdb":
				shutil.copy(os.path.join(workhost_dir, "ppTOPS", Tp), './')

		# Determine if the forcefield belong to the selected group - amber, charmm, gromos and opls
		if tff[0:4].lower() == Ligsff.lower() or tff[0:5].lower() == Ligsff.lower() or tff[0:6].lower() == Ligsff.lower():
			printNote("The forecfield in topol.top match the preferred forcefield group")
			time.sleep(5)
			print('\n')
		else:
			printNote("The forecfield in topol.top did not match the preferred forcefield group selected at setup")
			printNote("You may have changed the forecfield at a point. Please crosscheck and if need be rerun")
			time.sleep(5)
			print('\n')
			
		if tff[0:5].lower() == 'amber':
			print(tff, "is an amber forcefield")

		elif tff[0:6].lower() == 'charmm':
			print(tff, "is a charmm forcefield")

		elif tff[0:6].lower() == 'gromos':
			print(tff, "is a gromos forcefield")
			
		elif tff[0:4].lower() == 'opls':
			print(tff, "is an opls forcefield")

		time.sleep(5)

		# Prepare a version of amber tleap generated topology for suitable use with Gromacs
		chkSoldir = os.listdir()
		if (tlpgro in chkSoldir and tlptop in chkSoldir):
			tlpgmxtopol = TLtopol(rftopfile, tlptop, tff)
			shutil.copy(tlpgmxtopol, 'xtlptopol.top')

		# Time to prepared solvated structure
		chkSoldir = os.listdir()
		if (rftopfile in chkSoldir and rfpdbfile in chkSoldir):
			print(rname, "Solvation in progress.....")

			shutil.copy(rftopfile, 'topol.top')
			grosolvated = solvation(rfpdbfile, 'topol.top')
			
			if not grosolvated == 'fsolvated.gro':
				print(rname, "Solvation unsuccessful")

			else:
				print(rname, "Solvation was successful")

		else:
			print("Required file(s) for " + rname + " solvation not found / generated. Solvation can not continue")
			printNote("All needed files for manual solvation will be gathered into gmxmds subfolder")

		if 'notipw' in os.listdir():
			try:
				os.rename('fsolvated.gro', 'ufsolvate.gro')
			except Exception as e:
				print("Something went wrong with error", e)

		os.chdir('../')

		time.sleep(10)

		######################################
		# GATHERING FILES FOR MD SIMULATIONS #
		######################################
		# Making gmxmds directory and populate it with needed data
		print("Gathering files needed for MDS run ...")

		gmxname = "gmxmds"
		os.mkdir(gmxname)
		os.chdir(gmxname)

		soldirlist = os.listdir(os.path.join(workhost_dir, 'Solvation'))
		for file in soldirlist:
			if Path(file).suffix == ".itp":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == "topol.top":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == "tlptopol.top":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == "fsolvated.gro":
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')
			elif file == 'tlpSolvated.gro' or file == 'tlpSolvated.top':
				shutil.copy(os.path.join(workhost_dir, 'Solvation', file), './')

		# Removing files not needed
		print("Removing files not needed for MDS from gmxmds directory")
		print("These files can still be found in Solvation directory")

		removal = ['tlptopol_mn.itp', 'ffnonbonded.itp']
		for rmf in removal:
			try:
				os.remove(rmf)
			except:
				pass

		time.sleep(5)

		if not grosolvated == "fsolvated.gro":
			print("gmxmds subfolder has been populated and ready manual solvation")
			time.sleep(5)

		else:
			printNote("##############################################################################")
			print("# The current posre.itp restrain all heavy atoms which include Backbone atoms")
			print("# You can generate your desired restrain file if this does not meet your need")
			print("# This should be named posre_udp. To use posre.itp, define -DPOSRE in .mdp files")
			print("# To use posre_udp.itp instead, define -DPOSRE_UDP in .mdp files")
			printNote("##############################################################################")

			# Generating new restraint file if needed
			response = input("To generate a new restraint interactively, type YES/y, otherwise the default will be used: ")
			if response.lower() == "yes" or response.lower() == "y":
				success = 0
				try:
					subprocess.run('gmx genrestr -f ' + grosolvated + ' -o posre_udp.itp', shell=True)
				except subprocess.CalledProcessError as e:
					print(e)
					print("Something went wrong, Trying again ....")			
					time.sleep(5)
					Err1 = os.system('gmx genrestr -f ' + grosolvated + ' -o posre_udp.itp')
					if not Err1 == 0:
						print("Generating restraint failed. Please try this manually")
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

			else:
				shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')

				insertUDP = "; Include water topology"
				insertdetails('topol.top', 'mt-topol2.itp', insertUDP)
				os.remove('mt-topol2.itp')

				printNote("No posre_udp.itp has been generated. If need be, generate it manually")

			print("gmxmds subfolder has been populated and ready for use for MD simulation of", rname)
			time.sleep(5)

		return 'fsolvated.gro'

	# Create working folder using the generated name
	os.mkdir(foldername)
	os.chdir(foldername)

	work_dir = os.path.join(GMX_MDS, foldername)
	print("pyGROMODS current workspace directory set to ", work_dir)
	
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

	time.sleep(5)

	if 'Peptides' in listfPDB:		
		print("pyGROMODS Peptides project directory set to ", work_dir)
	
		pepname = "Peptide_"
		count = 0
		for pep in PEPfiles:
			# Create host directory for each pair of receptor and ligand
			count += 1
			pep_dir = pepname + str(count)
			os.mkdir(pep_dir)
			os.chdir(pep_dir)
		
			workhost_dir = os.path.join(work_dir, pep_dir)
			print("Current project host directory set to ", workhost_dir)

			####################################################
			# GENERATE PROTEIN STRUCTURE FROM PEPTIDE SEQUENCE #
			####################################################

			rname = 'peptide' + str(count)
			pname = rname + '.pdb'

			os.mkdir('pepPDB')
			os.chdir('pepPDB')

			# Copy needed files into the directory
			shutil.copy(os.path.join(fPDB, 'Peptides', pep), './')
			shutil.copy(os.path.join(PLD, 'peptopdbtleap.in'), './')
			shutil.copy(pep, 'pep.in')

			peppdb_dir = os.path.join(workhost_dir, 'pepPDB')
			print("Peptide structure data directory set to ", peppdb_dir)

			printNote("Convertion of Peptide sequence into pdb structure file...")

			# Run amber tleap to convert to pdb
			try:
				subprocess.run('tleap -s -f peptopdbtleap.in', shell=True)
			except subprocess.CalledProcessError as e:
				print(e)
				print("Something went wrong, Trying again ....")			
				time.sleep(5)
				Err2 = os.system('tleap -s -f peptopdbtleap.in')
				if not Err2 == 0:
					printWarning("Peptide preparation failed. The process will be aborted")
					printNote("Check the above error and make corrections and rerun. Read README.md for help")
					raise Exception("Process Aborted. Make necessary corrections and restart")

			shutil.copy('sequence.pdb', '../')
			os.chdir('../')
			os.rename('sequence.pdb', pname)

			shutil.copy(os.path.join(fPDB, 'Peptides', pep), './')
			os.rename(pep, 'ppmore.in')
			shutil.copy(os.path.join(PLD, 'peptoptleap.in'), './')
			tleapfile = "peptoptleap.in"

			pepsolvated = pptopsolgen(pname, rname, tleapfile)

			os.chdir('../')
			for fitem in os.listdir():
				if not os.path.isdir(fitem):
					os.remove(fitem)

			os.chdir('../')
		os.chdir('../')

	if 'Proteins' in listfPDB:
		print("pyGROMODS Proteins project directory set to ", work_dir)

		proname = "Protein_"
		count = 0
		for pro in PROfiles:
			# Create host directory for each protein to be setup for MDS
			count += 1
			pro_dir = proname + str(count)
			os.mkdir(pro_dir)
			os.chdir(pro_dir)

			workhost_dir = os.path.join(work_dir, pro_dir)
			print("Current project host directory set to ", workhost_dir)

			rrname = 'protein' + str(count)
			prname = rrname + '.pdb'
			shutil.copy(os.path.join(fPDB, 'Proteins', pro), './')
			shutil.copy(pro, prname)
			shutil.copy(os.path.join(PLD, 'protoptleap.in'), './')
			tleapfile = "protoptleap.in"

			protsolvated = pptopsolgen(prname, rrname, tleapfile)

			os.chdir('../')
			for fitem in os.listdir():
				if not os.path.isdir(fitem):
					os.remove(fitem)

			os.chdir('../') 
		os.chdir('../')

	print('\n')
	printNote("Setup with PPmore route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS")
