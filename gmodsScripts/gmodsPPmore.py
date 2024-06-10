#!/usr/bin/env python

"""
    pyGROMODS-v2024.01 Release

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

from gmodsScripts.gmodsHelpers import receptopol, complexgen, solvation, insertdetails, printWarning, printNote, tinput, defaults1, defaults2, gmxmdsFChecks

from gmodsScripts.gmodsTLptopol import TLtopol

def PPmore(appDIR, gmxDIR, fdefaults):
	print('\n')
	# Set some environment variables
	scriptDIR = appDIR
	GMX_MDS = gmxDIR

	print(f"User working directory set to: {GMX_MDS}")

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
			defaults = ["select", "select", "triclinic", 0.1, 60, "A"]
		else:
			defaults[5] = "C"

	response = tinput("To adjust further the selected default values of -d, -bt and timeout type YES/y: ", defaults[4], "n")
	if response.lower() == "yes" or response.lower() == "y":
		defaults[2], defaults[3], defaults[4] = defaults1() 

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
		print(f"Your preselected forcefield group is {Ligsff}")

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

	def pptopsolgen(pname, rname, tleapfile, rundefaults):
		# Get needed variables as global
		global Ligsff
		global tff
		global PLD
		global fPDB

		f = open("pperror.txt", "a")
		udefaults = rundefaults

		# Set workhost directory similar to the main function
		workhost_dir = Path.cwd()

		#######################################################
		# GENERATING GROMACS AND ALTERNATIVE AMBER TOPOLOGIES #
		#######################################################

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
				print(f" The forcefield in yout topol file - {tff} - does not match the preselected group: {Ligsff}")
				print("It is advisable to check pperror.txt file for details")
				printNote("PLEASE NOTE - If you choose to continue: ")
				print(f"****A). Your forcefield group will be changed to match {tff}")
				print("****B). By default, any actions related to former forcefield will be ignored")

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
						print(f"{tff} does not match any known forcefield group. Please check")
						printNote("This may happen if you used a self created or modified forcefield. As such standard naming convention for forcefield should be used. E.g. Amber group of forcefields starts with amber, Gromos with gromos, etc. OR it may happen if generation of topol.top fails.")
						printWarning("It is strongly recommended to abort the process, check uploaded file for correctness, and try again. Check README.md file for some troubleshooting tips")
						printNote("To abort, Type YES/y. To continue anyway, press ENTER")
						response = tinput("Response: ", udefaults[4], "y")
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

		if udefaults[5] == "B":
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

			# We will now lock these defaults by changing mode to C
			udefaults[5] = "C"
		else:
			udefaults[1] = defaults2(rftopfile)
			if udefaults[1] == "none":
				print("No water model was detected for your system")
			else:
				print(f"Your water model as contained in topol.top file is {udefaults[1]}")

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

			response = tinput("Type YES/y to continue. Otherwise press ENTER to abort: ", defaults[4], "y")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Process aborted. Make necessary corrections and rerun")

		else:
			ppmore_gro, ppmore_top = complexgen(tleapfile)

			if not (ppmore_gro == 'Complex.gro' and ppmore_top == 'Complex.top'):
				printWarning("tleap could not successfully generate alternative topologies")
				print("Please check, and if need be abort and attempt to fix the error")

				response = tinput("Type YES/y to continue. Otherwise press ENTER to abort: ", defaults[4], "y")
				if not (response.lower() == "yes" or response.lower() == "y"):
					raise Exception("Process aborted. Make necessary corrections and rerun")

			else:
				try:
					os.rename(ppmore_gro, tlpgro)
					os.rename(ppmore_top, tlptop)
				except Exception as e:
					print(f"One or more files could not be found. Renaming failed with error {e}")

		for file in ppatopsdir:
			if file == 'tlpSolvated.gro':
				os.rename('tlpSolvated.gro', tlpsolvgro)
			elif file == 'tlpSolvated.top':
				os.rename('tlpSolvated.top', tlpsolvtop)

		os.chdir('../')

		####################################
		# SOLVATION OF PROTEIN STRUCTURE   #
		####################################
		# Creating solvation directory and populating it with required files
		print("Generating Solvated Structure ...")
		time.sleep(2)

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
			printNote("The forecfield in topol.top did not match the preferred forcefield group")
			printNote("Please crosscheck and if need be rerun")
			time.sleep(5)
			print('\n')
			
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
			shutil.copy(tlpgmxtopol, 'xtlptopol.top')

		# Time to prepared solvated structure
		selwater = udefaults[1]
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

		######################################
		# GATHERING FILES FOR MD SIMULATIONS #
		######################################
		# Making gmxmds directory and populate it with needed data
		print("Gathering files needed for MDS run ...")
		time.sleep(2)

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
		print("Removing files not needed for MDS from gmxmds directory")
		print("These files can still be found in Solvation directory")

		gmxmdsrequired = ['fsolvated.gro', 'ufsolvate.gro', 'tlpSolvated.gro', 'tlpSolvated.top', 'tlptopol.top', 'topol.top']
		for itpf in os.listdir():
			if Path(itpf).suffix == ".itp":
				gmxmdsrequired.append(itpf)
		 
		for rmf in os.listdir():
			if not rmf in gmxmdsrequired: 
				os.remove(rmf)

		# Getting Include files ready and up-to-date
		neededIncludeFiles = gmxmdsFChecks(os.listdir())
		notNeededIncludeFiles = []

		for includefile in os.listdir():
			if Path(includefile).suffix == ".itp":
				if not includefile in neededIncludeFiles:
					if not includefile in notNeededIncludeFiles:
						notNeededIncludeFiles.append(includefile)
			else:
				pass
		
		os.mkdir("not4mds")
		for notfile in notNeededIncludeFiles:
			shutil.move(notfile, os.path.join(workhost_dir, 'gmxmds', 'not4mds'))
		
		# If required, new restraint file can now be generated
		listgmxmds = os.listdir()

		if "fsolvated.gro" in listgmxmds:
			grosolvated = "fsolvated.gro"

			printNote("##############################################################################")
			print("# The current posre.itp restrain all heavy atoms which include Backbone atoms")
			print("# You can generate your desired restrain file if this does not meet your need")
			print("# This should be named posre_udp. To use posre.itp, define -DPOSRE in .mdp files")
			print("# To use posre_udp.itp instead, define -DPOSRE_UDP in .mdp files")
			printNote("##############################################################################")

			# Generating new restraint file if needed
			response = tinput("To generate a new restraint interactively, type YES/y, or press ENTER: ", defaults[4], "n")
			if response.lower() == "yes" or response.lower() == "y":
				success = 0
				try:
					subprocess.run(['gmx', 'genrestr', '-f', grosolvated, '-o', 'posre_udp.itp'], stderr=subprocess.STDOUT, stdout=f, check=True, text=True)
				except subprocess.SubprocessError as e:
					print(f"Generating restraint failed with error {e}. Please try this manually")
					success = 1

				if success == 0:
					printNote("You have generated posre_udp.itp. To use it in MDS, define -DPOSRES_UDP in .mdp files. OR if already define -DPOSRE, back up posre.itp and rename posre_udp.itp to posre.itp")
				else:
					printNote("No posre_udp.itp has been generated. If need be, generate it manually")

			else:
				printNote("No posre_udp.itp has been generated. If need be, generate it manually")

			shutil.copy(os.path.join(scriptDIR, 'gmodsScripts', 'mt-topol2.itp'), './')
			insertUDP = "; Include water topology"
			insertdetails('topol.top', 'mt-topol2.itp', insertUDP)
			os.remove('mt-topol2.itp')

			print(f"gmxmds subfolder has been populated and ready for use for MD simulation of {rname}")
			f.close()
			return grosolvated, udefaults

		elif "ufsolvate.gro" in listgmxmds:
			grosolvated = "ufsolvate.gro"
			print("Manual solvation and/or generation of restraint is required if necessary")
			print("gmxmds subfolder has been populated and ready for use")
			f.close()
			return grosolvated, udefaults

		elif "tlpSolvated.gro" in listgmxmds:
			grosolvated = "tlpSolvated.gro"
			print("Manual generation of restraint for tleap generated solvated complex is required, if necessary")
			print("gmxmds subfolder has been populated and ready for use")
			f.close()
			return grosolvated, udefaults

		else:
			grosolvated = "none"
			print("No solvated or unsolvated complex was detected in gmxmds subfolder. Please check")
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
		print(f"pyGROMODS Peptides project directory set to {work_dir}")
	
		pepname = "Peptide_"
		count = 0
		for pep in PEPfiles:
			# Create host directory for each pair of receptor and ligand
			count += 1
			pep_dir = pepname + str(count)
			os.mkdir(pep_dir)
			os.chdir(pep_dir)
		
			workhost_dir = os.path.join(work_dir, pep_dir)
			print(f"Current project host directory set to {workhost_dir}")

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

			pepsolvated, udefaults = pptopsolgen(pname, rname, tleapfile, defaults)
			defaults = udefaults

			os.chdir('../')
			for fitem in os.listdir():
				if not os.path.isdir(fitem):
					os.remove(fitem)

			os.chdir('../')

	if 'Proteins' in listfPDB:
		print(f"pyGROMODS Proteins project directory set to {work_dir}")

		proname = "Protein_"
		count = 0
		for pro in PROfiles:
			# Create host directory for each protein to be setup for MDS
			count += 1
			pro_dir = proname + str(count)
			os.mkdir(pro_dir)
			os.chdir(pro_dir)

			workhost_dir = os.path.join(work_dir, pro_dir)
			print(f"Current project host directory set to {workhost_dir}")

			rrname = 'protein' + str(count)
			prname = rrname + '.pdb'
			shutil.copy(os.path.join(fPDB, 'Proteins', pro), './')
			shutil.copy(pro, prname)
			shutil.copy(os.path.join(PLD, 'protoptleap.in'), './')
			tleapfile = "protoptleap.in"

			protsolvated, udefaults = pptopsolgen(prname, rrname, tleapfile, defaults)
			defaults = udefaults

			os.chdir('../')
			for fitem in os.listdir():
				if not os.path.isdir(fitem):
					os.remove(fitem)

			os.chdir('../') 

	print('\n')
	print("PLEASE NOTE:")
	print("The files in folder 'not4mds' of 'gmxmds' are considered not necessary for MDS. Please check")

	print('\n')
	printNote("Setup with PPmore route completed. Please analyse the contents of 'check', 'check1' and/or 'check2' files and their backup versions in 'Solvation' folder before proceeding with MDS")
