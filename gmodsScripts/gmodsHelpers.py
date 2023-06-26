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
import glob
import subprocess
import time
import shutil
from pathlib import Path
from colored import Fore, Back, Style
from tkinter import Tk, filedialog
from inputimeout import inputimeout, TimeoutOccurred
from pytimedinput import timedInput


def tinput(message, timeout, default):
	def timer2(message, timeout, default):
		try:
			userinput = inputimeout(prompt=message, timeout=timeout)
		except TimeoutOccurred:
			userinput = default
		return userinput

	def timer1(message, timeout, default):
		userinput, texpired = timedInput(message, timeout=timeout)
		if(texpired):
			userinput = default
			return userinput
		else:
			return userinput

	try:
		Response = timer1(message, timeout, default)
	except:
		Response = timer2(message, timeout, default)

	return Response


def printWarning(message):
	try:
		print(Fore.WHITE + Back.RED + message + Style.RESET)
	except:
		print(message)


def printNote(message):
	try:
		print(Fore.PURPLE_4B + Back.LIGHT_YELLOW + message + Style.RESET)
	except:
		print(message)


def select_folder(title):
	root = Tk()
	root.withdraw()
	root.attributes('-topmost', True)
	sFolder = ""

	while True:
		sFolder = filedialog.askdirectory(title=title, initialdir=".")
		root.destroy()
		os.path.normpath(sFolder)
		if os.path.isdir(sFolder):
			break
		elif os.path.isfile(sFolder):
			print("You must select folder, not a file")
			continue
		else:
			print("The path you specify does not exist")
			continue

	return sFolder


def gmxtop():
	""" Scan for and obtain the list of available forcefields based on install version of gromacs """
	# Get the absolute path to the forcefield directory and get the awailabe ff
	gmxtopdir = " "
	gmxtopff = []

	def topdir():
		title = "Select Gromacs tops Directory"
		try:
			gmxTDir = os.path.join(os.environ.get('GMXDATA'), 'top')
		except:
			print("Autodetection of gromacs forcefield directory failed")
			print("Use the folder explorer to searh. Usually /path-to-gromacs/share/gromacs/top")
			while True:
				gmxTDir = select_folder(title)
				if not os.path.isdir(gmxTDir):
					print("Directory you specified does not exist. Please check and try again")
				else:
					break
		return gmxTDir

	while True:
		gmxtopdir = topdir()
		listtopdir = os.listdir(gmxtopdir)
		for fdir in listtopdir:
			if Path(fdir).suffix == ".ff":
				gmxtopff.append(Path(fdir).stem)
		
		if not len(gmxtopff) > 0:
			print("Selected directory is not a valid gromacs toplogy directory or it's empty")
			response = input("To continue anyway, type YES/y. Press ENTER to search again: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				continue
			else:
				break
		else:
			break

	return gmxtopff, gmxtopdir


def ligtopol(ligand, name):
	""" Generating parameters for ligands """
	# Setup some variables
	print('\n')
	liggro = name + ".gro"
	ligtop = name + ".top"
	ligmol = name + ".mol2"
	ligprm = name + ".prmtop"
	liginp = name + ".inpcrd"
	ligfrc = name + ".frcmod"
	ligpdb = name + "new" + ".pdb"
	ligfi = "pdb"
	flgro = "f" + liggro
	fltop = "f" + ligtop

	postCheck = 0
	f = open('ligerror.txt', 'a')

    # Run Antechamber with gaff2 and parmchk2 to generate ligand(s) parameters
	print("Running antechamber ...")
	if Path(ligand).suffix == ".mol2":
		ligfi = "mol2"

	lcmd1 = ['antechamber', '-i', ligand, '-fi', ligfi, '-o', 'LIG.mol2', '-fo', 'mol2', '-c', 'bcc', '-ek', 'scfconv=1.d-8', '-at', 'gaff2', '-pf', 'y']
	lcmd2 = ['parmchk2', '-i', 'LIG.mol2', '-f', 'mol2', '-o', 'LIG.frcmod']

	try:
		subprocess.run(lcmd1, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
		subprocess.run(lcmd2, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning("Something went wrong with the above error")			
		print(f"{ligand} preparation process may fail")
		postCheck += 1

    # Run Amber tleap with appropriate tleap source file
	print("Running tleap ...")
	try:
		subprocess.run(['tleap', '-s', '-f', 'ligtoptleap.in'], check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning(f"Something went wrong with the above error message. {ligand} preparation process may most likely fail")
		postCheck += 1	

	print("Renaming generated files ...")
	try:
		shutil.move(glob.glob("LIGx.*")[0], "LIG.prmtop")
		shutil.move(glob.glob("LIGy.*")[0], "LIG.inpcrd")
		shutil.move(glob.glob("LIGz_new.*")[0], "LIGnew.pdb")
	except:
		printWarning("Renaming of some needed files failed. This may affect your work")
		postCheck += 1

    # Run acpype to convert antechamber/tleap files to Gromacs or directly generate topology file
	lcmd3 = []
	if not postCheck > 0:
		print("Running acpype ...")
		lcmd3 = ['acpype', '-p', 'LIG.prmtop', '-x', 'LIG.inpcrd']
	else:
		print("Attempting to use acpype for topology generation ...")
		lcmd3 = ['acpype', '-i', ligand, '-b', name]
		
	try:
		subprocess.run(lcmd3, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning(f"Something went wrong with above error, {ligand} preparation process will most likely fail")
		print("Please check the content of ligerror.txt file")
		f.close()
		return flgro, fltop

	for folder in os.listdir():
		if os.path.isdir(folder) and not folder == "usertops":
			for facpype in os.listdir(folder):
				if Path(facpype).suffix == ".itp" or Path(facpype).suffix == ".top" or Path(facpype).suffix == ".gro":
					shutil.copy(os.path.join(Path.cwd(), folder, facpype), './')

			shutil.rmtree(folder)

	try:
		shutil.move(glob.glob("*_GMX.gro")[0], liggro)
		shutil.move(glob.glob("*_GMX.top")[0], ligtop)
		shutil.move("LIG.mol2", ligmol)
		shutil.move("LIG.prmtop", ligprm)
		shutil.move("LIG.inpcrd", liginp)
		shutil.move("LIG.frcmod", ligfrc)
		shutil.move("LIGnew.pdb", ligpdb)
	except:
		printWarning("Renaming of some needed files failed. This may affect your work")

	f.close()
	return liggro, ligtop


def receptopol(receptor, name, selff, selwater):
	""" Generating receptor topology and gro files """
	recpdb = name + '.pdb'
	rectop = name + '.top'
	rcmd = ['gmx', 'pdb2gmx', '-f', receptor, '-p', rectop, '-o', recpdb, '-ff', selff, '-water', selwater, '-ignh']

	try:
		subprocess.run(rcmd, check=True, stderr=subprocess.STDOUT, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		print(f"{receptor}, preparation failed with above error")
		response = input("To proceed interactively, type YES/y, otherwise, press ENTER to abort: ")
		if not (response.lower() == "yes" or response.lower() == "y"):
			raise Exception("Process Aborted. Make necessary corrections and restart")
		else:
			rcmd = ['gmx', 'pdb2gmx', '-f', receptor, '-p', rectop, '-o', recpdb, '-ff', 'select', '-water', 'select', '-ignh']
			subprocess.run(rcmd, check=True, stderr=subprocess.STDOUT, text=True)

	return rectop, recpdb, 'posre.itp'


def topolsplit(LFtop, name, lineindex):
	""" Spliting ligands topology into atomtypes and moleculetype """
	file1 = open(LFtop, "r")
	readline = file1.readlines()
	LIGS_at = name + "_at.itp"
	LIGS_mt = name + "_mt.itp"
	file2 = open(LIGS_at, "+a")
	file3 = open(LIGS_mt, "+a")
	x = lineindex['atomtypes']
	while x < lineindex['moleculetype']:
		try:
			file2.write(readline[x])
		except IndexError:
			pass
		x += 1

	y = lineindex['moleculetype']
	while y < lineindex['system']:
		try:
			file3.write(readline[y])
		except IndexError:
			pass
		y += 1

	file1.close()
	file2.close()
	file3.close()

	return LIGS_at, LIGS_mt


def indexoflines(LFtop):
	""" Determining selected line index of Gromacs compatible topology files """
	file1 = open(LFtop, "r")
	readline = file1.readlines()
	lineindex = ["x", "x", "x"]
	n = 0
	for line in readline:
		linelist = line.split()
		if "atomtypes" in linelist:
			lineindex[0] = n
			n += 1
		elif "moleculetype" in linelist:
			lineindex[1] = n
			n += 1
		elif "system" in linelist:
			lineindex[2] = n
			n += 1
		else:
			n += 1
	file1.close()

	Idx = 0
	while Idx < len(lineindex):
		if not str(lineindex[Idx]).isnumeric() == True:
			lineindex[Idx] = n + 1
			Idx += 1
		else:
			Idx += 1

	return {'atomtypes': lineindex[0], 'moleculetype': lineindex[1], 'system': lineindex[2]}


def complexgen(tleapfile):
	""" Generating needed complex structure(s) using amber tleap """
    # Run Amber tleap with appropriate tleap source file
	print("Running amber tleap ...")
	f = open('cplxerror.txt', 'a')
	try:
		subprocess.run(['tleap', '-s', '-f', tleapfile], check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning("The above errors detected and complex generation process failed. Please check 'cplxerror.txt' file and make corrections")
		printNote("However, attempt will be made to use alternative approach")
		f.close()
		return 'fComplex.gro', 'fComplex.top'

	try:
		shutil.move(glob.glob("complex1*")[0], "complex.inpcrd")
	except:
		pass

	# Convert Amber tleap generated prmtop and inpcrd files to Gromacs compatible format using acpype
	print("Running acpype ...")
	try:
		subprocess.run(['acpype', '-p', 'complex.prmtop', '-x', 'complex.inpcrd'], check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
	except subprocess.SubprocessError as e:
		print(e)
		printWarning("Preparing Complex or protein structure failed with above error")
		printNote("It is advisable to check 'cplxerror.txt' file, make corrections where necessary, and rerun")
		printNote("However, attempt will be made to use alternative approach")
		f.close()
		return 'fcomplex.gro', 'fcomplex.top'

	for folder1 in os.listdir():
		if os.path.isdir(folder1) and not (folder1 == "ligfrcmod" or folder1 == "ligmol" or folder1 == "ligpdb"):
			for fcomp in os.listdir(folder1):
				if Path(fcomp).suffix == ".itp" or Path(fcomp).suffix == ".top" or Path(fcomp).suffix == ".gro":
					shutil.copy(os.path.join(Path.cwd(), folder1, fcomp), './')

			shutil.rmtree(folder1)

	try:
		shutil.move('complex_GMX.gro', 'Complex.gro')
		shutil.move('complex_GMX.top', 'Complex.top')
	except Exception as e:
		print("One or more files could not be found. Renaming failed with error", e)
		printNote("Please check 'cplxerror.txt' file for details")
		f.close()
		return 'complex_GMX.gro', 'complex_GMX.top'

	checkfile = os.listdir()
	if ('tlpSolvated.prmtop' in checkfile and 'tlpSolvated.inpcrd' in checkfile):
		print("Running acpype conversion for amber generated solvated complex ...")
		try:
			subprocess.run(['acpype', '-p', 'tlpSolvated.prmtop', '-x', 'tlpSolvated.inpcrd'], check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
		except subprocess.SubprocessError as e:
			print(e)
			printWarning("Preparing solvated structure failed with above error. Check 'cplxerror.txt' for details")
			printNote("This may not affect your setup. If you need it, do it manually")
			f.close()
			return 'Complex.gro', 'Complex.top'

		for folder2 in os.listdir():
			if os.path.isdir(folder2) and not (folder2 == "ligfrcmod" or folder2 == "ligmol" or folder2 == "ligpdb"):
				for fambsol in os.listdir(folder2):
					if Path(fambsol).suffix == ".itp" or Path(fambsol).suffix == ".top" or Path(fambsol).suffix == ".gro":
						shutil.copy(os.path.join(Path.cwd(), folder2, fambsol), './')

				shutil.rmtree(folder2)

		try:
			shutil.move('tlpSolvated_GMX.gro', 'tlpSolvated.gro')
			shutil.move('tlpSolvated_GMX.top', 'tlpSolvated.top')
		except Exception as e:
			print("One or more files could not be found. Renaming failed with error", e)
			print("You may want to check 'cplxerror.txt' file for details of the error")

	f.close()
	return 'Complex.gro', 'Complex.top'


def pdbcatogro():
	""" Alternative generation of needed complex structure(s) using direct concatenation of relevant files """
	printNote("Exploring alternative approach to complex generation")
	print("Generating backup Protein - Ligand Complex ...")

	# Open a new catComplex.pdb file and write receptor pdb file into it
	complex = open("conComplex.pdb", "+a")
	R = open("nReceptor.pdb", "r")
	Rlines = R.readlines()
	for rline in Rlines:
		if not (rline.split()[0] == "TER" or rline.split()[0] == "END" or rline.split()[0] == "ENDMDL" or rline.split()[0] == "CONECT" or rline.split()[0] == "MASTER"):
			complex.write(rline)
	R.close()

	# Get ligand structure subdirectory to read pdb files and write into Complex.pdb
	curworkdir=Path.cwd()
	listpdb = os.listdir(os.path.join(curworkdir, 'ligpdb'))
	for Lpdb in listpdb:
		Lp = os.path.join(curworkdir, 'ligpdb', Lpdb)
		if Path(Lp).suffix == ".pdb":
			L = open(Lp, "r")
			Llines = L.readlines()
			for Lline in Llines:
				if Lline.split()[0] == "HETATM" or Lline.split()[0] == "HET" or Lline.split()[0][:3] == "HET":
					complex.write(Lline)
				elif Lline.split()[0] == "ATOM":
					replace = Lline.replace(Lline.split()[0], "HETATM", 1)
					complex.write(replace)
			L.close()

	# Append the end of pdb characters - MASTER, TER and END
	complex.write('TER')
	complex.write('\n')
	complex.write('END')
	complex.write('\n')
	complex.close()

	# Convert the pdb to gro using editconf
	f = open('cplxBKerror.txt', 'a')
	try:
		subprocess.run(['gmx', 'editconf', '-f', 'conComplex.pdb', '-o', 'conComplex.gro'], check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
	except subprocess.SubprocessError as e:
		printNote("Generating backup structure failed. Check cplxBKerror.txt file for details")
		f.close()
		return 'fconComplex.gro', 'fconComplex.pdb'
	
	f.close()
	print("Backup complex was successfully generated")
	return 'conComplex.gro', 'conComplex.pdb'


def defaults1():
	""" If needed, some default values to be used henceforth will be generated """
	# Define isfloat() function to use
	def isfloat(num):
		try:
			float(num)
			return True
		except ValueError:
			return False

	# Set needed variables
	change_bt = "triclinic"
	change_d = 0.1
	timedout = 60

	printNote("Only one of triclinic, cubic, dodecahedron or octahedron is acceptable for -bt")
	printNote("Only digit/float values are acceptable for -d and timeout")

	while True:
		new_timedout = input("Change timeout to: ")
		new_d = input("Change -d to: ")
		new_bt = input("Change -bt to: ")
		if new_d == " " or new_bt == " " or new_timedout == " " or new_d == "" or new_bt == "" or new_timedout == "":
			print("You must supply digit value for -d and timeout, and either retain triclinic or change -bt to cubic, dodecahedron or octahedron to proceed")
			continue
		elif not (isfloat(new_d) == True and new_timedout.isdigit() == True and new_bt.isalpha() == True):
			print("You must supply a float value for -d, digit for timeout, and only one of triclinic, cubic, dodecahedron or octahedron is acceptable for -bt to proceed")
			continue
		elif not (new_bt.lower() == "triclinic" or new_bt.lower() == "dodecahedron" or new_bt.lower() == "octahedron" or new_bt.lower() == "cubic"):
			print("Only one of triclinic, cubic, dodecahedron or octahedron is acceptable for -bt")
			continue
		else:
			break

	change_bt = new_bt
	change_d = float(new_d)
	timedout = int(new_timedout)

	return change_bt, change_d, timedout  


def defaults2(topol):
	""" If needed, some default values will be extracted from topol file to update defaults list """
	# Set needed variables
	select_water = "none"

	# Check topol file to identify the water model used
	Topen = open(topol, "r")
	Tread = Topen.readlines()
	n = 0
	nx = 0
	for tx in Tread:
		if 'Include' in tx.split() and 'water' in tx.split() and 'topology' in tx.split(): 
			nx = n + 1
			break
		else:
			n += 1

	if nx > 0:
		Topen.seek(0)
		Twateritp = Tread[nx].split()[1].split('"')[1].split('/')[1]
		Twater = Twateritp.split('.')[0]
		select_water = Twater
		return select_water
	else:
		return select_water


def solvation(grocomplex, topol, selwater, selbt, seld):
	""" Solvate the generated complex in preparation for MDS """
	# Create needed backup for error checkup when necessary
	for fchk in os.listdir():
		if fchk == 'checklog.txt':
			bk = '#'
			while True:
				newchk = bk + fchk + bk
				if newchk in os.listdir():
					bk = bk + '#'
					continue
				else:
					os.rename(fchk, newchk)
					break

	# Set needed variables
	select_topol = topol
	select_cs = "spc216"
	Twater = selwater

	change_d = str(seld)
	change_bt = selbt
	f = open('checklog.txt', 'a')

	if not Twater == "none":
		if Twater == "tip4p" or Twater == "tip4pew":
			select_cs = "tip4p"
		elif Twater == "tip5p" or Twater == "tip5pe":
			select_cs = "tip5p"
		else:
			select_cs = "spc216"

		printNote("If successful, the appropriateness of the generated fsolvated.gro should be checked")
		time.sleep(5)

		# Now set the box with supplied values -d and -bt and solvate
		print("Setting the box with the default parameters ...")
		solcmd1 = ['gmx', 'editconf', '-f', grocomplex, '-o', 'complex_new', '-d', change_d, '-bt', change_bt]
		try:
			subprocess.run(solcmd1, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
		except subprocess.SubprocessError as e:
			printWarning("gmx editconf failed. Solvation may not be successful")

		if not 'ions.mdp' in os.listdir():
			open('ions.mdp', 'w').close()

		# Setting up solvation and generating needed files
		print("Checking solvation with topol.top ...")
		rq = 0
		solcmd2 = ['gmx', 'solvate', '-cp', 'complex_new', '-cs', select_cs, '-p', select_topol, '-o', 'solvated']
		solcmd3 = ['gmx', 'grompp', '-f', 'ions.mdp', '-p', select_topol, '-c', 'solvated.gro', '-o', 'ions.tpr', '-maxwarn', '40']

		while True:
			try:
				subprocess.run(solcmd2, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
				subprocess.run(solcmd3, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
			except subprocess.SubprocessError as e:
				printWarning("Something appears to be wrong. Checking ....")

			# check to be sure required file, ions.tpr, was generated
			if not 'ions.tpr' in os.listdir():
				rq += 1
				printWarning("Warning...")
				print(f"{select_topol} failed to generate required file(s). For details, check the file named {f}")

				if rq > 1:
					f.close()
					return grocomplex

				elif not 'tlptopol.top' in os.listdir():
					print("tlptopol.top was not created. We can't try solvation with it")
					time.sleep(5)
					f.close()
					return grocomplex

				else:
					select_topol = "tlptopol.top"
					print("Trying solvation with", select_topol, "...")
					time.sleep(5)
					continue

			else:
				if select_topol == "topol.top":
					printNote("Checks Successful. Solvation can now proceed with topol.top!")
					printNote("The alternative tlptopol.top will be included in the gmxmds folder")
					time.sleep(5)
					break

				else:
					printNote("Checks successful. Solvation can now proceed with tlptopol.top!")
					print("We will now backup topol.top as topol.bk and rename tlptopol.top to topol.top")
					shutil.move('topol.top', 'topol.bk')
					shutil.move('tlptopol.top', 'topol.top')
					select_topol = "topol.top"
					time.sleep(5)
					break

		# Now add ions to nutralize the solvated complex
		print("Adding ions to neutralized the solvated complex ...")
		try:
			solcmd4 = f'printf "SOL" | gmx genion -s ions.tpr -p {select_topol} -pname NA -nname CL -neutral -conc 0.15 -o fsolvated.gro'
			subprocess.run(solcmd4, shell=True, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
		except subprocess.SubprocessError as e:
			print(f"Something went wrong with above error: {e}")
			printWarning(f"Please check file {f} for details")

		f.close()
		return 'fsolvated.gro'

	else:
		printNote("No water model could be detected for this setup. Therefore, no solvation can be performed")
		printNote("As such 'ufsolvate.gro' will be returned instead of 'fsolvated.gro' file")
		solcmd5 = ['gmx', 'editconf', '-f', grocomplex, '-o', 'ufsolvate.gro']
		try:
			subprocess.run(solcmd5, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
			f.close()
			return 'ufsolvate.gro'
		except subprocess.SubprocessError as e:
			print(e)
			printWarning("Something went wrong with above error. Please check")
			f.close()
			return grocomplex


def insertdetails(file1, file2, identifier):
	""" Inserting new details into any part of already created file """
	file1open = open(file1, "r+")
	file1readlines = file1open.readlines()
	file2open = open(file2, "r")
	file2read = file2open.read()

	countid = 0
	foundID = "NO"
	for LK in file1readlines:
		if LK.split() == identifier.split():
			foundID = "YES"
			file1readlines.insert(countid, file2read)
			file1open.seek(0)
			file1open.writelines(file1readlines)
			break
		else:
			countid += 1
	
	if foundID == "NO":
		alternativeID = " "
		if identifier == "[ moleculetype ]":
			printWarning("Your topol.top file lack [ moleculetype ] header. Please check and make correction if needed")
			alternativeID = "; Include forcefield parameters"
			printNote("Trying alternative Ligand atomtypes insertion point ...")
			time.sleep(5)

		elif identifier == "; Include water topology":
			printWarning("Your topol.top file lack the '; Include water topology' subheading. This is expected if you did not select any water for use with your system. Please check and make correction if needed")
			alternativeID = "[ system ]"
			printNote("Trying alternative Ligand moleculetype insertion point ...")		
			time.sleep(5)

		else:
			printNote("identifier for insertion point for needed information not found in your topol.top.file")
			printWarning("The setup will most likely fail. Check and make needed corrections")
			alternativeID = identifier
			time.sleep(5)

		file1open.seek(0)
		countid = 0
		for nLK in file1readlines:
			if nLK.split() == alternativeID.split():
				foundID = "YES"
				if alternativeID == "; Include forcefield parameters":
					countid += 2
				file1readlines.insert(countid, file2read)
				file1open.seek(0)
				file1open.writelines(file1readlines)
				break
			else:
				countid += 1

		if foundID == "NO":
			printNote("identifier for insertion point for needed information not found in your topol.top.file")
			printWarning("The setup will most likely fail. Check and make needed corrections")

	file1open.close()
	file2open.close()


def cleanup():
	""" Cleaning up opencl headers added dueing MDS """
	filesclh = glob.glob('*.clh')
	for clh in filesclh:
		try:
			os.remove(clh)
		except:
			pass

	filesh = glob.glob('*.h')
	for h in filesh:
		try:
			os.remove(h)
		except:
			pass

	filescl = glob.glob('*.cl')
	for cl in filescl:
		try:
			os.remove(cl)
		except:
			pass
