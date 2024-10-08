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
import glob
import subprocess
import time
import shutil
from pathlib import Path
try:
	from colored import Fore, Back, Style
except ImportError:
	try:
		from colored import fore, back, style
	except ImportError:
		pass
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
		try:
			print(fore.WHITE + back.RED + message + style.RESET)
		except:
			print(message)


def printNote(message):
	try:
		print(Fore.PURPLE_4B + Back.LIGHT_YELLOW + message + Style.RESET)
	except:
		try:
			print(fore.PURPLE_4B + back.LIGHT_YELLOW + message + style.RESET)
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
	attempt = 0

	def topdir():
		title = "Select GROMACS Topology Directory"
		try:
			gmxTDir = os.path.join(os.environ.get('GMXDATA'), 'top')
		except:
			printWarning("Autodetection of GROMACS forcefield directory failed")
			print("Tyr searching with folder explorer [Hint: /path-to-gromacs/share/gromacs/top]")
			while True:
				ndir = 0
				gmxTDir = select_folder(title)
				if not os.path.isdir(gmxTDir):
					print("Directory you specified does not exist. Please check and try again")
					ndir += 1
					if ndir > 3:
						printWarning("Too many attempts. Process will proceed but may fail")
						gmxTDir = "unknown"
						break
				else:
					break
		return gmxTDir

	while True:
		gmxtopdir = topdir()
		if gmxtopdir == "unknown":
			break
		listtopdir = os.listdir(gmxtopdir)
		for fdir in listtopdir:
			if Path(fdir).suffix == ".ff":
				gmxtopff.append(Path(fdir).stem)
		
		if not len(gmxtopff) > 0:
			print("Selected directory is not a valid GROMACS topology directory or it's empty")
			if not attempt > 3:
				response = input("To continue anyway, type YES/y. Press ENTER to search again: ")
				if not (response.lower() == "yes" or response.lower() == "y"):
					attempt += 1
					continue
				else:
					break
			else:
				printWarning("Too many attempts. Process will proceed but may fail")
				break
		else:
			break

	return gmxtopff, gmxtopdir


def ligtopol(ligand, name):
	""" Generating parameters for ligands """
	# Setup some variables
	print('\n')
	print(f'Generating topology for {ligand}')
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

	try:
		shutil.move(glob.glob("LIGx.*")[0], "LIG.prmtop")
		shutil.move(glob.glob("LIGy.*")[0], "LIG.inpcrd")
		shutil.move(glob.glob("LIGz_new.*")[0], "LIGnew.pdb")
	except:
		printNote("Renaming of some needed files failed. This may affect your work")
		postCheck += 1

    # Run acpype to convert antechamber/tleap files to Gromacs or directly generate topology file
	lcmd3 = []
	if not postCheck > 0:
		print("Running acpype ...")
		lcmd3 = ['acpype', '-p', 'LIG.prmtop', '-x', 'LIG.inpcrd']
	else:
		print(f"{ligand} topology generation failed. Attempting alternative approach ...")
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
		printNote("Renaming of some needed files failed. This may affect your work")

	f.close()
	return liggro, ligtop


def receptopol(receptor, name, selff, selwater):
	""" Generating receptor topology and gro files """
	f = open('recerror.txt', 'a')

	if selwater == "tip4pew" or selwater == "tip5pe":
		printNote(f"Your preferred {selwater} water model will have to be selelcted interactively")
		selwater = "select"

	recpdb = name + '.pdb'
	rectop = name + '.top'
	rcmd = ['gmx', 'pdb2gmx', '-f', receptor, '-p', rectop, '-o', recpdb, '-ff', selff, '-water', selwater, '-ignh']

	if selff == "select" or selwater == "select":
		print('\n')
		while True:
			try:
				subprocess.run(rcmd, check=True, stderr=subprocess.STDOUT, text=True)
			except subprocess.SubprocessError as e:
				print(e)
				print(f"{receptor}, preparation failed with above error")
				printNote("Retry with a different forcefield group is RECOMMENDED")
				response = input("Type YES/y to retry, OR press ENTER to abort: ")
				if not (response.lower() == "yes" or response.lower() == "y"):
					raise Exception("Process Aborted. Make necessary corrections and restart")
				else:
					continue
			else:
				break

	else:
		print('\n')
		print("Running pdb2gmx....")
		try:
			subprocess.run(rcmd, check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
		except subprocess.SubprocessError as e:
			print(e)
			print(f"{receptor}, preparation failed with above error. Please check error file in Receptor folder")
			response = input("To proceed interactively, type YES/y, otherwise, press ENTER to abort: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				raise Exception("Process Aborted. Make necessary corrections and restart")
			else:
				printNote("Interactvie mode: Using a different forcefield group is RECOMMENDED")
				while True:
					rcmd = ['gmx', 'pdb2gmx', '-f', receptor, '-p', rectop, '-o', recpdb, '-ff', 'select', '-water', 'select', '-ignh']
					try:
						subprocess.run(rcmd, check=True, stderr=subprocess.STDOUT, text=True)
					except subprocess.SubprocessError as e:
						print(e)
						print(f"{receptor}, preparation failed with above error")
						printNote("Retry with a different forcefield group is RECOMMENDED")
						response = input("Type YES/y to retry, OR press ENTER to abort: ")
						if not (response.lower() == "yes" or response.lower() == "y"):
							raise Exception("Process Aborted. Make necessary corrections and restart")
						else:
							continue
					else:
						break

	f.close()
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
		printWarning("Complex generation process failed. Please check abve error message and 'cplxerror.txt' file")
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
		printWarning("Process encountered error. Please check above error message and 'cplxerror.txt' file")
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
		printWarning("One or more files could not be found. Please check 'cplxerror.txt' file for details")
		f.close()
		return 'complex_GMX.gro', 'complex_GMX.top'

	checkfile = os.listdir()
	if ('tlpSolvated.prmtop' in checkfile and 'tlpSolvated.inpcrd' in checkfile):
		print("Running acpype conversion for amber generated solvated complex ...")
		try:
			subprocess.run(['acpype', '-p', 'tlpSolvated.prmtop', '-x', 'tlpSolvated.inpcrd'], check=True, stderr=subprocess.STDOUT, stdout=f, text=True)
		except subprocess.SubprocessError as e:
			print(e)
			printWarning("Preparing solvated structure failed. Check 'cplxerror.txt' for details")
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
			pass

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
					break

				elif select_topol == "tlptopol.top":
					printNote("Checks successful. Solvation can now proceed with tlptopol.top!")
					print("Performing backup and renaming of relevant files...")
					shutil.move('topol.top', 'topol.bk')
					shutil.move('posre.itp', 'posre.bk')
					shutil.move('tlptopol.top', 'topol.top')
					shutil.move('posre_tlptopol.itp', 'posre.itp')
					select_topol = "topol.top"
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
		if identifier == "[ moleculetype ]" or identifier == "[ atomtypes ]":
			printWarning(f"Your topol.top file lack '{identifier}' header. Please check and make correction if needed")
			alternativeID = "; Include forcefield parameters"
			printNote(f"Trying {alternativeID} as alternative insertion point ...")
			time.sleep(5)

		elif identifier == "; Include water topology":
			printWarning(f"Your topol.top file lack the '{identifier}' subheading.")
			print("This is expected if you did not select any water for use with your system. Please check and make correction if needed")
			alternativeID = "[ system ]"
			printNote(f"Trying {alternativeID} as alternative insertion point ...")		
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
					countid += 3
				file1readlines.insert(countid, file2read)
				file1open.seek(0)
				file1open.writelines(file1readlines)
				break
			else:
				countid += 1

		if foundID == "NO":
			printNote("No needed identifiers for insertion of needed information was found in your topol.top.file")
			printWarning("The setup will most likely fail. Check and make needed corrections")
		else:
			print("Insertion using alternative identifier appears successful")

	file1open.close()
	file2open.close()


def udrestraint(structureFile):
	print('\n')
	print("To generate restraint on multiple groups or subsection of a group, an index file is needed")
	print("Please check GROMACS help on how to generate index file with your preferred selections")
	print('\n')
	response = input("Type YES/y to generate an index file now. Otherwise press ENTER to continue: ")
	if response.lower() == "yes" or response.lower() == "y":
		printNote("Generation of index file in progress ...")
		print('\n')
		printNote("After selection, to exit the make_ndx function and save, type q or quit and press ENTER")
		print('\n')
		try:
			subprocess.run(['gmx make_ndx -f ' + structureFile + ' -o indexfile.ndx'], shell=True, stderr=subprocess.STDOUT, check=True, text=True)
		except subprocess.SubprocessError as e:
			print(e)
			printWarning("Something went wrong with above error. Please check")	
	print('\n')

	ndxfile = []
	for ndxf in os.listdir():
		if Path(ndxf).suffix == ".ndx":
			ndxfile.append(ndxf)
	
	if len(ndxfile) == 0:
		printNote("No index file was generated or found")
	elif len(ndxfile) == 1:
		if not ndxfile[0] == "indexfile.ndx":
			os.rename(ndxfile[0], "indexfile.ndx")
	elif len(ndxfile) > 1:
		if not "indexfile.ndx" in os.listdir():
			printNote("More than one index files were found. Index file will be ignored")
		else:
			printNote("More than one index files were found. The file named 'indexfile.ndx' will be used")

	success = 0
	if not "indexfile.ndx" in os.listdir():
		try:
			subprocess.run(['gmx genrestr -f ' + structureFile + ' -o posre_udp.itp'], shell=True, stderr=subprocess.STDOUT, check=True, text=True)
		except subprocess.SubprocessError as e:
			print(e)
			printWarning("Something went wrong with the above error message. Please check")			
			time.sleep(5)
			success = 1

	elif "indexfile.ndx" in os.listdir():
		try:
			subprocess.run(['gmx genrestr -f ' + structureFile + ' -n indexfile.ndx -o posre_udp.itp'], shell=True, stderr=subprocess.STDOUT, check=True, text=True)
		except subprocess.SubprocessError as e:
			print(e)
			printWarning("Something went wrong with the above error message. Please check")			
			time.sleep(5)
			success = 1

	if success == 0:
		printNote("Generation of user defined restraint file {posre_udp.itp} was successful")
		return "posre_udp.itp"
	else:
		printNote("Generation of user defined restraint file was unsuccessful. Try manually")
		return structureFile


def gmxmdsFChecks(listgmxmds):
	""" Determining selected line index of Gromacs compatible topology files """
	filestoretain = []
	for topolfile in listgmxmds:
		if Path(topolfile).suffix == ".top":
			fileT = open(topolfile, "r")
			readfileT = fileT.readlines()
	
			for lineT in readfileT:
				lineTlist = lineT.split()
				if "#include" in lineTlist and not lineTlist[1].split('"')[1] in filestoretain:
					filestoretain.append(lineTlist[1].split('"')[1])
			fileT.close()

	return filestoretain


def gmxmdsFClean(gfile):
	# Removing unnecessary spaces in the file
	gfcount = 0
	with open(gfile, "r") as gfr:
		gfreadlines = gfr.readlines()

	with open(gfile, "w") as gfw:
		for gline in gfreadlines:
			glc = gfcount + 1
			if gline == gfreadlines[-1]:
				gfw.write(gline)
				gfcount += 1
			elif not (gline.split() == [] and gfreadlines[glc].split() == []):
				gfw.write(gline)
				gfcount += 1
			else:
				gfcount += 1

	return gfile


def gmxmdsFEChecks(runsdefaults):
	print('\n')
	printWarning("It's appears the process failed to produced all the needed MDS input files")
	printNote("It's strongly advisable to check, as the error may affect subsequent MDS files generation")
	print("For more information, Please check the relevant error files as listed below: ")
	print(">>>>>>>> ligerror.txt in Ligand Directory or Subdirectories")
	print(">>>>>>>> recerror.txt in Receptor Directory")
	print(">>>>>>>> cplxerror.txt or cplxBKerror.txt in Complext Directory")
	print(">>>>>>>> checklog.txt in Solvation Directory")

	print('\n')
	printNote("If the error has to do the default values, we can attempt to make corrections, and re-generate MDS input files")
	print("Type YES/y to attempt regeneration, or press ENTER to continue")
	response = input("Response is: ")
	if not (response.lower() == "yes" or response.lower() == "y"):
		printNote("You have choosen not to attempt correction and restart generation")
		printNote("It is strongly recommended to abort the process and make necessary corrections")
		confirm = input("Type YES/y to abort or press ENTER to continue anyway: ")
		if confirm.lower() == "yes" or confirm.lower() == "y":
			return "abort", runsdefaults
		else:
			return "continue", runsdefaults
	else:
		print('\n')
		printNote("You have chosen to attempt correction and restart generation")
		ndefaults = []

		printNote("PLEASE NOTE YOUR CURRENT DEFAULT VALUES: ")
		if runsdefaults[5] == "A":
			print(f"Default mode of processing is Interactive")
		else:
			print(f"Default mode of processing is Noninteractive")
		print(f"Default forcefield is {runsdefaults[0]}")
		print(f"Default water model is {runsdefaults[1]}")
		print(f"Default editconf -bt option is {runsdefaults[2]}")
		print(f"Default editconf -d option is {runsdefaults[3]}")
		print(f"Default timeout for input() request is {runsdefaults[4]}")
		print('\n')

		print("Forcefield and water model will be selected Interactively first")
		print("Then you will have a choice to revert to Noninteractive or maintain Interactive")
		ndefaults.append("select")
		ndefaults.append("select")

		print("Type YES/y to adjust the default values of -d, -bt and timeout as may be required")
		print("Press ENTER to maintain the above default values of -d, -bt and timeout")
		adjust = input("Response is: ")
		if adjust.lower() == "yes" or adjust.lower() == "y":
			def2, def3, def4 = defaults1()
			
		ndefaults.append(def2)
		ndefaults.append(def3)
		ndefaults.append(def4)
		if runsdefaults[5] == "A":
			ndefaults.append("A")
		else:
			ndefaults.append("B")

		return "regenerate", ndefaults

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
