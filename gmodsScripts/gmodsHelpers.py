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
import glob
import subprocess
import time
from pathlib import Path
import shutil
import random
import string
import math
from colored import fore, back, style
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
	print(fore.WHITE + back.RED + message + style.RESET)


def printNote(message):
	print(fore.PURPLE_4B + back.LIGHT_YELLOW + message + style.RESET)


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
	title = "Select Gromacs tops Directory"
	# watermodel = ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p', 'tips3p']

	try:
		gmxtopdir = os.path.join(os.environ.get('GMXDATA'), 'top')
	except:
		print("We are unable to autodetect your forcefield directory")
		print("Please specify the directory using the open folder explorer")
		print("Usually /path-to-gromacs/share/gromacs/top")
		while True:
			gmxtopdir = select_folder(title)
			if not os.path.isdir(gmxtopdir):
				print("Directory you specified does not exist. Please check and try again")
			else:
				break

	listtopdir = os.listdir(gmxtopdir)
	for fdir in listtopdir:
		if Path(fdir).suffix == ".ff":
			gmxtopff.append(Path(fdir).stem)

	return gmxtopff, gmxtopdir


def ligtopol(ligand, name):
	""" Generating parameters for ligands """
	# Setup some variables
	liggro = name + ".gro"
	ligtop = name + ".top"
	ligmol = name + ".mol2"
	ligprm = name + ".prmtop"
	liginp = name + ".inpcrd"
	ligfrc = name + ".frcmod"
	ligpdb = name + "new" + ".pdb"

    # Run Antechamber with gaff2 and parmchk2 to generate ligand(s) parameters
	try:
		subprocess.run('antechamber -i ' + ligand + ' -fi pdb -o LIG.mol2 -fo mol2 -c bcc -ek scfconv=1.d-8 -at gaff2 -pf y', shell=True, check=True)
		subprocess.run('parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod', shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong with the above error, Trying again ....")			
		time.sleep(5)
		LEr1a = os.system('antechamber -i ' + ligand + ' -fi pdb -o LIG.mol2 -fo mol2 -c bcc -ek scfconv=1.d-8 -at gaff2 -pf y')
		LEr1b = os.system('parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod')
		if not (LEr1a == 0 and LEr1b == 0):
			print(ligand, "preparation error persist. The process may fail")

	time.sleep(5)

    # Run Amber tleap with appropriate tleap source file
	try:
		subprocess.run('tleap -s -f ligtoptleap.in', shell=True, check=True)
		subprocess.run('mv LIGx.* LIG.prmtop', shell=True, check=True)
		subprocess.run('mv LIGy.* LIG.inpcrd', shell=True, check=True)
		subprocess.run('mv LIGz_new.* LIGnew.pdb', shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong with the above error message, Trying again ....")			
		time.sleep(5)
		LEr2a = os.system('tleap -s -f ligtoptleap.in')
		LEr2b = os.system('mv LIGx.* LIG.prmtop')
		LEr2c = os.system('mv LIGy.* LIG.inpcrd')
		LEr2d = os.system('mv LIGz_new.* LIGnew.pdb')
		if not (LEr2a == 0 and LEr2b == 0 and LEr2c == 0 and LEr2d == 0):
			print(ligand, "preparation with tleap failed. The process may most likely fail")
			flgro = "f" + liggro
			fltop = "f" + ligtop
			return flgro, fltop

	time.sleep(5)

    # Run acpype to convert antechamber/tleap files to Gromacs
	try:
		subprocess.run('acpype -p LIG.prmtop -x LIG.inpcrd -a amber2', shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong with above error, Trying again ....")			
		time.sleep(5)
		LEr3a = os.system('acpype.py -p LIG.prmtop -x LIG.inpcrd -a amber2')
		if not LEr3a == 0:
			print(ligand, "preparation with acpype failed. Please check, make corrections and rerun. The process will most likely fail")
			flgro = "f" + liggro
			fltop = "f" + ligtop
			return flgro, fltop

	try:
		subprocess.run('mv *_GMX.gro ' + liggro, shell=True, check=True)
		subprocess.run('mv *_GMX.top ' + ligtop, shell=True, check=True)
		subprocess.run('mv LIG.mol2 ' + ligmol, shell=True, check=True)
		subprocess.run('mv LIG.prmtop ' + ligprm, shell=True, check=True)
		subprocess.run('mv LIG.inpcrd ' + liginp, shell=True, check=True)
		subprocess.run('mv LIG.frcmod ' + ligfrc, shell=True, check=True)
		subprocess.run('mv LIGnew.pdb ' + ligpdb, shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong with renaming one or more needed files, Trying again ....")			
		time.sleep(5)
		LEr4a = os.system('mv *_GMX.gro ' + liggro)
		LEr4b = os.system('mv *_GMX.top ' + ligtop)
		LEr4c = os.system('mv LIG.mol2 ' + ligmol)
		LEr4d = os.system('mv LIG.prmtop ' + ligprm)
		LEr4e = os.system('mv LIG.inpcrd ' + liginp)
		LEr4f = os.system('mv LIG.frcmod ' + ligfrc)
		LEr4g = os.system('mv LIGnew.pdb ' + ligpdb)
		if not (LEr4a == 0 and LEr4b == 0 and LEr4c == 0 and LEr4d == 0 and LEr4e == 0 and LEr4f == 0 and LEr4g == 0):
			printNote("renaming of some needed files failed. This may affect your work")
			flgro = "f" + liggro
			fltop = "f" + ligtop
			return flgro, fltop

	return liggro, ligtop


def receptopol(receptor, name, selff, selwater):
	""" Generating receptor topology and gro files """
	recpdb = name + ".pdb"
	rectop = name + ".top"

	try:
		subprocess.run('gmx pdb2gmx -f ' + receptor + ' -p ' + rectop + ' -o ' + recpdb  + ' -ff ' + selff + ' -water ' + selwater + ' -ignh', shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong. Check above error message. Trying again ....")			
		time.sleep(5)
		REr1a = os.system('gmx pdb2gmx -f ' + receptor + ' -p ' + rectop + ' -o ' + recpdb + ' -ff ' + selff + ' -water ' + selwater + ' -ignh')
		if not REr1a == 0:
			print(receptor, "preparation failed. The process can't continue. Make corrections and rerun")
			raise Exception("Process Aborted. Make necessary corrections and restart")

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
	try:
		subprocess.run('tleap -s -f ' + tleapfile, shell=True, check=True)
		subprocess.run('mv complex1* complex.inpcrd', shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong. Check above error message. Trying again ....")			
		time.sleep(5)
		CEr1a = os.system('tleap -s -f ' + tleapfile)
		CEr1b = os.system('mv complex1* complex.inpcrd')
		if not (CEr1a == 0 and CEr1b == 0):
			printWarning("Errors detected and complex generation process failed. Please check and make corrections")
			printNote("However, attempt will be made to use alternative approach")
			return 'fComplex.gro', 'fComplex.top'

	time.sleep(5)

	# Convert Amber tleap generated prmtop and inpcrd files to Gromacs compatible format using acpype
	try:
		subprocess.run('acpype -p complex.prmtop -x complex.inpcrd -a amber2', shell=True, check=True)
	except subprocess.CalledProcessError as e:
		print(e)
		printWarning("Something went wrong, Trying again ....")			
		time.sleep(5)
		CEr1c = os.system('acpype -p complex.prmtop -x complex.inpcrd -a amber2')
		if not CEr1c == 0:
			printWarning("Preparing Complex or protein structure with acpype failed")
			printNote("It is advisable to check, make corrections where necessary, and rerun")
			printNote("However, attempt will be made to use alternative approach")
			return 'fcomplex_GMX.gro', 'fcomplex_GMX.top'

	try:
		os.rename('complex_GMX.gro', 'Complex.gro')
		os.rename('complex_GMX.top', 'Complex.top')
	except Exception as e:
		print("One or more files could not be found. Renaming failed with error", e)
		return 'complex_GMX.gro', 'complex_GMX.top'

	checkfile = os.listdir()
	if ('tlpSolvated.prmtop' in checkfile and 'tlpSolvated.inpcrd' in checkfile):
		try:
			subprocess.run('acpype -p tlpSolvated.prmtop -x tlpSolvated.inpcrd -a amber2', shell=True, check=True)
		except subprocess.CalledProcessError as e:
			print(e)
			printWarning("Something went wrong, Trying again ....")			
			time.sleep(5)
			CEr1d = os.system('acpype -p tlpSolvated.prmtop -x tlpSolvated.inpcrd -a amber2')
			if not CEr1d == 0:
				printWarning("Preparing solvated structure with acpype failed")
				printNote("This may not affect your setup. If you need it, do it manually")
				return 'Complex.gro', 'Complex.top'

		try:
			os.rename('tlpSolvated_GMX.gro', 'tlpSolvated.gro')
			os.rename('tlpSolvated_GMX.top', 'tlpSolvated.top')
		except Exception as e:
			print("One or more files could not be found. Renaming failed with error", e)

	return 'Complex.gro', 'Complex.top'


def pdbcatogro():
	""" Alternative generation of needed complex structure(s) using direct concatenation of relevant files """

	# Open a new catComplex.pdb file and write receptor pdb file into it
	complex = open("conComplex.pdb", "+a")
	R = open("nReceptor.pdb", "r")
	Rlines = R.readlines()
	for rline in Rlines:
		if not (rline.split()[0] == "TER" or rline.split()[0] == "END" or rline.split()[0] == "ENDMDL"):
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
				if not (Lline.split()[0] == "TER" or Lline.split()[0] == "END" or Lline.split()[0] == "ENDMDL"):
					complex.write(Lline)
			L.close()

	# Append the end of pdb characters - TER and END
	complex.write('TER')
	complex.write('\n')
	complex.write('END')
	complex.write('\n')
	complex.close()

	# Convert the pdb to gro using editconf
	PEr1a = os.system('gmx editconf -f conComplex.pdb -o conComplex.gro')
	if not PEr1a == 0:
		printNote("pdb conversion to gro to generated alternative complex/protein structure failed")
		printNote("This may not affect your setup")
		return 'fconComplex.gro', 'fconComplex.pdb'
	
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
	select_topol = topol
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
		if fchk == 'check' or fchk == 'check1' or fchk == 'check2':
			bk = '#'
			while True:
				newchk = bk + fchk + bk
				if newchk in os.listdir():
					bk = bk + '#'
				else:
					os.rename(fchk, newchk)
					break

	# Set needed variables
	select_topol = topol
	select_cs = "spc216"
	Twater = selwater

	change_d = seld
	change_bt = selbt

	if not Twater == "none":
		if Twater == "tip4p" or Twater == "tip4pew":
			select_cs = "tip4p"
		elif Twater == "tip5p" or Twater == "tip5pe":
			select_cs = "tip5p"
		else:
			select_cs = "spc216"

		printNote("If successful, It's important to view and check the appropriateness of the generated fsolvated.gro")
		time.sleep(5)

		# Now set the box with supplied values -d and -bt and solvate
		try:
			subprocess.run('gmx editconf -f ' + grocomplex + ' -o complex_new -d ' + str(change_d) + ' -bt ' + change_bt, shell=True)
		except subprocess.CalledProcessError as e:
			solErr1 = os.system('gmx editconf -f ' + grocomplex + ' -o complex_new -d ' + str(change_d) + ' -bt ' + change_bt)
			if not solErr1 == 0:
				printWarning("gmx editconf failed. Solvation may not be successful")

		if not 'ions.mdp' in os.listdir():
			try:
				subprocess.run('touch ions.mdp', shell=True)
			except subprocess.CalledProcessError as e:
				solErr2 = os.system('touch ions.mdp')
				if not solErr2 == 0:
					printWarning("Unable to generate needed file. Solvation may not be successful")

		# Setting up solvation and generating needed files
		rq = 0
		while True:
			try:
				subprocess.run('gmx solvate -cp complex_new -cs ' + select_cs + ' -p ' + select_topol + ' -o solvated', shell=True)
				subprocess.run('gmx grompp -f ions.mdp -p ' + select_topol + ' -c solvated.gro -o ions.tpr -maxwarn 40 2>> check', shell=True)
			except subprocess.CalledProcessError as e:
				solErr3 = os.system('gmx solvate -cp complex_new -cs ' + select_cs + ' -p ' + select_topol + ' -o solvated')
				solErr4 = os.system('gmx grompp -f ions.mdp -p ' + select_topol + ' -c solvated.gro -o ions.tpr -maxwarn 40 2>> check')
				if not (solErr3 == 0 and solErr4 == 0):
					printWarning("Something appears to be wrong. Checking ....")

			# check to be sure required file, ions.tpr, was generated
			if not 'ions.tpr' in os.listdir():
				rq += 1
				rqCheck = 'check' + str(rq)
				try:
					os.rename('check', rqCheck)
				except:
					pass
				printWarning("Warning...")
				print(select_topol, "failed to generate required file(s). For details, check the file named", rqCheck)

				if rq > 1:
					return grocomplex

				elif not 'tlptopol.top' in os.listdir():
					print("tlptopol.top was not created. We can't try solvation with it")
					time.sleep(5)
					return grocomplex

				else:
					select_topol = "tlptopol.top"
					print("Trying solvation with", select_topol, "...")
					time.sleep(5)
					continue

			else:
				if select_topol == "topol.top":
					printNote("Solvation with topol.top was successful. Congratulations!")
					printNote("The alternative tlptopol.top will be included in the gmxmds folder")
					time.sleep(5)
					break

				else:
					printNote("It apprears all is well using tlptopol.top. Please check")
					print("We will now backup topol.top as topol.bk and rename tlptopol.top to topol.top")
					os.rename('topol.top', 'topol.bk')
					os.rename('tlptopol.top', 'topol.top')
					select_topol = "topol.top"
					time.sleep(5)
					break

		# Now add ions to nutralize the solvated complex
		try:
			subprocess.run('printf "SOL" | gmx genion -s ions.tpr -p ' + select_topol + ' -pname NA -nname CL -neutral -conc 0.15 -o fsolvated.gro', shell=True, check = True)
		except subprocess.CalledProcessError as e:
			print(e)
			printWarning("Something went wrong, Trying again ....")			
			time.sleep(5)
			SEr1a = os.system('printf "SOL" | gmx genion -s ions.tpr -p ' + select_topol + ' -pname NA -nname CL -neutral -conc 0.15 -o fsolvated.gro')
			if not SEr1a == 0:
				printWarning("Something went wrong with nutralizing solvated complex or protein. Do this manually")

		return 'fsolvated.gro'

	else:
		printNote("No water model could be detected for this setup. Therefore, no solvation can be performed")
		printNote("As such 'ufsolvate.gro' will be returned instead of 'fsolvated.gro' file")
		os.system('gmx editconf -f ' + grocomplex + ' -o ufsolvate.gro')
		return 'ufsolvate.gro'

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
