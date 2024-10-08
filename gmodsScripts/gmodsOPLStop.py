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
import glob

from gmodsScripts.gmodsHelpers import indexoflines, printWarning, printNote, gmxtop

def OPLStop(LFtop, ff, name):
	mergename = "opls" + name + ".top"
	oplstopopen = open(mergename, "+a")

	# Insert the header information into the opls[name].top file
	topheaders = ["; This file was generated by modifying tleap/acpype generated ligand topology", "; by replacing generated atom types with opls compatible atom types", "; as contained in the opls ffnonbonded.itp forcefield in gromacs.", "; Where opls compatible atom type(s) is/are not found, the tleap/acpype generated atom types will be retained"]

	for header in topheaders:
		oplstopopen.write(header)
		oplstopopen.write('\n')
	oplstopopen.write('\n')

	time.sleep(5)

	# We shall check for and remove duplicates in atomtypes between gmx standard and amber/tleap generated 
	if not ff[0:4] == 'opls':
		print(ff, "is not an opls forcefield. Generation cannot continue")
		return LFtop

	else:
		print("Compatible forcefiled detected:", ff)

	# Get the absolute path to the forcefield directory and copy ffnonbonded.itp file
	gmxtopdir = " "
	topffdir = " "
	gmxtopff = []

	nT = 0
	while True:
		nT += 1
		if nT > 3:
			print("You have exceeded maximum trying attempts")
			printWarning("Checked version of tlptopol cannot be generated")
			print("The platform will use the unchecked original file")
			time.sleep(5)
			return LFtop

		gmxtopff, topffdir = gmxtop()
		gmxtopdir = os.path.join(topffdir, ff)
		
		if not (len(gmxtopff) > 0 or Path(ff).stem in gmxtopff):
			print(f"The specified forcefield, {ff}, is missing from the selected directory")
			response = input("To continue anyway, type YES/y. Otherwise, press ENTER: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				print("Trying again ...")
				continue
			else:
				printWarning("Checked version of tlptopol cannot be generated")
				print("The platform will use the unchecked original file")
				time.sleep(5)
				return LFtop
			
		else:
			break

	lsgmxtopdir = os.listdir(gmxtopdir)
	for tp in lsgmxtopdir:
		if tp == "ffnonbonded.itp":
			shutil.copy(os.path.join(gmxtopdir, tp), './')

	# Get the ffnonbonded.itp opls atomtypes range to consider for use inside new file
	ffb = open("ffnonbonded.itp", "r")
	ffbreadlines = ffb.readlines()
	fn = 0
	for fb in ffbreadlines:
		if 'opls_128' in fb.split():
			break
		else:
			fn += 1
	ffb.seek(0)

	# Check the LFtop file to replace atomtypes and moleculetype atoms with opls types where possible
	lftindex = indexoflines(LFtop)
	atlp = int(lftindex['atomtypes'])
	mtlp = int(lftindex['moleculetype'])

	lftopopen = open(LFtop, "r")
	lftopreadlines = lftopopen.readlines()

	btlp = 0
	for bd in lftopreadlines:
		if 'bonds' in bd.split() and bd.split()[0] == '[' and bd.split()[1] == 'bonds':
			break
		else:
			btlp += 1

	lftopopen.seek(0)
	for aline in lftopreadlines[atlp:mtlp]:
		if not (aline.split() == [] or aline.split()[0][0] == ';' or aline.split()[0][0] == '['):
			atom = aline.split()[0]
			ffb.seek(0)
			afound = 0
			for fline in ffbreadlines[fn:]:
				if fline.split()[1] == atom or fline.split()[1].lower() == atom or fline.split()[1].upper() == atom.upper():
					afound += 1
			if not afound > 0:
				oplstopopen.write(aline)
		else:
			oplstopopen.write(aline)
	oplstopopen.write('\n')

	lftopopen.seek(0)
	for bline in lftopreadlines[mtlp:btlp]:
		if not (bline.split() == [] or bline.split()[0][0] == ';' or bline.split()[0][0] == '['):
			matom = bline.split()[1]
			ffb.seek(0)
			replace = " "
			bfound = 0
			for fline in ffbreadlines[fn:]:
				if fline.split()[1] == matom or fline.split()[1].lower() == matom or fline.split()[1] == matom.upper() and 'opls' in fline.split()[0].split('_'):
					replace = bline.replace(matom, fline.split()[0], 1)
					bfound += 1
					break
			if not bfound > 0:
				oplstopopen.write(bline)
			else:
				oplstopopen.write(replace)
		else:
			oplstopopen.write(bline)
	oplstopopen.write('\n')

	lftopopen.seek(0)
	for cline in lftopreadlines[btlp:]:
		oplstopopen.write(cline)

	ffb.close()
	lftopopen.close()
	oplstopopen.close()

	# Check to be sure oplstopol.top has been successfully generated
	checktlp = os.listdir()
	if not mergename in checktlp:
		printWarning("Something went wrong. Generating oplstopol.top file was not successful")
		return LFtop
	else:
		checkindex = indexoflines(mergename)
		a = int(checkindex['atomtypes'])
		b = int(checkindex['moleculetype'])
		c = int(checkindex['system'])
		checkopen = open(mergename, "r")
		checkreadlines = checkopen.readlines()
		if not ('atomtypes' in checkreadlines[a].split() and 'moleculetype' in checkreadlines[b].split() and 'system' in checkreadlines[c].split()):
			return LFtop
		else:
			printNote("opls compatible ligand topology has been generated successfully")
			return mergename


def OPLSacpype(LFtop, lig, name):
	# Create a temporary working folder and copy needed files into it
	f = open("oplserror.txt", "a")

	os.mkdir("ACPYPEDIR")
	shutil.copy(lig, os.path.join(Path.cwd(), "ACPYPEDIR"))
	os.chdir("ACPYPEDIR")

	# Run acpype to generate opls compatible topoloy file 
	try:
		subprocess.run(['acpype', '-i', lig, '-b', name], stderr=subprocess.STDOUT, stdout=f, check=True, text=True)
	except subprocess.SubprocessError as e:
		print(f"acpype failed with error {e}\n. Ligand preparation failed")
		return LFtop

	acpypefolder = []
	for folder in os.listdir():
		if os.path.isdir(folder):
			acpypefolder.append(folder)
			os.chdir(folder)
			for Afile in os.listdir():
				shutil.copy(Afile, '../')
			os.chdir('../')
			shutil.rmtree(folder)

	if not len(acpypefolder) > 0:
		print("No acpype generated folder was detected. This process may fail")

	try:
		shutil.move(glob.glob("*_GMX_OPLS.top")[0], "oplstop.top")		
		shutil.move(glob.glob("*_GMX_OPLS.itp")[0], "oplstop.itp")		
	except:
		printWarning("Some files failed renaming. Trying again")
		try:
			shutil.move(glob.glob("*_GMX.top")[0], "oplstop.top")		
			shutil.move(glob.glob("*_GMX.itp")[0], "oplstop.itp")		
		except:
			print("Error persist with renaming. Trying alternative")
			return LFtop

	mergename = "opls" + name + ".top"
	mergeopls = open(mergename, "+a")
	oplstopopen = open("oplstop.top", "r")
	oplstopread = oplstopopen.readlines()
	oplsitpopen = open("oplstop.itp", "r")
	oplsitpread = oplsitpopen.readlines()
	
	q = 0
	for oline in oplstopread:
		if not oline == [] and "system" in oline.split():
			break
		else:
			q += 1

	for itpline in oplsitpread:
		mergeopls.write(itpline)
	mergeopls.write('\n')

	oplstopopen.seek(0)
	while q < len(oplstopread):
		mergeopls.write(oplstopread[q])
		q += 1

	oplstopopen.close()
	oplsitpopen.close()
	mergeopls.close()

	shutil.copy(mergename, '../')
	os.chdir('../')

	# Check to be sure opls compatible ligand topology has been successfully generated
	checktlp = os.listdir()
	mergename = "opls" + name + ".top"
	if not mergename in checktlp:
		printWarning("Something went wrong. Generating oplstopol.top file was not successful")
		return LFtop
	else:
		checkindex = indexoflines(mergename)
		a = int(checkindex['atomtypes'])
		b = int(checkindex['moleculetype'])
		c = int(checkindex['system'])
		checkopen = open(mergename, "r")
		checkreadlines = checkopen.readlines()
		if not ('atomtypes' in checkreadlines[a].split() and 'moleculetype' in checkreadlines[b].split() and 'system' in checkreadlines[c].split()):
			return LFtop
		else:
			printNote("opls compatible ligand topology has been generated successfully")
			return mergename
