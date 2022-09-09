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

from gmodsScripts.gmodsHelpers import topolsplit, indexoflines, printWarning, printNote, select_folder, gmxtop

def checktopsimilarities(mtsavedir):
	""" To run this script you must navigate to the directory containing topology files to be compared """
	printNote("PLEASE TAKE NOTE OF THE FOLLOWING:")
	printNote("Only one of the Ligand(s) that shares similar or identical [ moleculetype ] name will be used while generating molecule type file. Others will be skipped. Where this occurs, please crosscheck the identified structure file(s), and if you want to include the file(s), ensure they have unique molecule name and rerun")

	time.sleep(10)

	# Get the list of files
	topoldir = os.listdir()

	file_mn = open("LIGS_mn.itp", "+a")

	ln = 1
	for file1 in topoldir:
		upload_name = "LIG" + str(ln)
		ligname = upload_name + ".pdb"
		f1_index = indexoflines(file1)
		x = f1_index['atomtypes']
		y = f1_index['moleculetype']
		z = f1_index['system']

		Aidentical = 0
		Bidentical = 0

		open1 = open(file1, "r")
		name1 = open1.readlines()

		for file2 in topoldir[0:ln]:
			f2_index = indexoflines(file2)
			a = f2_index['atomtypes']
			b = f2_index['moleculetype']
			c = f2_index['system']

			open2 = open(file2, "r")
			name2 = open2.readlines()

			if name1[y+2].split() == name2[b+2].split():
				Aidentical += 1

			try:
				if name1[z+1].split() == name2[c+1].split():
					Aidentical += 1
			except IndexError:
				if name1[len(name1)-1].split() == name2[len(name2)-1].split():
					Aidentical += 1

			count = 0
			if len(name1) > len(name2):
				for n in range(len(name2)):
					if name1[n].split() == name2[n].split():
						count += 1
			else:
				for n in range(len(name1)):
					if name1[n].split() == name2[n].split():
						count += 1

			percent_similarity = count * 100 / len(name1)
			if percent_similarity > 21:
				Bidentical += 1

			open2.close()

		if Aidentical > 2 and Bidentical > 1:
			print(ligname, "shares similar [ moleculetype ] name and some residues with a ligand already included in topology file. This file will be skipped")
			time.sleep(5)
			lastline = len(name1) - 1
			while lastline < len(name1):
				if name1[lastline].split() == "":
					lastline += -1
					continue
				else:
					file_mn.write(name1[lastline])
					break

			open1.close()
			ln += 1
			continue

		elif Aidentical > 2 and Bidentical <= 1:
			print(ligname, "[ moleculetype ] name is identitical with a ligand already included in topology file. This file will be skipped")
			time.sleep(5)
			lastline = len(name1) - 1
			while lastline < len(name1):
				if name1[lastline].split() == "":
					lastline += -1
					continue
				else:
					file_mn.write(name1[lastline])
					break

			open1.close()
			ln += 1
			continue

		else:
			if Bidentical > 1 and Aidentical <= 2:
				print(ligname, "shares some residues with one or more ligands. This may not affect your work. You may however wish to check to be doubling sure")
				time.sleep(5)

			lastline = len(name1) - 1
			while lastline < len(name1):
				if name1[lastline].split() == "":
					lastline += -1
					continue
				else:
					file_mn.write(name1[lastline])
					break

			ligand_at, ligand_mt = topolsplit(file1, upload_name, f1_index)

			shutil.move(ligand_mt, mtsavedir)

		open1.close()
		ln += 1

	file_mn.close()
	return 'LIGS_mn.itp'


def Checkligtop(ligtop, ff):
	# Set some needed variables
	tlpindex = indexoflines(ligtop)
	atlp = int(tlpindex['atomtypes'])
	mtlp = int(tlpindex['moleculetype'])

	topol = open(ligtop, "r")
	topolreadlines = topol.readlines()

	# Get the absolute path to the forcefield directory and copy ffnonbonded.itp file
	gmxtopdir = " "
	topffdir = " "
	gmxtopff = []

	nT = 0
	while True:
		nT += 1
		if nT > 3:
			break

		gmxtopff, topffdir = gmxtop()
		gmxtopdir = os.path.join(topffdir, ff)
		
		if not len(gmxtopff) > 0:
			print(f"The specified directory, {gmxtopdir}, is not a valid Gromacs forcefield directory")
			printNote("Trying again ...")
			continue
			
		elif not Path(ff).stem in gmxtopff:
			print(f"The specified forcefield, {ff}, is missing in the selected/detected directory")
			printNote("You might have selected a directory with an incomplete list of forcefields")
			printNote("Trying again ...")
			continue
			
		elif not os.path.isdir(gmxtopdir):
			print(gmxtopdir, "that was autodetected, is not a valid forcefield directory")
			printNote("Trying again ...")
			continue
			
		else:
			print("Your topology directory is", gmxtopdir)
			break

	time.sleep(5)

	lsgmxtopdir = os.listdir(gmxtopdir)
	for tp in lsgmxtopdir:
		if tp == "ffnonbonded.itp":
			shutil.copy(os.path.join(gmxtopdir, tp), './')

	# Get the atomtypes present in ffnonbonded.itp as list
	listffb = []

	ffb = open("ffnonbonded.itp", "r")
	ffbreadlines = ffb.readlines()

	for fb in ffbreadlines:
		try:
			if not (fb.split() == [] or fb.split()[0] == "[" or fb.split()[0] == ";" or fb.split()[0][0] == '#' or fb.split()[0][0] == ';' or fb.split()[0] in listffb):
				listffb.append(fb.split()[0])
		except IndexError:
			pass

	ffb.seek(0)

	# Check the ligand topology file (generated or uploaded) to remove duplicate atomtypes
	print("Checking for and removing duplicate atomtypes from ligand topology ...")
	time.sleep(5)
	mtlp = mtlp - 1
	topol.seek(0)
	ntopol = open("chkLIG.top", "+a")
	for aline in topolreadlines[0:atlp]:
		ntopol.write(aline)

	topol.seek(0)
	for bline in topolreadlines[atlp:mtlp]:
		try:
			if not (bline.split() == [] or bline.split()[0] in listffb):
				ntopol.write(bline)
		except IndexError:
			pass
	ntopol.write('\n')

	topol.seek(0)

	for cline in topolreadlines[mtlp:]:
		ntopol.write(cline)

	ffb.close()
	topol.close()
	ntopol.close()

	# Check to be sure new ligand topology has been successfully generated
	checktlp = os.listdir()
	if not "chkLIG.top" in checktlp:
		printWarning("Something went wrong. Generating checked version of ligand topology was not successful")
		print("The platform will use the unchecked original file")
		time.sleep(5)
		return ligtop
	else:
		printNote("checked version of ligand topology has been generated successfully")
		time.sleep(5)
		return 'chkLIG.top'
