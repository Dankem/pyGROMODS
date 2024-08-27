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
import time
from pathlib import Path
import shutil

from gmodsScripts.gmodsHelpers import topolsplit, indexoflines, printWarning, printNote, gmxtop

def checktopsimilarities(mtsavedir):
	""" To run this script you must navigate to the directory containing topology files to be compared """
	print('\n')
	print("Similarity checks in progress ...")
	time.sleep(3)

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
			print(f"{ligname} shares similar Moleculetypes & Recidues with existing topology. Skipped")
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
			print(f"{ligname} [ moleculetype ] name is identitical to exisitng topology file. Skipped")
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
				print(f"{ligname} shares some residues with one or more ligands. Included, but you may wish to check")
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
			print("You have exceeded maximum trying attempts")
			printWarning("Checked version of ligand topology cannot be generated")
			print("The platform will use the unchecked original file")
			time.sleep(5)
			return ligtop

		gmxtopff, topffdir = gmxtop()
		gmxtopdir = os.path.join(topffdir, ff)
		
		if not (len(gmxtopff) > 0 or Path(ff).stem in gmxtopff):
			print(f"The specified forcefield, {ff}, is missing in the supplied directory")
			response = input("To continue anyway, type YES/y. Otherwise, press ENTER: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				print("Trying again ...")
				continue
			else:
				printWarning("Checked version of ligand topology can not be generated")
				print("The platform will use the unchecked original file")
				time.sleep(5)
				return ligtop
			
		else:
			break

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
		printWarning("Generating checked version of ligand topology was not successful")
		print("The platform will use the unchecked original file")
		time.sleep(5)
		return ligtop
	else:
		printNote(f"checked version of {ligtop} has been generated successfully")
		time.sleep(5)
		return 'chkLIG.top'
