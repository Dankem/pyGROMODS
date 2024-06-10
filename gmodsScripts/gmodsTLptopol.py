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
from pathlib import Path
import time
import shutil

from gmodsScripts.gmodsHelpers import insertdetails, indexoflines, printWarning, printNote, tinput, gmxtop

def TLtopol(xrecfile, tlpcomtop, tff):
	print('\n')
	tlpcwdir = Path.cwd()
	tltopol = "tlptopol.top"
	tltopolopen = open(tltopol, "+a")

	# Insert the header information into the tlptopol.top file
	tltopolheaders = ["; This file was generated by modifying tleap/acpype generated topology file on one hand", "; And by adding selected pdb2gmx generated details on the other hands", "; The modification was done by removing duplicate atomtypes and/or replacing", "; them with equivalent compatible atoms as found in ffnonbonded.itp file of the selected forcefield", "; It was created to be used as alternative to pdb2gmx generated topologies when necessary", "; If used successfully it should have been renamed from tlptopol.top to topol.top for subsequent MDS"]

	for header in tltopolheaders:
		tltopolopen.write(header)
		tltopolopen.write('\n')
	tltopolopen.write('\n')

	# Insert the forcefiled parameter header from pdb2gmx topol.top file into tlptopol.top file
	xrecindex = indexoflines(xrecfile)
	tsindex = xrecindex['system']
	xrecfileopen = open(xrecfile, "r")
	xrecreadlines = xrecfileopen.readlines()

	at_id = 0
	for xL in xrecreadlines:
		xLat = xL.split()
		if ('Include' in xLat and 'forcefield' in xLat and 'parameters' in xLat):
			tltopolopen.write(xL)
			at_id += 1
			tltopolopen.write(xrecreadlines[at_id])
			tltopolopen.write('\n')
			break
		else:
			at_id += 1

	# Insert the content of tleap/acpype generated topol file from atomtypes into tlptopol.top
	tlpcindex = indexoflines(tlpcomtop)
	tlpcopen = open(tlpcomtop, "r")
	tlpcreadlines = tlpcopen.readlines()
	taindex = tlpcindex['atomtypes']

	for tline in tlpcreadlines:
		if not tline in tlpcreadlines[0:int(taindex)]:
			tltopolopen.write(tline)

	tltopolopen.close()
	tlpcopen.close()

	# Create a new file and populate it with relevant details found at the end of pdb2gmx topol file
	tlmnfile = "tlptopol_mn.itp"
	tlmnfileopen = open(tlmnfile, "+a")

	xrecfileopen.seek(0)
	xn = 0
	for xline in xrecreadlines:
		xlinelist = xline.split()
		if "POSRES" in xlinelist:
			tlmnfileopen.write('\n')
			xnid = xn - 1
			while xnid < int(tsindex):
				tlmnfileopen.write(xrecreadlines[xnid])
				xnid += 1
			tlmnfileopen.write('\n')
			break
		else:
			xn += 1

	tlmnfileopen.close()
	xrecfileopen.close()

	# Insert the details into the tlptopol.top file
	insertL = "[ system ]"
	insertdetails(tltopol, tlmnfile, insertL)

	# We shall check for and remove duplicates in atomtypes between gmx standard and amber/tleap generated 
	# Determine the appropriate forcefield directory selected at run time

	tlpindex = indexoflines(tltopol)
	atlp = int(tlpindex['atomtypes'])
	mtlp = int(tlpindex['moleculetype'])

	topol = open(tltopol, "r")
	topolreadlines = topol.readlines()
	nf = 1
	ff = ""
	for tline in topolreadlines:
		if ('Include' in tline.split() and 'forcefield' in tline.split() and 'parameters' in tline.split()):
			ln1 = topolreadlines[nf].split()
			ln2 = ln1[1].split('"')[1].split('/')
			for fd in ln2:
				if Path(fd).suffix == ".ff":
					ff = fd
					break
		else:
			nf += 1

	if not ff == tff:
		print('\n')
		print(f"The detected {ff} does not match what was detected earlier {tff}")
		print(f"To use the recommended forcefield, '{tff}', type YES/y")
		print("To continue with currently detected forcefield {", ff, "}, press ENTER")
		response = tinput("Response: ", 30, "y")
		if (response.lower() == "yes" or response.lower() == "y"):
			ff = tff
			print(f"The forcefield directory has been changed to {ff}")
		
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
			print('\n')
			return tlpcomtop

		gmxtopff, topffdir = gmxtop()
		gmxtopdir = os.path.join(topffdir, ff)
		
		if not (len(gmxtopff) > 0 or Path(ff).stem in gmxtopff):
			print(f"The specified forcefield, {ff}, is missing in the selected directory")
			response = input("To continue anyway, type YES/y. Otherwise, press ENTER: ")
			if not (response.lower() == "yes" or response.lower() == "y"):
				print("Trying again ...")
				continue
			else:
				printWarning("Checked version of tlptopol cannot be generated")
				print("The platform will use the unchecked original file")
				time.sleep(5)
				print('\n')
				return tlpcomtop
			
		else:
			break

	lsgmxtopdir = os.listdir(gmxtopdir)
	for tp in lsgmxtopdir:
		if tp == "ffnonbonded.itp":
			shutil.copy(os.path.join(gmxtopdir, tp), './')

	# Get the atomtypes present in ffnonbonded.itp as list
	listffb = []
	fn = 0

	ffb = open("ffnonbonded.itp", "r")
	ffbreadlines = ffb.readlines()

	if ff[0:4].lower() == "opls":
		for fb in ffbreadlines:
			if 'opls_128' in fb.split():
				break
			else:
				fn += 1
		ffb.seek(0)

		for fb in ffbreadlines[fn:]:
			try:
				if not (fb.split() == [] or fb.split()[0] == "[" or fb.split()[0] == ";" or fb.split()[0][0] == '#' or fb.split()[0][0] == ';' or fb.split()[1] in listffb):
					if fb.split()[0].split('_')[0].lower() == 'opls':
						listffb.append(fb.split()[1])
			except IndexError:
				pass
	else:
		for fb in ffbreadlines[fn:]:
			try:
				if not (fb.split() == [] or fb.split()[0] == "[" or fb.split()[0] == ";" or fb.split()[0][0] == '#' or fb.split()[0][0] == ';' or fb.split()[0] in listffb):
					listffb.append(fb.split()[0])
			except IndexError:
				pass

	ffb.seek(0)

	# Check the tlptopol.top file to remove duplicate by comparing with ffnonbonded.itp list of atomtypes
	print("Checking for and removing duplicate atomtypes from tlptopol.top...")
	time.sleep(5)
	mtlp = mtlp - 1
	topol.seek(0)
	ntopol = open("ntlptopol.top", "+a")
	for aline in topolreadlines[0:atlp]:
		ntopol.write(aline)

	topol.seek(0)
	for bline in topolreadlines[atlp:mtlp]:
		try:
			if not (bline.split() == [] or bline.split()[0].lower() in listffb or bline.split()[0].upper() in listffb):
				ntopol.write(bline)
		except IndexError:
			pass
	ntopol.write('\n')

	topol.seek(0)

	if ff[0:4].lower() == "opls":
		btlp = 0
		for bd in topolreadlines:
			if 'bonds' in bd.split() and bd.split()[0] == '[' and bd.split()[1] == 'bonds':
				break
			else:
				btlp += 1

		topol.seek(0)
		for bdline in topolreadlines[mtlp:btlp]:
			if not (bdline.split() == [] or bdline.split()[0][0] == ';' or bdline.split()[0][0] == '['):
				matom = bdline.split()[1]
				ffb.seek(0)
				replace = " "
				bfound = 0
				for fline in ffbreadlines[fn:]:
					if fline.split()[1] == matom or fline.split()[1].lower() == matom or fline.split()[1] == matom.upper() and 'opls' in fline.split()[0].split('_'):
						replace = bdline.replace(matom, fline.split()[0], 1)
						bfound += 1
						break
				if not bfound > 0:
					ntopol.write(bdline)
				else:
					ntopol.write(replace)
			else:
				ntopol.write(bdline)
		ntopol.write('\n')

		topol.seek(0)
		for cline in topolreadlines[btlp:]:
			ntopol.write(cline)

	else:	
		for cline in topolreadlines[mtlp:]:
			ntopol.write(cline)

	ffb.close()
	topol.close()
	ntopol.close()

	# Now we shall backup the tlptopol.top and rename ntlptopol.top as tlptopol.top
	os.rename('tlptopol.top', '##tlptopol##')
	os.rename('ntlptopol.top', 'tlptopol.top')

	# Check to be sure tlptopol.top has been successfully generated
	print('\n')
	checktlp = os.listdir()
	if not "tlptopol.top" in checktlp:
		printWarning("Generating tlptopol.top file was not successful")
		print("If you need it later, check and correct any error, then rerun")
		return tlpcomtop
	else:
		printNote("tlptopol.top has been generated successfully")
		print("PLEASE NOTE:")
		print("**** Two topology files have been generated - topol.top and tlptopol.top")
		print("**** By default, topol.top will be used, while tlptopol serves as backup")
		time.sleep(5)

		return 'tlptopol.top'
