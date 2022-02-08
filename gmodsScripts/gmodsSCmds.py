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

from gmodsScripts.gmodsHelpers import printWarning, printNote, cleanup

def SCmds():
    # Set global activation variables to include opncl headers in work subdirectories
    opencl_activate = 0

    # Set some environment variables
    scriptDIR = os.path.abspath(os.path.dirname(sys.argv[0]))

    GMX_MDS=Path.cwd()
    print("Current Working Directory set to ", GMX_MDS)

    # Check existence of required folders and moledular dynamic parameter (.mdp) files
    mdshomedir = os.path.join(GMX_MDS, 'MDSHOME')
    if not os.path.isdir(mdshomedir):
        raise Exception("MDSHOME folder missing. Run MDS setup and populate the folder")

    opencldir = os.path.join(mdshomedir, 'opencl')
    if not os.path.isdir(opencldir):
        opencl_activate = 0
    else:
        openclfiles = os.listdir(opencldir)
        if openclfiles == 0:
            opencl_activate = 0
        else:
            opencl_activate += 1

    gmxmdsdir = os.path.join(mdshomedir, 'mdsgmx')
    if not os.path.isdir(gmxmdsdir):
        raise Exception("mdsgmx folder missing. Run MDS setup and populate the folder")

    fMDP = os.path.join(scriptDIR, 'MDP')
    if not os.path.isdir(fMDP):
        raise Exception("MDP folder missing. Run setup and populate the folder")

    mdpfilename = ['minzsd.mdp', 'minzcg.mdp', 'equnvt.mdp', 'equnpt.mdp', 'equpmd.mdp', 'pmds.mdp']
    MDPfiles = os.listdir(fMDP)
    if len(MDPfiles) == 0:
        printWarning("No .mdp file found. Please make sure all files are correctly uploaded")
        printNote("Samples of .mdp files are contain in the sample folder of your Selected Route generated work directory. You can modify these files as needed")
        raise Exception("Needed mdp file(s) missing. Please upload needed files")

    nmdp = 0
    for file in MDPfiles:
        if not Path(file).suffix == ".mdp":
            nmdp += 1

    if not nmdp == 0:
        raise Exception(nmdp, "mdp file(s) lack the required .mdp extension. Please check and rename accordingly")

    nmdpn = 0
    for file in MDPfiles:
        if file not in mdpfilename:
            nmdpn += 1

    if not nmdpn == 0:
        raise Exception(nmdpn, "mdp file(s) missing. Please make sure all files are correctly uploaded")

    print("found needed .mdp files")

    # Perform initial checks of the required files in MDSHOME directory
    mdsfiles = os.listdir(gmxmdsdir)

    if 'posre.itp' not in mdsfiles:
        printWarning("posre.itp is missing. If defined -DPOSRE in .mdp file, MDS will fail without it.")
        response = input("To continue without it, type YES/y or press ENTER to abort: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("posre.itp is missing, please upload from mdsetup page or follow manual setup")

    if 'posre_udp.itp' not in mdsfiles:
        printWarning("posre_udp.itp, a user generated restraint file, is missing. If defined -DPOSRE_UDP in .mdp file, MDS will fail without it.")
        response = input("To continue without it, type YES/y or press ENTER to abort: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("posre_udp.itp is missing, please upload from mdsetup page")

    if 'posre.itp' in mdsfiles and 'posre_udp.itp' in mdsfiles:
        printNote("Both posre.itp and posre_udp.itp restraint files are present")
        printNote("If defined -DPOSRE in mdp file(s), posre.itp will be used")
        printNote("To use posre_udp.itp instead, type YES/y to backup posre.itp and rename posre_udp.itp to posre.itp")
        response = input("Response YES/y, otherwise press ENTER to continue: ")
        if response.lower() == "yes" or response.lower() == "y":
            os.rename(os.path.join(gmxmdsdir, 'posre.itp'), os.path.join(gmxmdsdir, 'posre.bk'))
            os.rename(os.path.join(gmxmdsdir,'posre_udp.itp'), os.path.join(gmxmdsdir, 'posre.itp'))
            print("\n")

        else:
            printNote("If defined -DPOSRE_UDP, posre_udp.itp will be used")
            printNote("To use posre.itp instead, type YES/y to backup posre_udp.itp and rename posre.itp to posre_udp.itp")
            response = input("Response YES/y, otherwise press ENTER to continue: ")
            if response.lower() == "yes" or response.lower() == "y":
                os.rename(os.path.join(gmxmdsdir, 'posre_udp.itp'), os.path.join(gmxmdsdir, 'posre_udp.bk'))
                os.rename(os.path.join(gmxmdsdir,'posre.itp'), os.path.join(gmxmdsdir, 'posre_udp.itp'))

    if 'topol.top' not in mdsfiles:
        printWarning("topol.top file is missing")
        printNote("Please include the file in your working folder - mdsgmx, or if need be rename your generated top file to topol.top and type YES/y to continue or press ENTER to abort")
        response = input("Response: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("Please upload necessary files from mdsetup page")

    if 'LIGS_at.itp' not in mdsfiles and 'LIGS_mt.itp' not in mdsfiles:
        printWarning("Incomplete files. Important Ligand(s) topology file(s) are missing")
        printNote("This may not affect the MDS if using tlptopol.top instead of topol.top. In that case tlptopol.top would have been renamed to topol.top. Or if you are working with protein and peptides without ligands")
        response = input("To continue without these files, type YES/y or press ENTER to abort: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("Please upload necessary files from mdsetup page")

    if 'fsolvated.gro' not in mdsfiles:
        printWarning("Incomplete files. The fsolvated.gro file is missing. Check the following:")
        printNote("1). You may need to rename manually generated solvated.gro file to fsolvated.gro")
        printNote("2). OR, if genrated, rename 'ufsolvate.gro' file to fsolvated.gro")
        printNote("Now, reupload your files and try again")
        raise Exception("The fsolvated.gro is missiing. Please upload necessary files from mdsetup page")

    print("MDSHOME folder has been created to host all MDS run. It location is ", mdshomedir)
    print("mdsgmx subfolder has been created inside MDSHOME. It location is ", gmxmdsdir)

    time.sleep(10)

    os.chdir(mdshomedir)

	# Generate unique id number for the project
    printNote("mdsgmx subfolder will now be renamed with added unique ID reflecting the current MDS run")
    printNote("PLEASE NOTE: Any former MDS run - failed or succeed - can be assessed either as backup folder or as renamed folder inside MDSHOME")

    while True:
        idnumber = random.randint(0, 9)
        idalpha1 = random.choice(string.ascii_letters)
        idalpha2 = random.choice(string.ascii_letters)
        ID = str(idnumber) + idalpha1 + idalpha2
        newWDIR = 'mdsgmx_' + ID
        if os.path.isdir(newWDIR):
            continue
        else:
            os.rename(gmxmdsdir, newWDIR)
            break

    print("Your Unique ID for MDS run is: ", ID)

    # Change working folder name by adding the unique ID
    mdswork_dir = os.path.join(mdshomedir, newWDIR)
    print("Current workspace directory set to ", mdswork_dir)

    os.chdir(mdswork_dir)

    printNote("A new directory will be created for each level of MDS at each step in the workflow for maximum organization")

    time.sleep(5)

    # User can set marwarn for grompp if needed
    maxwarn = 0
    print("The default value for -maxwarn is 0. You may wish to increase this as needed")
    print("To change this value, type YES/y. Otherwise, press ENTER to continue with the default")
    response = input("Response: ")
    if response.lower() == "yes" or response.lower() == "y":
        while True:
            maxwarn = input("Change maxwarn to: ")
            if maxwarn == " " or maxwarn == "" or maxwarn.isdigit() == False:
                print("You must supply digit value for maxwarn to continue")
                continue
            else:
                print("maxwarn has been changed to: ", maxwarn)
                break

    #################################
    # ENERGY MINIMIZATION 1: STEEP  #
    #################################
    print("Starting Minimization for Solvated Complex...")

    try:
        os.mkdir("mdps")
    except FileExistsError:
        shutil.rmtree("mdps")
        os.mkdir("mdps")

    os.chdir("mdps")
    for file in MDPfiles:
        shutil.copy(os.path.join(fMDP, file), './')

    os.chdir('../')

    try:
        os.mkdir("EM")
    except FileExistsError:
        shutil.rmtree("EM")
        os.mkdir("EM")

    os.chdir("EM")
    if opencl_activate > 0:
        oclfiles = os.listdir(opencldir)
        for ocl in oclfiles:
            shutil.copy(os.path.join(opencldir, ocl), './')

    # Iterative calls to grompp and mdrun to run the simulations
    try:
        subprocess.run('gmx grompp -f ../mdps/minzsd.mdp -c ../fsolvated.gro -p ../topol.top -o minzsd.tpr -maxwarn ' + str(maxwarn), shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Something went wrong with above error. Trying again...")		
        Err1a = os.system('gmx grompp -f ../mdps/minzsd.mdp -c ../fsolvated.gro -p ../topol.top -o minzsd.tpr -maxwarn ' + str(maxwarn))
        if not Err1a == 0:
            print("Generating needed files failed. Please check error message")

    try:
        subprocess.run('gmx mdrun -deffnm minzsd', shell=True)
    except subprocess.CalledProcessError as e:
        Err1b = os.system('gmx mdrun -deffnm minzsd')
        if not Err1b == 0:
            print("Steepest Descent Minimization failed. Please check error message")

    if not ('minzsd.tpr' in os.listdir() or 'minzsd.gro' in os.listdir()):
        printWarning("Steepest Descent Minimization failed. Please check error message")
    else:
        printNote("Steepest Descent Minimization Completed Siccessfully")

    time.sleep(5)

    try:
        subprocess.run('gmx grompp -f ../mdps/minzcg.mdp -c minzsd.gro -p ../topol.top -o minzcg.tpr -maxwarn ' + str(maxwarn), shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Something went wrong with above error. Trying again...")		
        Err1c = os.system('gmx grompp -f ../mdps/minzcg.mdp -c minzsd.gro -p ../topol.top -o minzcg.tpr -maxwarn ' + str(maxwarn))
        if not Err1c == 0:
            print("Generating needed files failed. Please check error message")

    try:
        subprocess.run('gmx mdrun -deffnm minzcg', shell=True)
    except subprocess.CalledProcessError as e:
        Err1d = os.system('gmx mdrun -deffnm minzcg')
        if not Err1d == 0:
            printWarning("Something was not right")

    if not ('minzcg.tpr' in os.listdir() or 'minzcg.gro' in os.listdir()):
        printWarning("Conjugate Gradient Minimization failed. Please check error message")
        printWarning("There were errors during Minimization Phases")
        printNote("It is advisable to abort the processs, check your files and errors and rerun")
        response = input("To abort type YES/y. Otherwise process will continue anyway ")
        if response.lower() == "yes" or response.lower() == "y":
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Conjugate Gradient Minimization Completed Siccessfully")
        printNote("Minimization Phases Completed Successfully")

    os.chdir('../')
    time.sleep(10)

    #####################
    # NVT EQUILIBRATION #
    #####################
    print("Starting Constant Volume Equilibration...")

    try:
        os.mkdir("NVT")
    except FileExistsError:
        shutil.rmtree("NVT")
        os.mkdir("NVT")

    os.chdir("NVT")
    if opencl_activate > 0:
        oclfiles = os.listdir(opencldir)
        for ocl in oclfiles:
            shutil.copy(os.path.join(opencldir, ocl), './')

    try:
        subprocess.run('gmx grompp -f ../mdps/equnvt.mdp -c ../EM/minzcg.gro -r ../EM/minzcg.gro -p ../topol.top -o equnvt.tpr -maxwarn ' + str(maxwarn), shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Something went wrong with above error. Trying again...")		
        Err2a = os.system('gmx grompp -f ../mdps/equnvt.mdp -c ../EM/minzcg.gro -r ../EM/minzcg.gro -p ../topol.top -o equnvt.tpr -maxwarn ' + str(maxwarn))
        if not Err2a == 0:
            print("Generating needed files for NVT failed. Please check error message")

    try:
        subprocess.run('gmx mdrun -deffnm equnvt', shell=True)
    except subprocess.CalledProcessError as e:
        Err2b = os.system('gmx mdrun -deffnm equnvt')
        if not Err2b == 0:
            printWarning("Something was not right")

    if not ('equnvt.tpr' in os.listdir() or 'equnvt.gro' in os.listdir()):
        printWarning("There were errors during Constant Volume Equilibration Phase")
        printNote("It is advisable to abort the processs, check your files and errors and rerun")
        response = input("To abort type YES/y. Otherwise process will continue anyway ")
        if response.lower() == "yes" or response.lower() == "y":
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Constant Volume Equilibration Successfully Completed")

    os.chdir('../')
    time.sleep(10)

    #####################
    # NPT EQUILIBRATION #
    #####################
    print("Starting Constant Pressure Equilibration...")

    try:
        os.mkdir("NPT")
    except FileExistsError:
        shutil.rmtree("NPT")
        os.mkdir("NPT")

    os.chdir("NPT")
    if opencl_activate > 0:
        oclfiles = os.listdir(opencldir)
        for ocl in oclfiles:
            shutil.copy(os.path.join(opencldir, ocl), './')

    try:
        subprocess.run('gmx grompp -f ../mdps/equnpt.mdp -c ../NVT/equnvt.gro -r ../NVT/equnvt.gro -p ../topol.top -t ../NVT/equnvt.cpt -o equnpt.tpr -maxwarn ' + str(maxwarn), shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Something went wrong with above error. Trying again...")		
        Err3a = os.system('gmx grompp -f ../mdps/equnpt.mdp -c ../NVT/equnvt.gro -r ../NVT/equnvt.gro -p ../topol.top -t ../NVT/equnvt.cpt -o equnpt.tpr -maxwarn ' + str(maxwarn))
        if not Err3a == 0:
            print("Generating needed files for NPT failed. Please check error message")

    try:
        subprocess.run('gmx mdrun -deffnm equnpt', shell=True)
    except subprocess.CalledProcessError as e:
        Err3b = os.system('gmx mdrun -deffnm equnpt')
        if not Err3b == 0:
            printWarning("Something was not right")

    if not ('equnpt.tpr' in os.listdir() or 'equnpt.gro' in os.listdir()):
        printWarning("There were errors during Constant Pressure Equilibration Phase")
        printNote("It is advisable to abort the processs, check your files and errors and rerun")
        response = input("To abort type YES/y. Otherwise process will continue anyway ")
        if response.lower() == "yes" or response.lower() == "y":
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Constant Pressure Equilibration Successfully Completed")

    os.chdir('../')
    time.sleep(10)

    ####################
    # EQ PRODUCTION MD #
    ####################
    print("Starting Equilibration Production MD simulation...")

    try:
        os.mkdir("EMD")
    except FileExistsError:
        shutil.rmtree("EMD")
        os.mkdir("EMD")

    os.chdir("EMD")
    if opencl_activate > 0:
        oclfiles = os.listdir(opencldir)
        for ocl in oclfiles:
            shutil.copy(os.path.join(opencldir, ocl), './')

    try:
        subprocess.run('gmx grompp -f ../mdps/equpmd.mdp -c ../NPT/equnpt.gro -r ../NPT/equnpt.gro -p ../topol.top -t ../NPT/equnpt.cpt -o equpmd.tpr -maxwarn ' + str(maxwarn), shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Something went wrong with above error. Trying again...")		
        Err4a = os.system('gmx grompp -f ../mdps/equpmd.mdp -c ../NPT/equnpt.gro -r ../NPT/equnpt.gro -p ../topol.top -t ../NPT/equnpt.cpt -o equpmd.tpr -maxwarn ' + str(maxwarn))
        if not Err4a == 0:
            print("Generating needed files for EQ Production MD failed. Please check error message")

    try:
        subprocess.run('gmx mdrun -deffnm equpmd', shell=True)
    except subprocess.CalledProcessError as e:
        Err4b = os.system('gmx mdrun -deffnm equpmd')
        if not Err4b == 0:
            printWarning("Something was not right")

    if not ('equpmd.tpr' in os.listdir() or 'equpmd.gro' in os.listdir()):
        printWarning("There were errors during Equilibration Production MD Phase")
        printNote("It is advisable to abort the processs, check your files and errors and rerun")
        response = input("To abort type YES/y. Otherwise process will continue anyway ")
        if response.lower() == "yes" or response.lower() == "y":
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Equilibration Production MD Successfully Completed")

    os.chdir('../')
    time.sleep(10)

    #################
    # PRODUCTION MD #
    #################
    print("Starting Production MD simulation...")

    try:
        os.mkdir("PMD")
    except FileExistsError:
        shutil.rmtree("PMD")
        os.mkdir("PMD")

    os.chdir("PMD")
    if opencl_activate > 0:
        oclfiles = os.listdir(opencldir)
        for ocl in oclfiles:
            shutil.copy(os.path.join(opencldir, ocl), './')

    try:
        subprocess.run('gmx grompp -f ../mdps/pmds.mdp -c ../EMD/equpmd.gro -r ../EMD/equpmd.gro -p ../topol.top -t ../EMD/equpmd.cpt -o pmds.tpr -maxwarn ' + str(maxwarn), shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Something went wrong with above error. Trying again...")		
        Err5a = os.system('gmx grompp -f ../mdps/pmds.mdp -c ../EMD/equpmd.gro -r ../EMD/equpmd.gro -p ../topol.top -t ../EMD/equpmd.cpt -o pmds.tpr -maxwarn ' + str(maxwarn))
        if not Err4a == 0:
            print("Generating needed files for Production MD failed. Please check error message")

    try:
        subprocess.run('gmx mdrun -deffnm pmds', shell=True)
    except subprocess.CalledProcessError as e:
        Err5b = os.system('gmx mdrun -deffnm pmds')
        if not Err5b == 0:
            printWarning("Something was not right")

    if not ('pmds.tpr' in os.listdir() or 'pmds.gro' in os.listdir()):
        printWarning("There were errors during Production MD Phase")
        printNote("Please check your files and errors and be sure all is in order. Otherwise rerun")
    else:
        printNote("Production MD Successfully Completed")

    os.chdir('../')
    time.sleep(5)

	##################
    # CLEAN UP FILES #
    ##################
    if opencl_activate > 0:
        printNote("Cleaning up the opencl headers...")
        time.sleep(5)

        os.chdir("EM")
        cleanup()

        os.chdir("../NVT")
        cleanup()

        os.chdir("../NPT")
        cleanup()

        os.chdir("../EMD")
        cleanup()

        os.chdir("../PMD")
        cleanup()

        os.chdir('../../')
        try:
            shutil.rmtree(opencldir)
        except:
            pass

        os.chdir('../')

    else:
        os.chdir('../../')

    time.sleep(10)

    # End
    print("Ending. Job Completed Successfully for MDS codenamed", newWDIR)
