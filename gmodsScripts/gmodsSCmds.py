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

from gmodsScripts.gmodsHelpers import printWarning, printNote, tinput, cleanup

scriptDIRa = os.path.abspath(os.path.dirname(sys.argv[0]))

def SCmds(appDIR, gmxDIR):
    print('\n')
    # Set some environment variables
    scriptDIR = appDIR
    GMX_MDS = gmxDIR

    print(f"User MDS working directory set to: {GMX_MDS}")

    # Set global activation variables to include opncl headers in work subdirectories
    opencl_activate = 0

    # Check existence of required folders and molecular dynamic parameter (.mdp) files
    mdshomedir = os.path.join(GMX_MDS, 'MDSHOME')
    if not os.path.isdir(mdshomedir):
        raise Exception("MDSHOME folder missing. Run MDS setup and populate the folder")

    opencldir = os.path.join(mdshomedir, 'opencl')
    if not os.path.isdir(opencldir):
        opencl_activate = 0
    else:
        openclfiles = os.listdir(opencldir)
        if len(openclfiles) == 0:
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
    print('\n')

    # Perform initial checks of the required files in MDSHOME directory
    mdsfiles = os.listdir(gmxmdsdir)

    if 'posre.itp' not in mdsfiles:
        printWarning("posre.itp is missing. If defined -DPOSRE in .mdp file, MDS will fail without it.")
        response = tinput("To continue without it, type YES/y or press ENTER to abort: ", 10, "y")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("posre.itp is missing, please upload from mdsetup page or follow manual setup")

    if 'posre_udp.itp' not in mdsfiles:
        printWarning("posre_udp.itp, a user generated restraint file, is missing. If defined -DPOSRE_UDP in .mdp file, MDS will fail without it.")
        response = tinput("To continue without it, type YES/y or press ENTER to abort: ", 10, "y")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("posre_udp.itp is missing, please upload from mdsetup page")

    if 'posre.itp' in mdsfiles and 'posre_udp.itp' in mdsfiles:
        printNote("Both posre.itp and posre_udp.itp restraint files are present")
        printNote("Only the one defined {-DPOSRE or -DPOSRE_UDP} in mdp file(s), will be used.")

    if 'topol.top' not in mdsfiles:
        printWarning("topol.top file is missing. You must include it in the working folder")
        print("If need be, you may also rename tlptopol.top to topol.top for use instead")
        raise Exception("topol.top file is missing. Please rerun and upload this file")

    if 'LIGS_at.itp' not in mdsfiles and 'LIGS_mt.itp' not in mdsfiles:
        printWarning("Important Ligand(s) topology file(s) are missing")
        printNote("Don't worry, if all molecule and atom types are captured in the topology file")
        response = tinput("To continue without these files, type YES/y or press ENTER to abort: ", 10, "y")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("Please upload necessary files from mdsetup page")

    if 'fsolvated.gro' not in mdsfiles:
        printWarning("The fsolvated.gro file is missing. Check and reupload")
        printNote("You may need to manually rename generated solvated.gro or ufsolvate.gro file to fsolvated.gro")
        raise Exception("The fsolvated.gro is missiing. Please upload necessary files from mdsetup page")

    print(f"MDSHOME folder has been created to host all MDS run. It location is {mdshomedir}")
    print(f"mdsgmx subfolder has been created inside MDSHOME. It location is {gmxmdsdir}")

    os.chdir(mdshomedir)

    print('\n')
	# Generate unique id number for the project
    printNote("mdsgmx subfolder will now be renamed with added unique ID reflecting the current MDS run")

    while True:
        idnumber = random.randint(0, 9)
        idalpha1 = random.choice(string.ascii_lowercase)
        idalpha2 = random.choice(string.ascii_lowercase)
        ID = str(idnumber) + idalpha1 + idalpha2
        newWDIR = 'mdsgmx_' + ID
        if os.path.isdir(newWDIR):
            continue
        else:
            os.rename(gmxmdsdir, newWDIR)
            break

    print(f"Your Unique ID for MDS run is: {ID}")

    # Change working folder name by adding the unique ID
    mdswork_dir = os.path.join(mdshomedir, newWDIR)
    print(f"Current workspace directory set to {mdswork_dir}")

    os.chdir(mdswork_dir)

    printNote("A new directory will be created for each level of MDS at each step in the workflow for maximum organization")

    # User can set marwarn for grompp if needed
    maxwarn = 0
    print("The default value for -maxwarn is 0. You may wish to increase this as needed")
    print("To change this value, type YES/y. Otherwise, press ENTER to continue with the default")
    response = tinput("Response: ", 30, "n")
    if response.lower() == "yes" or response.lower() == "y":
        while True:
            maxwarn = input("Change maxwarn to: ")
            if maxwarn == " " or maxwarn == "" or maxwarn.isdigit() == False:
                print("You must supply digit value for maxwarn to continue")
                continue
            else:
                print(f"maxwarn has been changed to: {maxwarn}")
                break

    # Copy mdp files into mdps sub-folder in the working directory
    try:
        os.mkdir("mdps")
    except FileExistsError:
        shutil.rmtree("mdps")
        os.mkdir("mdps")

    os.chdir("mdps")
    for file in MDPfiles:
        shutil.copy(os.path.join(fMDP, file), './')

    os.chdir('../')
    print('\n')

    #################################
    # ENERGY MINIMIZATION 1: STEEP  #
    #################################
    print("Executing: Minimization for Solvated Complex...")

    try:
        os.mkdir("EM")
    except FileExistsError:
        shutil.rmtree("EM")
        os.mkdir("EM")

    os.chdir("EM")
    fem = open("emerror.txt", "a")

    if opencl_activate > 0:
        oclfiles = os.listdir(opencldir)
        for ocl in oclfiles:
            shutil.copy(os.path.join(opencldir, ocl), './')

    # Iterative calls to grompp and mdrun to run the simulations
    cmdsd = ['gmx', 'grompp', '-f', '../mdps/minzsd.mdp', '-c', '../fsolvated.gro', '-p', '../topol.top', '-o', 'minzsd.tpr', '-maxwarn', str(maxwarn)]
    try:
        subprocess.run(cmdsd, stderr=subprocess.STDOUT, stdout=fem, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'minzsd'], stderr=subprocess.STDOUT, stdout=fem, check=True, text=True)
    except subprocess.SubprocessError as e:
        print(f"Steepest Descent Minimization encountered error {e}")

    if not ('minzsd.tpr' in os.listdir() or 'minzsd.gro' in os.listdir()):
        printWarning("Steepest Descent Minimization failed. Please check emerror.txt file")
    else:
        printNote("Steepest Descent Minimization Completed Successfully")

    cmdcg = ['gmx', 'grompp', '-f', '../mdps/minzcg.mdp', '-c', 'minzsd.gro', '-p', '../topol.top', '-o', 'minzcg.tpr', '-maxwarn', str(maxwarn)]
    try:
        subprocess.run(cmdcg, stderr=subprocess.STDOUT, stdout=fem, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'minzcg'], stderr=subprocess.STDOUT, stdout=fem, check=True, text=True)
    except subprocess.SubprocessError as e:
        print(f"Something went wrong with error {e}")
        printWarning("Generating needed files failed. Please check error message")

    if not ('minzcg.tpr' in os.listdir() or 'minzcg.gro' in os.listdir()):
        printWarning("Conjugate Gradient Minimization failed. Please check error message")
        printNote("It is advisable to abort the processs, check emerror.txt file")
        response = tinput("To abort type YES/y. Otherwise press ENTER to continue anyway ", 10, "y")
        if response.lower() == "yes" or response.lower() == "y":
            fem.close()
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Conjugate Gradient Minimization Completed Successfully")
        printNote("Minimization Phases Completed Successfully")

    fem.close()
    os.chdir('../')
    print('\n')

    #####################
    # NVT EQUILIBRATION #
    #####################
    print("Executing: Constant Volume Equilibration...")

    try:
        os.mkdir("NVT")
    except FileExistsError:
        shutil.rmtree("NVT")
        os.mkdir("NVT")

    os.chdir("NVT")
    fnvt = open("nvterror.txt", "a")
    if opencl_activate > 0:
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

    cmdnvt = ['gmx', 'grompp', '-f', '../mdps/equnvt.mdp', '-c', '../EM/minzcg.gro', '-r', '../EM/minzcg.gro', '-p', '../topol.top', '-o', 'equnvt.tpr', '-maxwarn', str(maxwarn)]

    try:
        subprocess.run(cmdnvt, stderr=subprocess.STDOUT, stdout=fnvt, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equnvt'], stderr=subprocess.STDOUT, stdout=fnvt, check=True, text=True)
    except subprocess.SubprocessError as e:
        print(e)
        printWarning("Something was not right with above error. Checking....")

    if not ('equnvt.tpr' in os.listdir() or 'equnvt.gro' in os.listdir()):
        printWarning("There were errors during Constant Volume Equilibration Phase")
        printNote("It is advisable to abort, and check nvterror.txt file for details")
        response = tinput("To abort type YES/y. Otherwise press ENTER to continue anyway: ", 10, "y")
        if response.lower() == "yes" or response.lower() == "y":
            fnvt.close()
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Constant Volume Equilibration Successfully Completed")
        fnvt.close()

    os.chdir('../')
    time.sleep(5)
    print('\n')

    #####################
    # NPT EQUILIBRATION #
    #####################
    print("Executing: Constant Pressure Equilibration...")

    try:
        os.mkdir("NPT")
    except FileExistsError:
        shutil.rmtree("NPT")
        os.mkdir("NPT")

    os.chdir("NPT")
    fnpt = open("npterror.txt", "a")
    if opencl_activate > 0:
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

    cmdnpt = ['gmx', 'grompp', '-f', '../mdps/equnpt.mdp', '-c', '../NVT/equnvt.gro', '-r', '../NVT/equnvt.gro', '-p', '../topol.top', '-t', '../NVT/equnvt.cpt', '-o', 'equnpt.tpr', '-maxwarn', str(maxwarn)]
    try:
        subprocess.run(cmdnpt, stderr=subprocess.STDOUT, stdout=fnpt, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equnpt'], stderr=subprocess.STDOUT, stdout=fnpt, check=True, text=True)
    except subprocess.SubprocessError as e:
        printWarning("Something was not right. Checking...")

    if not ('equnpt.tpr' in os.listdir() or 'equnpt.gro' in os.listdir()):
        printWarning("There were errors during Constant Pressure Equilibration Phase")
        printNote("It is advisable to abort. Please check npterror.txt file for details")
        response = tinput("To abort type YES/y. Otherwise press ENTER to continue anyway: ", 10, "y")
        if response.lower() == "yes" or response.lower() == "y":
            fnpt.close()
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Constant Pressure Equilibration Successfully Completed")
        fnpt.close()

    os.chdir('../')
    time.sleep(5)
    print('\n')

    ####################
    # EQ PRODUCTION MD #
    ####################
    print("Executing: Equilibration Production MD simulation...")

    try:
        os.mkdir("EMD")
    except FileExistsError:
        shutil.rmtree("EMD")
        os.mkdir("EMD")

    os.chdir("EMD")
    feq = open("eqerror.txt", "a")
    if opencl_activate > 0:
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

    cmdeq = ['gmx', 'grompp', '-f', '../mdps/equpmd.mdp', '-c', '../NPT/equnpt.gro', '-r', '../NPT/equnpt.gro', '-p', '../topol.top', '-t', '../NPT/equnpt.cpt', '-o', 'equpmd.tpr', '-maxwarn', str(maxwarn)]
    try:
        subprocess.run(cmdeq, stderr=subprocess.STDOUT, stdout=feq, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equpmd'], stderr=subprocess.STDOUT, stdout=feq, check=True, text=True)
    except subprocess.SubprocessError as e:
        printWarning("Something was not right. Checking...")

    if not ('equpmd.tpr' in os.listdir() or 'equpmd.gro' in os.listdir()):
        printWarning("There were errors during Equilibration Production MD Phase")
        printNote("It is advisable to abort, and check eqerror.txt file for details")
        response = tinput("To abort type YES/y. Otherwise press ENTER to continue anyway ", 10, "y")
        if response.lower() == "yes" or response.lower() == "y":
            feq.close()
            raise Exception("Process Aborted. Make necessary corrections and restart")
    else:
        printNote("Equilibration Production MD Successfully Completed")
        feq.close()

    os.chdir('../')
    time.sleep(5)
    print('\n')

    #################
    # PRODUCTION MD #
    #################
    print("Executing: Production MD simulation...")

    try:
        os.mkdir("PMD")
    except FileExistsError:
        shutil.rmtree("PMD")
        os.mkdir("PMD")

    os.chdir("PMD")
    fpmd = open("pmderror.txt", "a")
    if opencl_activate > 0:
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

    cmdpmd = ['gmx', 'grompp', '-f', '../mdps/pmds.mdp', '-c', '../EMD/equpmd.gro', '-r', '../EMD/equpmd.gro', '-p', '../topol.top', '-t', '../EMD/equpmd.cpt', '-o', 'pmds.tpr', '-maxwarn', str(maxwarn)]

    try:
        subprocess.run(cmdpmd, stderr=subprocess.STDOUT, stdout=fpmd, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'pmds'], stderr=subprocess.STDOUT, stdout=fpmd, check=True, text=True)
    except subprocess.SubprocessError as e:
        printWarning("Something was not right. Checking...")

    if not ('pmds.tpr' in os.listdir() or 'pmds.gro' in os.listdir()):
        printWarning("There were errors during Production MD Phase. Check pmderror.txt file for details")
    else:
        printNote("Production MD Successfully Completed")

    fpmd.close()
    os.chdir('../')
    time.sleep(5)
    print('\n')

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

    time.sleep(5)

    # End
    print(f"Ending. Job Completed Successfully for MDS codenamed {newWDIR}")
