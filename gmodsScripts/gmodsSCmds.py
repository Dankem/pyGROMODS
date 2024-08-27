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
import random
import string
from datetime import datetime
from tkinter import Tk, messagebox

from gmodsScripts.gmodsHelpers import printWarning, printNote, tinput, cleanup

scriptDIRa = os.path.abspath(os.path.dirname(sys.argv[0]))

def SCmds(appDIR, gmxDIR):
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

    # Send out popup message to capture user's attention to a need for final setup
    messages = "LET'S GET THINGS READY FOR SMOOTH MDS RUN. YOU NEED TO ANSWER SOME QUESTIONS ON THE TERMINAL"

    crp = Tk()
    crp.title("MDS FINAL SETUP ATTENTION!!!")
    crp.withdraw()
    crp.attributes("-topmost", True)
    crp.geometry("250x750")

    checkmessage = messagebox.askquestion("MDS SETUP MESSAGE:", f"{messages} \n\n ARE YOU READY FOR FINAL SETUP FOR MDS?")
    crp.destroy()

    if not (checkmessage.lower() == "yes" or checkmessage.lower() == "y"):
        raise Exception("Operation Interrupted by User. Make necessary corrections and restart")
    else:
        printNote("Please respond to questions as appropriate")
    time.sleep(5)
    print('\n')

    # Checking needed files for correctness 
    mdsfiles = os.listdir(gmxmdsdir)
    udp = 0

    if not ('posre.itp' in mdsfiles and 'posre_udp.itp' in mdsfiles):
        printNote("No restraint file was found. If defined -DP0SRE OR -DPOSRE_UDP in .mdp file, you must include posre.itp or posre_udp.itp respectively in the working folder")
        response = input("Type YES/y to continue without retraint file or press ENTER to abort")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("Operation Interrupted by User. Make necessary corrections and restart")
        else:
            print("User choose to continue without restriant file")
            print('\n')

    if 'posre_udp.itp' in mdsfiles and 'posre.itp' in mdsfiles:
        printNote("Two restraint files were found. Only the one defined in .mdp file will be used")
        print("If defined -DPOSRES in .mdp file, posre.itp will be used")
        print("If defined -DPOSRES_UDP in .mdp file, posre_udp.itp will be used")        
        print("However, assuming default -DPOSRES in .mdp file, You can choose to use posre_udp.itp instead of posre.itp")
        response = input("Type YES/y to use posre_udp.itp or press ENTER to continue with posre.itp: ")
        if response.lower() == "yes" or response.lower() == "y":
            print("Backing up posre.itp and renaming posre_udp.itp ...")
            os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'), os.path.join(mdshomedir, 'mdsgmx', '#posre.itp#'))
            os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre_udp.itp'), os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'))
            print("User defined restriant file will be used")
            print("Please note that any associated index file must be included as well")
            udp += 1
        else:
            print("User choose to continue with posre.itp file")
        print('\n')

    if not 'posre_udp.itp' in mdsfiles and 'posre.itp' in mdsfiles:
        print("We found the default posre.itp file and will be used")
        print('\n')

    if 'posre_udp.itp' in mdsfiles and not 'posre.itp' in mdsfiles:
        print("Only user defined restriant file was found. To use it, you must define -DPOSRES_UDP in .mdp file")
        print("However, if already defined -DPOSRE in .mdp file, the file must be renamed to posre.itp")
        response = input("Type YES to rename. Otherwise, press ENTER to continue: ")
        if response.lower() == "yes" or response.lower() == "y":
            os.rename('posre_udp.itp', 'posre.itp')
            print("User defined restriant file will be used")
            print("Please note that any associated index file must be included as well")
            udp += 1
        print('\n')

    if 'topol.top' not in mdsfiles and 'tlptopol.top' not in mdsfiles:
        printWarning("topology file is missing. You must include either topol.top or tlptopol.top in the working folder")
        raise Exception("topol.top file is missing. Please rerun and upload this file")

    if 'topol.top' in mdsfiles and not 'tlptopol.top' in mdsfiles:
        print("We found only the default topol.top and will be used")

    if 'topol.top' not in mdsfiles and 'tlptopol.top' in mdsfiles:
        print("We found only tlptopol.top in the working directory. This will be renamed to topol.top for use in MDS")
        os.rename(os.path.join(mdshomedir, 'mdsgmx', 'tlptopol.top'), os.path.join(mdshomedir, 'mdsgmx', 'topol.top'))
        
        if not udp > 0 and 'posre_tlptopol.itp' in mdsfiles:
            if not 'posre.itp' in mdsfiles:
                print("Renaming posre_tlptopol.itp to posre.itp")
                os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre_tlptopol.itp'), os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'))
            else:
                print("Backing up posre.itp and renaming posre_tlptopol.itp to posre.itp")
                os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'), os.path.join(mdshomedir, 'mdsgmx', '#posre.itp#'))
                os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre_tlptopol.itp'), os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'))
        elif udp > 0:
            print("User defined restraint file will be used")
            print("Please note that any associated index file must be included as well")
        print('\n')

    if 'topol.top' in mdsfiles and 'tlptopol.top' in mdsfiles:
        print("We found topol.top and tlptopol.top in the working directory. By default, topol.top will be used")
        response = input("To use tlptopol.top instead, type YES/y or press ENTER to continue: ")
        if response.lower() == "yes" or response.lower() == "y":
            print("Backing up topol.top and renaming tlptopol.top ...")
            os.rename(os.path.join(mdshomedir, 'mdsgmx', 'topol.top'), os.path.join(mdshomedir, 'mdsgmx', '#topol.top#'))
            os.rename(os.path.join(mdshomedir, 'mdsgmx', 'tlptopol.top'), os.path.join(mdshomedir, 'mdsgmx', 'topol.top'))

            if not udp > 0 and 'posre_tlptopol.itp' in mdsfiles:
                if not 'posre.itp' in mdsfiles:
                    print("Renaming posre_tlptopol.itp to posre.itp")
                    os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre_tlptopol.itp'), os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'))
                else:
                    print("Backing up posre.itp and renaming posre_tlptopol.itp to posre.itp")
                    os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'), os.path.join(mdshomedir, 'mdsgmx', '#posre.itp#'))
                    os.rename(os.path.join(mdshomedir, 'mdsgmx', 'posre_tlptopol.itp'), os.path.join(mdshomedir, 'mdsgmx', 'posre.itp'))
            elif udp > 0:
                print("User defined restraint file will be used")
                print("Please note that any associated index file must be included as well")
        else:
            print("The default topol.top will be used. You may try tlptopol.top if the process fails")
        print('\n')

    if 'indexfile.ndx' in mdsfiles:
        printNote("Found index file. It will be used by default, except otherwise instructed")
        print("Index file used to generate restraint file must be used together with the file")
        print("If not needed, we have to delete it from this directory")
        print("Type YES/y to confirm the inclusion of index file, OR press ENTER to delete it")
        response = input("Response: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            confirm = input("ARE YOU SURE YOU WANT TO DELETE THE INDEX FILE? (YES/n): ")
            if confirm.lower() == "yes" or confirm.lower() == "y":
                os.remove('indexfile')

    if 'LIGS_at.itp' not in mdsfiles and 'LIGS_mt.itp' not in mdsfiles:
        printNote("Important Ligand(s) topology file(s) are missing. This may not be important")
        printNote("However, Check that all molecule and atom types are captured in the topology file")
        response = input("To continue without these files, type YES/y or press ENTER to abort: ")
        if not (response.lower() == "yes" or response.lower() == "y"):
            raise Exception("Please upload necessary files from mdsetup page")
        print('\n')

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
    print('\n')

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
                print(f"maxwarn has been changed to: {maxwarn}")
                break
    else:
        print(f"maxwarn for this run is: {maxwarn}")

    # Copy mdp files into mdps sub-folder in the working directory
    print('\n')
    print("Copying needed .mdp files to the working directory ...")
    try:
        os.mkdir("mdps")
    except FileExistsError:
        shutil.rmtree("mdps")
        os.mkdir("mdps")

    os.chdir("mdps")
    for file in MDPfiles:
        shutil.copy(os.path.join(fMDP, file), './')
    os.chdir('../')

    printNote("Setup Completed. To continue with MDS, type YES/y or press ENTER to abort")
    response = input("Response: ")
    if not (response.lower() == "yes" or response.lower() == "y"):
        raise Exception("Operation Interrupted by User. Make necessary corrections and restart")
    else:
        # Set the starting date and time
        mdstime = datetime.now()
        Tstart = mdstime.strftime("%B %d, %Y %H:%M:%S")
        print(f'MDS Process begins: {Tstart}')

    time.sleep(5)
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

    if not ('minzsd.tpr' in os.listdir() and 'minzsd.gro' in os.listdir()):
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

    if not ('minzcg.tpr' in os.listdir() and 'minzcg.gro' in os.listdir()):
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

    cmdnvt = []
    if not 'indexfile.ndx' in mdsfiles:
        cmdnvt = ['gmx', 'grompp', '-f', '../mdps/equnvt.mdp', '-c', '../EM/minzcg.gro', '-r', '../EM/minzcg.gro', '-p', '../topol.top', '-o', 'equnvt.tpr', '-maxwarn', str(maxwarn)]
    else:
        cmdnvt = ['gmx', 'grompp', '-f', '../mdps/equnvt.mdp', '-c', '../EM/minzcg.gro', '-r', '../EM/minzcg.gro', '-n', '../indexfile.ndx', '-p', '../topol.top', '-o', 'equnvt.tpr', '-maxwarn', str(maxwarn)]

    try:
        subprocess.run(cmdnvt, stderr=subprocess.STDOUT, stdout=fnvt, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equnvt'], stderr=subprocess.STDOUT, stdout=fnvt, check=True, text=True)
    except subprocess.SubprocessError as e:
        print(e)
        printWarning("Something was not right with above error. Checking....")

    if not ('equnvt.tpr' in os.listdir() and 'equnvt.gro' in os.listdir()):
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

    cmdnpt = []
    if not 'indexfile.ndx' in mdsfiles:
        cmdnpt = ['gmx', 'grompp', '-f', '../mdps/equnpt.mdp', '-c', '../NVT/equnvt.gro', '-r', '../NVT/equnvt.gro', '-p', '../topol.top', '-t', '../NVT/equnvt.cpt', '-o', 'equnpt.tpr', '-maxwarn', str(maxwarn)]
    else:
        cmdnpt = ['gmx', 'grompp', '-f', '../mdps/equnpt.mdp', '-c', '../NVT/equnvt.gro', '-r', '../NVT/equnvt.gro', '-n', '../indexfile.ndx', '-p', '../topol.top', '-t', '../NVT/equnvt.cpt', '-o', 'equnpt.tpr', '-maxwarn', str(maxwarn)]

    try:
        subprocess.run(cmdnpt, stderr=subprocess.STDOUT, stdout=fnpt, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equnpt'], stderr=subprocess.STDOUT, stdout=fnpt, check=True, text=True)
    except subprocess.SubprocessError as e:
        printWarning("Something was not right. Checking...")

    if not ('equnpt.tpr' in os.listdir() and 'equnpt.gro' in os.listdir()):
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

    cmdeq = []
    if not 'indexfile.ndx' in mdsfiles:
        cmdeq = ['gmx', 'grompp', '-f', '../mdps/equpmd.mdp', '-c', '../NPT/equnpt.gro', '-r', '../NPT/equnpt.gro', '-p', '../topol.top', '-t', '../NPT/equnpt.cpt', '-o', 'equpmd.tpr', '-maxwarn', str(maxwarn)]
    else:
        cmdeq = ['gmx', 'grompp', '-f', '../mdps/equpmd.mdp', '-c', '../NPT/equnpt.gro', '-r', '../NPT/equnpt.gro', '-n', '../indexfile.ndx', '-p', '../topol.top', '-t', '../NPT/equnpt.cpt', '-o', 'equpmd.tpr', '-maxwarn', str(maxwarn)]

    try:
        subprocess.run(cmdeq, stderr=subprocess.STDOUT, stdout=feq, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'equpmd'], stderr=subprocess.STDOUT, stdout=feq, check=True, text=True)
    except subprocess.SubprocessError as e:
        printWarning("Something was not right. Checking...")

    if not ('equpmd.tpr' in os.listdir() and 'equpmd.gro' in os.listdir()):
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

    cmdpmd = []
    if not 'indexfile.ndx' in mdsfiles:
        cmdpmd = ['gmx', 'grompp', '-f', '../mdps/pmds.mdp', '-c', '../EMD/equpmd.gro', '-r', '../EMD/equpmd.gro', '-p', '../topol.top', '-t', '../EMD/equpmd.cpt', '-o', 'pmds.tpr', '-maxwarn', str(maxwarn)]
    else:
        cmdpmd = ['gmx', 'grompp', '-f', '../mdps/pmds.mdp', '-c', '../EMD/equpmd.gro', '-r', '../EMD/equpmd.gro', '-n', '../indexfile.ndx', '-p', '../topol.top', '-t', '../EMD/equpmd.cpt', '-o', 'pmds.tpr', '-maxwarn', str(maxwarn)]

    try:
        subprocess.run(cmdpmd, stderr=subprocess.STDOUT, stdout=fpmd, check=True, text=True)
        subprocess.run(['gmx', 'mdrun', '-deffnm', 'pmds'], stderr=subprocess.STDOUT, stdout=fpmd, check=True, text=True)
    except subprocess.SubprocessError as e:
        printWarning("Something was not right. Checking...")

    if not ('pmds.tpr' in os.listdir() and 'pmds.gro' in os.listdir()):
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

    time.sleep(2)

    # End
    print(f"Ending. Job Completed Successfully for MDS codenamed {newWDIR}")

    # Set the ending date and time
    mdstime2 = datetime.now()
    Tends = mdstime2.strftime("%B %d, %Y %H:%M:%S")
    print(f'MDS Process ends: {Tends}')
