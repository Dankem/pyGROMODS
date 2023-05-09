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
import subprocess
from pathlib import Path
import time
import shutil

from gmodsScripts.gmodsHelpers import printWarning, printNote, tinput, cleanup

def nvtSCmds():
    print('\n')
    # Set global activation variables to include opncl headers in work subdirectories
    opencl_activate = 0

    # Check existence of opencl files if required
    if not os.path.isdir('../opencl'):
        opencl_activate = 0
    else:
        openclfiles = os.listdir('../opencl')
        if len(openclfiles) == 0:
            opencl_activate = 0
        else:
            opencl_activate += 1

    # User can set marwarn for grompp if needed
    maxwarn = 0
    print("The default value for -maxwarn is 0. To change the value: ")
    print("Type YES/y. Otherwise, press ENTER to continue with the default")
    response = tinput("Response: ", 20, "n")
    if response.lower() == "yes" or response.lower() == "y":
        while True:
            maxwarn = input("Change maxwarn to: ")
            if maxwarn == " " or maxwarn == "" or maxwarn.isdigit() == False:
                print("You must supply digit value for maxwarn to continue")
                continue
            else:
                print(f"maxwarn has been changed to: {maxwarn}")
                break

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
            shutil.rmtree("opencl")
        except:
            pass

        os.chdir('../')

    else:
        os.chdir('../../')

    time.sleep(5)

    # End
    printNote("Ending. Continuation MDS Job Completed Successfully")


def nptSCmds():
    # Set global activation variables to include opncl headers in work subdirectories
    opencl_activate = 0

    # Check existence of opencl files if required
    if not os.path.isdir('../opencl'):
        opencl_activate = 0
    else:
        openclfiles = os.listdir('../opencl')
        if len(openclfiles) == 0:
            opencl_activate = 0
        else:
            opencl_activate += 1

    # User can set marwarn for grompp if needed
    maxwarn = 0
    print("The default value for -maxwarn is 0. You may wish to increase this as needed")
    print("To change this value, type YES/y. Otherwise, press ENTER to continue with the default")
    response = tinput("Response: ", 20, "n")
    if response.lower() == "yes" or response.lower() == "y":
        while True:
            maxwarn = input("Change maxwarn to: ")
            if maxwarn == " " or maxwarn == "" or maxwarn.isdigit() == False:
                print("You must supply digit value for maxwarn to continue")
                continue
            else:
                print(f"maxwarn has been changed to: {maxwarn}")
                break

    time.sleep(5)

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
            shutil.rmtree("opencl")
        except:
            pass

        os.chdir('../')

    else:
        os.chdir('../../')

    time.sleep(5)

    # End
    printNote("Ending. Continuation MDS Job Completed Successfully")


def emdSCmds():
    # Set global activation variables to include opncl headers in work subdirectories
    opencl_activate = 0

    # Check existence of opencl files if required
    if not os.path.isdir('../opencl'):
        opencl_activate = 0
    else:
        openclfiles = os.listdir('../opencl')
        if len(openclfiles) == 0:
            opencl_activate = 0
        else:
            opencl_activate += 1

    # User can set marwarn for grompp if needed
    maxwarn = 0
    print("The default value for -maxwarn is 0. You may wish to increase this as needed")
    print("To change this value, type YES/y. Otherwise, press ENTER to continue with the default")
    response = tinput("Response: ", 20, "n")
    if response.lower() == "yes" or response.lower() == "y":
        while True:
            maxwarn = input("Change maxwarn to: ")
            if maxwarn == " " or maxwarn == "" or maxwarn.isdigit() == False:
                print("You must supply digit value for maxwarn to continue")
                continue
            else:
                print("maxwarn has been changed to: ", maxwarn)
                break

    time.sleep(5)

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
            shutil.rmtree("opencl")
        except:
            pass

        os.chdir('../')

    else:
        os.chdir('../../')

    time.sleep(5)

    # End
    printNote("Ending. Continuation MDS Job Completed Successfully")


def pmdSCmds():
    # Set global activation variables to include opncl headers in work subdirectories
    opencl_activate = 0

    # Check existence of opencl files if required
    if not os.path.isdir('../opencl'):
        opencl_activate = 0
    else:
        openclfiles = os.listdir('../opencl')
        if len(openclfiles) == 0:
            opencl_activate = 0
        else:
            opencl_activate += 1

    # User can set marwarn for grompp if needed
    maxwarn = 0
    print("The default value for -maxwarn is 0. You may wish to increase this as needed")
    print("To change this value, type YES/y. Otherwise, press ENTER to continue with the default")
    response = tinput("Response: ", 20, "n")
    if response.lower() == "yes" or response.lower() == "y":
        while True:
            maxwarn = input("Change maxwarn to: ")
            if maxwarn == " " or maxwarn == "" or maxwarn.isdigit() == False:
                print("You must supply digit value for maxwarn to continue")
                continue
            else:
                print(f"maxwarn has been changed to: {maxwarn}")
                break

    time.sleep(5)

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
            shutil.rmtree("opencl")
        except:
            pass

        os.chdir('../')

    else:
        os.chdir('../../')

    time.sleep(5)

    # End
    printNote("Ending. Continuation MDS Job Completed Successfully")
