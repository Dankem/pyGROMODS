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

from tkinter import Tk, filedialog
from inputimeout import inputimeout, TimeoutOccurred
from pytimedinput import timedInput

from gmodsScripts.gmodsHelpers import printWarning, printNote, select_folder, tinput

def nvtSCmds():
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
    response = tinput("Response: ", 30, "n")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
        response = tinput("To abort type YES/y. Otherwise process will continue anyway ", 10, "y")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
        response = tinput("To abort type YES/y. Otherwise process will continue anyway ", 10, "y")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
        response = tinput("To abort type YES/y. Otherwise process will continue anyway ", 10, "y")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
    response = tinput("Response: ", 10, "n")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
        response = tinput("To abort type YES/y. Otherwise process will continue anyway ", 10, "y")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
        response = tinput("To abort type YES/y. Otherwise process will continue anyway ", 10, "y")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
    response = tinput("Response: ", 10, "n")
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
    print("Starting Equilibration Production MD simulation...")

    try:
        os.mkdir("EMD")
    except FileExistsError:
        shutil.rmtree("EMD")
        os.mkdir("EMD")

    os.chdir("EMD")
    if opencl_activate > 0:
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
        response = tinput("To abort type YES/y. Otherwise process will continue anyway ", 10, "y")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
    response = tinput("Response: ", 10, "n")
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
        oclfiles = os.listdir('../../opencl')
        for ocl in oclfiles:
            os.system('cp ../../opencl/' + ocl + ' ./')

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
    printNote("Ending. Continuation MDS Job Completed Successfully")
