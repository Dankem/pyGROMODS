#!/usr/bin/env python
"""
    This code is released under GNU General Public License V3.

          <<<  NO WARRANTY AT ALL!!!  >>>

	Daniyan, Oluwatoyin Michael, B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) Biochemistry
    Department of Pharmacology, Faculty of Pharmacy
    Obafemi Awolowo Universiy, Ile-Ife, Nigeria.
    >>http://www.oauife.edu.ng<<

    mdaniyan@oauife.edu.ng; toyinpharm@gmail.com

    This code checks for the installation of needed packages
"""
import sys

if sys.version_info < (3, 5):
    raise Exception("Python 3.5 or a more recent version is required.")

import os
import string
import platform
from os.path import join

def reqPackages_Check():
    check_results = []
    check_files = ['gmx', 'antechamber', 'parmchk2', 'tleap', 'acpype', 'acpype.py', 'gmx.exe', 'antechamber.exe', 'parmchk2.exe', 'teleap.exe']
    availableDrives = ['%s:' % d for d in string.ascii_uppercase if os.path.exists('%s:' % d)]

    def check_Windows(availableDrives, check_files):
        windows_checks = []

        for i in range(len(availableDrives)):
            for root, dirs, files in os.walk(availableDrives[i]):
                for ckfile in check_files:
                    if ckfile in files and ckfile not in windows_checks:
                        windows_checks.append(ckfile)

        for i in range(len(availableDrives)):
            driveFolder = join(availableDrives[i], "/")
            availableFolders = os.listdir(driveFolder)
            tocheck = ['Program Files', 'Program Files (x86)', 'ProgramData']
            for f in availableFolders:
                if f == 'gromacs' or f == 'amber' or f == 'acpype' or f[0:7] == 'gromacs' or f[0:5] == 'amber' or f[0:6] == 'acpype':
                    targetFolder = join(driveFolder, f)
                    for root, dirs, files in os.walk(targetFolder):
                        for ckfileB in check_files:
                            if ckfileB in files and ckfileB not in windows_checks:
                                windows_checks.append(ckfileB)

                elif f in tocheck:
                    newFolder = join(driveFolder, f)
                    listnewFolder = os.listdir(newFolder)
                    for j in listnewFolder:
                        if j == 'gromacs' or j == 'amber' or j == 'acpype' or j[0:7] == 'gromacs' or j[0:5] == 'amber' or j[0:6] == 'acpype':
                            targetFolder = join(newFolder, j)
                            for root, dirs, files in os.walk(targetFolder):
                                for ckfileB in check_files:
                                    if ckfileB in files and ckfileB not in windows_checks:
                                        windows_checks.append(ckfileB)
        return windows_checks

    def check_Linux(availableDrives, check_files):
        linux_checks = []

        for i in range(len(availableDrives)):
            for root, dirs, files in os.walk(availableDrives[i]):
                for ckfile in check_files:
                    if ckfile in files and ckfile not in linux_checks:
                        linux_checks.append(ckfile)

        try:
            import whichcraft as wch
        except Exception as e:
            print(f"whichcraft module is not installed or not found. Detect error {e}")
            return linux_checks

        for file in check_files:
            found = wch.which(file)
            if found is not None and file not in linux_checks:
                linux_checks.append(file)

        return linux_checks

    if platform.system().lower() == "windows":
        check_results = check_Windows(availableDrives, check_files)

    elif platform.system().lower() == "linux" or os.name == "posix":
        print("Assuming Linux or Posix Platform")
        check_results = check_Linux(availableDrives, check_files)

    else:
        print("We could not check your system at this time")
        print("Please make sure the following Packages are installed and discoverable:") 
        print("gromacs")
        print("antechamber")
        print("parmchk2")
        print("tleap")
        print("acpype")

    if not len(check_results) > 0:
        return "Unable to Check"

    passed_Checks = []
    failed_Checks = []

    if "gmx" in check_results or "gmx.exe" in check_results:
        passed_Checks.append("gromacs")
    else:
        failed_Checks.append("gromacs")

    if "antechamber" in check_results or "antechamber.exe" in check_results:
        passed_Checks.append("antechamber")
    else:
        failed_Checks.append("antechamber")

    if "parmchk2" in check_results or "parmchk2.exe" in check_results:
        passed_Checks.append("parmchk2")
    else:
        failed_Checks.append("parmchk2")

    if "tleap" in check_results or "tetleap.exe" in check_results:
        passed_Checks.append("tleap")
    else:
        failed_Checks.append("tleap")

    if "acpype" in check_results or "acpype.py" in check_results:
        passed_Checks.append("acpype")
    else:
        failed_Checks.append("acpype")

    if len(passed_Checks) > 0:
        print("The Installation of the following Package(s) was/were detected")
        for pp in passed_Checks:
            print(pp)

    if len(failed_Checks) > 0:
        print("The Installation of the following Package(s) could not be detected")
        for ff in failed_Checks:
            print(ff)
        print("Please ensure you install the package(s) and/or make the package discoverable")
        print("On Windows: Install in 'C:/', 'Program Files', 'Program Files (x86)', or 'ProgramData'")
        print("On Linux: Install where you have write access and provide the path in .bashrc file")

    return f"Successfully Checked: {len(passed_Checks)} Passed AND {len(failed_Checks)} Failed"
