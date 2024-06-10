#!/usr/bin/env python
"""
    pyGROMODS-v2024.01 Release
    
            <<<  NO WARRANTY AT ALL!!!  >>>
    
    This Script can be run from any where and can be used to update dependencies of python 
    packages as contained in the requirements.txt file. It can also generate requirements.txt 
    file and produce an updated version of the file. It can be placed anywhere, but preferably,
    in a place where it can be easily accessed, such as in linux, you can set path 
    to it  for ease of use across board
"""
import sys

if sys.version_info < (3, 5):
    raise Exception("Python 3.5 or a more recent version is required.")

import os
import shutil
import subprocess
from tkinter import Tk, filedialog

# Select the desired python package directory
def main():
    root = Tk()
    root.withdraw()
    root.attributes('-topmost', True)

    uppDIR = filedialog.askdirectory(title="Select Python Package Directory", initialdir=".")

    root.destroy()
    os.path.normpath(uppDIR)

    # Create working folder inside the current directory
    try:
        os.mkdir("getUpdated")
    except FileExistsError:
        try:
            shutil.rmtree("getUpdated")
        except:
            os.remove("getUpdated")
        finally:
            pass
        os.mkdir("getUpdated")

    os.chdir("getUpdated")

    # Check for the presence of requirements.txt file in the project folder
    if not "requirements.txt" in os.listdir(uppDIR):
        print("Missing requirements.txt file detected")
        print("Trying to generate new requirements.txt file ....")

        # Check for required pipreqs python module
        try:
            subprocess.run("pip freeze | grep pipreqs > chkreq", shell=True, stderr=subprocess.STDOUT, check=True, text=True)
        except subprocess.SubprocessError as e:
            if not "chkreq" in os.listdir():
                print("Unable to detect installation of required python module, pipreqs")

        if not os.path.getsize("chkreq") > 0:
            print("Getting ready to install required package: pipreqs")
            response = input("Type YES/y to install or press ENTER to abort")
            if not (response.lower() == "yes" or response.lower() == "y"):
                sys.exit()
            else:
                try:
                    subprocess.run("pip install pipreqs", shell=True, stderr=subprocess.STDOUT, check=True, text=True)
                except subprocess.SubprocessError as e:
                    print(f"Installation of pipreqs failed with error {e}")
                    print("Please include requirements.txt file in your package directory")
                    sys.exit()

        # Generate the requirements.txt file
        try:
            subprocess.run("pipreqs " + uppDIR, shell=True, stderr=subprocess.STDOUT, check=True, text=True)
        except subprocess.SubprocessError as e:
            print(f"Generating requirements.txt file failed with error {e}")
            print("Please include requirements.txt file in your package directory")
            sys.exit()

        # Re-check for the presence of requirements.txt file
        if not "requirements.txt" in os.listdir(uppDIR):
            print("Could not locate the generated requiremnts.txt file")
            sys.exit()
        else:
            print("Generation of requirements.txt file was successful")

    # Copy the requirements.txt file into the getUpdated directory
    outdated = "requirements.txt"
    shutil.copy(os.path.join(uppDIR, outdated), './')

    # Open and parse the outdated file for one by one pip upgrading
    file = open(outdated, "r")
    readline = file.readlines()
    newreqs = "newreqs.txt"
    openreqs = open(newreqs, "+a")

    failed = []
    for line in readline:
        package = ""
        if line.count("=") == 2:
            package = line.strip().replace("=", ">", 1)
        elif line.count("=") == 1 and line.count(">") == 1:
            if line.strip().index("=") > line.strip().index(">"):
                package = line.strip()
            else:
                package = line.strip().replace(">", "=").replace("=", ">", 1)
        else:
            pass
        
        try:
            print(f"Running: Upgrading {package} ....")
            subprocess.run("pip install --upgrade --upgrade-strategy=only-if-needed " + package, shell=True, stderr=subprocess.STDOUT, check=True, text=True)
        except subprocess.SubprocessError as e:
            print(f"Installation of {package} instalation failed {e}")
            failed.append(package)

        try:
            package = line.strip().split("==")[0]
            openreqs.write(package)
            openreqs.write('\n')
        except:
            pass

    openreqs.close()
    file.close()

    if len(failed) > 0:
        print("The following package failed upgrading. Check the individually generated log file by version number for details")
        for f in failed:
            print(f)

    # Generate new requirements.txt file following upgrade
    print("Generating new requirements.txt file ...")
    os.remove(os.path.join(uppDIR, outdated))
    os.rename(outdated, "oldreqs.txt")

    try:
        subprocess.run("pipreqs " + uppDIR, shell=True, stderr=subprocess.STDOUT, check=True, text=True)
    except:
        try:
            subprocess.run("pip freeze | grep -f newreqs.txt > required.txt", shell=True, stderr=subprocess.STDOUT, check=True, text=True)
            shutil.copy("required.txt", os.path.join(uppDIR, "requirements.txt"))        
        except:
            print("Generation of new requirements.txt file failed")
            print("Restoring the old file ...")
            shutil.copy("oldreqs.txt", os.path.join(uppDIR, "requirements.txt"))

    shutil.copy("oldreqs.txt", os.path.join(uppDIR, "oldreqs.txt"))
    os.chdir(uppDIR)

    # Updating the new requirement.txt file
    print("Updating the new requirement.txt file...")

    os.rename("requirements.txt", "required.txt")

    fileB = open("required.txt", "r")
    readBline = fileB.readlines()
    latestreqs = "requirements.txt"
    openlatestreqs = open(latestreqs, "+a")

    for lineB in readBline:
        pkgupdate = ""
        if lineB.count("=") == 2:
            pkgupdate = lineB.strip().replace("=", ">", 1)
        elif lineB.count("=") == 1 and lineB.count(">") == 1:
            if lineB.strip().index("=") > lineB.strip().index(">"):
                pkgupdate = lineB.strip()
            else:
                pkgupdate = lineB.strip().replace(">", "=").replace("=", ">", 1)
        else:
            pass
        
        openlatestreqs.write(pkgupdate)
        openlatestreqs.write('\n')

    openlatestreqs.close()
    fileB.close()

    try:
        os.rename("oldreqs.txt", "reqs-backup.txt")
    except FileExistsError:
        os.remove("reqs-backup.txt")
        os.rename("oldreqs.txt", "reqs-backup.txt")

    os.remove("required.txt")
    
    print("Please compare 'reqs-backup.txt' with 'requirements.txt' and update the later with any missing module")
    print("IF INSTALLATION FAILS USING THIS PACKAGE. TRY THE FOLLOWING:")
    print("1). Navigate to the target package folder that contains 'requirements.txt' file")
    print("2). Try this command: pip install --upgrade --upgrade-strategy=only-if-needed -r requirements.txt")
    print("3). Check on-screen output to confirm success or otherwise")

if __name__ == "__main__":
    main()