# pyGROMODS
#### Video Demo:  https://youtu.be/Pt0WMhvGkv8
#### Description:

## Requirements: 	
			Python 3 or higher
			Antechamber (from AmberTools preferably)
			OpenBabel (strongly recommended for use with acpype)
			acpype (latest version recommended with all its requirements)
			Gromacs (Compulsory)
			flask to be installed via pip
			flaskwebgui to be installed via pip
			Check requirements.txt file for list of required python modules

## Operating Systems:
			It is expected to work on all platforms that can successfully install all the requirements above. However, the scripts were only tested on:
				Windows 10 (PowerShell with miniconda3)
				Windows 10 Windows Subsystem Linux (running Ubuntu 20.04.2 LTS)
				Ubuntu 20.04.2 LTS
				MSYS2 and Cygwin (GUI can be launched by copy and paste flask link on browser)

## DISCLAIMER
		This code is released under the MIT License.

          <<<  PLEASE READ LICENSE.md FILE!!!  >>>

		You can copy, reuse, distribute and make amendments with prior permission of the Author

## APPRECIATION:
		This work was inspired by:
			CS50 online training of which the initial code serves as part of the final project

		My Special Appreciation therefore goes to ** Prof. David J. Malan ** (malan@harvard.eduand) and His Team of Experts, for making this course available for free to millions of persons like me. Thank You.

## AKNOWLEDGEMENT:
		The Author aknowledge the use of freely available python packages imported (as contained in the requirements.txt file) into this module.

## AUTHORS DETAILS:
		Daniyan, Oluwatoyin Michael, 
		(B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) in Biochemistry)
		Department of Pharmacology, Faculty of Pharmacy, 
		Obafemi Awolowo Universiy, Ile-Ife, Nigeria.
		mdaniyan@oauife.edu.ng; toyinpharm@gmail.com

*** PLEASE Read this README.md file and also follow instructions on the GUI and/or Terminal ***

## INTRODUCTION:
	This project was inspired by a desire to make computational drug discovery less strenous and easy to apply for many scientists in biological and biomedical sciences who are not conversant with or have limited knowledge of command line interface, that are often necessary while making use of available computational resources.

	One important aspect of study biological systems computationally is Molecular Dynamic Simulation. From proper ligand preparation to setting up of receptor - ligand complex for simulation. The pyGROMODS platform, here presented, provides opportunity for more insterested scientists, who are hitherto take aback by lack of interest in often strenous computer skills, to take advantage of enormous possibilities in applying computational approaches to research. pyGROMODS is a Graphical User Interface, making it easy to use by all. All needed files can be easily uploaded and needed files for MDS will be generated.


## INSTALLATION:
	pyGROMODS can be obtained by direct download from https://github.com/Dankem/pyGROMODS or by using Git or checkout with SVN. The downloaded pyGROMODS must be installed in a location where the user has write access. This is becuase folders and files will be generated in the process. In linux environment, you can set the PATH to the pyGROMODS in the .bashrc file. Setting the PATH is the recommended approach for installation. It allows user to have working directory in any other place on the system without tampering with the scripts. After setting the PATH and source the .bashrc file, just type gmodsapp.py on the terminal to launch the platform. However, if PATH can not be set, you can navigate to the directory containing pyGROMODS folder. To launch the application, type python pyGROMODS/gmodsapp.py or python3 pyGROMODS/gmodsapp.py. You will be required to choose a working directory before the platform GUI interface will launch. DON'T CHOOSE INSTALLATION PATH AS YOUR WORKING DIRECTORY.All folders and file will be generated and saved in the working directory. 

	Also, since this package is just a gui interface for setting up and running molecular dynamic simulation with gromacs, all dependent programs (gromacs, ambertools antechamber, acpype, etc) are correctly installed and their executables can be found in the path. This package works was tested with various versions of gromacs, ambertools and acpype (from 2014 to 2022). Therefore, it is advisable to install the latest version of these program for use with pyGROMODS.


## USAGE INSTRUCTIONS:
	*** PLEASE NOTE: ***
	This is not a docking platform. It is therefore essential that where applicable, all dockings would have been done using your preferred docking software. Also, this platform does not curently provide for MDS analysis using gromacs. However, the author may consider including some commonly used analysis in the near future. Also, there are other freely available analysis tools that can be used 

	### **PREPARING NEEDED FILES - RECEPTOR, LIGAND AND/OR PEPTIDES:**
	***For Ligand:*** Remove the ligand(s) from the complex, careful to maintain its docked orientation and position, then, Prepare the ligand(s) by adding hydrogen and checking that the charge is zero, then perform an intial manual check of the ligand using Antechamber, and where necessary, correct any errors.

 	***For Receptor/Proteins:*** While maintaning it in its docked orientation, prepare the protein by removing all hydrogen and all connectivity using your preferred software. Then save the protein separately. Errors in downloaded protein structures can be corrected using SwissPDB viewer, pdb4amber, or tleap depending on the issue at hand. Check the suitability of protein manually using gmx pdb2gmx and tleap commands. Correct any other errors where necessary.

	***For User Supplied Ligand(s) Topology(ies):*** Depending on your choice of preferred forcefiled (amber, charmm, gromos or opls), you may want to generate corresponding ligand topologies from other available web servers. You will be provided options to upload these files along with your ligands. All user suppied topology files must have have at least the following subheadings [ moleculetype ], [ atomtypes ] and [ system ] to be successfully used by the platform </p>

	***FIRST STEP: SELECT THE DESIRED ROUTE***

	*** FOR RECEPTOR - LIGAND COMPLEXES: ***
	This platform provides three routes to follow in generating solvated receptor - ligand complex. The route to use depend on your needs and are summarized as follows:

		***RLmulti Route:***
		This script will generates Receptor - multi_Ligands complex. This is good for a protein bound with more than one ligand or bound to a major ligand and minor moeities like Zn, PO4, SO4, etc. These minor moeities should have been shown to be important for ligand interactions with the receptor by docking. If they are not essential for your work, please remove them in the process of protein preparation. If important, they must be prepared as ligands. This script will produce only one complex having all the supplied ligand interacting with the receptor. Only one receptor can be supplied to this route, but multiple ligands are allowed.

		***RLmany Route:***
		This route will generates separate Receptor - Ligand complex between the receptor and each of the ligand you supplied. This route is important where multiple ligands have been screened by docking against a target protein. The route is not suitable where minor moeities are important, requiring multi-complex preparation. Only one receptor can be supplied to this route, but multiple ligands are allowed.

		***RLsingle Route:***
		This route generates Receptor - Ligand complex between pairs of protein and ligand. This route accepts more than one receptors and ligands but must be uploaded in matching pair of protein and ligand.

	*** FOR PROTEINS AND PEPTIDES: ***
		***PPmore Route:***
		This route will prepare peptides and proteins and generate topology files suitable for molecular dynamic simulation with gromacs. To use this route, the uploaded structure must not have any none standard animo acids residues. As such, protein - protein complexes can be uploaded as protein here. If the peptide has been converted to protein structure (.pdb), it should be uploaded as protein.

	***THE SECOND STEP: UPLOAD THE NECESSARY FILES***
	Having selected the desired route, your preferred forcefield group (amber, charmm, gromos or opls) and weither or not to upload user generated ligand topology(ies), proceed to upload menu to upload all the necessary files. 

	***THE THIRD STEP: RUN THE SELECTED SETUP ROUTE***
	Now, a new working directory will be created for each job using the unique ID number and your supplied name. If no name was supplied, the default name identifying each setup route will be used. In each directory at each step in the workflow sub-directories will be created for maximum job organization.

	Please check sample folder in your working directory for samples of .mdp and peptide files for use with this platform. Each peptide must be prepared in the format presented in the samples. You can upload as many peptide files as possible when using the platform.

	***THE FOURTH STEP: CHECK WITH SHORT MDS USING SUPPLIED .mdp FILES***
	Thereafter, you can check the suitability of the generated solvated complex with short MDS run using the supplied sample MDP files via the MDS page

	***THESE STEPS ARE FURTHER SUMMARIZED BELOW***

	***To Setup Receptor - Ligand(s) Complex:***
		1). Read the QGuides in the GUI interface
		2). Click on SetRoute menu to choose which route to run in generating the complex
		3). Immediately you set the route other menu buttons will become available
		4). Select your preferred forcefield group and weither or not to upload ligand topology file(s) or in case of PPmore route, you are to select either Peptide or Protein
		5). Click Uploads menu to uploads either your ligand(s) and receptor files (and if required, ligand(s) topology(ies)) or Peptides or Proteins
		6). Go to Generate menu to run the selected route
			*** PLEASE NOTE: ***
			The platform provides for two modes for generating MDS input files. They are Interactive or Noninteractive modes. Interactive allows the users to actively interact with the system, make informed decision by taking time to answer questions from input request and make changes as may be needed. This interactive approach is recommended for new users, and will be very useful when comparing suitability of different forcefields, box types, box sizes, etc before taking final decision on default parameters to use. RLmany and RLsingle can come handy when troubleshooting with interactive mode. Experienced users or users with well defined default parameters can work with noninteractive mode. In both cases, other default parameters (timeout, editconf -d and -bt) should be set. The default values you see are just ok for noninteractive mode, although, you can extend the timeout, to give you ample time to study the platform. Meanwhile, Noninteractive mode will move faster by reducing the timeout.	
		7). If succesfful, a message will be displayed
			*** PLEASE NOTE: ***
				a). Both posre.itp and posre_udp.itp (if generated at run time) restraint files are present.
				b). To use any, it must be defined in .mdp file. Else, only the one defined will be used.
				c). Alternatively, you may backup one and rename other to correspond to either -DPOSRE (for posre.itp) or 
				-DPOSRE_UDP (for posre_udp.itp) defined in .mdp
		8). ***PLEASE NOTE: There are other onscreen instruction that may require your attention as scripts runs on terminal. Stay on the Terminal while script runs***

	***To Perform Molecular Dynamic Simulation:***
		1). Click on MDS and select FRESH or CONTINUATION. FRESH will required upload of necessary files as indicated. All needed files for MDS are saved in gmxmds folder during preparatory stage above. The "samples" subfolder inside the working directory also contains samples of mdp files for you to edit and upload. Alternatively, if you generated your MDS input files from another source, please follow Manual setup approach below.
		2). Make sure you have your .mdp (molecular dynamic parameter) files ready. Files must be uploaded for minimization, equilibration and production runs. While user may not be interested in all the stages defined for mdp, none can be left empty. As such, upload mdp files most suitable for each defined stages. For instance, user can upload same minimization methods, steepest descent or conjugate gradient or different one, for both, as long as such methods are compatible with gromacs. This also applies to equilibration.
		3). Click upload after selecting the files, you will be redirected to MDSrun page where you can choose to activate opencl. Activating openCL is only neccessary if you need it or you are having problem compiling openCL headers.
		4). Click MDSrun and MDS simulation will start. Folders and files will be generated
		5). NOTE: MDS menu is available at startup of the GUI. You can use it for already generated complex.
		6). If succesfful, a message will be displayed
		7). PLEASE NOTE: There are other onscreen instruction that may require your attention as scripts runs on terminal. Stay on the Terminal while script runs

	*** PLEASE NOTE: Uploading needed files from MDS menu is the recommended approach. Nevertheless, you can follow guidelines below to do it manually ***

	***To Perform Molecular Dynamic Simulation By Manual Setup of needed files:***
		1). Create 'MDP' folder anywhere and populate it with files with renamed as follows:
				a. minzsd.mdp --- for Steepest Descent Minimization
				b. minzcg.mdp --- for Conjugate Gradient Minimization
				c. equnvt.mdp --- for Constant Volume Equilibration
				d. equnpt.mdp --- for Constant Pressure Equilibration
				e. equpmd.mdp --- for Equilibration Production MD
				f. pmds.mdp --- for Production MD
				You can check Samples of .mdp files in the 'sample' folder inside your setup work directory. These can be adjusted to suite your purpose.
				PLEASE NOTE: All these five files must be provided. However, you do not need to follow the suggested stages as provided. If you don't need any of the stages, just duplicate the one you need to make up the number, but your files must be named as given. For example, you may be interested in 10,000 steps of steepest descent minimization only. All you need do is to put the same settings with 5000 steps each in minzsd.mdp and minzcg.mdp files.
		2). Create 'MDSHOME' folder inside the working directory, if not present
		3). Backup any former 'mdsgmx' folder inside MDSHOME if not having unique ID attached
		4). Create a new 'gmxmds' folder anywhere and populate with the following files labelled as such:
				a. topol.top
				b. LIGS_at.itp --- for Ligand atomtypes, if needed
				c. LIGS_mt.itp --- for Ligand moleculetypes, if needed
				d. posre.itp and/or posre_udp.itp --- for restraint file, if needed
				e. fsolvated.gro --- for final solvated structure file
		5). Go to the GUI platform, click MDS menu, select fresh, and upload necessary files as prepared. Following succcessful upload, you will be redirected to MDSrun page.
		6). HAVING PROBLEM COMPILING OPENCL AT MDS RUN TIME, Try the following:
				a. Create a reuseable OpenCL folder
				b. Copy all files found inside 'ewald', 'gpu_utils', 'mdlib', and 'pbcutil' subdirectories of "\gromacs\share\gromacs\opencl\gromacs" into OpenCL folder
				c. Check and edit all files in OpenCL folder, making sure especially that the path to '#include ..." point to the this directory.
				d. Copy any other #include files found somewhere else (e.g. /gromacs/include/gromacs/utility) into OpenCL folder
				e. On the MDSrun page, click ACTIVATE to activate and redirect to page where you can upload these files for use
				f. Alternatively, though not recommended, you can create 'opencl' folder inside 'MDSHOME' directory and copy content of OpenCL folder into it
		7). On the platform, MDSrun page, click MDSrun to run the MD Simulation.


## EDITABLE FILES:
	You can edit the files contained in the following folders as needed. However, before editing any of the files, backup these original files into a safe place. They are:
	***FILES IN gmodsTSF FOLDERS***
		1). This folder contains some tleap source files. You can edit to change e.g. forcefields, water etc. BUT THE NAME OF THE FILES AS GIVEN MUST BE MAINTAINED
		2). The supplied source files take care of generating topologies and/or structure files for ligands, non-standard residues and peptides, as well as generating receptor - multi ligands complex for up to ten ligands
		3). Copies of these files have been copied into 'samples/tleap' subfolder in your work directory. You can edit these copies and use it to replace the content of gmodsTSF folder
		4). However, for most usage, editing these files will not be necessary

	***FILES IN sample SUBFOLDER IN YOUR pyGROMODS WORKSPACE DIRECTORY***
		A 'samples' folder is created in your workplace directory. It contains samples of ".mdp" files in 'mdps' subfolder and ".in" files for peptides in 'peptides' subfolders. The '.mdp' files can be adjusted as needed and thereafter uploaded via MDS menu on the platform. Remember, while moving, the files need to be renamed as indicated above.


### UPGRADING / INSTALLING PYTHON MODULES DEPENDENCIES :
	The pyGROMODS folder contained requirements.txt file. "pip install -r requirements.txt --upgrade" can be used to install all the python dependencies. However, the package comes with a python script, "uppreqs.py" that can be called to do the job and produce updated version of requirements.txt file for futher use. It can also be used to create requirements.txt file in case it's missing. Fortunately, this sript is multi-purpose, as it can be used to update dependencies of any python packages. This script can be called from anywhere.


### TROUBLESHOOTINGS:
	When faced with failure to generate the required solvated structure using topol.top and/or tlptopol.top, You can try the following:
		1). Look into error files 'check' and 'check2' and their backups versions and attempt to fix the errors and rerun
		2). Try manual solvation with files gathered into gmxmds folder
		3). If the above fail to resolve the error(s), try the following
				a). Uncomment the three rows above the last row of relevant tleap source file
				b). Rerun the setup to generate tleap solvated files - tlpSolvated.gro and tlpSolvated.top
				c). Open these files in any text editor and change WAT to SOL in both files
				d). Try running Gromacs MDS with the updated files
