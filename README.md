# pyGROMODS
Python based Graphical User Interface for generating input files for and running molecular dynamic simulation.

#### Video Demo:  NA
#### Description:

### Requirements: 	Python 3 or higher
						Antechamber (from AmberTools preferably)
            OpenBabel (strongly recommended for use with acpype)
						acpype (latest version recommended with all its requirements)
						Gromacs (Compulsory)
						flask to be installed via pip
						flaskwebgui and/or pyfladesk to be installed via pip

### Operating Systems:
			It is expected to work on all platforms that can successfully install all the requirements above. However, the scripts were only tested on:
						Windows 10 (PowerShell with miniconda3)
						Windows 10 Windows Subsystem Linux (running Ubuntu 20.04.2 LTS)
						Ubuntu 20.04.2 LTS
						MSYS2 and Cygwin (GUI can be launched by copy and paste flask link on browser)

### DISCLAIMER
		This code is released under the MIT License.

          <<<  PLEASE READ LICENSE.md FILE!!!  >>>

		You can copy, reuse, distribute and make amendments with prior permission of the Author

### APPRECIATION:
		This work was inspired by:
			CS50 online training of which the initial code serves as part of the final project

		My Special Appreciation therefore goes to ** Prof. David J. Malan ** (malan@harvard.eduand) and His Team of Experts, for making this course available for free to millions of persons like me. Thank You.

### AUTHORS DETAILS:
		Daniyan, Oluwatoyin Michael, B.Pharm. M.Sc. (Pharmacology) and Ph.D. (Science) Biochemistry
		Department of Pharmacology, Faculty of Pharmacy
		Obafemi Awolowo Universiy, Ile-Ife, Nigeria.
		>>http://www.oauife.edu.ng<<
		mdaniyan@oauife.edu.ng; toyinpharm@gmail.com

*** PLEASE Read this README.md file and also follow instructions on the GUI and/or Terminal ***

### INTRODUCTION:
This project was inspired by a desire to make computational drug discovery less strenous and easy to apply for many scientists in biological and biomedical sciences who are not conversant with computer language, that are often necessary while making use of available computational resources.

One important aspect of study biological systems computationally is Molecular Dynamic Simulation. From proper ligand preparation to setting up of receptor - ligand complex for simulation. The pyGROMODS platform here presented provide opportunity for more insterested scientists, who are hitherto take aback by lack of interest in often strenous computer skills, to take advantage of enormous possibilities in applying computational approaches to research. pyGROMODS is a Graphical User Interface, making it easy to use by all. All needed files can be easily uploaded and files will be generated

### INSTALLATION:
The pyGROMODS must be installed in a location where the user has write access. This is becuase folders and files will be generated in the process. To install, download the package from github (https://github.com/Dankem/pyGROMODS) and place the folder and its content in location you have write access. In linux environment, you can set the PATH to the pyGROMODS in the .bashrc file. Setting the PATH is the recommended approach for installation. It allows user to have working directory in any other place on the system without tampering with the scripts. After setting the PATH and source the .bashrc file, just type gmodsapp.py on the terminal inside any new working directory, the platform GUI interface will launch.

However, if PATH can not be set, create a working directory where you have write access and easy location e.g, on Windows Desktop or in Linux, /home/User/ or usr/local folder. Install pyGROMODS in this folder. To launch the application, navigate to the working directory from the terminal and type python pyGROMODS/gmodsapp.py or python3 pyGROMODS/gmodsapp.py. All folders and file will be generated and saved in the working directory. PLEASE AVOID USING pyGROMODS folder as working directory.

### USAGE INSTRUCTIONS:
	*** PLEASE NOTE: ***
	This is not a docking platform. It is therefore essential that where applicable, all dockings would have been done using your preferred docking software.

	*** PREPARING NEEDED FILES - RECEPTOR, LIGAND AND/OR PEPTIDES:: ***
	***For Ligand:*** Remove the ligand(s) from the complex, careful to maintain its docked orientation and position, then, Prepare the ligand(s) by adding hydrogen and checking that the charge is zero, then perform an intial manual check of the ligand using Antechamber, and where necessary, correct any errors.

 	***For Receptor/Proteins:*** While maintaning it in its docked orientation, prepare the protein by removing all hydrogen and all connectivity using your preferred software. Then save the protein separately. Errors in downloaded protein structures can be corrected using SwissPDB viewer, pdb4amber, or tleap depending on the issue at hand. Check the suitability of protein manually using gmx pdb2gmx and tleap commands. Correct any other errors where necessary.

	***For User Supplied Ligand(s) Topology(ies):*** Depending on your choice of preferred forcefiled (amber, charmm, gromos or opls), you may want to generate corresponding ligand topologies from other available web servers. You will be provided options to upload these files along with your ligands. All user suppied topology files must have have at least the following headings [ moleculetype ], [ atomtypes ] and [ system ] to be successfully used by the platform </p>

	***FIRST STEP IS TO SELECT THE DESIRED ROUTE***

	*** FOR RECEPTOR - LIGAND COMPLEXES: ***
	This platform provides three routes to follow in generating solvated receptor - ligand complex. The route to use depend on your needs and are summarized as follows:

		***RLmulti Route:***
		This script will generates Receptor - multi_Ligands complex. This is good for a protein bound with more than one ligand or bound to a major ligand and minor ones like Zn, PO4, SO4, etc. These minor moeities should	have been shown to be important for ligand interactions with the receptor by docking. This script will produce only
		one complex having all the supplied ligand interacting together with the receptor. Only one receptor can be supplied to this route, but multiple ligands are allowed.

		***RLmany Route:***
		This route will generates separate Receptor - Ligand complex between the receptor and each of the ligand you supplied. This route is important where multiple ligands have been screened by docking against a target protein. Only one receptor can be supplied to this route, but multiple ligands are allowed.

		***RLsingle Route:***
		This route generates Receptor - Ligand complex between pairs of protein and ligand. This route accepts more than one receptors and ligands but must be uploaded in matching pair of protein and ligand.

	*** FOR PROTEINS AND PEPTIDES: ***
		***PPmore Route:***
		This route will prepare peptides and proteins and generate topology files suitable for molecular dynamic simulation with gromacs. To use this route, the uploaded structure must not have any none standard animo acids residues. As such, protein - protein complexes can be uploaded as protein here. If the peptide has been
		converted to protein structure (.pdb), it should be uploaded as protein.

	***THE SECOND STEP IS TO UPLOAD THE NECESSARY FILES***

	***THE THIRD STEP IS TO RUN THE SELECTED SETUP ROUTE***

	Now, a new working directory will be created for each job using the unique ID number and your supplied name. In each directory at each step in the workflow sub-directories will be created for maximum job organization.

	Thereafter, you can check the suitability of the generated solvated structures and topology files with short MDS run using the supplied sample MDP files	via the platform MDS page.

	Please check sample folder in your working directory for examples of .mdp and peptide files for use with this platform. Each peptide must be prepared in the format presented in the samples. You can upload as many peptide files as possible when using the platform.

	***THE FOURTH STEP IS TO CHECK WITH SHORT MDS USING SUPPLIED .mdp FILES***

	Thereafter, you can check the suitability of the generated solvated complex with short MDS run using the supplied sample MDP files via the MDS page

	***THESE STEPS ARE FURTHER SUMMARIZED BELOW***

	***To Setup Receptor - Ligand(s) Complex:***
		1). Read the QGuides in the GUI interface
		2). Click on SetRoute menu to choose which route to run in generating the complex
		3). Immediately you set the route other menu buttons will become available
		4). Select your preferred forcefield group and weither or not upload ligand topology file(s) or in case of PPmore route, you are to select either Peptide or Protein
		5). Click Uploads menu to uploads either your ligand(s) and receptor files (and if required, ligand(s) topology(ies)) or Peptides or Proteins
		6). Go to Generate menu to run the selected route
		7). If succesfful, a message will be displayed
				*** PLEASE NOTE: ***
					a). Both posre.itp and posre_udp.itp (if generated at run time) restraint files are present.
					b). To use any, it must be defined in .mdp file. Else, only the one defined will be used.
					c). Alternatively, you may backup one and rename other to correspond to either -DPOSRE (for posre.itp) or -DPOSRE_UDP (for posre_udp.itp) defined in .mdp
		8). ***PLEASE NOTE: There are other onscreen instruction that may require your attention as scripts runs on terminal. Stay on the Terminal while script runs***

	***To Perform Molecular Dynamic Simulation:***
		1). Click on MDS and upload necessary files as indicated. All needed files for MDS are saved in gmxmds folder during generation of complex
		2). Make sure you have your .mdp (molecular dynamic parameter) files ready. You can use the samples provided in 'sample' folder in your setup workspace directory
		3). Click upload after selecting the files, you will be redirected to MDSrun page
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
				You can check Samples of .mdp files in the 'sample' folder inside your setup work directory. These can be adjusted to suite your purpose
		2). Create 'MDSHOME' folder inside the working directory, if not present
		3). Backup any former 'mdsgmx' folder inside MDSHOME if not having unique ID 	 attached
		4). Create a new 'gmxmds' folder anywhere and populate with the following files labelled as such:
				a. topol.top
				b. LIGS_at.itp --- for Ligand atomtypes, if neeeded
				c. LIGS_mt.itp --- for Ligand moleculetypes, if needed
				d. posre.itp and/or posre_udp.itp --- for restraint file, if needed
				e. fsolvated.gro --- for final solvated structure file
		5). Go to the GUI platform, click MDS menu, upload necessary files as prepared. Following succcessful upload, you will be redirected to MDSrun page.
		6). HAVING PROBLEM COMPILING OPENCL AT MDS RUN TIME, Try the following:
				a. Create a reuseable OpenCL folder
				b. Copy all files found inside 'ewald', 'gpu_utils', 'mdlib', and 'pbcutil' subdirectories of "\gromacs\share\gromacs\opencl\gromacs" into OpenCL folder
				c. Check and edit all files in OpenCL folder, making sure especially that #include point to the files in this directory.
				d. Copy any other #include files found somewhere else (e.g. /gromacs/include/gromacs/utility) into OpenCL folder
				e. On the MDSrun page, click ACTIVATE to activate and redirect to page where you can upload these files for use
				f. Alternatively, though not recommended, you can create 'opencl' folder inside 'MDSHOME' directory and copy content of OpenCL folder into it
		7). On the platform, MDSrun page, click MDSrun to run the MD Simulation.

### EDITABLE FILES:
	You can edit the files contained in the following folders as needed. However, before editing any of the files, backup these original files into a safe place. They are:
	***FILES IN gmodsTSF FOLDERS***
		1). This folder contains some tleap source files. You can edit to change e.g. forcefields, water etc. BUT THE NAME OF THE FILES AS GIVEN MUST BE MAINTAINED
		2). The supplied source files take care of generating topologies and/or structure files for ligands, non-standard residues and peptides, as well as generating receptor - multi ligands complex for up to ten ligands
		3). Copies of these files have been copied into 'samples/tleap' subfolder in your work directory. You can edit these copies and use it to replace the content of gmodsTSF folder
		4). However, for most usage, editing these files will not be necessary

	***FILES IN sample SUBFOLDER IN YOUR pyGROMODS WORKSPACE DIRECTORY***
		A 'samples' folder is created in your workplace directory. It contains samples of ".mdp" files in 'mdps' subfolder and ".in" files for peptides in 'peptides' subfolders. The '.mdp' files can be adjusted as needed and thereafter uploaded via MDS menu on the platform. Remember, while moving, the files need to be renamed as indicated above.

### TROUBLESHOOTINGS:
	When faced with failure to generate the required solvated structure using topol.top and/or tlptopol.top, You can try the following:
		1). Look into error files 'check' and 'check2' and their backups versions and attempt to fix the errors and rerun
		2). Try manual solvation with files gathered into gmxmds folder
		3). If the above fail to resolve the error(s), try the following
				a). Uncomment the three rows above the last row of relevant tleap source file"
				b). Rerun the setup to generate tleap solvated files - tlpSolvated.gro and tlpSolvated.top"
				c). Open these files in any text editor and change WAT to SOL in both files")
				d). Try running Gromacs MDS with the updated files")
