# FInd_Surrounding_Residues
Python Script to find the surrounding residues around "selection" throughout MD simulation trajectory
#
This project started on Sat Oct 15 2015 as part of a "hack day" at the CECAM workshop for "Analysing simulation data".

Requirements:
1) Python: Scientific Python packages

2) MDAnalysis software should be installed

Application:

Counts the % occurance of surrounding residues within a specified cutoff distance from the "selection"

	:Arguments:
            *grofile*
                string of the name of the GROMACS gro file
		'Input File name of gro file (Example: docked_ami_100ns_nowater.gro).'

            *xtcfile*
                trajectory file
		'Input File name of gromacs trajectory file (Example: docked_ami_100ns_nowater.xtc).')
                
            *find around*
                String of the selection atoms around which you want to find residues. Specfy the selections in double quoates (" ")
		(Example: "resname UNK" or "protein" or "resid 458").')
            
            *Selection to find in*
                String of the selection of atoms in which you want to find residues. Specfy the selections in double quoates (" ")
		 (Example: "resname POPC" or "protein" or "resid 1-458").')

             *cutoff*
                cutoff distance from the ligand that is taken into account in angstroms

             *interval*
                Input iteger to specify the Interval between reading frames from MD trajectory

	     *output*
	        'Name of OUTPUT DATA file', default='occupancy.dat')
        
	:Returns:
            *Residues and % Occupancy*
                a dictionary of residues and how many times (%) during the simulation they are within a cutoff value from the ligand


Command:

	python FInd_Surrounding_Residues.py -i traj_at.gro -t traj_at.xtc -fa "resname POPC" -fin "protein" -d 5 -step 10 -o occupancy.dat



