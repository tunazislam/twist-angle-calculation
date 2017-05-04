1) Use putty or any ssh to login either "atria.cs.odu.edu" or "sirius.cs.odu.edu"
2) keep necessary folders in a directory.

	For Image we have this folder "Image_sampling_skip8"
	For True case we have three folders named "true skip_0", "true skip_1", "true skip_2"

3) go to that directory writing "cd" command
   
    cd "type the path  here"
	
   For image each protein folder has two different folders named "good" and "bad":
   For example : For good case of 1d1t protein image
   cd "/home/tislam/Spring2017/Graphs in Bioinformatics/Project/class project/Final Project submission_01038325/Image_sampling_skip8/1dtd/good"

   
   For True case: 
   For example : For true case of 1d1t protein 
   cd "/home/tislam/Spring2017/Graphs in Bioinformatics/Project/class project/Final Project submission_01038325/true skip_0/TrueCase/1dtd"

4) For image we used several input files according to the number of strands 

	For example: s1.pdb_data8.pdb, s2.pdb_data8.pdb,.............................(As we skip 8 points for image)
	
   For true case we used several input files according to the number of strands and skipping points 

	For example: s1.pdb, s2.pdb,............................. (if we skip nothing)
				 s1.pdb_data1.pdb, s2.pdb_data1.pdb,.............................(if we skip 1 point)
				 s1.pdb_data2.pdb, s2.pdb_data2.pdb,.............................(if we skip 2 points)
   
5) To run the code 

	python filename.py
	
	For image :
	 for example: python twist_image_0_30.py 
	For true case :
	 for example: python twist_true_1dtd.py 

6) Output file will be generated in current directory as a name of "twist_angles.txt"
	Inside output files you can see the details of twist angles.
	for example: twist angle between two neighbouring strands, avg. twist angles between two neighbouring strands,
	as well as avg. twist angles of 2 longest strands.