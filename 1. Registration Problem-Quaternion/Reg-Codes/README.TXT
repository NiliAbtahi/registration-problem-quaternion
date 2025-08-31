
   Registration Algorithm in Robotics using quaternion representation of SO(3), special orthogonal group.

   Created by: Nili Abtahi
   Contact: Nili Abtahi (abtahi.niloofar.me@gmail.com)
   -------------------------------------------------------------------------------------------------------
   How to use:
   1. Open Matlab/Octave software and navigate it to the folder containing the m-files.
   2. In Matlab/Octave command window write:
   
            main   firstFile.txt   secondFile.txt
   where the files have a format similar to those in test folder.
   
   3. You can also directly call the data_registration function in command window by 
           data_registration( a_i_3D_points , A_i_3D_points )
   where the arguments of position vectors of 3D points in two coordinate systems.
   
   ------------------------------------------------------------------------------------------------------
   
   Three exmaples provided in test folder:
   
   main   ./test/data1-firstFrame.txt   ./test/data1-secondFrame.txt
   main   ./test/data2-firstFrame.txt   ./test/data2-secondFrame.txt
   main   ./test/data3-firstFrame.txt   ./test/data3-secondFrame.txt
   

