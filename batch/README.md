These scripts retrieve the binary files for a given LOAD and FPGA board, and decode them to produce a ROOT tree for further processing.  You can also run the other executables in this process, eg. adcPlot and laser_exposure.

To submit jobs to the farm via swif, run this command:

swif_decodeMAPMT.py my_workflow load_10_DIRC_075_076 DIRC_075_201_196_197_2018_03_16_16_20

where the 3 arguments are: 

my_workflow
load_name 
fpga_name

where my_workflow can be any string that uniquely labels your work.  load_name = the directory on /cache/halld/detectors/DIRC/MAPMT_test/ which contains the files you want to process and fpga_name is directory for a specific FPGA board withing load_name. 

Some additional commands that are helpful with swif:

swif list    			   = check on the status of your workflows
swif cancel my_workflow -delete    = cancel my_workflow and delete all jobs
