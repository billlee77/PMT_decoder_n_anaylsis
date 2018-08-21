#!/usr/bin/env python

##########################################################################################################################
#
# SWIF DOCUMENTATION:
# https://scicomp.jlab.org/docs/swif
# https://scicomp.jlab.org/docs/swif-cli
# https://scicomp.jlab.org/help/swif/add-job.txt #consider phase!
#
##########################################################################################################################

from optparse import OptionParser
import os.path
import os
import sys
import re
import subprocess
import glob

#################################################### GLOBAL VARIABLES ####################################################

# DEBUG
VERBOSE    = True

# PROJECT INFO
PROJECT    = "gluex"          # http://scicomp.jlab.org/scicomp/#/projects
TRACK      = "analysis"		   # https://scicomp.jlab.org/docs/batch_job_tracks

# RESOURCES
NCORES     = "1"               # Number of CPU cores
DISK       = "10GB"            # Max Disk usage
RAM        = "2GB"             # Max RAM usage
TIMELIMIT  = "1440minutes"     # Max walltime
OS         = "centos7"         # Specify CentOS7 machines

# SOURCE DATA INFORMATION
DATA_SOURCE_TYPE      = "file" #"mss" for tape files, "file" for disk files
DATA_SOURCE_BASE_DIR  = "/volatile/halld/home/jrsteven/2018-mapmt/decodeMAPMT/"

# OUTPUT DATA LOCATION
DATA_OUTPUT_BASE_DIR    = "/volatile/halld/home/%s/2018-mapmt/fitMAPMT/"%(os.environ['USER'])   ## CHANGE IF YOU WANT TO

# JOB EXECUTION
SCRIPTFILE        = "/work/halld2/home/jrsteven/2018-mapmt/PMT_decoder_n_anaylsis/batch/script_fitMAPMT.sh"
ENVFILE           = "/work/halld2/home/jrsteven/2017-spring-ana/setup.csh"

####################################################### FIND FILES #######################################################

def find_files(DATA_SOURCE_DIR):

	# CHANGE TO THE DIRECTORY CONTAINING THE INPUT FILES
	current_dir = os.getcwd()
        if not os.path.isdir(DATA_SOURCE_DIR):
                return []
	os.chdir(DATA_SOURCE_DIR)

	# SEARCH FOR THE FILES
	file_signature = "adc_plots*.hist.root"
	file_list = glob.glob(file_signature)
	if(VERBOSE == True):
		print "size of file_list is " + str(len(file_list))

	# CHANGE BACK TO THE PREVIOUS DIRECTORY
	os.chdir(current_dir)
	return file_list

######################################################## ADD JOB #########################################################

def add_job(WORKFLOW, DATA_SOURCE_DIR, FILENAME, LOAD, FPGA, FILENO):

	# PREPARE NAMES
	STUBNAME = LOAD[0:7] + "_" + FPGA[0:8] + "_" + FILENO
	JOBNAME = WORKFLOW + "_" + STUBNAME

	# CREATE ADD-JOB COMMAND
	# job
	add_command = "swif add-job -workflow " + WORKFLOW + " -name " + JOBNAME
	# project/track
	add_command += " -project " + PROJECT + " -track " + TRACK
	# resources
	add_command += " -cores " + NCORES + " -disk " + DISK + " -ram " + RAM + " -time " + TIMELIMIT + " -os " + OS
	# inputs
	add_command += " -input " + FILENAME + " " + DATA_SOURCE_TYPE + ":" + DATA_SOURCE_DIR + "/" + FILENAME
	# stdout
	add_command += " -stdout " + DATA_OUTPUT_BASE_DIR + "/" + LOAD + "/log/" + "stdout." + STUBNAME + ".out"
	# stderr
	add_command += " -stderr " + DATA_OUTPUT_BASE_DIR + "/" + LOAD + "/log/" + "stderr." + STUBNAME + ".err"
	# command
	add_command += " " + SCRIPTFILE + " " + ENVFILE + " " + FILENAME + " " + DATA_OUTPUT_BASE_DIR  + " " + LOAD + " " + FPGA + " " + FILENO

	if(VERBOSE == True):
		print "job add command is \n" + str(add_command)

	# ADD JOB
	status = subprocess.call(add_command.split(" "))


########################################################## MAIN ##########################################################
	
def main(argv):
	parser_usage = "swif_fitMAPMT.py workflow load_name fpga_name"
	parser = OptionParser(usage = parser_usage)
	(options, args) = parser.parse_args(argv)

	if(len(args) != 3):
		parser.print_help()
		return

	# GET ARGUMENTS
	WORKFLOW = args[0]
	LOADNAME = args[1]
	FPGANAME = args[2]

	# CREATE WORKFLOW IF IT DOESN'T ALREADY EXIST
	WORKFLOW_LIST = subprocess.check_output(["swif", "list"])
	if WORKFLOW not in WORKFLOW_LIST:
	    status = subprocess.call(["swif", "create", "-workflow", WORKFLOW])

	# FIND/ADD JOBS
	DATA_SOURCE_DIR = DATA_SOURCE_BASE_DIR + "/" + LOADNAME + "/" + FPGANAME
	print DATA_SOURCE_DIR
	if(os.path.exists(DATA_SOURCE_DIR)):
		file_list = find_files(DATA_SOURCE_DIR)
	else:
		exit
	if(len(file_list) == 0):
		exit

	# Add jobs to workflow
	for FILENAME in file_list:
		FILENO = FILENAME[10:16]

		# only fit one configuration for now
		if int(FILENO) == 1:
			print FILENAME
			print FILENO

			add_job(WORKFLOW, DATA_SOURCE_DIR, FILENAME, LOADNAME, FPGANAME, FILENO)

if __name__ == "__main__":
   main(sys.argv[1:])

