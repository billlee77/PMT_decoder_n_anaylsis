#!/bin/csh -f

# SET INPUTS
setenv ENVIRONMENT $1
setenv INPUTFILE $2
setenv OUTDIR $3
setenv LOAD $4
setenv FPGA $5
setenv FILENO $6

# PRINT INPUTS
echo "ENVIRONMENT       = $ENVIRONMENT"
echo "INPUTFILE         = $INPUTFILE"
echo "OUTDIR            = $OUTDIR"
echo "LOAD              = $LOAD"
echo "FPGA              = $FPGA"
echo "FILENO            = $FILENO"

# ENVIRONMENT
source $ENVIRONMENT
echo pwd = $PWD
printenv

# COPY INPUT FILE TO WORKING DIRECTORY
# This step is necessary since the cache files will be created as soft links in the current directory, and we want to avoid large I/O processes.
# We first copy the input file to the current directory, then remove the link.
ls -al
cp $INPUTFILE ./tmp_file
rm -f $INPUTFILE
mv tmp_file $INPUTFILE
ls -al

# RUN FITTER
root -b -q $INPUTFILE

ls -al

# RETURN CODE
set RETURN_CODE = $?
echo "Return Code = " $RETURN_CODE
if ($RETURN_CODE != 0) then
	exit $RETURN_CODE
endif

# save output histograms
mkdir -p -m 775 ${OUTDIR}/${LOAD}/${FPGA}/
if (-e fit_plots.hist.root) then
        cp -v fit_plots.hist.root ${OUTDIR}/${LOAD}/${FPGA}/fit_plots_${FILENO}.hist.root
        chmod 664 ${OUTDIR}/${LOAD}/${FPGA}/fit_plots_${FILENO}.hist.root
endif

