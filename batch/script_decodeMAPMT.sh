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

# RUN DECODER
/work/halld2/home/jrsteven/2018-mapmt/PMT_decoder_n_anaylsis/decodeFE $INPUTFILE
/work/halld2/home/jrsteven/2018-mapmt/PMT_decoder_n_anaylsis/adcPlot run_${FILENO}.bin.root
/work/halld2/home/jrsteven/2018-mapmt/PMT_decoder_n_anaylsis/laserExposure run_${FILENO}.bin.root

ls -al

# RETURN CODE
set RETURN_CODE = $?
echo "Return Code = " $RETURN_CODE
if ($RETURN_CODE != 0) then
	exit $RETURN_CODE
endif

# save output tree
mkdir -p -m 775 ${OUTDIR}/${LOAD}/${FPGA}/
if (-e run_${FILENO}.bin.root) then
	cp -v run_${FILENO}.bin.root ${OUTDIR}/${LOAD}/${FPGA}/run_${FILENO}.bin.root
	chmod 664 ${OUTDIR}/${LOAD}/${FPGA}/run_${FILENO}.bin.root
endif

# save output histograms
mkdir -p -m 775 ${OUTDIR}/${LOAD}/${FPGA}/
if (-e adc_plots.hist.root) then
        cp -v adc_plots.hist.root ${OUTDIR}/${LOAD}/${FPGA}/adc_plots_${FILENO}.hist.root
        chmod 664 ${OUTDIR}/${LOAD}/${FPGA}/adc_plots_${FILENO}.hist.root
endif

# save output images
mkdir -p -m 775 ${OUTDIR}/plots
if (-e adc_plots.pdf) then
        cp -v adc_plots.pdf ${OUTDIR}/plots/adc_plots_${LOAD}_${FPGA}_${FILENO}.pdf
        chmod 664 ${OUTDIR}/plots/adc_plots_${LOAD}_${FPGA}_${FILENO}.pdf
endif

mkdir -p -m 775 ${OUTDIR}/plots
if (-e laser_exposure.pdf) then
        cp -v laser_exposure.pdf ${OUTDIR}/plots/laser_exposure_${LOAD}_${FPGA}_${FILENO}.pdf
        chmod 664 ${OUTDIR}/plots/laser_exposure_${LOAD}_${FPGA}_${FILENO}.pdf
endif

