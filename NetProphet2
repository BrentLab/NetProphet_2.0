#!/bin/bash
#SBATCH -D ./
#SBATCH -o LOG/NetProphet2_masterlog_%A.out
#SBATCH -e LOG/NetProphet2_masterlog_%A.err

usage() {
cat << EOF
NetProphet2:
	NetProphet2 is a “data light” algorithm for mapping transcription factor regulatory network.

Usage:
	[Parallel Run]
		sbatch [slurm-option] NetProphet2 -f <config_file>
	[Serial Run]
		./NetProphet2 -s -f <config_file>
Options:
	-s				No value required.Use serial process if SLURM is not available. No value required.
	-f <config_json>		Custom config file.
	--mail-type=<type>		[Slurm-options] Types are NONE, BEGIN, END, FAIL, REQUEUE, ALL. Use "," as delimiter.
	--mail-user=<your_email>	[Slurm-options] Email address to recieve updates.
EOF
}

checkFileParam(){
	if [ -z $2 ]; then
		echo "Option: -${1} is not defined!"
		usage
		exit 1
	fi
	if [ ! -f $2 ]; then
		echo "ERROR: The argument of option -${1} cannot be found: ${2}"
		usage
		exit 1
	fi
}

USE_SERIAL=false
CONFIG_FILE=

while getopts “:hsf:” OPTION
do
	case $OPTION in
		h)
			usage
			exit 0
			;;
		s)
			USE_SERIAL=true
			;;
		f)
			CONFIG_FILE=$OPTARG
			checkFileParam ${OPTION} ${OPTARG}
			;;
		?)
			usage
			exit
			;;
	esac
done

if [[ -z $CONFIG_FILE ]]; then
	echo "ERROR: You must specify arguments for options: -f"
	usage
	exit 1
fi

if $USE_SERIAL; then
	echo -e "\e[1mRunning NetProphet2 in serial fashion.\e[0m"
	snakemake --cores 2 --latency-wait 30 --nolock --configfile $CONFIG_FILE -s run_serial.Snakefile all
else
	echo -e "\e[1mRunning NetProphet2 in parallel fashion.\e[0m"
	module load R/3.2.1
	module load openmpi/1.8.8
	module load python/2.7.15
	module load py-numpy/1.16.2-python-2.7.15
	module load py-scipy/1.2.1-python-2.7.15
	snakemake --cores 2 --latency-wait 172800 --nolock --configfile $CONFIG_FILE -s run_parallel.Snakefile all
	# snakemake --cores 6 --latency-wait 172800 --nolock --configfile $CONFIG_FILE -s run_parallel_split.Snakefile all
fi
