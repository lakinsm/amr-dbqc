#!/bin/bash

## This script analyses the sensitivity of database clustering and multiple alignment
## based on varying levels of sequence identity with cd-hit.  It is being used in
## the development of the Hidden Markov Model for the Microbial Ecology Group.

## Author: Steven Lakin
## Email: Steven.Lakin@colostate.edu

##########
## Vars ##
##########
infile=""
outbase=""
annot=""
check=0


#############
## Methods ##
#############
display_help() {
    echo "
        Usage: sensitivity.sh [options]

            -h | --help		Display this menu

        Input Options

	    -i | --infile STR	Input database fasta file			
            -o | --outbase STR	Base name for output file
            -a | --annotations STR	Input annotation csv file
            -t | --threads INT	Number of threads to use

    "
}


##########
## Main ##
##########
while [[ "${1+defined}" ]]; do
    case "$1" in
		-i | --infile )
			infile=$2
			shift 2
			;;
		-a | --annotations )
			annot=$2
			shift 2
			;;
        -o | --outbase )
            outbase=$2
            shift 2
            ;;
        -h | --help )
            display_help
            exit 0
            ;;
        -t | --threads )
            threads=$2
            shift 2
            ;;
        --)  # End of options
            shift 1
            break
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
        *)  # No more options
            break
            ;;
    esac
done


for sim in $(seq 0.8 0.02 1.0); do
	outfile=${outbase}_${sim}
    if (( $(echo "0.80 <= $sim && $sim < 0.85" | bc -l) )); then
		cd-hit-est -i ${infile} -o ${outfile} -c ${sim} -g 1 -p 1 -d 10000 -T $threads -M 16000 -n 5 
	elif (( $(echo "0.85 <= $sim && $sim < 0.88" | bc -l) )); then
		cd-hit-est -i ${infile} -o ${outfile} -c ${sim} -g 1 -p 1 -d 10000 -T $threads -M 16000 -n 6 
	elif (( $(echo "0.88 <= $sim && $sim < 0.9" | bc -l) )); then
		cd-hit-est -i ${infile} -o ${outfile} -c ${sim} -g 1 -p 1 -d 10000 -T $threads -M 16000 -n 7 
	elif (( $(echo "0.9 <= $sim && $sim < 0.92" | bc -l) )); then
		cd-hit-est -i ${infile} -o ${outfile} -c ${sim} -g 1 -p 1 -d 10000 -T $threads -M 16000 -n 8 
	elif (( $(echo "0.92 <= $sim && $sim < 0.95" | bc -l) )); then
		cd-hit-est -i ${infile} -o ${outfile} -c ${sim} -g 1 -p 1 -d 10000 -T $threads -M 16000 -n 9 
	elif (( $(echo "0.95 <= $sim" | bc -l) )); then
		cd-hit-est -i ${infile} -o ${outfile} -c ${sim} -g 1 -p 1 -d 10000 -T $threads -M 16000 -n 10 
	fi
	if [ $check == 1 ]; then
		annot=${prevfile}_annotations.csv
	fi
	check=1
	./cdhit_impute_annotations.py ${outfile}.clstr ${annot} ${outfile}_annotations.csv
	prevfile=${outfile}
done













exit 0

