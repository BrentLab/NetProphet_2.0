#!/bin/bash
LOG_FILE=$1
EXPECTATION=$2
FLAG_FILE=$3

exp_cnt=$( cat $EXPECTATION | sed '/^\s*$/d' | wc -l )

while true; do
	cnt=$( cat $LOG_FILE | sed '/^\s*$/d' | wc -l )
	if (( "$cnt" < "$exp_cnt" )); then
		sleep 5
	else
		touch $FLAG_FILE
		echo "DONE"
		break
	fi
done
