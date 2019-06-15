#!/bin/bash
C_THRD_FLAG=false
C_MTRX_FLAG=false
NO_CHECK_MODE=false
LOG_FLAG=false
AVG_FLAG=false
PRINT_FLAG="false"
PROGRESSBAR_FLAG=false
TIME_PAR_MODE=""
NO_CHECK_MODE=""
E_ARG=55	#usage error
WRONG_CNT=0
SUM_TIME_INIT=0
SUM_TIME_LUDE=0
SUM_TIME_CHCK=0

# Usage info
show_help()
{
	echo "Usage: [-f ./filename] [-H] [-r repeat] <size|file> <threads>
-f	./filename	Executable file name.	
-r	repeat		Repeat execution.
-p			Print result one by one for repeated work.
-a			Calculate average time.
-b			Print progress bar. #not implemented
-l			Log file. #have to use with -a option.
-c			Continuously increase thread count 1 to 16 (matrix size must be less than 1000)
-C			Continuously increase matrix size 10 to 2000 (thread count must be over 10)
[-n]			No check mode.
-H			Help
		
For Example: ./run.sh -f ./lu -r 1000 3000 10
"
exit ${E_ARG}
}

while getopts "cClnatdpHr:f:" OPT
do
	case $OPT in
		H) 
			show_help	;;
		r)
			re="${OPTARG}"	;;
		f)
			executefile="${OPTARG}" ;;
		p)
			PRINT_FLAG="true"
			;;
		b)
			PROGRESSBAR_FLAG=true
			;;
		a)
			AVG_FLAG=true
			;;
		t)
			TIME_PAR_MODE="-t"
			TIME_PAR_FLAG=true
			;;
		n)
			NO_CHECK_MODE="-n" ;;
		\?)
			echo "Invalid option: -$OPTARG. Use -H flag for help."
			exit ${E_ARG}
			;;
	esac
done

#OPTIND shift
shift $(( $OPTIND - ))

arg1=$1			#matrix size : n
arg2=$2			#random seed
arg3=$3			#thread num
arg4=$4			#print matrix


MATRIX_SIZE=${arg1}
RANDOM_SEED=${arg2:-1}
THREADS_NUM=${arg3:-16} 
DEBUG_FLAG=${arg4:-0}

if [ $# = 0 ]
then
	echo "Invalid argument: Use -H flag for help."
	exit $E_ARG
fi

#make use Makefil
make

#anounce
echo "[PROGRAM INFO]
 -------------------------------------------------------------------------------
| program name	| size or file	| threads num	| repeat count	| par time	|
| ${executefile}		| ${MATRIX_SIZE}		| ${THREADS_NUM}threads	| ${re:-1}times	| "false"		|
 -------------------------------------------------------------------------------"

if ${PRINT_FLAG}; then
	echo "
init    lude    check"
fi

RE_CNT=${re:-1}
while [ ${re:-1} -ge 1 ]; do

	${executefile:-./lu} ${MATRIX_SIZE} ${RANDOM_SEED} ${THREADS_NUM} ${DEBUG_FLAG} ${NO_CHECK_MODE}>&1
#	${executefile:-./lu} ${MATRIX_SIZE} ${RANDOM_SEED} ${THREADS_NUM} ${DEBUG_FLAG} ${NO_CHECK_MODE}>temp
#	exec 3<>temp

	if ${PRINT_FLAG}; then
		echo <&3
	fi

	if ${AVG_FLAG}
	then
		read a1<&3
		#a1 to string array
		stringarr=($a1)
		time_init=${stringarr[0]}
		time_lude=${stringarr[1]}
		time_chck=${stringarr[2]}

		SUM_TIME_INIT=`echo $SUM_TIME_INIT + $time_init | bc`
		SUM_TIME_LUDE="$( bc <<<"$SUM_TIME_LUDE + $time_lude" )"
		SUM_TIME_CHCK="$( bc <<<"$SUM_TIME_CHCK + $time_chck" )"
	fi

#	len=${#a1}
#	correct=${a1:$len-1:$len-0}
#	if [ $correct != 0 ]
#	then
#		((WRONG_CNT++))
#	fi

	re=$(( ${re}-1 ))
done


#It is correct when wrong count is 0.
if [ ${WRONG_CNT} == 0 ]
then
	echo "
 -----------------------------------------------
|	      C  O  R  R  E  C  T		|
 -----------------------------------------------"
else
	echo "
 -----------------------------------------------
|	W R O N G ${WRONG_CNT}times	|
 -----------------------------------------------"
fi

if ${AVG_FLAG}
then
	AVG_TIME_INIT="$( bc <<<"scale=5;$SUM_TIME_INIT/$RE_CNT" )"
	AVG_TIME_LUDE="$( bc <<<"scale=5;$SUM_TIME_LUDE/$RE_CNT" )"
	AVG_TIME_CHCK="$( bc <<<"scale=5;$SUM_TIME_CHCK/$RE_CNT" )"
	printf "| INIT: %.5f  LUDE: %.5f  CHCK: %.5f\t|\n" "${AVG_TIME_INIT}" "${AVG_TIME_LUDE}" "${AVG_TIME_CHCK}"
	echo " -----------------------------------------------"
	echo "INIT: ${AVG_TIME_INIT}	LUDE: ${AVG_TIME_LUDE}	CHCK: ${AVG_TIME_CHCK}			SIZE: ${MATRIX_SIZE}		TRD: ${THREADS_NUM}		REP: ${RE_CNT}">>log.txt
fi
