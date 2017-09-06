#!/usr/bin/env bash

HELP=NO
############################################################################################################
if [[ $# -eq 0 ]]; then
  HELP=YES
fi

PAIRS=""
while [[ $# -gt 0 ]]
do
    key="$1"
    case ${key} in
        -h|--help)
        HELP=YES
        shift
        ;;

        *) # Last arg as session name
        WORK_DIR="$2"
        PAIRS="${PAIRS};${key}=${WORK_DIR}"
        shift
        shift
        ;;
    esac
done

if [ ${HELP} == "YES" ]; then
    echo "For each given SRX_ID_i and WORK_DIR_i pair this script downloads all *.sra from
rsync://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRXnnn/SRXnnnnn
to given directory.

Usage:
  geo_rsync.sh SRX_ID_1 WORK_DIR_1 [SRX_ID_2 WORK_DIR_2 ..]"
    exit 0
fi

IFS=';' list=(${PAIRS})
for item in "${list[@]}"; do
    if [[ ! -z ${item} ]]; then
        WORK_DIR="${item#*=}"
        SRXID="${item%=*}"

        URL="rsync://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/${SRXID:0:6}/${SRXID}"
        echo "Downloading: ${URL} to ${WORK_DIR}"
        # Options:
        #   -a, --archive            archive mode; equals -rlptgoD (no -H,-A,-X)
        #   -z, --compress           compress file data during the transfer
        #   -v, --verbose            increase verbosity
        #   -i,  --itemize-changes   output a change-summary for all updates
        #   -t, --times              preserve modification times
        #   —partial (or -P = —progress —partial): enable partial transmission,
        #                            but not for local fs
        #   -r, --recursive          recurse into directories
        #   -h, --human-readable     output numbers in a human-readable format
        rsync -azvvit --partial-dir=.rsync-partial --human-readable --progress ${URL} ${WORK_DIR}
    fi
done