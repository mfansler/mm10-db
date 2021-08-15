#!/usr/bin/env bash

#########################
# VERIFY RUN CONDITIONS #
#########################

if [ "$1" == "" ]; then
    echo "[ERROR] No input file provided!" >&2
    exit 1
fi

gtf="$1"

if [ ! -f $gtf ]; then
    echo "[ERROR] File not found: $FILE" >&2
    exit 1
fi

# Temp file to track transcripts with missing ends
txids=`mktemp ${TMPDIR}/removeGENCODEmRNA_end_NF.XXXXXX`
if [ $? -ne 0 ]; then
    echo "[ERROR] $0: Can't create temp file, exiting..." >&2
    exit 1
fi

############################
# FIRST PASS: FIND BAD TXS #
############################

# Stream file into AWK
# Output transcript IDs
if [[ $gtf =~ \.gz$ ]]; then
    gzip -cd $gtf
else
    cat $gtf
fi | gawk '
BEGIN {
    FS="\t";
    txs=0;
    print "[INFO] Processing transcripts..." > "/dev/stderr";
}
{
    if ($3 ~ /transcript/ && $9 ~ /tag "mRNA_end_NF"/) {
        txs++;
        match($9, /transcript_id "([^;]+)";/, txid);
	print txid[1];
    }
}
END {
    print "[INFO] Done.\n[INFO] Removing " txs " transcripts." > "/dev/stderr";
}' >> $txids

##########################################
# SECOND PASS: FILTER ASSOCIATED ENTRIES #
##########################################

# Pass saved transcript IDs file as first argument (reads it first)
# Stream full GTF into AWK as second argument, processing on the fly
# Output will go to STDOUT
if [[ $gtf =~ \.gz$ ]]; then
    gzip -cd $gtf
else
    cat $gtf
fi | gawk '
BEGIN {
    FS="\t"
    print "[INFO] Filtering entries..." > "/dev/stderr";
}
(NR==FNR) {  # Load transcript_ids to be removed
    bad_ids[$0]=FNR;
}
(NR!=FNR) {  # Process full file
    if ($9 ~ /transcript_id/) {
        match($9, /transcript_id "([^;]+)";/, txid);
        if (!(txid[1] in bad_ids)) {
     	    print $0;
        }
    } else { 
        print $0;
    }
}
END {
    print "[INFO] Filtering complete." > "/dev/stderr";
}' $txids - 
