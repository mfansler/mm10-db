#!/usr/bin/gawk -f
#' Extract TX -> GENE map from GTF, including chromosome
BEGIN {
    FS  = "\t"
    OFS = "\t"
}
{
    if ($3 ~ /transcript/) {
	match($9, /gene_id "(ENSMUSG[0-9]+\.[0-9]+)";.*transcript_id "([^"]+\.[0-9]+)"/, ids);
	match($9, /gene_name "([^;]+)";/, symbol);
	print ids[2], ids[1], symbol[1], $1;
    }
}
