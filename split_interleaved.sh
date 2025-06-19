#!/bin/bash

# --- Configuration ---
# Set your input gzipped FASTQ file here
INPUT_FASTQ="SRR10488341.fastq.gz"

# --- Main script ---
if [ -z "$INPUT_FASTQ" ]; then
    echo "Error: INPUT_FASTQ variable is not set. Please specify your input file."
    exit 1
fi

if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input file '$INPUT_FASTQ' not found."
    exit 1
fi

echo "Splitting '$INPUT_FASTQ' into forward_reads.fastq and reverse_reads.fastq..."

# Decompress the gzipped file and pipe its content to awk
gunzip -c "$INPUT_FASTQ" | awk -f /dev/stdin << 'AWK_SCRIPT_END'
# Each FASTQ read is 4 lines: ID, Sequence, +, Quality
# We're interested in the ID line (line 1 of each 4-line block)
if (NR % 4 == 1) {
    # Check if the ID ends with /1 (forward read)
    if ($1 ~ /\/1$/) {
        print > "forward_reads.fastq"
        read_type = "forward"
    }
    # Check if the ID ends with /2 (reverse read)
    else if ($1 ~ /\/2$/) {
        print > "reverse_reads.fastq"
        read_type = "reverse"
    } else {
        # Fallback for other naming conventions or errors
        # You might need to adjust the regex for your specific IDs
        print "Warning: Unrecognized read ID format: " $1 > "/dev/stderr"
        # As a fallback, you could try to assume interleaving based on NR
        # if ( (NR-1)/4 % 2 == 0) { print > "forward_reads.fastq"; read_type = "forward" }
        # else { print > "reverse_reads.fastq"; read_type = "reverse" }
    }
}
# For lines 2, 3, 4, print to the file determined by the ID line
else {
    if (read_type == "forward") {
        print > "forward_reads.fastq"
    } else if (read_type == "reverse") {
        print > "reverse_reads.fastq"
    }
}
AWK_SCRIPT_END

echo "Splitting complete."