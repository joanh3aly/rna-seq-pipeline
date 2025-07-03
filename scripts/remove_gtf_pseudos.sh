

awk 'BEGIN {OFS="\t"; FS="\t"} {
    line = $0;

    if ($3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS") {
        attrs_str = $9;
        gsub(/"/, "", attrs_str);
        split(attrs_str, attrs_array, "; ");

        gene_name = "";
        gene_type = ""; # Kept for potential future use or safety rule
        is_pkd1p_gene = 0;

        for (i=1; i<=length(attrs_array); i++) {
            attr = attrs_array[i];
            if (attr ~ /gene_name/) {
                gene_name = substr(attr, index(attr, "gene_name ") + 10);
                # THIS IS THE CRUCIAL TRIM:
                gsub(/^[ \t]+|[ \t]+$/, "", gene_name); # Remove leading/trailing spaces/tabs
            }
            if (attr ~ /gene_type/) {
                gene_type = substr(attr, index(attr, "gene_type ") + 10);
                # Trim gene_type too, just in case
                gsub(/^[ \t]+|[ \t]+$/, "", gene_type);
            }
        }

        # DEBUGGING: Print extracted gene_name with delimiters to expose ANY whitespace
        # This output goes to standard error (your terminal), not the output file
        if (gene_name != "") { # Only print if gene_name was found
            print "DEBUG: Extracted gene_name: ['" gene_name "'] (length " length(gene_name) ") for line starting with: " substr($0, 1, 80) "..." > "/dev/stderr";
        }


        if (gene_name ~ /^PKD1P/) {
            print "DEBUG: MATCHED /^PKD1P/ for gene_name: " gene_name > "/dev/stderr"; # Indicate a match
            is_pkd1p_gene = 1;
        }

        # Safety rule: Keep main PKD1 protein_coding gene
        if (gene_name == "PKD1" && gene_type == "protein_coding") {
            print "DEBUG: Keeping PKD1 protein_coding gene: " gene_name " (protein_coding)" > "/dev/stderr";
            is_pkd1p_gene = 0;
        }

        if (is_pkd1p_gene == 1) {
            next; # Skip this line
        }
    }
    print line; # Print non-filtered lines
}' gencode.v48.basic.annotation.gtf > gencode.v48.basic.annotation_no_PKDPs_2.gtf

# awk 'BEGIN {OFS="\t"; FS="\t"} {
#     line = $0;

#     if ($3 == "gene" || $3 == "transcript" || $3 == "exon" || $3 == "CDS") {
#         attrs_str = $9;
#         gsub(/"/, "", attrs_str);
#         split(attrs_str, attrs_array, "; ");

#         gene_name = "";
#         gene_type = "";
#         transcript_type = "";
#         is_pseudogene = 0;

#         for (i=1; i<=length(attrs_array); i++) {
#             attr = attrs_array[i];
#             if (attr ~ /gene_name/) {
#                 gene_name = substr(attr, index(attr, "gene_name ") + 10);
#             } else if (attr ~ /gene_type/) {
#                 gene_type = substr(attr, index(attr, "gene_type ") + 10);
#             } else if (attr ~ /transcript_type/) {
#                 transcript_type = substr(attr, index(attr, "transcript_type ") + 16);
#             }
#         }

#         if (gene_name ~ /^PKD1P/ && (gene_type == "pseudogene" || transcript_type == "lncRNA")) {
#             is_pseudogene = 1;
#         }

#         if (gene_name == "PKD1" && gene_type == "protein_coding") {
#             is_pseudogene = 0;
#         }

#         if (is_pseudogene == 1) {
#             next;
#         }
#     }
#     print line;
# }' gencode.v48.basic.annotation.gtf > gencode.v48.basic.annotation_no_PKDPs.gtf