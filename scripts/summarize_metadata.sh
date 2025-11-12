#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 comprehensive_metadata.tsv"
    exit 1
fi

META_FILE="$1"

##############################
# 0. Make paired-only subset #
##############################
# We assume:
#   col2  = study_accession
#   col16 = library_layout
#   col17 = read_count
#
# We keep:
#   - header
#   - rows where column16 (library_layout) is PAIRED (case-insensitive)

echo "Generating comprehensive_metadata_paired.tsv ..."

awk -F'\t' '
NR==1 { print; next }
{
    layout = toupper($16);
    if (layout == "PAIRED")
        print;
}
' "$META_FILE" > comprehensive_metadata_paired.tsv


##############################################
# All further analysis will use ONLY paired  #
# i.e. comprehensive_metadata_paired.tsv     #
##############################################

PAIRED_FILE="comprehensive_metadata_paired.tsv"

### 1. Generate per-study read depth stats (paired only)
### ------------------------------------------------
# Output: per_study_readcount_stats.tsv
# Columns:
# Study  N_samples  Mean_Reads  Median_Reads  SD_Reads  Min_Reads  Max_Reads  CV
#
# We do robust numeric handling to avoid bogus zeros.

awk -F'\t' '
function human(x){
    if (x>=1e9) return sprintf("%.2fG", x/1e9);
    else if (x>=1e6) return sprintf("%.2fM", x/1e6);
    else if (x>=1e3) return sprintf("%.2fK", x/1e3);
    else return x;
}

NR>1 && $17 ~ /^[0-9]+$/ {
    study = $2;
    i = ++count[study];
    vals[study, i] = $17 + 0;   # numeric
}

END {
    print "Study\tN_samples\tMean_Reads\tMedian_Reads\tSD_Reads\tMin_Reads\tMax_Reads\tCV" > "per_study_readcount_stats.tsv.tmp";

    for (s in count) {
        n = count[s];

        # collect and compute basic stats
        sum = 0;
        for (i=1; i<=n; i++) {
            v = vals[s, i];
            arr[i] = v;
            sum += v;
        }

        mean = sum / n;

        # standard deviation (population SD)
        sd_sum = 0;
        for (i=1; i<=n; i++) {
            v = vals[s, i];
            sd_sum += (v - mean) * (v - mean);
        }
        sd = sqrt(sd_sum / n);

        # sort for median/min/max
        asort(arr);  # arr now sorted ascending, 1..n
        if (n % 2) {
            med = arr[(n+1)/2];
        } else {
            med = (arr[n/2] + arr[n/2 + 1]) / 2;
        }
        minv = arr[1];
        maxv = arr[n];

        cv = (mean > 0 ? sd/mean : 0);

        printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%.3f\n", \
            s, n, human(mean), human(med), human(sd), human(minv), human(maxv), cv \
            >> "per_study_readcount_stats.tsv.tmp";
    }
}
' "$PAIRED_FILE"

# sort by N_samples desc, keep header first
{
    head -n1 per_study_readcount_stats.tsv.tmp
    tail -n +2 per_study_readcount_stats.tsv.tmp | sort -k2,2nr
} > per_study_readcount_stats.tsv

rm per_study_readcount_stats.tsv.tmp


### 2. Generate per-study paired/single layout counts
### ------------------------------------------------
# Output: per_study_layout_counts.tsv
#
# BUT: we are now operating on ONLY PAIRED rows.
# So Single_Count should always be 0.
#
# Still, we’ll compute both for clarity / sanity check.

awk -F'\t' '
NR>1 {
    study = $2;
    layout = toupper($16);

    if (layout == "PAIRED")
        paired[study]++;

    if (layout == "SINGLE" || layout == "SINGLE-END" || layout == "SINGLEEND")
        single[study]++;
}

END {
    print "Study\tPaired_Count\tSingle_Count\tTotal_Libraries" > "per_study_layout_counts.tsv.tmp";

    # union of all studies that appeared at all
    for (s in paired) seen[s]=1;
    for (s in single) seen[s]=1;

    for (s in seen) {
        p  = (s in paired ? paired[s] : 0);
        si = (s in single ? single[s] : 0);
        total = p + si;
        printf "%s\t%d\t%d\t%d\n", s, p, si, total >> "per_study_layout_counts.tsv.tmp";
    }
}
' "$PAIRED_FILE"

# sort by Total_Libraries desc, keep header first
{
    head -n1 per_study_layout_counts.tsv.tmp
    tail -n +2 per_study_layout_counts.tsv.tmp | sort -k4,4nr
} > per_study_layout_counts.tsv

rm per_study_layout_counts.tsv.tmp


### 3. Join both tables into one summary
### ------------------------------------------------
# Output: per_study_summary.tsv
# Columns:
# Study  N_samples  Mean_Reads  Median_Reads  SD_Reads  Min_Reads  Max_Reads  CV  Paired_Count  Single_Count  Total_Libraries
#
# We use only PAIRed data (PAIRED_FILE), so this is already paired-only view.

# Prep, strip headers, sort by Study for join
tail -n +2 per_study_readcount_stats.tsv | sort -k1,1 > _stats_nh.tsv
tail -n +2 per_study_layout_counts.tsv | sort -k1,1 > _layout_nh.tsv

# Join on Study (field 1). If some study somehow has stats but no layout row (shouldn't
# happen, but just in case), fill missing with 0.
join -t $'\t' -a1 -e 0 \
     -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2,2.3,2.4 \
     _stats_nh.tsv _layout_nh.tsv \
     | sort -k2,2nr \
     > _summary_body.tsv

echo -e "Study\tN_samples\tMean_Reads\tMedian_Reads\tSD_Reads\tMin_Reads\tMax_Reads\tCV\tPaired_Count\tSingle_Count\tTotal_Libraries" \
    > per_study_summary.tsv
cat _summary_body.tsv >> per_study_summary.tsv

rm _stats_nh.tsv _layout_nh.tsv _summary_body.tsv

echo "Done."
echo "Created:"
echo " - comprehensive_metadata_paired.tsv   (filtered input: PAIRED only)"
echo " - per_study_readcount_stats.tsv       (paired only)"
echo " - per_study_layout_counts.tsv         (paired only view; Single_Count should be 0)"
echo " - per_study_summary.tsv               (paired only summary)"