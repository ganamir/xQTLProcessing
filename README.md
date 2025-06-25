# xQTLProcessing

Theorized pipeline:
1. Trim files:
````
#!/bin/bash

DIR="/mnt/d/xQTL_2025_Data/23152-05-06162025_145113-457684232/BaseSpace_CLI_2025-06-16_17_37_50Z->
# Iterate over all R1 fastq.gz files in the directory
for R1 in "$DIR"/*_R1_001.fastq.gz; do
    # Extract base filename without path and suffix
    BASE=$(basename "$R1")
    BASE=${BASE%%_R1_001.fastq.gz}

    # Construct full paths to R1 and R2
    R1_FILE="$DIR/${BASE}_R1_001.fastq.gz"
    R2_FILE="$DIR/${BASE}_R2_001.fastq.gz"

    # Check if both files exist
    if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
        echo "Trimming: $BASE"
        trim_galore --paired --length 40 --max_n 1 -q 20 -j 8 -o "/mnt/d/xQTL_2025_Data/xQTL_Data>    else
        echo "Missing file(s) for $BASE"
    fi
done
````

1. bwa-mem > generate bam files (include raw FASTQ for each of the pooled samples + pB DSPR pop) <<< make.bams.txt - ref >>>
````
# Directory containing .gz files
input_dir="/mnt/d/xQTL_2025_Data/xQTL_Data_Real_Samples"
output_dir="$input_dir/bams"

mkdir -p "$output_dir"

# Reference genome
reference_genome="/mnt/d/xQTL_2025_Data/ref/dm6.fa"

# Number of threads for BWA
threads=50

# Iterate over all R1 files in the specified directory
for file1 in "$input_dir"/*_R1_*_val_1.fq.gz; do
    # Construct the expected R2 filename
    file2="${file1/_R1_/_R2_}"
    file2="${file2/_val_1/_val_2}"
    echo "$file1"
    echo "$file2"

    # Check if the R2 file exists
    if [ -f "$file2" ]; then
        # Output filename
        output_file1=$(basename "${file1%_R1_*.fq.gz}_temp_aligned.bam")
        output_file2=$(basename "${file1%_R1_*.fq.gz}_aligned.bam")

        bwa mem "$reference_genome" -M -t "$threads" -v 3 "$file1" "$file2" | samtools view -bS - > "$output_dir/$outpu>        samtools sort "$output_dir/$output_file1" -o "$output_dir/$output_file2"
        samtools index "$output_dir/$output_file2"
    else
        echo "Corresponding R2 file for $file1 not found."
    fi
done
````

2. bcftoolsmpileup > combined all of the pooled samples (control and treatment) as well as founder samples <<< call.SNPs.txt - ref >>>
````
#!/bin/bash

# Reference genome
ref="/mnt/d/xQTL_2025_Data/Practice_data/ref/dm6.fa"

# Chromosome list
declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")

# Directory with BAM files
bam_dir="/mnt/d/xQTL_2025_Data/Practice_data/TrimOut/bwaAligned"

# Output directory
outdir="accuracyhap"
mkdir -p "$outdir"

# Collect all BAMs into a list
bam_files=("$bam_dir"/*.bam)

# Loop through each chromosome
for mychr in "${chrs[@]}"; do
    echo "Processing all BAMs for chromosome $mychr"

    bcf_out="$outdir/combined.${mychr}.bcf"
    echo "Writing BCF to $bcf_out"

    # Run bcftools mpileup + call for all BAMs at once
    bcftools mpileup -v 3 --threads 50 -I -d 1000 -r "$mychr" -a "FORMAT/AD,FORMAT/DP" -f "$ref" "${bam_files[@]}" |
    bcftools call -mv -Ob -o "$bcf_out"

    # Frequency table
    bcftools query -e 'GT="./." || QUAL<20 || FORMAT/DP<20 || FORMAT/DP>200' \
      -f '%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' "$bcf_out" |
    grep -v '\.' |
    perl scripts/accuracy.freqtab.pl > "$outdir/combined.temp.${mychr}.txt"

    # Count table
    bcftools query -e 'GT="./." || QUAL<20 || FORMAT/DP<20 || FORMAT/DP>200' \
      -f '%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' "$bcf_out" |
    grep -v '\.' |
    perl scripts/accuracy.counttab.pl > "$outdir/combined.temp.count.${mychr}.txt"

done

````

3. Assemble a SNP table: <<< call.SNPs.txt - ref >>>
````
#####   SNP tables....
###  concatenate the chromosomes...and build SNP tables

# real data
echo -ne "CHROM\tPOS\tfreq_A1\tfreq_A2\tfreq_A3\tfreq_A4\tfreq_A5\tfreq_A6\tfreq_A7\tfreq_AB8\tfreq_B1\tfreq_B2\tfreq_B3\tfreq_B4\tfreq_B5\tfreq_B6\tfreq_B7\t" > July24/SNP.accuracy.freq.txt
echo -ne "freq_A.T.100\tfreq_B.T.100\tfreq_C.T.100\tfreq_D.T.100\tfreq_A.C.100\tfreq_B.C.100\tfreq_C.C.100\tfreq_D.C.100\tfreq_A.T.25\tfreq_B.T.25\tfreq_C.T.25\tfreq_D.T.25\tfreq_A.C.25\tfreq_B.C.25\tfreq_C.C.25\tfreq_D.C.25\tfreq_A.T.10\tfreq_B.T.10\tfreq_C.T.10\tfreq_D.T.10\tfreq_A.C.10\tfreq_B.C.10\tfreq_C.C.10\tfreq_D.C.10\n" >> July24/SNP.accuracy.freq.txt
cat July24/temp.chrX.txt >> July24/SNP.accuracy.freq.txt
cat July24/temp.chr2L.txt >> July24/SNP.accuracy.freq.txt
cat July24/temp.chr2R.txt >> July24/SNP.accuracy.freq.txt
cat July24/temp.chr3L.txt >> July24/SNP.accuracy.freq.txt
cat July24/temp.chr3R.txt >> July24/SNP.accuracy.freq.txt

echo -ne "CHROM\tPOS\tN_A1\tN_A2\tN_A3\tN_A4\tN_A5\tN_A6\tN_A7\tN_AB8\tN_B1\tN_B2\tN_B3\tN_B4\tN_B5\tN_B6\tN_B7\t" > July24/SNP.accuracy.N.txt
echo -ne "N_A.T.100\tN_B.T.100\tN_C.T.100\tN_D.T.100\tN_A.C.100\tN_B.C.100\tN_C.C.100\tN_D.C.100\tN_A.T.25\tN_B.T.25\tN_C.T.25\tN_D.T.25\tN_A.C.25\tN_B.C.25\tN_C.C.25\tN_D.C.25\tN_A.T.10\tN_B.T.10\tN_C.T.10\tN_D.T.10\tN_A.C.10\tN_B.C.10\tN_C.C.10\tN_D.C.10\n" >> July24/SNP.accuracy.N.txt
cat July24/temp.count.chrX.txt >> July24/SNP.accuracy.N.txt
cat July24/temp.count.chr2L.txt >> July24/SNP.accuracy.N.txt
cat July24/temp.count.chr2R.txt >> July24/SNP.accuracy.N.txt
cat July24/temp.count.chr3L.txt >> July24/SNP.accuracy.N.txt
cat July24/temp.count.chr3R.txt >> July24/SNP.accuracy.N.txt

rm July24/temp.count.chr*
rm July24/temp.chr*
````

4. Haplotype calling: <<< call.haplotypes.txt - ref >>>
````
#######################
#### call haplotypes...
#######################

# haplotype call real data & pool data
qsub -t 1 scripts/haplotyper.justhapfreq.sh helperfile/justhaps.samples.txt July24 July24/SNP.accuracy.freq.txt helperfile/founders.txt
qsub -t 2-24 scripts/haplotyper.justhapfreq.sh helperfile/justhaps.samples.txt July24 July24/SNP.accuracy.freq.txt helperfile/founders.txt
qsub -t 1-12 scripts/haplotyper.justhapfreq.sh helperfile/haplotype.samples.pool.experiment.txt accuracyhap accuracyhap/SNP.accuracy.freq.txt helperfile/founders.txt

# merge real data & pool data
#read
cat July24/B.T.T_hap_freq.txt | head -n 1 > July24/allhaps.txt
awk FNR-1 July24/*_hap_freq.txt >> July24/allhaps.txt
cat July24/allhaps.txt | gzip -c > allhaps.200kb.txt.gz

#pool
cat accuracyhap/A_hap_freq.txt | head -n 1 > accuracyhap/accuracyhap.txt
awk FNR-1 accuracyhap/*_hap_freq.txt >> accuracyhap/accuracyhap.txt
cat accuracyhap/accuracyhap.txt | gzip -c > accuracyhap.200kb.txt.gz

###  testing
pool="A.T.100"
OutName="A.T.F"
folder="newhap"
SNPtable="SNP.freq.txt"
foundernames="founders.txt"

library(limSolve)
source("scripts/haplotyper.justhapfreq.code.R")
runscan(pool, OutName, folder, SNPtable, foundername)
### end testing

######################
## done calling haplotypes
######################

allhaps.txt has the following columns
1.  chr
2.	pos
3.	pool name
4.	founder name
5.	adjusted haplotype frequency estimate
###############
````

