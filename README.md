# xQTLProcessing

Theorized pipeline:

1. bwa-mem > generate bam files (include raw FASTQ for each of the pooled samples + pB DSPR pop) <<< make.bams.txt - ref >>>
````
bwa mem -t 8 -M $ref raw/${R1} raw/${R2} | samtools view -bS - > $dir1/$shortname.temp.bam
samtools sort $dir1/$shortname.temp.bam -o $dir1/$shortname.bam
samtools index $dir1/$shortname.bam
````

2. bcftoolsmpileup > combined all of the pooled samples (control and treatment) as well as founder samples <<< call.SNPs.txt - ref >>>
````
module load samtools/1.9
module load bcftools/1.9
ref="/share/adl/tdlong/DSPR/DNAseq/ref/dm6.fa"
declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SGE_TASK_ID - 1]}
# files.txt is a list of bam files to combine
# it should include all the pooled samples AND the founders as BAMs!
bcftools mpileup -I -d 1000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b helperfile/files.txt | bcftools call -mv -Ob > July24/calls.$mychr.bcf  
# frequencies @ SNPs by sample
bcftools query -e'GT ="./."'  -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' July24/calls.$mychr.bcf | grep -v '\.' | perl scripts/accuracy.freqtab.pl >July24/temp.$mychr.txt
# counts @ SNPs by sample
bcftools query -e'GT ="./."'  -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' July24/calls.$mychr.bcf | grep -v '\.' | perl scripts/accuracy.counttab.pl >July24/temp.count.$mychr.txt
####

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

