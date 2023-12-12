# WGS

## Raw data
- sequencing was performed in Manchester (?) and then transferred to Sanger
- they provided a download script "download_fastqs.py" to automate the transfer

```bash
# downloaded files
KH1_S21_R1_001.fastq.gz
KH1_S21_R2_001.fastq.gz
KH2_S22_R1_001.fastq.gz
KH2_S22_R2_001.fastq.gz
KH3_S23_R1_001.fastq.gz
KH3_S23_R2_001.fastq.gz
KH4_S24_R1_001.fastq.gz
KH4_S24_R2_001.fastq.gz

```
- sample names to experiment names 
    - control_kh1
    - sample_d1_kh2
    - sample_d2_kh3
    - sample_d4_kh4

- to prepare for mapping, made a sample manifest

```bash
# sample manifest - "sample_manifext.txt"
ID,R1,R2
control_kh1,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH1_S21_R1_001.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH1_S21_R2_001.fastq.gz
sample_d1_kh2,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH2_S22_R1_001.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH2_S22_R2_001.fastq.gz
sample_d2_kh3,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH3_S23_R1_001.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH3_S23_R2_001.fastq.gz
sample_d4_kh4,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH4_S24_R1_001.fastq.gz,/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/KH4_S24_R2_001.fastq.gz


```


## Reference
- want to map to both the T. muris genome and the transgenesis plasmid
- ideally, want to find reads that map in both the plasmid and the genome, which will show intergration sites
- however, evidence of sequencing coverage on the plasmid genome alone is pretty good evidence of integration

```bash
# muris genome - from WBPS v18
ln -s ../trichuris_muris.PRJEB126.WBPS18.genomic.fa

# transgenesis plasmid genome from Kelly
transgenesis_plasmid.fa

cat trichuris_muris.PRJEB126.WBPS18.genomic.fa transgenesis_plasmid.fa > tm_genome_plus_plasmid.fa

```

## mapping 
```bash

module load mapping-helminth/v1.0.8

mapping-helminth --input sample_manifext.txt --reference trichuris_muris.PRJEB126.WBPS18.genomic.fa --outdir tmuris_transgenesis_wgs_mapping

# once completed, ran multiqc to collect mapping stats
multiqc .
```
![](../04_analysis/tmuris_transgenesis_multiqc_report.html)


## Coverage
```bash
# collect the bam files
cd /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS

OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/COVERAGE
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'

cd /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS/WGS/COVERAGE

# calculate coverage of the bams
bsub.py 10 cov_stats "~sd21/bash_scripts/run_cov_stats 10000"


# extract mean coverage of the chromosome and plasmid, and calculate the ratio of the two
echo -e "name\tchromosome_cov\tplasmid_cov\tratio" > chromosome_to_plasmid_cov.stat
for i in *.chr.cov; do
     name=${i%.chr.cov};
     chromosome=$(grep "TMUE_LG1\|TMUE_LG2\|TMUE_LG3" ${i} | datamash mean 5 );
     plasmid=$(grep "plasmid" ${i} | datamash mean 5);
     ratio=$(echo "scale=3; ${plasmid}/${chromosome}" | bc);
     echo -e "${name}\t${chromosome}\t${plasmid}\t${ratio}";
done >> chromosome_to_plasmid_cov.stat

```
| name          | chromosome_cov  | plasmid_cov | ratio |
|---------------|-----------------|-------------|-------|
| control_kh1   | 341.187         | 0.486406    | .001  |
| sample_d1_kh2 | 292.493         | 13.3445     | .045  |
| sample_d2_kh3 | 321.883         | 27.6059     | .085  |
| sample_d4_kh4 | 325.63766666667 | 49.4403     | .151  |


# reads with both plasmid and worm sequence
- reads that span the plasmid and the worm will inform integration sites
- could be R1 in plasmid and R2 in worm, and vice versa 
- note that because the integration is random, and the sequence was performed on DNA from a pool, integration sites might be extremely low coverage - even single reads
- also note - the hits between WGS and amplicon experiements might not be concordant - the amplicon sequencing is only sampling a fraction of integration sites that meet the criteria - ie have restriction sites near to integration sites

```bash
# extract lines in bams with both plasmid and TMUE in the line - should pick up paired reads, one in the worm and one in plasmid

for i in *.bam; do
    samtools view ${i} | grep "plasmid" | grep "TMUE" > ${i%bam}.plasmid-worm.hits;     
    done

wc -l *plasmid-worm.hits
#    0 control_kh1..plasmid-worm.hits
#   19 sample_d1_kh2..plasmid-worm.hits
#   24 sample_d2_kh3..plasmid-worm.hits
#   31 sample_d4_kh4..plasmid-worm.hits
```
- very few hits in total - perhaps expected
- none in the control, which is good
- most hits not in the chromosomes, which is unfortunate
- most hits are split to some degree, which might be expected if spanning a integration site, but also might be non-specific


## extract plasmid reads for visualisation
```bash

for i in *bam; do 
    samtools view --bam -o ${i%.bam}.plasmid.bam ${i} plasmid;
    samtools index ${i%.bam}.plasmid.bam;
    done

```
- screenshot over coverage from jBrowse
[](../04_analysis/tmuris_transgenesis_wgs_plasmid-coverage.png)