# Transgenesis of Trichuris muris

- in collaboration with Richard Grencis' group at the University fo Manchester


- Background:
    - protocol and informatics processing: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4829050/
    - sequencing data: https://www.ebi.ac.uk/ena/browser/view/PRJEB53152

- Samples
    - 6881STDY12870367	sample1 - D1	ERR10822048
    - 6881STDY12870368	sample2 - L1	ERR10822049
    - 6881STDY12870369	sample3 -Tm	ERR10822050
    - 6881STDY12870370	sample4 - t	ERR10822051
    - 6881STDY12870371	sample5 - D2	ERR10822052
    - 6881STDY12870372	sample6 - D4	ERR10822053

- where:
    - S1: Plasmid transfected L1
    - S2: Plasmid transfected cells (mammalian) – positive control for integration
    - S3: Trichuris DNA – negative control
    - S4: 2nd plasmid transfected ells (mammalian) – positive control for integration
    - S5: 2nd plasmid transfected L1
    - S6: 3rd plasmid transfected L1
    - Samples 1, 5,6 are the test samples, all done on different days and DNA taken form lots of worms (couldn’t pick out singly positive worms) and all samples positive for plasmid by PCR. As there are multiple worms and this was done multiple times it is likely that there will be multiple integration sites into the genome (fingers crossed the process worked) as lentivirus is randomly integrated.

### get data and references
```bash
# download data from ENA using wget
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/048/ERR10822048/ERR10822048_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/048/ERR10822048/ERR10822048_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/049/ERR10822049/ERR10822049_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/049/ERR10822049/ERR10822049_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/050/ERR10822050/ERR10822050_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/050/ERR10822050/ERR10822050_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/051/ERR10822051/ERR10822051_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/051/ERR10822051/ERR10822051_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/052/ERR10822052/ERR10822052_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/052/ERR10822052/ERR10822052_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/053/ERR10822053/ERR10822053_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR108/053/ERR10822053/ERR10822053_2.fastq.gz


# get trichuris muris genome assembly from wormbase parasite
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/trichuris_muris/PRJEB126/trichuris_muris.PRJEB126.WBPS18.genomic.fa.gz
gunzip trichuris_muris.PRJEB126.WBPS18.genomic.fa.gz

```


### mapping reads to the reference genome
```bash
# map samples to the reference genome
ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//g' | awk '{print $1,$1}' OFS="\t" > lanes_samples.list

./run_mapping.sh

```

tgtgactctggtaactagagatccctc


zcat ERR10822048_1.fastq.gz |  sed -n '2~4p' | head -n 1000000 | sort | uniq -c | sort -n | tail




@ERR10822048.1 MS3_45020:1:1114:8721:5555/1
CGATGTGAGATCCCTCAGACCCTTTTAG-TCAGTGTGGAAAATCTCTAGCAGGG-CCCTAAAAACATGACAAGAGGGAGGCGACACCTCGGCCGAACAAGCATGGTTAGGTAGGACCTTAACAACGAACTCGTTCAATTTTTCGCCTGTCTT
CGATGTGAGATCCCTCAGACCCTTTTAG-TCAGTGTGGAAAATCTCTAGCAGGG-CCCCTAGGCTGTTGAGGTAACACATTAGGTGGTTTAGGACGCAGTCTACGGTTAGTCCCTTAAGCGGAGGCCCTATAGTGAGTCGTATTACAGATCG
CGATGTGAGATCCCTCAGACCCTTTTAG-TCAGTGTGGAAAATCTCTAGCAGGG-CCCCTAGGCTGTTGAGGTATCACATTTGGTGGTTTAGGACGCAATATACGGTTAGTCCCTTTAGCGGAGGCCCTATAGTGAGTCGTATTACTAATCG
CGATGTGAGATCCCTCAGACCCTTTTAG-TCAGTGTGGAAAATCTCTAGCAGGG-CCCCTAGGCTGTTGACTTAACACGTTAGGTGTTATAGGACGCAATCTACGGTTAGTCCCTTAGGCGGGACCCCTCTCGTGAGTCTTATTACAGAGCG
CGATGTGAGATCCCTCAGACCCTTTTAG-TCAGTGTGGAAAATCTCTAGCAGGG-GCGCCCGGACGTGGACTTAAAATCTAAAGTGTTACAGGACGCACTCTCTGGTCAGAAGCCTCGGCGTGCTCAATCTCATACGGCAATTTACAAGGCG
+


@ERR10822048.1 MS3_45020:1:1114:8721:5555/2
GTAATACGACTCACTATAGGGC-CTCCGCTTAAGGGACTAA-AGACATCCCAAAAGTTGAACTAGAGCGATGTCAAATCCTCGCCCAACCCTGCGGATTCCTCCGACCTGTCTTCTACTGTCTGTCATGCTTTAAATTTTCTGCAATAGCTT
GTAATACGACTCACTATAGGGC-CTCCGCTTAAGGGACTAA-CCGTAGACTGCGTCCTAAACCACCTAATGTGTTACCTCAACAGCCTAGGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCGAGATCG
GTAATACGACTCACTATAGGGC-CTCCGCTTAAGGGACTAA-CCGTAGATTGCGCCCTCAACCAACTAGTGTGTTACCTCAACCACCCCGGGATCCTGATAGAGATTTTCAGCACTAACTCAAAGAGTCTGAGGGACCTCGCATCGAGATCG
GTAATACGACTCACTATAGGGC-CTCCGCTTAAGGGACTAA-CCGTAGACTTCGCCCTCAGCCACCTAGTGTGTTACCTCACCAACCCCGAGATCCTGCTAGCGATTTTCCGCACTGGCTAATAGGTTCTGGTGGTTCTCGCATCATTATCG
GTAATACGACTCACTATAGGGC-CTCCGCTTAAGGGACTAA-CCGTAGATTGCGTCCTAAAACACCTAATGTGTTACCTCAATAGCCTAGAGGTCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGTTCTCACATCAAGATCG



samtools view ERR10822048.bam | cut -f10 | head
GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC-TAACCGTAGACTGTGTCCTAAACCACCTAATGTGTTACCTCTATAGCCCAGAGGCCCTGCTAGAGATTTTCCACAATGACTTAAAGGGCCTGAGGGATCCCACATCAAGATCG
GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC-TAACCGTAGACTGCGCCCTAAACCACCTAATGTGTTACCTCAACAGCCTAGGGGCCCTGCTAGAGATTTTCAACACTGACTAAAAGGGTTTGAGGGATCTCACATCGTTATCG
GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC-TAACCGTAGACTGTGTCCTAAACCACCTAATGTGTTACCTCAACAGCCTAGGGGCCCTGCTAGAGACTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCGAGATCG
GTAATACGACTCACTATAGCGCCTCCGCTTAAGGGAC-TAACCGTAGACTGCTTCCTAAACCACCTAATGTGTTACCTCTACACCCCAGGAGCCCTGCTAGAGATTTTGCACACTGACTAAAATGGTCAGAAGGATCTCACATCTAGATCG
CGATCTGTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC-TAACCGTAGACTGCGTCCTAAACCACCAAATGTGTTACCTCAACAGCCTAGGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCG
CGATCTGTAATACCACTCGCTATAGGGCCTACGCTTAAGGGAC-TAACCGTAGACTGCGACCTAAACCACCAAATGTGTTACCTCAACAGCCTAGGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCG
CGATCTGTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC-TAACCGTAGCCTGCGTCCTAAACCACCTAATGTGTTACCTCAACAGCCTAGGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCG
CGATCTGTAATACGACTCACTATAGGGGCTCCGGTTAAGGGAC-TAACCGTAGACTGTGACCTAAACCACCTAATGTGTTACCTAAAAAGCACAGGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCG
CGATCTGTAATACGACGCACTATAGGGCAGCCACTTAAGGGAC-TAACCGTGGACTGCGTCCTAAACCACCTAATGTGTTACCTCAGAAGCCCACGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCG
CGATCTGTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC-TAACCGTAGACTGCGTCATAAACCACCTAATGTGTTACCTCAACAGCCAAGGGGCCCTGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCACATCG




# creating sequnence logo of first 20 bp of trichuris sequennce
samtools view ERR10822048.bam | cut -f10  | sort | uniq | grep "GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC" | sed 's/^.*GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC//g' | sed 's/TGCTAGAGATTTTCC.*$//g' | cut -c-20 | sort | uniq -c | sort -n | awk '{if(length($2)==20) print $2}'





zcat ERR10822048_1.fastq.gz | grep -c "@ERR10822048" | head
>> 1583841 reads 

zcat ERR10822048_1.fastq.gz | grep -c "CTCCGCTTAAGGGACTAA" | head
>> 0 reads 

zcat ERR10822048_2.fastq.gz | grep -c "CTCCGCTTAAGGGACTAA" | head
>> 757690 reads with linker and RD site - linker only found in read 2 - means the assay is directional when sequneced

zcat ERR10822048_2.fastq.gz | grep -c "CTCCGCTTAAGGGAC" | head
>> 1472806 reads - more reads when RD site missing


zcat ERR10822048_1.fastq.gz | grep -c "TTAGTCCCTTAAGCGGAGG" | head
>> 360020 reads with linker in RC with RD - only found in read 1, not in read 2

zcat ERR10822048_1.fastq.gz | grep -c "GTCCCTTAAGCGGAGG" | head
>> 956483  reads with linker in RC but with no RD - only found in read 1, not in read 2

So there is a good proportion of reads with no restriction site 



zcat ERR10822048_2.fastq.gz | grep "CTCCGCTTAAGGGACTAA" | grep "TGCTAGAGATTTTCC" | wc -l
>> 459873 reads with both linker and LTR sequences


zcat ERR10822048_1.fastq.gz | grep "GACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGG" | wc -l
>> 354000

zcat ERR10822048_1.fastq.gz | grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGGAGGCCC" | wc -l
>> 337669


zcat ERR10822048_2.fastq.gz | grep "GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGAC" | grep "CTGACTAAAAGGGTCTGAGGGATCTC" | wc -l
>> 880055 reads with both linker, LTR , but no RD site

zcat ERR10822048_2.fastq.gz | grep "GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGACTAA" | grep "CTGACTAAAAGGGTCTGAGGGATCTC" | wc -l
>> 319454 reads with both linker, LTR and RD site

zcat ERR10822048_1.fastq.gz | grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "GTCCCTTAAGCGGAGGCCCTATAGTGAGTCGTATTAC" | wc -l
>> 864882 reads with both linker, LTR , but no RD site

zcat ERR10822048_1.fastq.gz | grep --colour "GAGATCCCTCAGACCCTTTTAGTCAG" | grep --colour "TTAGTCCCTTAAGCGGAGGCCCTATAGTGAGTCGTATTAC" | wc -l
>> 316316 reads with both linker, LTR and RD site


for i in *_1.fastq.gz; do echo $i ; zcat $i | grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGGAGGCCC" | wc -l; done
ERR10822048_1.fastq.gz
337669
ERR10822049_1.fastq.gz
136833
ERR10822050_1.fastq.gz
261272
ERR10822051_1.fastq.gz
33863
ERR10822052_1.fastq.gz
217562
ERR10822053_1.fastq.gz
231912

for i in *_2.fastq.gz; do echo $i ; zcat $i | grep "GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGACTAA" | grep "CTGACTAAAAGGGTCTGAGGGATCTC"  | wc -l; done
ERR10822048_2.fastq.gz
319454
ERR10822049_2.fastq.gz
113456
ERR10822050_2.fastq.gz
191552
ERR10822051_2.fastq.gz
30504
ERR10822052_2.fastq.gz
164190
ERR10822053_2.fastq.gz
191632


for i in *_1.fastq.gz; do 
    echo $i ; zcat $i |\
     grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGGAGGCCC" |\
     sed 's/^.*GAGATCCCTCAGACCCTTTTAGTCAG//g' | sed 's/TTAGTCCCTTAAGCGGAGGCCC.*$//g' |\
     sort | uniq -c | wc -l; done


ERR10822053_1.fastq.gz

for i in ERR10822053_1.fastq.gz; do 
    echo $i ; zcat $i |\
     grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGGAGGCCC" |\
     sed 's/^.*GAGATCCCTCAGACCCTTTTAGTCAG//g' | sed 's/TTAGTCCCTTAAGCGGAGGCCC.*$//g' |\
     sort | uniq -c ; done



Sequence preference
for i in ERR10822052_1.fastq.gz; do      echo $i ; zcat $i |     grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGGAGGCCC" |     sed 's/^.*GAGATCCCTCAGACCCTTTTAGTCAG//g' | sed 's/TTAGTCCCTTAAGCGGAGGCCC.*$//g' |     sort | uniq -c ; done | sort -n

ERR10822048_1.fastq.gz
TGTGGAAAATCTCTAGCAGGGCCC

ERR10822049_1.fastq.gz
TGTGGAAAATCTCTAGCA

ERR10822050_1.fastq.gz
TGTGGAAAATCTCTAGCAGGGCCC

ERR10822051_1.fastq.gz
TGTGGAAAATCTCTAGCA

ERR10822052_1.fastq.gz
TGTGGAAAATCTCTAGCAGGGCCC

ERR10822053_1.fastq.gz
TGTGGAAAATCTCTAGCAGGGCCC


### length of insert between LTR and linker primers
```bash
for i in *_1.fastq.gz; do 
    zcat $i |  grep "GAGATCCCTCAGACCCTTTTAGTCAG" | grep "TTAGTCCCTTAAGCGGAGGCCC" | sed 's/^.*GAGATCCCTCAGACCCTTTTAGTCAG//g' | sed 's/TTAGTCCCTTAAGCGGAGGCCC.*$//g' | awk -v i=${i%_1.fastq.gz} '{if(length>=0) print i,"R1",length }' OFS="\t" ; 
    done > length.histo

for i in *_2.fastq.gz; do 
    zcat $i | grep "GGGCCTCCGCTTAAGGGACTAA" | grep "CTGACTAAAAGGGTCTGAGGGATCTC"  | sed 's/^.*GGGCCTCCGCTTAAGGGACTAA//g' | sed 's/CTGACTAAAAGGGTCTGAGGGATCTC.*$//g' | awk -v i=${i%_2.fastq.gz} '{if(length>=0) print i,"R2",length }' OFS="\t" ; 
    done >> length.histo


```R
library(tidyverse)

data <- read.table("length.histo", header=F)

ggplot(data, aes(V3)) + 
    geom_histogram(binwidth=1) +
    facet_grid(V1~V2, scales="free_y") + 
    labs(x="Distance between LTR and linker sites (bp)")

```



### scaffolds with exact hits to LTR and linkers
```bash
samtools view ERR10822050.bam | grep --color=always "TTAGTCCCTTAAGCGGAGGCCC\|GGGCCTCCGCTTAAGGGACTAA" | grep --color=always  "GAGATCCCTCAGACCCTTTTAGTCAG\|CTGACTAAAAGGGTCTGAGGGATCTC" | cut -f3 | sort | uniq -c
```



### Identifying reads with linker sequences
- select reads that have both LTR and adaptor sequences, trim them, and then map them
- filtered for good mapping qualities, unique mapping
- selecting sites with at least 100 reads

```bash
module load cutadapt/1.15--py36_0

cd /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/TRANSGENESIS

for i in ERR10822048 \
ERR10822049 \
ERR10822050 \
ERR10822051 \
ERR10822052 \
ERR10822053; do 

cutadapt --quiet --trimmed-only -g GAGATCCCTCAGACCCTTTTAGTCAG ${i}_1.fastq.gz | cutadapt  --quiet --trimmed-only -a TTAGTCCCTTAAGCGGAGGCCCTATAGTGAGTCGTATTAC - > ${i}.trimmed.fastq

cutadapt --quiet --trimmed-only -g GTAATACGACTCACTATAGGGCCTCCGCTTAAGGGACTAA ${i}_2.fastq.gz | cutadapt --quiet --trimmed-only -a CTGACTAAAAGGGTCTGAGGGATCTC - >> ${i}.trimmed.fastq

minimap2 -x sr trichuris_muris.PRJEB126.WBPS18.genomic.fa ${i}.trimmed.fastq  > ${i}.tm.minimap.paf
minimap2 -a -x sr trichuris_muris.PRJEB126.WBPS18.genomic.fa ${i}.trimmed.fastq | samtools sort -o ${i}.tm.minimap.bam --write-index -

minimap2 -x sr GRCh38_latest_genomic.fna ${i}.trimmed.fastq  > ${i}.hs.minimap.paf
minimap2 -a -x sr GRCh38_latest_genomic.fna ${i}.trimmed.fastq | samtools sort -o ${i}.hs.minimap.bam --write-index -


awk '{if($12>=50) print }' ${i}.tm.minimap.paf | cut -f6,8 | sort | uniq -c | awk '{if($1>=100) print}' OFS="\t" | sort -k2,2 -k3,3n > ${i}.insertion.coords
awk '{if($12>=50) print }' ${i}.hs.minimap.paf | cut -f6,8 | sort | uniq -c | awk '{if($1>=100) print}' | sort -k2,2 -k3,3n >> ${i}.insertion.coords;
done


# generate a filtered set of paf data for plotting.
for i in *tm.minimap.paf; do 
    awk '{if($12>=50) print}' > ${i%paf}filtered.paf; 
    done
```


samtools faidx trichuris_muris.PRJEB126.WBPS18.genomic.fa
cat trichuris_muris.PRJEB126.WBPS18.genomic.fa.fai | grep "TMUE_LG" | sort | cut -f1,2 > tmuris.genome

bedtools makewindows -g tmuris.genome -w 10000 > tmuris.genome.10k.bed

for i in *.tm.minimap.filtered.paf; do
    cat ${i} | cut -f6,8,9 > ${i%.paf}.bed
    bedtools coverage -a tmuris.genome.10k.bed -b ${i%.paf}.bed > ${i%.paf}.10k.counts
    awk -v name=${i%.tm.minimap.filtered.paf} '{print name, $0}' OFS="\t" ${i%.paf}.10k.counts > tmp; mv tmp ${i%.paf}.10k.counts ; 
done

cat *.10k.counts > data.10k.counts

```R
library(tidyverse)

data <- read.table("data.10k.counts", header=F)

ggplot(data, aes(V3, log(V5))) + geom_point() + facet_grid(V1~V2)

```



@ERR10822050.155617 MS3_45020:1:1101:16003:2069/1
CGATGTGAGATCCCTCAGACCCTTTTAGTCAG TGTGGAAAATCTCTAGC AGGGCCCCTGTTGGGAATCACTGGAGTTGAGGCTGATTTCATAACTGAATTACCGGTGAGTGCGAATCCAACTGAAATGATTCACTTTAGTCCCTTAAGCG


ERR10822050.155617
AGTGAATCATTTCAGTTGGATTCGCACTCACCGGTAATTCAGTTATGAAATCAGCCTCAACTCCAGTGATTCCCAACAGGGGCCCT GCTAGAGATTTTCCACA






