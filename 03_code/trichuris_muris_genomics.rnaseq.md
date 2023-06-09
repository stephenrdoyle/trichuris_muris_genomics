# trichuris_muris_genomics.rnaseq

### stephen doyle

- code used to explore bulk RNAseq data from Maria Duque's lab comparing T. muris development in mice and in organoids



## Download reference from WormBase ParaSite
```bash
# working directory
cd /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/RNAseq/REF

# genome
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/trichuris_muris/PRJEB126/trichuris_muris.PRJEB126.WBPS17.genomic.fa.gz

# annotation
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS17/species/trichuris_muris/PRJEB126/trichuris_muris.PRJEB126.WBPS17.annotations.gff3.gz

gunzip trichuris_muris.PRJEB126.WBPS17.genomic.fa.gz
gunzip trichuris_muris.PRJEB126.WBPS17.annotations.gff3.gz
```



## Get the raw RNAseq data
- Maria gave me some metadata of samples all from the lane 46825_1
- need to download these from iRods, and then convert from cram to fastq format

```bash
# raw data directory
cd /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/RNAseq/RAW_DATA

# download data from iRods as cram files
for i in {1..36}; do 
    iget /seq/illumina/runs/46/46825/lane1/plex${i}/46825_1#${i}.cram; 
    done

# convert crams to fastq files
bsub.py --queue long 10 cram2fq "./run_cram2fastq"

# once finished, remove crams
rm *cram
```


- where "run_cram2fastq" is :
```bash
#!/bin/bash

for i in *cram; do \
	samtools view -ub --threads 7 ${i} -T /lustre/scratch117/core/sciops_repository/references/Trichuris_muris/V5_250416/all/fasta/Trichuris_muris_V5_250416.fa |\
	samtools sort -n - |\
	samtools fastq --threads 7 -0 ${i%.cram}.fastq.gz -1 ${i%.cram}_1.fastq.gz -2 ${i%.cram}_2.fastq.gz - ;\
done
```




## Run Kallisto

```bash 
# working directory
cd /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/RNAseq/MAPPING

# get reference data
ln -s ../REF/trichuris_muris.PRJEB126.WBPS17.genomic.fa REF.fa
ln -s ../REF/trichuris_muris.PRJEB126.WBPS17.annotations.gff3 ANNOTATION.gff3

# make a transcripts fasta
gffread -x TRANSCRIPTS.fa -g REF.fa ANNOTATION.gff3


# load kallisto
module load kallisto/0.46.2--h4f7b962_1

# index the transcripts
kallisto index --index TRANSCRIPTS.ixd TRANSCRIPTS.fa

# run kallisto
while read LANE SAMPLE_ID; do \
kallisto quant \
     --bias \
     --index TRANSCRIPTS.ixd \
     --output-dir kallisto_${SAMPLE_ID}_out \
     --bootstrap-samples 100 \
     --threads 7 \
     --fusion \
     ../RAW_DATA/${LANE}_1.fastq.gz ../RAW_DATA/${}_2.fastq.gz;
done < /nfs/users/nfs_s/sd21/lustre_link/trichuris_muris/RNAseq/lanes_samples.list

mkdir KALLISTO_MAPPED_SAMPLES
mv kallisto_* KALLISTO_MAPPED_SAMPLES/


```


```bash

# extract TPMs per sample
for i in ` ls -1d *out `; do 
  echo $i > ${i}.tpm ; cat ${i}/abundance.tsv | cut -f5 | sed '1d' >> ${i}.tpm; 
  done

# generate a "transcripts" list, taken from the TRANSCRIPTS.fa file
#echo "ID" > transcripts.list; grep ">" ../TRANSCRIPTS.fa | cut -f1 -d  " " | sed 's/>//g' >> transcripts.list
# due to Apollo giving long unique codes, the transcript IDs are obscure. Here is the fix
#awk '$3=="mRNA" {print $9}' ../ANNOTATION.gff3 | cut -f3,5 -d";" | sed -e 's/ID=//g' -e 's/;Name=/\t/g' > mRNA_IDtoNAME_conversion.txt

#while read ID NAME; do sed -i "s/${ID}/${NAME}/g" transcripts.list; done < mRNA_IDtoNAME_conversion.txt &

# ALTERNATE WAY, direct from the annotaiton
echo "ID" > transcripts.list; grep ">" ../TRANSCRIPTS.fa | cut -f1 -d" " | sed -e 's/>//g' >> transcripts.list



# make a data frame containing all TMP values from all samples
cp transcripts.list tmp; ls -1v *_out.tpm | while read SAMPLE; do paste tmp $SAMPLE > tmp2; mv tmp2 tmp; done

mv tmp kallisto_allsamples.tpm.table

sed -i -e 's/Transcript://g' -e 's/kallisto_//g' -e 's/_out//g' kallisto_allsamples.tpm.table
```


```R
# load libraries
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(gplots)

# load data
data<-read.table("kallisto_allsamples.tpm.table", header=T, row.names=1, sep="\t" )

# set a TPM cutoff,
data<-(data > 1) * (data - 1) + 1
data<-log10(data)
data<-as.matrix(data)

# fix infinite values to NA
is.na(data) <- sapply(data, is.infinite)

# calculate variance per row
#https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
RowVar <- function(x, ...) {
  rowSums(na.rm=TRUE,(x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
var<-as.matrix(RowVar(data))
var<-as.data.frame(var)
var <- var[order(var), ,drop = FALSE]
var_filter <- tail(var,500)   # set number of genes here
var_filter <- rownames_to_column(var_filter)

# bring all data together
data<-as.data.frame(data)
data<-rownames_to_column(data)

data_filtered <- dplyr::semi_join(data, var_filter, by = "rowname")
data_filtered <- column_to_rownames(data_filtered,'rowname')
data_filtered<-as.matrix(data_filtered)
#var<-as.matrix(RowVar(data))
#data<-cbind(data, variance = var )


# dendrogram for genes only
heatmap.2(data_filtered,trace="none",na.color="grey",labRow=F,dendrogram='row',Colv=FALSE,col= colorRampPalette(brewer.pal(8, "Blues"))(25), margins=c(12,8))

# dendrogram for both
heatmap.2(data_filtered,trace="none",na.color="grey",labRow=F,dendrogram='both',col= colorRampPalette(brewer.pal(8, "Blues"))(25), margins=c(12,8))


```




```R
#Create a matrix from our table of counts
pca_matrix <- read.table("kallisto_allsamples.tpm.table", header=T) %>%
#Add one to all TPM values
mutate(across(where(is.numeric),~ .x + 1)) %>%
#Log-transform the TPM+1 values
mutate(across(where(is.numeric),log2)) %>%
#Transpose the counts table so the "gene" column becomes the rownames
column_to_rownames("ID") %>%
#Coerce into a matrix
as.matrix() %>%
#Transpose the matrix so that rows=samples and columns=variables
t()
#Check for #N/A and non-numeric values
all(is.finite(pca_matrix))

#Perform the PCA
sample_pca <- prcomp(pca_matrix)

#Check output of prcomp
pca_matrix[1:10, 1:5]

#Class of the object
class(sample_pca)

#Structure of the object
str(sample_pca)

# "sdev" contains the standard deviation explained by each PC, so if we square it we get the eigenvalues (or explained variance)
# "rotation" contains the variable loadings for each PC, which define the eigenvectors
# "x" contains the PC scores, i.e. the data projected on the new PC axis
# "center" in this case contains the mean of each gene, which was subtracted from each value
# "scale" contains the value FALSE because we did not scale the data by the standard deviation

# eigenvalues squared
pc_eigenvalues <- sample_pca$sdev^2

#Variable loadings
pc_loadings <- sample_pca$rotation

ncol(pc_scores)
length(pc_eigenvalues)
ncol(pc_loadings)

#PC scores
pc_scores <- sample_pca$x

#Variance explained by PCs
#Convert matrix to tibble object in order to make a plot in ggplot2
# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

#Print the result
pc_eigenvalues

#Produce a Scree Plot to show the fraction of the total variance explained by each PC
pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x
pc_scores <- pc_scores %>% 
#Convert to a tibble retaining the sample names as a new column
as_tibble(rownames = "sample")

#Print the result
pc_scores

#Load sample info i.e. metadata 
sample_info <- read_tsv("~/Desktop/STAR_stringtie/sample_info.txt")

#Extract the loadings from prcomp
pca_plot1 <- sample_pca$x %>% 
#Convert to a tibble retaining the sample names as a new column
as_tibble(rownames = "sample") %>% 
#Join with "sample_info" table
#Create the plot
ggplot(aes(x=PC1, y=PC2)) +
geom_point(aes(),size=5) +
#geom_text_repel(aes(), point.padding=0.8, size=3, box.padding = 0.3, min.segment.length=1, max.overlaps=Inf) +
theme (panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  legend.key = element_blank(),
  axis.line = element_line(colour = "gray")) +
scale_alpha(guide = 'none') +
scale_colour_viridis_d()

#Print and save the plot
pca_plot1
```