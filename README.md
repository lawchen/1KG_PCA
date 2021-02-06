# README for 1000 Genomes Phase III Project Data

Documented by Lawrence Chen

Start Date: June 16, 2019

The data collected and commands here follow [Kevin Blighe's tutorial on producing a PCA plot for 1000 Genomes (1KG) Phase III](https://www.biostars.org/p/335605/)

## Final Output

1KG whole genome (MAF filtered, LD-pruned) = Merged.bed/.bim/.fam

PCA results = plink.eigenval/.eigenvec

PCA Plot = 1KG_PCA_plot.pdf

## Command Log

1. Download the files as VCF.gz (and tab-indices)

```
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ;
for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done
```

2. Download 1000 Genomes PED file

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;
```

3. Download the GRCh37 / hg19 reference genome

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;
gunzip -c human_g1k_v37.fasta.gz > human_g1k_v37.fasta ; # (I slightly modified Blighe's command here)
```

4. Convert the 1000 Genomes files to BCF

```
for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
    ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
    bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' |
    bcftools norm -Ob --rm-dup both \
    > ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    
    bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done
```

5. Convert the BCF files to PLINK format

```
for chr in {1..22}; do
    plink --noweb --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed \
    --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done
```

6. Exclude variants not on the coding strand

Not performed because Blighe suggest it's not applicable.

7. Prune variants from each chromosome

```
mkdir Pruned ;

for chr in {1..22}; do
    plink --noweb --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --maf 0.10 --indep 50 5 1.5 \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;

    plink --noweb --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.prune.in --make-bed \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done
```

8. Get a list of all PLINK files

```
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;

sed -i 's/.bim//g' ForMerge.list ;
```

9. Merge all projects into a single PLINK file

```
plink --merge-list ForMerge.list --out Merge ;
```

10. Perform PCA

```
plink --bfile Merge --pca
```

11. Generate plots in R

```
R

options(scipen=100, digits=3)

# read in the eigenvectors, produced in PLINK
eigenvec <- data.frame(read.table("plink.eigenvec", header=FALSE, skip=0, sep=" "))
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste("Principal Component ", c(1:20), sep="")

# read in the PED data
PED <- data.frame(read.table("20130606_g1k.ped", header=TRUE, skip=0, sep="\t"))
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
all(PED$Individual.ID == rownames(eigenvec)) == TRUE
[1] TRUE

# set colours
#BiocManager::install("RColorBrewer")
require("RColorBrewer")

# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
        "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
        "CLM","MXL","PEL","PUR",
        "CDX","CHB","CHS","JPT","KHV",
        "CEU","FIN","GBR","IBS","TSI",
        "BEB","GIH","ITU","PJL","STU"))

col <- colorRampPalette(c(
        "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
        "forestgreen","forestgreen","forestgreen","forestgreen",
        "grey","grey","grey","grey","grey",
        "royalblue","royalblue","royalblue","royalblue","royalblue",
        "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

# generate PCA bi-plots
project.pca <- eigenvec
summary(project.pca)


pdf("1KG_PCA_plot.pdf", width=16, height=8) # I slightly modified Blighe's code here to adjust the image size and the font size
par(mar=c(5,5,5,5), cex=1, cex.main=3, cex.axis=1.75, cex.lab=1.75, mfrow=c(1,2))

plot(project.pca[,1], project.pca[,2], type="n", main="A", adj=0.5, xlab="First component", ylab="Second component")
points(project.pca[,1], project.pca[,2], col=col, pch=20, cex=1.25)
legend("bottomright", bty="n", cex=1, title="", c("African","Hispanic","East-Asian","Caucasian","South Asian"), fill=c("yellow","forestgreen","grey","royalblue","black"))

plot(project.pca[,1], project.pca[,3], type="n", main="B", adj=0.5, xlab="First component", ylab="Third component")
points(project.pca[,1], project.pca[,3], col=col, pch=20, cex=1.25)

dev.off()
```
