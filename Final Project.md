# Abiageal Keegan EE282 Final Project

## *Investigating hiPSC-CMs through bulk RNA-seq data*

#### I. Rationale

My ultimate goal is to investigate genes that are up-regulated or down-regulated in hiPSC-CMs under oxygen-deprived conditions. I will utilize these patterns in gene expression to determine if the treatment of bone marrow or cardiac stromal cell extracellular vesicles can rescue phenotypes observed in cardiomyocytes under oxygen-deprived hypoxic conditions.

#### II. Introduction

Heart disease is the leading cause of death in the United States, accounting for nearly 1 in every 5 deaths in 2022 and costing over 250 billion dollars from 2019 to 2020 (CDC, 2024). While cardiomyopathy, chronic heart disease, represents a majority of cardiac-related fatalities and can be treated through medications and lifestyle changes, myocardial infarctions are often silent and thus are difficult to study (NIH, 2022). Hypoxia is a condition that underlies myocardial infarction when human body tissues do not receive enough oxygen (Cleveland Clinic, 2022). This may cause the heart to fail at delivering oxygen-rich blood to the rest of the body, wreaking havoc on the body as a whole. 

Recent studies have utilized human-induced pluripotent stem cell-derived cardiomyocytes (hiPSC-CMs) to investigate heart conditions in vitro. In order to analyze the impact of hypoxia on gene expression relating to heart failure, I will analyze a dataset of hIPSC-CMs under both normal and hypoxic conditions.

#### III. Methods

I began with a dataset of raw Illumina bulk RNA-seq Fastq files from the dataset, "Human iPSC-derived cardiomyocytes under hypoxia/normoxia, treated with bone marrow or cardiac stromal cell extracellular vesicles" from NCBI. The dataset contained 12 different runs, 3 per condition. The four conditions were normal healthy CMs (positive control, H), CMs under hypoxia (negative control, D), CMs under hypoxia treated with cardiac stromal cell extracellular vesicles (treatment 1, EV), and bone marrow treated hypoxic CMs (EV2). I obtained 1 read from normal CMs and 1 read from hypoxia CMs from the ncbi dataset. I first began by installing the SRA toolkit to extract the datasets from ncbi as they are large in size. I then added the directory to my path to utilize the software properly.

``` bash
wget wget https://ftp.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xvzf sratoolkit.current-centos_linux64.tar.gz
cd sratoolkit.current-centos_linux64/bin
export PATH=$PATH:$(pwd)/bin
```

I then was able to utilize the built-in prefetch command to obtain 1 dataset per condition (H1, D1)

``` bash
prefetch SRX25745531 #H1
prefetch SRX25745534 #D1
```

I then renamed all of these directories to match the conditions provided in the abstract (H1, D1). Now I must convert each of these files to the fastq format for gene alignment and gene expression analysis purposes.

``` bash
fasterq-dump SRR30284852.sra #H1
fasterq-dump SRR30284849.sra #D1
```

Now that I have obtained all of the fasta files, I need to download a reference genome in order to investigate gene expression in the future. I downloaded the commonly used human reference genome GrCh38 from NCBI. I included the fasta file and the genome annotation file for gene expression analysis purposes. I downloaded the zipped file to my directory. I opted to use the bowtie2 package to align my genomes. To being this process, I built a bowtie2 reference index using GrCh38 fasta.

``` bash
bowtie2-build GRCh38.fasta GRCh38_index
```

Now, I have obtained all of the fastq files and the reference genome. I will use the tool bowtie2 to align my genome with the sample sequences to result in a SAM file. Here is an example of the code for the D1 sample.

``` bash
bowtie2 -x GRCh38_index \
>   -1 /data/homezvol3/keegana/FinalProject/D1/SRR30284849_1.fastq \
>   -2 /data/homezvol3/keegana/FinalProject/D1/SRR30284849_2.fastq \
>   -S /data/homezvol3/keegana/FinalProject/D1/D1.sam
```

Now that I have the SAM files for each condition, I will utilize the package samtools to convert each SAM file to a BAM format and then sort it in order to reveal trends in gene expression. Here is an example of the code for the D1 sample

``` bash
samtools view -bS D1.sam > D1.bam
samtools sort -o D1_sorted.bam D1.bam 
```

Once I obtained the sorted bam files, I opted to use the tool htseq that I installed via mamba. Here is an example of the code for D1.

``` bash
htseq-count -f bam -r pos -s no -i gene_id /data/homezvol3/keegana/FinalProject/D1/D1_sorted.bam  /data/homezvol3/keegana/FinalProject/genomic.gtf > D1counts.txt
```

This allowed me to obtain a count matrix for each condition with a list of the human genes and counts for each condition. I sorted it by reverese numerical order of the counts to highlight the genes in each condition. I then exported these matrices to R to perform data analysis.

In order to analyze my results, I first wanted to observe the patterns of gene expression amongst each of the samples. To do this, I began by merging the count matrices from H1 and D1 that I derived from my previous analysis. I log transformed the feature counts in order to make the dataset easy to work with.

``` R
D1sort <- read.table("/Users/abiagealkeegan/Downloads/D1sortcounts.txt")
H1sort <- read.table("/Users/abiagealkeegan/Downloads/H1sortcounts.txt")
colnames(D1sort) <- c("Gene, "CountD")
colnames(H1sort) <- c("Gene", "CountH")
mergedset <- merge(H1sort, D1sort, by= "Gene", all= TRUE)
mergedset$LogCountH <-log10(mergedset$CountH +1)
mergedset$LogCountD <-log10(mergedset$CountD +1)
```

I then performed a density analysis to observe the patterns in gene expression

``` R
plot(density(mergedset$LogCountH, na.rm = TRUE), 
+      main = "Density Plot of LogCountH and LogCountD", 
+      xlab = "Log Counts", 
+      col = "blue", lwd = 2)
lines(density(mergedset$LogCountD, na.rm = TRUE), col = "red", lwd = 2)
legend("topright", legend = c("LogCountH", "LogCountD"), 
+        col = c("blue", "red"), lwd = 2)
```

I then made a histogram of the feature expression value for each dataset (H1 and D1)

``` R
hist(mergedset$LogCountH, main = "Histogram of LogCountH", xlab = "LogCountH", col = "blue", breaks = 20)
hist(mergedset$LogCountD, main = "Histogram of LogCountD", xlab = "LogCountD", col = "red", breaks = 20)
```

I then identified the top 5 most highly expressed genes in each group

``` R
top5D1 <- mergedset %>% arrange(desc(LogCountD)) %>% slice(1:5)
top5H1 <- mergedset %>% arrange(desc(LogCountH)) %>% slice(1:5)
```

I then wanted to see the top 5 genes with the biggest change from Hypoxia to normal conditions. I did so by calculating the difference between the two count columns, identifying the largest 5 differences, and plotting them.

``` R
mergedset$AbsChange <- abs(mergedset$LogCountH - mergedset$LogCountD)
mergedset_sorted <- mergedset[order(-mergedset$AbsChange), ]
top_genes <- head(mergedset_sorted, 5)
barplot(top_genes$AbsChange, 
        names.arg = top_genes$Gene,  # Gene names from the second column
        las = 2, 
        col = "blue", 
        main = "Top 5 Genes with Largest Change", 
        xlab = "Genes", 
        ylab = "Absolute Change (|LogCountH - LogCountD|)")
```

#### IV. Results

Through my analysis, I obsereved that the pattern of gene expression was much more broad in cells under oxygen-stress, hypoxia, as compared to healthy normal conditions. Genes in H1 were expressed at lower levels as comapred to genes in the D1 hypoxia conditions.

![](images/Density Plot.png)

A closer look at the data reveals that a majority of genes in the D1 condition are actually lowly expressed, but there are still genes expressed at higher levels.

![](images/Histogram D.png)

The histogram of the H1 feature counts confirms the density plot in that all the genes are expressed at lower levels.

![](images/Histogram H.png)

The genes that were most highly expressed in D1 as compared to H1 were RNR2, RNR1, FLNC, TTN, and HSPA5. RNR2 and RNR1 are ribosmal cluster genes. FLNC is an actin filament cytoskeleton protein. TTN encodes for a striated muscle protein associated with myopathy. HSPA5 is a heat-shock protein found in the endoplasmic reticulum. All information was obtained from *Gene Cards: The Human Gene Database.*

![](images/Top5 big change.png)

#### V. Discussion 

The gene expression analysis revealed that cells treated with hypoxia had broader gene expression patterns and more highly expressed genes. However, there are a few ways these results could be skewed. Firstly, the reference genome used is the widely utilized human GRCh38 genome from the ncbi database. While this has been proven to be highly accurate, it may not be the proper reference for this context. My experience in differentiating human-induced pluripotent stem cells has revealed that cardiomyocytes derived from these cells tend to stay fetal in nature. During development, there are increases in gene expression of many regulatory genes and proteins that aid in the development of specific tissues. These genes tend to be highly expressed at certain points in development and decrease in expression over time. Previous studies in the same context, human cardiomyocytes, have been able to identify major differences in gene expression patterns at single-cell resolution between the fetal and adult heart (Mehdiabadi et.al, 2022). This could have been a confounding factor in this analysis. Despite this possible confounding variable, there are still a few key takeaways from these results.

The two most differentialy expressed genes in the D1 Oxygen-stressed hypoxia hIPSC-CMs were both ribosomal proteins, RNR1 and RNR2. These proteins are ribonuclease reductases and have been shown to catalyze the conversion of nucleotides to deoxynucleotides in vitro (Greene et.al, 2020). This process is essential for DNA repair and replication. Although hypoxia does not directly influence or generate genome instability or DNA damage, hypoxia is in fact associated with high mutation rates, over-replication, and a decreased cellular capacity for DNA repair pathways (Kaplan et.al, 2020). Thus, the upregulation of these ribonuclease reductase genes may hint towards a decreased capacity of these cardiomyocytes to respond to DNA damage induced by hypoxia.

The differential gene analysis also revealed the upregulation of two sarcomeric structural proteins: FLNC and TTN. FLNC is an actin-cytoskeletal protein encoding for filamin C, a protein associated with mechanical stabilization and intracellular signaling in striated muscle (Leber et.al, 2024). This gene also encodes for a protein involved in damage response similar to RNR1 and RNR2; filamin C can respond to increased workload induced damage in both myofibrils and cardiomyocytes (Leber, et.al, 2024). Hypoxia has been shown to induce muscle damage at the cellular level; recent studies have shown oxygen-deprivation can alter the calcium gradient in muscle cells and induce damage (Allen et.al, 2024). Additionally, creatine kinase levels, indicative in muscle damage, were elevated in a study in which male athletes were placed in hypoxic conditions (Chen et. al, 2022). The TTN gene encodes for titin, one of the largest proteins in the human body. Titin is another sarcomeric protein that aids in myofibril organization and cell signaling (Akhtar et.al 2020). While there is no published studies directly linking titin and hypoxia, the DNA damage inflicted by hypoxia may allow for the recruitment of titin to repair damage at the sarcomeres. Overall, the upregulation of these two protiens suggests a structural damage response to oxygen-deprivation.

Another gene of note is the HSPA5 gene upregulated in hiPSC-CMs under hypoxia. While this gene appears to be substantially upregulated in hypoxic cells, the correlation between the gene and condition is unclear. HSPA5 is a member of the heat shock protein family and an endoplasmic reticulum chaperone (Rehati et. al, 2023). The protein has been shown to regulate lipid metabolism and has been associated with nonalcoholic fatty liver disease (Rehati et.al, 2023). Other than the aforementioned paper, there are no studies that connect HSPA5 with striated muscle, cardiomyocytes, or heart disease.

Overall, the upregulation of RNR1, RNR2, FLNC, and TTN connect hypoxia to intracellular stress. Future projects could be conducted to adress hypoxia at the 3D tissue level in cardiomyocytes. Additionally, future analysis could be conducted on the presented dataset to look at the replicates for each condition and gather a more accurate picture of the data. However, this analysis is still extremely informative and confirms the sarcomeric damage induced by a lack of oxygen.

#### VI. References

US Centers for Disease Control and Prevention. 2024. Heart Disease Facts. <https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR12815092&display=metadata> (Acccessed October 16, 2024)

National Heart, Lung, and Blood Institute. 2022. Cardiomyopathy treatment. <https://www.nhlbi.nih.gov/health/cardiomyopathy/treatment> (Acccessed October 16, 2024)

Cleveland Clinic. 2022. Hypoxia. <https://my.clevelandclinic.org/health/diseases/23063-hypoxia> (Acccessed October 16, 2024)

Mehdiabadi NR, Boon Sim C, Phipson B, Kalathur RKR, Sun Y, Vivien CJ, Ter Huurne M, Piers AT, Hudson JE, Oshlack A, Weintraub RG, Konstantinov IE, Palpant NJ, Elliott DA, Porrello ER. Defining the Fetal Gene Program at Single-Cell Resolution in Pediatric Dilated Cardiomyopathy. Circulation. 2022 Oct 4;146(14):1105-1108. doi: 10.1161/CIRCULATIONAHA.121.057763. Epub 2022 Oct 3. PMID: 36191067; PMCID: PMC9528943.

Greene BL, Kang G, Cui C, Bennati M, Nocera DG, Drennan CL, Stubbe J. Ribonucleotide Reductases: Structure, Chemistry, and Metabolism Suggest New Therapeutic Targets. Annu Rev Biochem. 2020 Jun 20;89:45-75. doi: 10.1146/annurev-biochem-013118-111843. PMID: 32569524; PMCID: PMC7316142.

Kaplan AR, Glazer PM. Impact of hypoxia on DNA repair and genome integrity. Mutagenesis. 2020 Feb 13;35(1):61-68. doi: 10.1093/mutage/gez019. PMID: 31282537; PMCID: PMC7317153.

Yvonne Leber, Avnika A. Ruparelia, Gregor Kirfel, Peter F.M. van der Ven, Bernd Hoffmann, Rudolf Merkel, Robert J. Bryson-Richardson, Dieter O. Fürst, Filamin C is a highly dynamic protein associated with fast repair of myofibrillar microdamage, *Human Molecular Genetics*, Volume 25, Issue 13, 1 July 2016, Pages 2776–2788.

Chen PW, Hsu CC, Lai LF, Chi CP, Yu SH. Effects of Hypoxia-Hyperoxia Preconditioning on Indicators of Muscle Damage After Acute Resistance Exercise in Male Athletes. Front Physiol. 2022 Apr 19;13:824210. doi: 10.3389/fphys.2022.824210. PMID: 35514339; PMCID: PMC9062696.

Allen DG, Whitehead NP, Yeung EW. Mechanisms of stretch-induced muscle damage in normal and dystrophic muscle: role of ionic changes. J Physiol. 2005 Sep 15;567(Pt 3):723-35. doi: 10.1113/jphysiol.2005.091694. Epub 2005 Jul 7. PMID: 16002444; PMCID: PMC1474216.

Akhtar MM, Lorenzini M, Cicerchia M, Ochoa JP, Hey TM, Sabater Molina M, Restrepo-Cordoba MA, Dal Ferro M, Stolfo D, Johnson R, Larrañaga-Moreira JM, Robles-Mezcua A, Rodriguez-Palomares JF, Casas G, Peña-Peña ML, Lopes LR, Gallego-Delgado M, Franaszczyk M, Laucey G, Rangel-Sousa D, Basurte M, Palomino-Doza J, Villacorta E, Bilinska Z, Limeres Freire J, Garcia Pinilla JM, Barriales-Villa R, Fatkin D, Sinagra G, Garcia-Pavia P, Gimeno JR, Mogensen J, Monserrat L, Elliott PM. Clinical Phenotypes and Prognosis of Dilated Cardiomyopathy Caused by Truncating Variants in the *TTN* Gene. Circ Heart Fail. 2020 Oct;13(10):e006832. doi: 10.1161/CIRCHEARTFAILURE.119.006832. Epub 2020 Sep 23. PMID: 32964742.

Rehati A, Abuduaini B, Liang Z, Chen D, He F. Identification of heat shock protein family A member 5 (HSPA5) targets involved in nonalcoholic fatty liver disease. Genes Immun. 2023 Jun;24(3):124-129. doi: 10.1038/s41435-023-00205-y. Epub 2023 May 8. PMID: 37156995; PMCID: PMC10266971.
