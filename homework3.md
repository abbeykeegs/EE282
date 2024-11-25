# Homework 3

### Abiageal Keegan

#### Part 1: Summarize Genome Assembly

In order to summarize the *Drosophila melanogaster* all chromosomes genome file, I will first begin by loading in the data and processing it via `faSize`.

``` bash
wget https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz
faSize -detailed dmel-all-chromosome-r6.60.fasta.gz > all-chromosome-processed.fasta.gz
```

Now I have the original file and the processed data. However, before I continue I want to verify the file integrity using a checksum. The flybase website utilizes a `md5sum`. I will compare the value listed on the flybase website text file and the value I get from running the following checksum code block.

``` bash
md5sum dmel-all-chromosome-r6.48.fasta.gz
```

I have now verfied that my file has not been tampered with and can continue processing the data as usual. My goal is to determine the total number of nucleotides, total number of Ns, and total number of sequences. In order to do this, I want to eliminate the sequence IDs. I will do this by sorting the file for ease of reading and then cutting the first column of sequence IDs.

``` bash
sort all-chromosome-processed.fasta.gz > all-chromosome-sorted.fasta.gz
cut -f 2 all-chromosome-sorted.fasta.gz > all-chromosome-cut.fasta.gz
```

Now I am left with a file that contains only the sequence sizes. This will allow me to determine the number of sequences and the total number of nucleotides. To determine the number of sequences, I will deliver a grep command to count each entry of the file.

``` bash
grep -c '' all-chromosome-cut.fasta.gz
```

Based on the result of running this command, **there are 1870 sequences in the file**. To determine the number of nucleotides, I want to sum up the length of all the sequences present. To do so, I will utilize the following amk command that will sum each entry of the first (and only) column.

``` bash
bioawk '{sum += $0} END {print sum}' all-chromosome-cut.fasta.gz
```

Based on the result of running this command, **there are 143726002 nucleotides in the file**. Now, I must determine the number of Ns in the file. In order to do this, I will have to go back to the original sequence file. This file contains the list of nucleotides. I will first utilize `zcat` to unzip the file and then tee it to `grep` to count the number of characters that are N in the file.

``` bash
zcat dmel-all-chromosome-r6.60.fasta.gz | grep -o 'N' | wc -l
```

Based on the result of running the command, **there are 1154851 Ns present in the sample**. This is about 0.8% of all the nucleotides.

#### Part 2: Summarize Genome Annotation File

I will begin by once again performing a file integrity check using a checksum. The flybase website utilizes a `md5sum`. I will compare the value listed on the flybase website text file and the value I get from running the following checksum code block.

``` bash
md5sum dmel-all-r6.60.gtf.gz 
```

Now that I have verified the file integrity, my goal is to identify the number of features of each type, sorted from the most common to the least common. I will do so using a bioawak command that will identify the third field of the file (gene), list the results, sorting the results, identifying unique entries and combining them, and then once again sorting the results in reverse order.

``` bash
bioawk -c gff '{print $3}' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr
```

**This results in the following answer:**

| Commonality Rank | Prevalence (Count) | Type        |
|------------------|--------------------|-------------|
| 1                | 190038             | exon        |
| 2                | 163253             | CDS         |
| 3                | 46806              | 5' UTR      |
| 4                | 33741              | 3' UTR      |
| 5                | 30888              | start codon |
| 6                | 30828              | stop codon  |
| 7                | 30802              | mRNA        |
| 8                | 17872              | gene        |
| 9                | 3059               | ncRNA       |
| 10               | 485                | miRNA       |
| 11               | 365                | pseudogene  |
| 12               | 312                | tRNA        |
| 13               | 270                | snoRNA      |
| 14               | 262                | pre miRNA   |
| 15               | 115                | rRNA        |
| 16               | 32                 | snRNA       |

My goal is to now identify the number of genes per chromosome arm. To do so, I will once again utilize a bioawk command. I will again identify the third field of the gtf file and specifically search for lines that begin with "gene". I will then fetch the first field of the gene entries that specifies their location on the chromosome arm. I will then once again sort and combine entries with the same location using the following command.

``` bash
bioawk -c gff '$3 == "gene" {print $1}' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr
```

**This results in the following answer:**

| Commonality Rank | Prevalence (Count) | Location |
|------------------|--------------------|----------|
| 1                | 4226               | 3R       |
| 2                | 3649               | 2R       |
| 3                | 3508               | 2L       |
| 4                | 3481               | 3L       |
| 5                | 2704               | X        |
| 6                | 114                | 4        |
| 7                | 113                | Y        |
