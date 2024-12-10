# Homework 4

### Abiageal Keegan

#### Part 1: Summarize Partitions of a Genome Assembly

I began by downloading the *Drosophila Melanogaster* genome as follows:

``` bash
wget https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz
```

Then, I used `faFilter` to partition the genome into 2 groups, sequences over 100 kb (large) and sequences less than 100 kb (small).

``` bash
faFilter -minSize=100000 dmel-all-chromosome-r6.60.fasta.gz lgfilter.fasta
faFilter -maxSize=100000 dmel-all-chromosome-r6.60.fasta.gz smfilter.fasta
```

I then was able to utilize the `faSize` to determine the Number of sequences, N's, and nucleotides in each partition by utilizing the following commands.

``` bash
faSize lgfilter.fasta
faSize smfilter.fasta
```

These commands resulted in the following summaries:

| Partition       | Nucleotides | N's    | Sequences |
|-----------------|-------------|--------|-----------|
| Large Sequences | 137547960   | 490385 | 7         |
| Small Sequences | 6178042     | 662593 | 1863      |

#### Part 2: Plotting the Genome Summaries

To plot the cummulative length using the `plotCDF` command, I first had to begin by creating a text file with the sorted lengths of each sequence in each partition. I utilized the following `bioawk` commands to do so:

``` bash
bioawk -c fastx '{ print length($seq) }' lgfilter.fasta | sort -nr > lgsort.txt
bioawk -c fastx '{ print length($seq) }' smfilter.fasta | sort -nr > smsort.txt
```

Then, I was able to use the `plotCDF` command as described below that resulted in the following plots

``` bash
plotCDF lgsort.txt lgCDF.png
plotCDF smsort.txt smCDF.png
```

##### Large Sequences Cumulative Length Distribution

![](homework4/lgCDF.png)

##### Small Sequences Cumulative Length Distribution

![](homework4/smCDF.png)

Then, I wanted to plot a histogram of the sequence length distribution for each partition. I accomplished this by using R in the terminal with the following commands resulting in the graphs below.

``` r
smlength <- scan("smsort.txt")
png("smlengths.png"),
hist(smlength, main = "Histogram of Sequence Lengths less than 100 kb", xlab = "Sequence Length in kb")
dev.off

lglength <- scan("lgsort.txt)
png("lglengths.png"),
hist(lglength, main = "Histogram of Sequence Lengths larger than 100 kb", xlab = "Sequence Length in kb")
dev.off()
```

![](homework4/smlengths.png)

![](homework4/lglengths.png)

To plot the GC content of each partition, I had to begin by determining the GC content of each sequence utilizing the following `bioawk` commands:

``` bash
bioawk -c fastx '{ print gc($seq) }' lgfilter.fasta > lggc.txt
bioawk -c fastx '{ print gc($seq) }' smfilter.fasta > smgc.txt
```

Then, I utilized R in the terminal and input the following commands to result in graphs of the GC distribution in each of the partitions.

``` r
smgc <- scan("smgc.txt")
png("smgc.png")
hist(smgc, main = "Histogram of GC in Sequence Lengths less than 100 kb", xlab = "GC content", ylab = "Frequency") 
dev.off()

lggc <- scan("lggc.txt")
png("lggc.png")
hist(lggc, main = "Histogram of GC in Sequence Lengths larger than 100 kb", xlab = "GC content", ylab = "Frequency") 
dev.off()
```

![](images/smgc.png)

![](images/lggc.png)

#### Part 3: Assembling a Genome

I began assembling the genome utilizing PacBio HiFi Reads by downloading the reads from the specified folder [`/data/class/ee282/public/ISO1_Hifi_AdaptorRem.40X.fasta.gz`](https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz). Then, assembled the reads using `hifiasm` and then converted the reads to the correct format with the following commands.

``` bash
hifiasm -o iso1_assembly -t 2 ISO1_Hifi_AdaptorRem.40X.fasta.gz
awk '/^S/{print ">"$2"\n"$3}' iso1_assembly.bp.p_ctg.gfa > assembly.fasta
```

#### Part 4: Assembly Assesment

To calculate the N50 of the assembled reads, I utilized bioawk to take the assembly, find the sequence lengths, sort the lengths, sum the lengths, dividing by 2, and getting the base at the 50% mark.

``` bash
bioawk -c fastx '{ print length($seq) }' assembly.fasta | sort -nr > lengths.txt
awk '{ sum += $1; if (sum >= total_length / 2) { print $1; exit } }' total_length=$(bioawk -c fastx '{ sum += length($seq) } END { print sum }' assembly.fasta)
```

This resulted in an N50 of 33. The scaffold N50 listed on the website is 25.3 Mb and the contig N50 is 21.5 Mb. All these values are in the same ballpark.

I unfortunately ran out of time to produce graphs, however, I would utilize the lengths.txt from the previous problem to plot the contiguity plot.

``` bash
plotCDF lengths.txt > contiquity.png
```
