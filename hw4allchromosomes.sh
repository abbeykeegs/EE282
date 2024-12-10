wget https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz
faFilter -minSize=100000 dmel-all-chromosome-r6.60.fasta.gz lgfilter.fasta
faSize lgfilter.fasta
faFilter -maxSize=100000 dmel-all-chromosome-r6.60.fasta.gz smfilter.fasta
faSize smfilter.fasta
bioawk -c fastx '{ print length($seq) }' lgfilter.fasta >lglengths.txt
bioawk -c fastx '{ print length($seq) }' smfilter.fasta >smlengths.txt
bioawk -c fastx '{ print gc($seq) }' lgfilter.fasta > lggc.txt
bioawk -c fastx '{ print gc($seq) }' smfilter.fasta > smgc.txt
bioawk -c fastx '{ print length($seq) }' lgfilter.fasta | sort -nr > lgsort.txt
plotCDF lgsort.txt lgCDF.png
bioawk -c fastx '{ print length($seq) }' smfilter.fasta | sort -nr > smsort.txt
plotCDF smsort.txt smCDF.png
