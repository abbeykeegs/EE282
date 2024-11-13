wget https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/fasta/dmel-all-chromosome-r6.60.fasta.gz
faSize -detailed dmel-all-chromosome-r6.60.fasta.gz > all-chromosome-processed.fasta.gz
sort all-chromosome-processed.fasta.gz > all-chromosome-sorted.fasta.gz
cut -f 2 all-chromosome-sorted.fasta.gz > all-chromosome-cut.fasta.gz
awk '{sum += $0} END {print sum}' all-chromosome-cut.fasta.gz
grep -c '' all-chromosome-cut.fasta.gz

wget https://ftp.flybase.net/releases/FB2024_05/dmel_r6.60/gtf/dmel-all-r6.60.gtf.gz
md5sum dmel-all-r6.60.gtf.gz
bioawk -c gff '{print $3}' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr
bioawk -c gff '$3 == "gene" {print $1}' dmel-all-r6.60.gtf.gz | sort | uniq -c | sort -nr