for i in *sam.bam
do
featureCounts $i -a /camp/lab/ulej/working/Rupert/genomes/hs/annotation/gencode.v29.annotation.gtf -T 8 -o $i\.rfp_counts
done