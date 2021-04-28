for i in *.bam
do
sbatch --wrap="featureCounts $i -a /camp/lab/ulej/working/Rupert/genomes/hs/annotation/gencode.v29.annotation.gtf -T 8 -o $i\.cds_counts"
done