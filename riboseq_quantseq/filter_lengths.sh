for i in *.bam 
do 
samtools view -h $i | awk 'length($10) >= 29 || length($10) <= 35 || $1 ~ /^@/' | samtools view -bS - > $i\filtered_29_35.bam
done
