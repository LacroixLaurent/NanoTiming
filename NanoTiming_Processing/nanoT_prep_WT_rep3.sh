### bash script to prepare for nanoT_parsing:

NanoTimingDir=YourChoice
Exp=WT_rep3

path2work="$NanoTimingDir"/"$Exp"

NREADS=50000
NCORES=24

cd "$path2work"/mega
samtools view -H mod_mappings.bam > header
samtools view -@ "$NCORES" mod_mappings.bam | split - mod_splitted_ -l "$NREADS" -a 3 -d --filter='cat header - | samtools view -@ "$NCORES" -b - > "$FILE".bam'
find "$path2work"/mega -name "mod_splitted_*.bam" > "$path2work"/bam_list.txt
