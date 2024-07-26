###bash script to export bam coverage to bed

NanoTimingDir=YourChoice
Exp=WT_rep3

path2work="$NanoTimingDir"/"$Exp"

filename="$path2work"/"$Exp"\_nanoT
samtools view -@ 24 -b "$path2work"/mega/mod_mappings.bam | bedtools bamtobed -i - | grep -v chrM > "$filename".bed
gzip -f "$filename".bed

