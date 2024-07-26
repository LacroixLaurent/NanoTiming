### bash script to prepare for nanoT_parsing:

NanoTimingDir=YourChoice
Exp=WT_rep3

path2work="$NanoTimingDir"/"$Exp"

rm -rf "$path2work"/Parsed
rm -rf "$path2work"/bam_list.txt
rm -rf "$path2work"/mega/header
rm -rf "$path2work"/mega/mod_splitted_*.bam

mkdir -p "$path2work"/script_done
mv "$path2work"/nanoT_* "$path2work"/script_done/

