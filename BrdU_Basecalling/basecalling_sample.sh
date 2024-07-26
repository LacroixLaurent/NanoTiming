# Sample of the basecalling procedure used with our BrdU trained model
# This requires the ONT megalodon and guppy software
# We used guppy v4.4.1 GPU and megalogon v2.2.9 with minimap2 v2.24.
# the configuration file should point to the modified base model
# see BrdU_Configutation4megalodon.cfg
# (the model file is provided as a compressed archive (model_adversial_bug.json.zip))
# a conda env was created using the yml file provided within miniconda3
conda env create -f meg2.2.9env.yml
# then within this env, we used megalodon with the following parameters
conda activate mega_2.2.9
WORKDIR=NanoTiming
EXP=WT_rep3
REF="$WORKDIR"/Reference_Genome/BT1_multiUra.fa

GUPPYPATH="~/src/ont-guppy_4.4.1"
CFG=BrdU_Configuration4megalodon.cfg
# the config file is placed in the "$GUPPYPATH"/data folder
INPUT="$WORKDIR"/data/"$EXP"/fast5
OUTPUT="$WORKDIR"/data/"$EXP"/mega

nohup megalodon "$INPUT" --guppy-server-path "$GUPPYPATH"/bin/guppy_basecall_server --outputs mod_mappings --reference "$REF"  --output-dir "$OUTPUT"  --guppy-config "$CFG" --disable-mod-calibration --overwrite --device cuda:0 --processes 20 --guppy-timeout 90 --mod-min-prob 0  &
# please notice than for experiment with many long reads, we had to increase the guppy-timeout to 600
# also notice that this computer was equiped with a NVidia RTX2080Ti GPU

# this generate a megalodon's style mod-mappings.bam file in the "$EXP"/mega folder
