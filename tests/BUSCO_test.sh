#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 16 --mem 24G -p short --out busco.%a.log -J BUSCO
# for augustus training
# set to a local dir to avoid permission issues and pollution in global
module unload miniconda3
module load busco
#export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
export AUGUSTUS_CONFIG_PATH=$(realpath augustus_config)

module load workspace/scratch

CPU=${SLURM_CPUS_ON_NODE}
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
  if [ -z $N ]; then
	  N=1
  fi
fi
if [ ! $CPU ]; then
     CPU=2
fi
export NUMEXPR_MAX_THREADS=$CPU

GENOMEFOLDER=genomes
EXT=fasta
LINEAGE=capnodiales_odb10
OUTFOLDER=BUSCO
SEED_SPECIES=anidulans
mkdir -p $OUTFOLDER
for ASM in $(ls */Rant.sorted.fasta | sed -n ${N}p)
do
    FOLDER=$(dirname $ASM)
    echo "running $FOLDER"
    if [ -d "$OUTFOLDER/${FOLDER}" ]; then
	    if [ ! -f $OUTFOLDER/${FOLDER}/short_summary.specific.$LINEAGE.${FOLDER}.txt ]; then
		    busco -c $CPU -o ${FOLDER} --out_path ${OUTFOLDER} --restart  -m genome -l $LINEAGE  --offline --augustus_species $SEED_SPECIES  --in $ASM --download_path $BUSCO_LINEAGES
	    fi
    else 
      busco -m genome -l $LINEAGE -c $CPU -o ${FOLDER} --out_path ${OUTFOLDER} --offline --augustus_species $SEED_SPECIES \
  		      --in $ASM --download_path $BUSCO_LINEAGES
    fi
done
