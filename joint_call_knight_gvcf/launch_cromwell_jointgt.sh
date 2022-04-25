#!/bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J launch_cromwell
#SBATCH --mem=16000

java -Dconfig.file=slurm.conf -jar \
/gpfs/ycga/project/lek/shared/tools/jars/cromwell-56.jar run \
Multisample_jointgt_GATK4.wdl \
-i Multisample_jointgt_GATK4_inputs.json \
-o cromwell.options
