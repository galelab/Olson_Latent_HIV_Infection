#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -t 5800-12
#SBATCH --partition HoldingPen
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out

#AUTHOR: Leanne Whitmore 
#DESCRIPTION: HIVsc_latency Alignment
#NOTE: change from 5th to 6th run is the HIV genome being used (HIV pNL4-3)

PATH=$PATH:/share/tools/10x/cellranger-4.0.0
export PATH

GENOME_DIR="./refdata-cellranger-hg19-HIVpNL4-3_HIVonly-4.0.0/"
PATH_FASTQ="./cell-Ranger-demultiplex/run01/HY7WYBGXC/outs/fastq_path/HY7WYBGXC/"
PATH_FASTQ_run02="./cell-Ranger-demultiplex/run02/HY5WTBGXC/outs/fastq_path/HY5WTBGXC/"
PATH_FASTQ_run03="./cell-Ranger-demultiplex/run03/HYC32BGXC/outs/fastq_path/HYC32BGXC/"
PATH_FASTQ_run04="./cell-Ranger-demultiplex/run04/HY75YBGXC/outs/fastq_path/HY75YBGXC/"
PATH_FASTQ_run05="./cell-Ranger-demultiplex/run05/HY5Y3BGXC/outs/fastq_path/HY5Y3BGXC/"
PATH_FASTQ_run06="./cell-Ranger-demultiplex/run06/HYC7VBGXF/outs/fastq_path/HYC7VBGXF/"
PATH_FASTQ_run07="./cell-Ranger-demultiplex/run07/H2HGMBGXG/outs/fastq_path/H2HGMBGXG/"
PATH_FASTQ_run08="./cell-Ranger-demultiplex/run08/HHTWLBGXF/outs/fastq_path/HHTWLBGXF/"

##PARAMETERS FOR cellranger count 
#--expect-cells - got 8,000 thefrom 10x sample processing folder Z:\New U Drive\Projects\P51 Projects\HIV and SIV\HIV01_scLatency\Sample Processing\HIV01_scLatency_tracking_sheet_20201002_es.xlsx
#--sample=Sample name (each sample needs its own cellranger count command)
#--id=name for each count run (should probably be different each run)
#--fastqs=path to fastq files (cellrange creates a link in its own folder named fastq_path to the fastq folders)
#--transcriptome=path to reference transcriptome 
#--jobmode=slurm (tells cellranger to work slurm)

samples=("$PATH_FASTQ"G00*)
samplesrun02=("$PATH_FASTQ_run02"G00*)
samplesrun03=("$PATH_FASTQ_run03"G00*)
samplesrun04=("$PATH_FASTQ_run04"G00*)
samplesrun05=("$PATH_FASTQ_run05"G00*)
samplesrun06=("$PATH_FASTQ_run06"G00*)
samplesrun07=("$PATH_FASTQ_run07"G00*)
samplesrun08=("$PATH_FASTQ_run08"G00*)
echo "Number of samples = ${#samples[*]} (should be 8)"
echo "Number of expected cell numbers = ${#expectedcells[*]} (should be 8)"
 
###G001    
tempsample="${samples[0]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun05[0]}"
sample_name2=${tempsample#$PATH_FASTQ_run05}
tempsample="${samplesrun06[0]}"
sample_name3=${tempsample#$PATH_FASTQ_run06}
tempsample="${samplesrun08[0]}"
sample_name4=${tempsample#$PATH_FASTQ_run08}
echo "Processing sample $sample_name, $sample_name2, $sample_name3 and $sample_name4"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[0]},${samplesrun05[0]},${samplesrun06[0]},${samplesrun08[0]} & 
wait
    
###G002
tempsample="${samples[1]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun05[1]}"
sample_name2=${tempsample#$PATH_FASTQ_run05}
tempsample="${samplesrun06[1]}"
sample_name3=${tempsample#$PATH_FASTQ_run06}
tempsample="${samplesrun08[1]}"
sample_name4=${tempsample#$PATH_FASTQ_run08}
echo "Processing sample $sample_name, $sample_name2, $sample_name3 and $sample_name4"

srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[1]},${samplesrun05[1]},${samplesrun06[1]},${samplesrun08[1]} &             
wait

###G003
tempsample="${samples[2]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun02[0]}"
sample_name2=${tempsample#$PATH_FASTQ_run02}
tempsample="${samplesrun04[0]}"
sample_name3=${tempsample#$PATH_FASTQ_run04}
tempsample="${samplesrun06[2]}"
sample_name4=${tempsample#$PATH_FASTQ_run06}
tempsample="${samplesrun08[2]}"
sample_name5=${tempsample#$PATH_FASTQ_run08}

echo "Processing sample $sample_name, $sample_name2, $sample_name3,  $sample_name4 and $sample_name5"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4,$sample_name5 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[2]},${samplesrun02[0]},${samplesrun04[0]},${samplesrun06[2]},${samplesrun08[2]} &     
wait

###G004
tempsample="${samples[3]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun03[0]}"
sample_name2=${tempsample#$PATH_FASTQ_run03}
tempsample="${samplesrun04[1]}"
sample_name3=${tempsample#$PATH_FASTQ_run04}
tempsample="${samplesrun08[3]}"
sample_name4=${tempsample#$PATH_FASTQ_run08}

echo "Processing sample $sample_name, $sample_name2, $sample_name3, $sample_name4"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4  --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[3]},${samplesrun03[0]},${samplesrun04[1]},${samplesrun08[3]} &
wait

###G005
tempsample="${samples[4]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun05[2]}"
sample_name2=${tempsample#$PATH_FASTQ_run05}
tempsample="${samplesrun07[0]}"
sample_name3=${tempsample#$PATH_FASTQ_run07}

echo "Processing sample $sample_name, $sample_name2, $sample_name3"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[4]},${samplesrun05[2]},${samplesrun07[0]} &
wait

###G006
tempsample="${samples[5]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun05[3]}"
sample_name2=${tempsample#$PATH_FASTQ_run05}
tempsample="${samplesrun07[1]}"
sample_name3=${tempsample#$PATH_FASTQ_run07}
tempsample="${samplesrun08[4]}"
sample_name4=${tempsample#$PATH_FASTQ_run08}

echo "Processing sample $sample_name, $sample_name2, $sample_name3 and $sample_name4"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[5]},${samplesrun05[3]},${samplesrun07[1]},${samplesrun08[4]} &
wait

###G007
tempsample="${samples[6]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun02[1]}"
sample_name2=${tempsample#$PATH_FASTQ_run02}
tempsample="${samplesrun04[2]}"
sample_name3=${tempsample#$PATH_FASTQ_run04}
tempsample="${samplesrun08[5]}"
sample_name4=${tempsample#$PATH_FASTQ_run08}

echo "Processing sample $sample_name, $sample_name2, $sample_name3 and$sample_name4"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[6]},${samplesrun02[1]},${samplesrun04[2]},${samplesrun08[5]} &
wait

###G008
tempsample="${samples[7]}"
sample_name=${tempsample#$PATH_FASTQ}
tempsample="${samplesrun03[1]}"
sample_name2=${tempsample#$PATH_FASTQ_run03}
tempsample="${samplesrun04[3]}"
sample_name3=${tempsample#$PATH_FASTQ_run04}
tempsample="${samplesrun05[4]}"
sample_name4=${tempsample#$PATH_FASTQ_run05}
tempsample="${samplesrun08[6]}"
sample_name5=${tempsample#$PATH_FASTQ_run08}

echo "Processing sample $sample_name, $sample_name2, $sample_name3, and $sample_name4, $sample_name5"
srun -c 2 cellranger count --sample=$sample_name,$sample_name2,$sample_name3,$sample_name4,$sample_name5 --id=$sample_name --jobmode=slurm --expect-cells=8000 --transcriptome=$GENOME_DIR --fastqs=${samples[7]},${samplesrun03[1]},${samplesrun04[3]},${samplesrun05[4]},${samplesrun08[6]} &
wait

mutt -s "Cellranger count done for HIV SC latency run 7" lwhitmo@uw.edu -c lwhitmo@uw.edu
