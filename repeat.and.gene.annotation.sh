#!/usr/bin/env bash
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AUX="${REPO_ROOT}/scripts"

#repeat sequence annotation
earlGrey -g path/PPA.final.genome.fna -s PPA_carnivora -o ./ -t 40 -r carnivora -c yes -d yes -e yes

#Ab initio prediction
#ANNEVO
cd path/PPA_carnivora
grep -vE 'Simple_repeat|Low_complexity' PPA_carnivora.filteredRepeats.bed > PPA_carnivora.filteredRepeats.noTR.bed
bedtools maskfasta -fi path/PPA.final.genome.fna -bed PPA_carnivora.filteredRepeats.noTR.bed -fo PPA_carnivora.softmasked.noTR.fasta -soft
python path/ANNEVO/annotation.py --genome path/PPA_carnivora/PPA_carnivora.softmasked.noTR.fasta --model_path path/ANNEVO/ANNEVO_model/ANNEVO_Mammalia2.pt --output PPA.ANNEVO.gff --threads 100 &> PPA_ANNEVO.log

#Tiberius
python path/Tiberius/tiberius.py --genome path/PPA_carnivora/PPA_carnivora.softmasked.noTR.fasta --out PPA.Tiberius.gtf --model_cfg path/Tiberius/model_cfg/mammalia_softmasking_v2.yaml --batch_size 8
perl "${AUX}/Tiberius_gtf2gff.pl" PPA.Tiberius.gtf > PPA.Tiberius.gff

#Augustus
augustus --species=human path/PPA_carnivora/PPA_carnivora.softmasked.noTR.fasta | perl "${AUX}/ConvertFormat_augustus.pl" - PPA.augustus.gff

#Homology-based prediction
for i in bbear cat cattle dog gp human pig sheep;
do
echo "java -Xmx150g -jar path/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=50 AnnotationFinalizer.r=NO pc=true o=true t=path/PPA.final.genome.fna GAF.a=\"iAA>=0.75\" i=$i a=path/homology_species/${i}.genome.gff g=path/homology_species/${i}.genome.fna outdir=path/$i GeMoMa.m=500000; perl \"${AUX}/ConvertFormat_GeMoMa.pl\" path/$i/final_annotation.gff"
done > run_GeMoMa.sh
echo "#!/usr/bin/env bash
#SBATCH --job-name GMM
#SBATCH -c 50 --mem 150G
#SBATCH --get-user-env
#SBATCH -o path/GMM.slurm.out
#SBATCH -e path/GMM.slurm.err
#SBATCH -D path
#SBATCH -N 1
#SBATCH --partition=bnode
bash run_GeMoMa.sh
" > GeMoMa_sbatch.sh
sbatch GeMoMa_sbatch.sh

#RNA-seq annotation
## preparing the cleaned reads using fastp with deafalut paramaters.
## The name format of read files: tissue_name_R1.fastq.gz, tissue_name_R2.fastq.gz.
perl "${AUX}/RNAseq_annot_pipe.pl" --fq_dir ./RNA_seq --exec_type submit --run_pasa yes --genome path/PPA_carnivora/PPA_carnivora.softmasked.noTR.fasta --pasa_path path/pasa/opt/pasa-2.5.3 --species PPA --max_intronlen 500000 --genome_guided --threads 40 --partition fast_debug

#Intergrating these results
##The EVidenceModeler was used to combine these gff3 files derived from three level evidences.
