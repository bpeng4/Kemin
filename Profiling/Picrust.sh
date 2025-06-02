#!/bin/bash
#SBATCH --time=168:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=5G       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=Picrust
#SBATCH --error=/work/benson/bpeng4/Kemin/Kemin_All/err/job.%J.err
#SBATCH --output=/work/benson/bpeng4/Kemin/Kemin_All/err/job.%J.output
#SBATCH --partition=benson,batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36


module load picrust2/2.4
# Use an extra cd to go to your work directory (not shown)
cd picrust2

#picrust2_pipeline.py -s seqs.fna -i raw_counts.tsv -p 2 -o results \
#                     --stratified --verbose

cd /work/benson/bpeng4/Kemin/Kemin_All/picrust2/results

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
                    
# Decompress tsv files
find . -name "*desc*" -exec pigz -dkv \{\} \;



cp ./picrust2/results/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv  ./picrust2/EC.tsv
cp ./picrust2/results/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv  ./picrust2/KO.tsv
cp ./picrust2/results/pathways_out/path_abun_unstrat_descrip.tsv  ./picrust2/path.tsv
