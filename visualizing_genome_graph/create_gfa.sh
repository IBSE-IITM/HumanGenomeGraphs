#PBS -S /bin/bash
#PBS -l nodes=cn4:ppn=1
#PBS -N gfa_creation
#PBS -J 1-322
#PBS -l mem=110gb
###PBS -o gfa_creation.out   # uncomment twice if you want all output and error to be in one file
###PBS -e gfa_creation.err

cd $PBS_O_WORKDIR/


#echo "$(date):The file used is $(sed -n "${PBS_ARRAY_INDEX} p" /apps/venkatesh/genome_graph/bash_scripts/CCMB_sample_names.txt)" >> test_log.txt


bed_file=/cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/gfa_files/10MB_subgraphs/10M_windows.bed

path=$(head -n ${PBS_ARRAY_INDEX} $bed_file | tail -n 1 | cut -f 1) # For all rows, column 1 has chr
start_position=$(head -n ${PBS_ARRAY_INDEX} $bed_file | tail -n 1 | cut -f 2)
let "start_position=start_position+1"
end_position=$(head -n ${PBS_ARRAY_INDEX} $bed_file | tail -n 1 | cut -f 3)
let "end_position=end_position+1"
vg find -x /cn4/data4/venkatesh_data4/1KGP_genome_graphs/complete_1KGP_graph/vg_manual_construct/1KGP_complete_graph.xg -p $path:$start_position-$end_position -E | vg view -v - --gfa > gfa_files/${path}_${start_position}_${end_position}.gfa