#! /bin/bash

export motif_model=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/glam2.txt
export tmp_dir=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tmp_dir/
#export n=100

genome_dir=/home/troyalty/Documents/projects/luxbox_motif/data/tara_MAGs/
#random_sequences=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tara/tmp_dir/random_sequences.fasta
#random_scores=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tara/tmp_dir/random_scores.tsv
cutoff_score=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/query_scores.tsv

#echo -e sequence'\t'start'\t'sequence'\t'stop'\t'strand'\t'bitscore'\t'signifiance_cutoff'\t'genome >$query_score
echo -e bitscore'\t'genome'\t'gc_content'\t'genome_size'\t'replicate > $cutoff_score
rm -fr $tmp_dir 2>/dev/null
mkdir $tmp_dir




function parallel_query_glam2scan {
	file=$1
	tmp=.tmp #2 line orginal sequences
	tmp_genome=$(basename $file | cut -f 1 -d '.')
	tmp_file=$tmp_dir$tmp_genome$tmp
	>$tmp_file
	seq_len=$(bioawk -c fastx 'BEGIN {seqlen=0}{seqlen+=length($seq)}END{print seqlen}' $file)
	gc=$(perl -lane 'unless (/^>/) { $l += length(); $gc++ while /[GC]/ig } END { print $gc/$l}' $file)
	for i in {1..1000}; do
		glam2scan -2 -t -n 1 n $motif_model <(fasta-shuffle-letters -copies 1 $file) |\
		sed '1,6d;8d' |\
		sed 's/\s\{1,\}/\t/g' |\
		cut -f 6 |\
		sed "s/$/\t$tmp_genome/" |\
		sed "s/$/\t$gc/" |\
		sed "s/$/\t$seq_len/" |\
		sed "s/$/\t$i/"  >> $tmp_file
	done
	cat $tmp_file
}

export -f parallel_query_glam2scan


parallel -j 30 --eta "parallel_query_glam2scan {}" ::: $genome_dir/* >> $cutoff_score
