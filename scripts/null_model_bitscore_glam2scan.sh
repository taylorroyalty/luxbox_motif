#! /bin/bash

motif_model=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/kher_alignment.txt
sequence_file=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/april_fasta_sequences_200.fasta
random_sequences=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/random_sequences.fasta
random_dir=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/tmp_dir/
random_score=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/random_scores.tsv
n=1000

echo -e sequence'\t'start'\t'sequence'\t'stop'\t'strand'\t'bitscore >$random_score
rm -fr $random_dir 2>/dev/null
mkdir $random_dir


fasta-shuffle-letters -line $n -copies $n $sequence_file > $random_sequences

while read -r line1; read -r line2; do

	echo $line1 > $random_dir/tmp
	echo $line2 >> $random_dir/tmp
	glam2scan -2 -t -n 1 n $motif_model $random_dir/tmp | sed '1,6d;8d' | sed 's/\s\{1,\}/\t/g' >> $random_score

done<$random_sequences
