#! /bin/bash

motif_model=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/glam2.txt
genome_dir=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tara/MAGs
tmp_dir=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tara/tmp_dir/
random_sequences=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tara/tmp_dir/random_sequences.fasta
query_score=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tara/query_scores.tsv
n=100

echo -e sequence'\t'start'\t'sequence'\t'stop'\t'strand'\t'bitscore >$query_score
rm -fr $tmp_dir 2>/dev/null
mkdir $tmp_dir

for f in genome_dir/TOBG*; do
	while read -r line1; read -r line2; do

		echo $line1 > $tmp_dir/tmp
		echo $line2 >> $tmp_dir/tmp
		fasta-shuffle-letters -line $n -copies $n $tmp_dir/random > $random_sequences
		glam2scan -2 -t -n 1 n $motif_model $tmp_dir/tmp | sed '1,6d;8d' | sed 's/\s\{1,\}/\t/g' >> $query_score

	done< <(perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' $f)
done
