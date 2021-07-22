#! /bin/bash

#Inputs
seed_sequences=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/kher.fasta
sequence_file=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/april_fasta_sequences_200.fasta
seq_len=200
n=100

#Subdirectories
motif_model=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model
alignments=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/alignment
glam_alignment=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/alignment/alignment
tmp_dir=/home/troyalty/Documents/projects/luxbox_motif/data/meme_suite/glam2/motif_model/tmp_dir/

#Files
random_sequences=random_sequences.fasta
random_score=random_scores.tsv
query_score=query_scores.tsv
new_sequences=new_sequeneces.fasta

#iterators
i=0
initial=$(echo _$i/)

#Configure subdirectories
rm -fr $motif_model 2>/dev/null
mkdir $motif_model
mkdir $alignments
mkdir $tmp_dir


################################################################################################################
#Initialize an motif alignment based on seed sequences
################################################################################################################
glam2 -Q -O $glam_alignment$initial -2 -z 2 -a 15 -b 20 -w 10 -r 10 -n 8000 -D 0.1 -E 2.0 -I 0.02 -J 1.0 n $seed_sequences
align_score_o=$(grep "Score:" $glam_alignment$initial/glam2.txt | head -n 1 | cut -f 2 -d ' ')
cp $seed_sequences $glam_alignment$initial$new_sequences

echo -e sequence'\t'start'\t'sequence'\t'stop'\t'strand'\t'bitscore > $glam_alignment$initial$random_score
>$glam_alignment$initial$query_score

fasta-shuffle-letters -line $seq_len -copies $n $sequence_file > $glam_alignment$initial$random_sequences

while read -r line1; read -r line2; do

	echo $line1 > $tmp_dir/tmp
	echo $line2 >> $tmp_dir/tmp
	glam2scan -2 -t -n 1 n $glam_alignment$initial/glam2.txt $tmp_dir/tmp | sed '1,6d;8d' | sed 's/\s\{1,\}/\t/g' >> $glam_alignment$initial$random_score

done<$glam_alignment$initial$random_sequences

bitscore_95=$(awk '{print $6}' $glam_alignment$initial$random_score | sed '1d' | sort -n | perl -e '$d=.95;@l=<>;print $l[int($d*$#l)]')

while read -r line1; read -r line2; do

	echo $line1 > $tmp_dir/tmp
	echo $line2 >> $tmp_dir/tmp
	glam2scan -2 -t -n 1 n $glam_alignment$initial/glam2.txt $tmp_dir/tmp | sed '1,6d;8d' | sed 's/\s\{1,\}/\t/g' | sed "s/$/\t$bitscore_95/" >> $glam_alignment$initial$query_score

done<$sequence_file


while read id start sequence stop strand bitscore; do

	echo \>$id >> $glam_alignment$initial$new_sequences
	echo $sequence | sed -e 's/\(.*\)/\U\1/' >> $glam_alignment$initial$new_sequences

done< <(cat $glam_alignment$initial$query_score | awk -v cutoff=$bitscore_95 '$6>cutoff')
################################################################################################################
################################################################################################################

################################################################################################################
#Dynamically modify motif using input sequences
################################################################################################################
check=0
while [ $check -lt 2 ]; do #continue iterations until alignment score does not improve by more than 5% for two iterations in a row
	it_o=$(echo _$i/)
	let "i=i+1"
	it=$(echo _$i/)
	glam2 -Q -O $glam_alignment$it -2 -z 2 -a 15 -b 20 -w 10 -r 10 -n 8000 -D 0.1 -E 2.0 -I 0.02 -J 1.0 n $glam_alignment$it_o$new_sequences
	
	align_score_n=$(grep "Score:" $glam_alignment$it/glam2.txt | head -n 1 | cut -f 2 -d ' ') #find alignment score
	d_score=$(echo "($align_score_n - $align_score_o)/$align_score_o" | bc -l) #calculate change in alignment score
	echo $align_score_o,$align_score_n,$d_score
	align_score_o=$align_score_n #set old alignment score to newly calculated alignment score
	
	if (( $(echo "$d_score <= 0.05" | bc) )); then #increase check if score did not improve by at least 5%
		let "check=check+1"
	else
		check=0
	fi
	
	cp $seed_sequences $glam_alignment$it$new_sequences

	echo -e sequence'\t'start'\t'sequence'\t'stop'\t'strand'\t'bitscore > $glam_alignment$it$random_score
	>$glam_alignment$it$query_score

	fasta-shuffle-letters -line $seq_len -copies $n $sequence_file > $glam_alignment$it$random_sequences

	while read -r line1; read -r line2; do

		echo $line1 > $tmp_dir/tmp
		echo $line2 >> $tmp_dir/tmp
		glam2scan -2 -t -n 1 n $glam_alignment$it/glam2.txt $tmp_dir/tmp | sed '1,6d;8d' | sed 's/\s\{1,\}/\t/g' >> $glam_alignment$it$random_score

	done<$glam_alignment$it$random_sequences

	bitscore_95=$(awk '{print $6}' $glam_alignment$it$random_score | sed '1d' | sort -n | perl -e '$d=.95;@l=<>;print $l[int($d*$#l)]')

	while read -r line1; read -r line2; do

		echo $line1 > $tmp_dir/tmp
		echo $line2 >> $tmp_dir/tmp
		glam2scan -2 -t -n 1 n $glam_alignment$it/glam2.txt $tmp_dir/tmp | sed '1,6d;8d' | sed 's/\s\{1,\}/\t/g' | sed "s/$/\t$bitscore_95/" >> $glam_alignment$it$query_score

	done<$sequence_file


	while read id start sequence stop strand bitscore; do

		echo \>$id >> $glam_alignment$it$new_sequences
		echo $sequence | sed -e 's/\(.*\)/\U\1/' >> $glam_alignment$it$new_sequences

	done< <(cat $glam_alignment$it$query_score | awk -v cutoff=$bitscore_95 '$6>cutoff')

done
################################################################################################################
################################################################################################################

let "it=it-1"
cp $glam_alignment$it/glam2.txt $motif_model
