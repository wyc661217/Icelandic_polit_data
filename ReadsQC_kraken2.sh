#This is for processing the icelandic polit data generated in 2 flowcells

##############################
#reads QC
##############################

##adpter residue sequences
echo ">adapter1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>adapter2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>residue1
AGATCGGAAGAG
>residue2
CTCCAGTCAC
>residue3
GAAAGAGTGT
>residue4
ATCTCGTATG
>residue5
GTGTAGATCT
" > adpater_list.fa

#the two data folder
 ../SCIENCE-SNM-Sekvenslab-FASTQData/221212_A00706_0681_BH3TMMDSX5_Xihan_WRGRX/IcelandicPilot-eDNALib029
 ../SCIENCE-SNM-Sekvenslab-FASTQData/221212_A00706_0680_AHYLLHDSX3_Xihan_ZRR4B/IcelandicPilot-eDNALib028

#adpter removal and concatinate pair-ended
for file in *_R1_001.fastq.gz
do

  file2=`echo "${file/_R1/_R2}"`
  bname=`echo $file | cut -d"_" -f1-3`
  echo -e "\n\nprocessing" $bname "\n\n"

  AdapterRemoval --file1 $file --file2 $file2 --basename $bname --threads 20 \
  --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
  --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
  --trimns --trimqualities --mm 3 --collapse --minalignmentlength 8 --gzip

  zcat $bname.collapsed.gz $bname.pair1.truncated.gz $bname.pair2.truncated.gz $bname.singleton.truncated.gz > $bname.adRemoved.fq
  rm $bname.collapsed.gz $bname.collapsed.truncated.gz $bname.discarded.gz $bname.pair1.truncated.gz $bname.pair2.truncated.gz $bname.singleton.truncated.gz

  fastp --in1 $bname.adRemoved.fq --out1 $bname.residue_cleaned.fq \
  --adapter_fasta adpater_list.fa --html $bname.fastp.html --json $bname.fastp.json \
  --dedup --dont_eval_duplication --low_complexity_filter -Q -l 30 -w 20 \
  --trim_poly_g --poly_g_min_len 5 --trim_poly_x --poly_x_min_len 5

  gzip $bname.residue_cleaned.fq
  rm $bname.adRemoved.fq

done &> adapter_removal_log.txt

#merge 4 lanes, duplicates and low-complexity removal
for file in *_L001.residue_cleaned.fq.gz
do

  bname=`echo $file | cut -d"_" -f1-2`
  echo -e "\n\nprocessing" $bname "\n\n"

  zcat ${bname}_*.residue_cleaned.fq.gz > $bname.L1234.fq
  rm ${bname}_*.residue_cleaned.fq.gz

  sga preprocess --dust-threshold=30 -m 30 $bname.L1234.fq -o $bname.prep.fq
  sga index --algorithm=ropebwt -t 20 --threads=40 $bname.prep.fq
  sga filter --threads=40 --no-kmer-check $bname.prep.fq -o $bname.cleaned.fq

  rm $bname.L1234.fq *.bwt *.discard.fa *.rbwt *.rsai *.sai $bname.prep.fq
  gzip $bname.cleaned.fq

done &> cleaning_log.txt



##############################
#kranken2 subsetting
##############################

#download NCBI nt (16th Nov 2022) and Refseq (release 214), construct kraken2 db on the UK-crop server
kraken2-build --download-taxonomy --use-ftp --db kraken2
kraken2-build --download-library nt --db kraken2 --threads 40 --no-masking --use-ftp
kraken2-build --add-to-library bacteria.fa --db kraken2 --threads 40 --no-masking
kraken2-build --add-to-library others_prok.fa --db kraken2 --threads 40 --no-masking
kraken2-build --add-to-library others_euka.fa --db kraken2 --threads 40 --no-masking
kraken2-build --add-to-library vert_mam.fa --db kraken2 --threads 40 --no-masking
kraken2-build --add-to-library invert.fa --db kraken2 --threads 40 --no-masking
kraken2-build --add-to-library vert_other.fa --db kraken2 --threads 40 --no-masking
kraken2-build --add-to-library complete.fa --db kraken2 --threads 40 --no-masking
kraken2-build --build --db kraken2 --kmer-len 28 --minimizer-len 27 --threads 40 --minimizer-spaces 6

#keaken2 claasification
kraken2 --paired --threads 40 -db $combined_db --confidence 1 --output output/ lib028/*.cleaned.fq.gz
kraken2 --paired --threads 40 -db $combined_db --confidence 1 --output output/ lib029/*.cleaned.fq.gz


#merge kraken2 classified reads of different sampels into one big fq
mkdir classified_read
for path in *gz_out
do
  bname=`basename $path | cut -d"-" -f4 | cut -d"_" -f1`
  bname=`echo "CGG3_"$bname`
  echo $bname
  mv $path/classified_out classified_read/$bname.fq
done

cat *.fq > lib28.fastq


#mapping
##the R script for subseting the bam headers
echo "args <- commandArgs(TRUE)
header_full <- args[1]
header_new <- args[2]
df1 = read.csv(header_new, header = F, stringsAsFactors = F)
df1 = unique(df1[,1])
df1 = paste("SN:",df1,sep="")
header = read.csv(header_full, header = F, stringsAsFactors = F, sep = "\t")
header = header[which(header[,2] %in% df1),]
write.table(header, "header_subset.2.txt", row.names = F, col.names = F, quote = F,sep = "\t")" > bam_header_subset.R

##mapping with slurm
###this solution that subset reads for bowtie2 mapping does not save much time, give up
sbatch s_i2.2_mapping.sh
--------------------
#!/bin/bash
#SBATCH --job-name=i2.2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=22
#SBATCH --exclude=dandycmpn01fl,dandycmpn02fl,dandycmpn03fl,dandycmpn04fl
#SBATCH --mem=250000
module load bowtie2/2.4.2
module load gcc/11.2.0
module load R/4.2.2
module load samtools/1.12

ID=4
cd /home/zsb202/project/icaland/kraken2/batch_2_conflim_1/classified_read
database="/home/zsb202/shared_db/NCBI_nt_refseq_20221206_ycw/db_list/list.$ID.txt"
for file in lib29_b2.fastq.gz
do

  bname=$(basename "$file" | cut -d. -f1)
  echo -e "\n\nprocessing" $bname "\n\n"
  mkdir $bname.$ID
  cd $bname.$ID

  cat $database | while read DB
  do

    bDB=$(basename $DB)
    bowtie2 --threads 22 -k 1000 -x $DB -U ../$file --no-unal -S ${bname}.$bDB.sam 2> ${bname}.$bDB.log.txt

    grep "@HD" ${bname}.$bDB.sam > header_subset.1.txt
    grep "@SQ" ${bname}.$bDB.sam > header_subset.2.txt
    awk '!seen[$0]++' header_subset.2.txt > header_subset.2.1.txt
    rm header_subset.2.txt
    mv header_subset.2.1.txt header_subset.2.txt
    grep "@PG" ${bname}.$bDB.sam > header_subset.3.txt
    grep '^@' -v ${bname}.$bDB.sam > alignment.txt
    rm ${bname}.$bDB.sam

    cat alignment.txt | cut -f3 | sort -u > header_new.txt
    Rscript /home/zsb202/script/bam_header_subset.R header_subset.2.txt header_new.txt
    rm header_new.txt
    cat header_subset.1.txt header_subset.2.txt header_subset.3.txt alignment.txt | samtools view -@ 30 -b -o ${bname}.$bDB.bam
    rm header_subset.1.txt header_subset.2.txt header_subset.3.txt alignment.txt

  done

  for txt in *.log.txt
  do
    echo "$txt" >> $bname.nt_refseq_log.txt
    cat $txt >> $bname.nt_refseq_log.txt
  done

  rm *.log.txt
  cd ..

done &>alngnment_log.i2.$ID.txt
--------------------


##############################
#here to play the kraken2 results
##############################
#this approach need to be refined

##generate taxa list in R
--------------------
#download the flora and fauna list from ni.is
#library(readr)
#ni_is_list <- read_table2("Desktop/Iceland/taxa_list/ni_is_list.txt",col_names = FALSE)
#ni_is_list = ni_is_list[!duplicated(ni_is_list[,1]),] #240 genus
#write.table(ni_is_list,"Desktop/niis_genus.txt",col.names = F,row.names = F,quote = F)

#download the whole panel list of iceland taxa from gbif
#GBIF_all_panel <- read.delim("~/Desktop/Iceland/taxa_list/GBIF_all_panel.tsv", quote="")
#GBIF_all_panel = GBIF_all_panel[c(which(GBIF_all_panel$taxonRank == "GENUS"),which(GBIF_all_panel$taxonRank == "SPECIES")),]
#GBIF_all_panel = GBIF_all_panel[!duplicated(GBIF_all_panel$genus),] #5309 genus
#GBIF_all_panel = GBIF_all_panel[,18]
#write.table(GBIF_all_panel,"Desktop/gbif_genus.txt",col.names = F,row.names = F,quote = F)
--------------------

##generate the taxaID list
cd /home/zsb202/project/icaland/kraken2/taxa_list
/home/zsb202/app/taxonkit/taxonkit name2taxid niis_genus.txt --data-dir /home/zsb202/shared_db/NCBI_nt_refseq_20221206_ycw/taxonomy -o niis_genus.taxaID.txt

re='^[0-9]+$'
cat /home/zsb202/project/icaland/kraken2/taxa_list/niis_genus.taxaID.txt | while read line
do
  tn=`echo $line | cut -d" " -f1`
  ti=`echo $line | cut -d" " -f2`
  if [[ $ti =~ $re ]] ; then
    echo $ti >> IDlist1.txt
  fi
done

cat IDlist1.txt | uniq > IDlist2.txt
sed 's/^/ /' IDlist2.txt > IDlist3.txt
sed 's/$/:/' IDlist3.txt > IDlist4.txt

mv IDlist2.txt niis_genus.taxaID_list.txt
mv IDlist4.txt niis_genus.taxaID_list_formatted.txt


##found that too less taxa were identified for plant, here add the species also into the list
cat niis_genus.txt ni_is_plant_list.txt > niis_genus_and_species_list.txt
/home/zsb202/app/taxonkit/taxonkit name2taxid niis_genus_and_species_list.txt \
--data-dir /home/zsb202/shared_db/NCBI_nt_refseq_20221206_ycw/taxonomy -o niis_genus_and_species_list.taxaID.txt
/home/zsb202/app/taxonkit/taxonkit name2taxid niis_genus_and_species_list.txt \
--data-dir /home/zsb202/shared_db/NCBI_nt_refseq_20221206_ycw/taxonomy -o niis_genus_and_species_list.with_taxaID.txt

re='^[0-9]+$'
cat /home/zsb202/project/icaland/kraken2/taxa_list/niis_genus_and_species_list.taxaID.txt | while read line
do
  tn=`echo $line | cut -d" " -f1`
  ti=`echo $line | cut -d" " -f2`
  if [[ $ti =~ $re ]] ; then
    echo $ti >> IDlist1.txt
  fi
done

cat IDlist1.txt | uniq > IDlist2.txt
sed 's/^/ /' IDlist2.txt > IDlist3.txt
sed 's/$/:/' IDlist3.txt > IDlist4.txt

mv IDlist2.txt niis_genus_and_species.taxaID_list.txt
mv IDlist4.txt niis_genus_and_species.taxaID_list_formatted.txt
rm IDlist1.txt IDlist3.txt


##extract kmer frequency
cd /home/zsb202/project/icaland/kraken2/lib028
mkdir controls
mv LibNTC_S96_L003.residue_cleaned.fq.gz_out LV7001884361-LibPTC22101201-c-LibPTC_S95.cleaned.fq.gz_out LV7001884367-ExrPTC22101201-c-ExrPTC_S47.cleaned.fq.gz_out LV7001884379-ExrNTC22101201-c-ExrNTC_S48.cleaned.fq.gz_out controls/
mv LV7001883389-LibPTC22101201-LibPTC_S95.cleaned.fq.gz_out LV7001883389-LibPTC22101201-LibPTC_S95.cleaned.fq.gz_out LV7001883395-ExrPTC22101201-ExrPTC_S47.cleaned.fq.gz_out LV7001883407-ExrNTC22101201-ExrNTC_S48.cleaned.fq.gz_out controls/

for path in *.gz_out
do
  bname=`basename $path | cut -d"-" -f4 | cut -d"_" -f1`
  bname=`echo "CGG3_"$bname`
  echo processing $bname now

  cd $path
  grep "^C" k2_output.txt > $bname.koutput.txt
  grep -on -f /home/zsb202/project/icaland/kraken2/taxa_list/niis_genus_and_species.taxaID_list_formatted.txt \
  $bname.koutput.txt | uniq | awk -F: '{count[$2]++} END {for (word in count) print word, count[word]}' > $bname.niis.kmer_freq.txt
  cd ..
done

##gzip files
for path in *.gz_out
do
  cd $path
  bname=`basename $path | cut -d"-" -f4 | cut -d"_" -f1`
  bname=`echo "CGG3_"$bname`
  echo processing $bname now

  rm k2_output.txt
  line=`wc -l *.koutput.txt`
  gzip *.koutput.txt
  mv classified_out $bname.fq
  gzip $bname.fq
  mv $bname.fq.gz ../classified_read/
  echo $bname,$line >> ../classified_read.count.txt
  cd ..
done

##sort files
cd /home/zsb202/project/icaland/kraken2/lib029
mkdir niis_restuls
for path in *gz_out
do
  mv $path/*.niis.kmer_freq.txt niis_restuls/
done


##############################
#reads qc statistics
##############################

##adpter removal and fastp
echo "sample,raw_reads,adpter_removed,fastp_input,fastp_after_filtering,fastp_passed_filters,fastp_too_short,fastp_low complexity" >QC_report1.csv
for file in *.settings
do
  bname=`basename $file | cut -d. -f1`
  echo $bname
  n1=`grep "Total number of read pairs:" $file | cut -d" " -f6`
  n2=`grep "Number of retained reads:" $file | cut -d" " -f5`
  n3=`grep "total reads:" $bname.fastp.html | head -n1 | cut -d">" -f5 | cut -d" " -f1`
  n4=`grep "total reads:" $bname.fastp.html | tail -n1 | cut -d">" -f5 | cut -d" " -f1`
  n5=`grep "reads passed filters:" $bname.fastp.html | cut -d">" -f5 | cut -d" " -f1`
  n6=`grep "reads too short:" $bname.fastp.html | cut -d"(" -f2 | cut -d")" -f1`
  n7=`grep "low complexity:" $bname.fastp.html | cut -d"(" -f2 | cut -d")" -f1`
  echo $bname,$n1,$n2,$n3,$n4,$n5,$n6,$n7 >>QC_report1.csv
done


##the final cleaned files
cd /home/zsb202/project/icaland/fq/lib029/cleaned_fq_lib028
/home/zsb202/app/seqkit/seqkit stats -a *.fq.gz > lib028.cleaned.statreport.txt
cd /home/zsb202/project/icaland/fq/lib029/cleaned_fq_lib029
/home/zsb202/app/seqkit/seqkit stats -a *.fq.gz > lib029.cleaned.statreport.txt

##read length distrubution
for file in *.gz
do
  zcat $file | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file.read_length.txt
done


##the duplicates between two flowcells
cd /home/zsb202/eDNA/iceland/lib028/cleaned_fq_lib028
cat 15.txt | while read file
do
  bname=`echo $file | cut -d"_" -f1 | cut -d"-" -f3,4`
  echo processing $bname now

  mkdir $bname
  zcat $file /home/zsb202/eDNA/iceland/lib029/cleaned_fq_lib029/*$bname* > $bname/$bname.combined.fq
  cd $bname

  sga preprocess --no-primer-check --dust-threshold=30 -m 30 $bname.combined.fq -o $bname.prep.fq
  /home/zsb202/data/soft/seqkit/seqkit stats -a $bname.prep.fq >> ../2flowcells_merged.complexfilter30.txt

  sga index --algorithm=ropebwt -t 45 --threads=40 $bname.prep.fq
  sga filter --threads=45 --no-kmer-check $bname.prep.fq -o $bname.cleaned.fq
  /home/zsb202/data/soft/seqkit/seqkit stats -a $bname.cleaned.fq >> ../2flowcells_merged.statreport.txt

  sga preprocess --no-primer-check --dust-threshold=30 -m 3 $bname.combined.fq -o $bname.prep1.fq
  /home/zsb202/data/soft/seqkit/seqkit stats -a $bname.prep1.fq >> ../2flowcells_merged.complexfilter3.txt

  cd ..
  rm -r $bname
done
