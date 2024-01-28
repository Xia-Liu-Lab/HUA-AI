#Fastp (trim and filter adaptors)
mkdir clean_data
for j in `cat checklist.csv`;do
  {
  cd /clean_data/
  mkdir ${j}
  fastp -i /rawdata/${j}/${j}_1.fastq.gz -o /clean_data/${j}/${j}_1.fastq.gz -I rawdata/${j}/${j}_2.fastq.gz -O /clean_data/${j}/${j}_2.fastq.gz --html /clean_data/${j}/${j}.html --json /clean_data/${j}/${j}.json
    }
done

#Bowtie2 (remove human contaminants)
mkdir unaligned_data
for j in `cat checklist.csv`;do
  {
  cd /unaligned_data/
  mkdir ${j}
  bowtie2 -p 15 -x /database/GRCh38/GRCh38 -1 /clean_data/${j}/${j}_1.fastq.gz -2 /clean_data/${j}/${j}_2.fastq.gz -S /unaligned_data/${j}/${j}_L1.sam --un-conc /unaligned_data/${j}/unaligned_${j}.fastq
  cd /unaligned_data/${j}
  rm ${j}_L1.sam
    gzip *
    cat unaligned_${j}.1.fastq.gz unaligned_${j}.2.fastq.gz > unaligned_${j}.fastq.gz
    }
done


#metaphlan (taxonomic profiling)
source activate humann
mkdir metaphlan
for j in `cat checklist.csv`;do
  {
  cd /metaphlan/
  metaphlan /unaligned_data/${j}/unaligned_${j}.fastq.gz --bowtie2out ${j}.bowtie2.bz2 --nproc 5 --input_type fastq -o ${j}.txt --bowtie2db /database/metaphlan_database/ -x mpa_v30_CHOCOPhlAn_201901
  rm ${j}.bowtie2.bz2
  }
done

#humann (functional profiling)
source activate humann
mkdir humann
for j in `cat  checklist.csv`;do
  {
  cd /humann/
  mkdir ${j}
  cd /unaligned_data/${j}/
  humann --input unaligned_${j}.fastq.gz --output /humann/${j} --metaphlan-options="--bowtie2db /database/metaphlan_database/"
  cd /humann/${j}
  humann_renorm_table --input ${j}_pathabundance.tsv --output ${j}_pathabundance_relab.tsv --units relab --special n
  rm -r ${j}_humann_temp/
  }
done

##KEGG
#megahit
mkdir megahit
for j in `cat checklist.csv`;do
  {
  cd /unaligned_data/${j}/
  megahit -1 unaligned_${j}.1.fastq.gz -2 unaligned_${j}.2.fastq.gz -o /megahit/${j}/
  }
done

#kofamscan
mkdir KEGG/translation
mkdir KEGG/kofamscan
for j in `cat checklist.csv`;do
  {
  prodigal -i /megahit/${j}/final.contigs.fa -a /KEGG/translation/${j}.fa
  cd /software/kofamscan/bin/kofam_scan-1.3.0/
  ./exec_annotation -f mapper -o /KEGG/kofamscan/${j}.txt /KEGG/translation/${j}.fa --tmp-dir=/KEGG/kofamscan/temp/${j}/

  rm -r /KEGG/kofamscan/temp/${j}/
  }
done


