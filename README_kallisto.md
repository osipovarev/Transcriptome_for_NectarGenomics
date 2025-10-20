###################################################
#                                                 #
###      Run Kallisto RNAseq analysis            ##
#                                                 #
###################################################



## Trim raw reads
### (NB: added ploy-A and poly-T in all_adapters file!)
```
for l in $(cut -f1 batch1_lib_info.files.tsv | grep -v LibraryID);
do
 d=$(grep $l batch1_lib_info.files.tsv | cut -f5);
 echo -e "trimmomatic PE -threads 16 $d/${d}_R1.fastq.gz $d/${d}_R2.fastq.gz fastq_files_batch1/${l}.trim_R1.fastq.gz fastq_files_batch1/unpaired.${l}.trim_R1.fastq.gz fastq_files_batch1/${l}.trim_R2.fastq.gz fastq_files_batch1/unpaired.${l}.trim_R2.fastq.gz ILLUMINACLIP:adapter_all.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
done > jobs.trimm.batch1
```


## Prepare annotation fasta files
```
genome_dir=/projects/hillerlab/genome/gbdb-HL/
anno_dir=/projects/project-eosipova/NectarivoryProject/Genome_annotation_2021/


for db in HLphyNov1 HLtriMol2 HLnymHol2 HLcalAnn5 HLapuApu1 HLtaeGut4;
do
 ANNOBED=geneNames.$db.final_anno.bed
 ANNOFA=geneNames.$db.final_anno.fa
 echo -e "bedtools getfasta -split -name -fi $genome_dir/$db/$db.fa -bed $anno_dir/$db/$ANNOBED -fo $anno_dir/$db/$ANNOFA";
done > jobs.getfasta
```


## Prepare kallisto index
```
anno_dir=/projects/project-eosipova/NectarivoryProject/Genome_annotation_2021/

for db in HLphyNov1 HLtriMol2 HLnymHol2 HLcalAnn5 HLtaeGut4 HLapuApu1;
do
 #ANNOFA=noFused.geneNames.$db.final_anno.fa
 ANNOFA=noFused.validated.$db.final_anno.fa;
 anno=$anno_dir/$db/$ANNOFA;
 index=noFused.kallisto_transcripts_${db}.idx;
 echo -e "kallisto index -i $index $anno";
done > jobs.kallisto.index
```


## Run Kallisto quant for batch1 
```
md Kallisto_quant_ALL/

for l in $(cut -f1 batch1_lib_info.files.tsv | grep -v LibraryID);
do
 db=$(grep $l batch1_lib_info.files.tsv | cut -f3);
 r1=fastq_files_batch1/$l.trim_R1.fastq.gz;
 r2=fastq_files_batch1/$l.trim_R2.fastq.gz;
 index=noFused.kallisto_transcripts_${db}.idx;
 echo -e "kallisto quant -i $index -o Kallisto_quant_ALL/$l -b 100 $r1 $r2";
done > jobs.kallisto.qunat.1

para make -memoryMb 200000 jobs.kallisto.qunat.1 jobs.kallisto.qunat.1
```

## Kallisto quant for batches 1 2 3 4
```
for i in {1..4};
do
  for l in $(cut -f1 batch${i}_lib_info.files.tsv | grep -v LibraryID);
  do
    db=$(grep $l batch${i}_lib_info.files.tsv | cut -f3);
    r1=$(grep $l batch${i}_lib_info.files.tsv | cut -f5);
    r2=$(grep $l batch${i}_lib_info.files.tsv | cut -f6);
    index=kallisto_transcripts_${db}.idx;
    echo -e "kallisto quant -i $index -o Kallisto_quant_ALL/$l -b 100 $r1 $r2";
  done | sed 's/,/ /g' > jobs.kallisto.qunat.$i;
done

para make -memoryMb 200000 jobs.kallisto.qunat jobs.kallisto.qunat.$i
```

the next step is: 
## [Differential gene expression analysis with DESeq2](https://github.com/osipovarev/Transcriptome_for_NectarGenomics/blob/main/README_deseq2.md)
