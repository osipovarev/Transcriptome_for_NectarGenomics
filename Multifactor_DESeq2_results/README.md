# Multifactorial DESeq2 analysis

### Set up
```
RSCRIPTS=/Users/osipova/Documents/Scripts/rnaseq_tools/

info=ALL_info.tsv
QUANTDIR=Kallisto_quant_ALL/
QUANTFILE=validated.abundance.tsv
QUANTPATH=/projects/project-eosipova/Bird_transcriptomics/

### which tissues will be considered as factors
tissues=liver,pectoralis


### set 1
sp1=Annas_hummingbird
sp2=common_swift

### set 2
sp1=New_Holland_honeyeater
sp2=zebra_finch

### set 3
sp1=rainbow_lorikeet
sp2=cockatiel
```

### run multifactorial DESeq2 analysis
```
Rscript $RSCRIPTS/deseq2_expression_analysis.R \
  -multi \
  -w $(pwd) \
  -d $QUANTDIR \
  -c same_len.gene_${QUANTFILE} \
  -m $info \
  -t $tissues \
  -s1 $sp1 \
  -s2 $sp2 \
  -od multi_deseq2_res.$sp1.tsv  \
  -on norm_counts.multi_deseq2_res.$sp1.tsv;
```


### run GO enrichment analysis
```
hg38_dict=~/Documents/LabDocs/NectarivoryProject/absrel/absrel_analysis_2024/galGal6_gene.hg38_gene_symbol.tsv
sp=

for d in up down; 
do
	echo $d;
	genes=$d.$sp.multi_deseq2_res.top_0.05.lst;
	renameToHLscaffolds.py -c 1 -a $genes -d  <(sed 's/\t/,/' $hg38_dict) > hg38.$genes; 
	goenrich_genelist.R -w $(pwd) -g hg38.$genes -u ../background_genes.all_tissues.txt -o goenrich.hg38.${genes%lst}tsv;
done
```
