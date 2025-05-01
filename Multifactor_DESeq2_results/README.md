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

### Run multifactorial DESeq2 analysis
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

### Fix headers
```
for f in $(ls *multi_deseq2_res.*.tsv);
do
 echo $f;
 sed 's/^baseMean/gene\tbaseMean/' $f > file;
 mv file $f;
done
```

### Get genes of interest:
1) select genes that are differentially expressed between tissues (liver-pectoralis) in CONTROL (ref)
2) out of those, select genes that have a different 'between-tissue-expression' pattern in NECTAR (interaction)
-> UP = genes that have significant & strong tissue effect in CONTROL, and significantly higher tissue-specific expression in NECTAR vs CONTROL
-> DOWN = genes that have significant & strong tissue effect in CONTROL, and significantly lower tissue-specific expression in NECTAR vs CONTROL
=> genes with tissues-specific expression that show a different pattern in nectar vs control (evolutionary divergent pattern)
```
for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	ref=ref.multi_deseq2_res.$sp.tsv; \
	awk '($7<0.01 && ($3<=-2 || $3>=2)){print $1}' $ref > sign.${ref%tsv}lst; \
done

for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	interaction=interaction.multi_deseq2_res.$sp.tsv; \
	genes=$sp.multi_deseq2_res.lst; \
	filter_annotation_with_list.py -b -c 1 -l sign.ref.multi_deseq2_res.in_any.lst -a <(awk '($7<0.05 && $3>1 ){print $1}' $interaction | grep -v baseMean ) > up.$genes; \
	filter_annotation_with_list.py -b -c 1 -l sign.ref.multi_deseq2_res.in_any.lst -a <(awk '($7<0.05 && $3<-1 ){print $1}' $interaction | grep -v baseMean ) > down.$genes; \
done
```


## Run GO enrichment analysis with ClusterProfiler
```
hg38_dict=~/Documents/LabDocs/NectarivoryProject/absrel/absrel_analysis_2024/galGal6_gene.hg38_gene_symbol.tsv

for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	for d in up down; \
	do \
		echo $d; \
		genes=$sp.multi_deseq2_res.lst; \
		renameToHLscaffolds.py -c 1 -a $d.$genes -d <(sed 's/\t/,/' $hg38_dict) > hg38.$d.$genes; \ 
		goenrich_genelist.R -w $(pwd) -g hg38.$d.$genes -u ../background_genes.all_tissues.txt -o goenrich.hg38.$d.${genes%lst}tsv; \
	done; \
done
```


### get representative GO terms
```
for sp in Annas_hummingbird New_Holland_honeyeater rainbow_lorikeet; \
do \
	for d in up down; \
	do \
		echo $d; \
		genes=$sp.multi_deseq2_res.lst; \
		find_representative_go.R -w $(pwd) -e goenrich.hg38.$d.${genes%lst}tsv -o represent_GO.goenrich.hg38.$d.${genes%lst}tsv; \
done
```


### get rank2 sets
```
for d in up down; \
do \
	for g in $(cut -f1 goenrich.hg38.down.*multi_deseq2_res.tsv | g -v ^ID | s | uniq -c | awk '$1>=2{print $2}'); \
	do \
		grep $g goenrich.hg38.$d.Annas_hummingbird.multi_deseq2_res.tsv ; \
	done | cut -f1,2 > 2_way_convergent_terms.multi_deseq2.${d}_pectoralis.tsv' \
done
```


### exclude children
```
GOOBO=/Users/osipova/Documents/LabDocs/GO_terms_genes/go.obo

for d in up down; \
do \
	f=rank2.goenrich.hg38.$d.multi_deseq2_res.tsv \
	golist=$(cut -f2 $f | tail +2 | tr '\n' ','); \

	for g in $(cut -f2 $f | tail +2); \
	do \
		get_go_children.py -f $GOOBO -go $g -l $golist ; \
	done | grep "has parents" | awk '{print $1}' > to_exclude_go.lst; \

	filter_annotation_with_list.py -b -c 2 -a $f -l to_exclude_go.lst > noChildren.$f; \
done
```

