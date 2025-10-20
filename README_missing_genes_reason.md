## GLUTs

### SLC2A7
not labled in galGal6;  <br>
the real one is: reg_10682	rna-XM_004945643.3, rna-XM_004945642.3,	rna-XM_025155583.1  <br>
labeled as SLC2A11L5 in validated dataset

### SLC2A8
not labled in galGal6;  <br>
the real one is: reg_5963	ENSGALT00000001251.5, rna-XM_015279297.2

### SLC2A14
not labled in galGal6;  <br>
the real one is: reg_2730	ENSGALT00000037225.4  

SLC2A8 missing HLcalAnn5 HLtaeGut4 HLphyNov1 - TOGA says one2zero in them.  <br>
SLC2A11 missing in HLtriMol2 HLnymHol2 - TOGA says one2zero in them.  <br>
SLC2A11L missing in HLnymHol2 - one2many  <br>
SLC2A12 missing in HLtriMol2 - one2many  <br>
SLC5A9 missing in HLphyNov1 - one2many  <br>


## glycolysis/gluconeogenesis/glycogen synthesis

### GCK 
missing HLcalAnn5 HLtriMol2 - one2zero;  <br>
annotated in HLtaeGut4 as reg_15241 - fixed now. 

### HK2
missing HLapuApu1 HLphyNov1  <br>
 in HLapuApu1:  <br>
 HK2_ENSGALT00000108207.1_mRNA14998, <br>
 HK2_ENSGALT00000108207.1_rna-XM_030290093.1.98, <br>
 HK2_ENSGALT00000108207.1_rna-XM_030290094.1.98<br>


 in HLphyNov1:  <br>
 NBCI annotated a very truncated iso (<50% TOGA iso) <br> => added it manually to one2ones: HK2_rna-XM_030290094.1.513, HK2_mRNA15514

HK3 missing HLtriMol2 - one2zero  <br>
GPI missing HLtriMol2 - one2many  <br>
PFKL missing HLtriMol2 - one2zero  <br>
PDHA1 missing HLtriMol2 - one2many  <br>
PDHB missing HLtaeGut4 - NBCI annotated a very truncated iso (<50% TOGA iso) => added it manually to one2ones: PDHB_ENSGALT00000048041.2_rna-NM_001245636.1  <br>
GAPDH missing HLcalAnn5 - one2zero :(  <br>
PFKM missing : HLcalAnn5, HLapuApu1, HLphyNov1 - one2zero; HLtaeGut4 - ona2many;  <br>


MLXIPL missing HLtaeGut4 - one2zero  <br>


## Fix mislabeled

### SLC2A7
```
for db in HLcalAnn5 HLapuApu1 HLphyNov1 HLtriMol2 HLnymHol2 HLtaeGut4; \
do 
	echo $db; \
	sed 's/reg_10682$/SLC2A7/g' $db.one2ones.toga.my_anno.dict_GENES.tsv > file.$db; \
	mv file.$db $db.one2ones.toga.my_anno.dict_GENES.tsv; \
done
```

### SLC2A8
```
for db in HLcalAnn5 HLapuApu1 HLphyNov1 HLtriMol2 HLnymHol2 HLtaeGut4; \
do \
	echo $db; \
	sed 's/reg_5963$/SLC2A8/g' $db.one2ones.toga.my_anno.dict_GENES.tsv > file.$db; \
	mv file.$db $db.one2ones.toga.my_anno.dict_GENES.tsv; \
done
```

### SLC2A14
```
for db in HLcalAnn5 HLapuApu1 HLphyNov1 HLtriMol2 HLnymHol2 HLtaeGut4; \
	do echo $db; \
	sed 's/reg_2730$/SLC2A14/g' $db.one2ones.toga.my_anno.dict_GENES.tsv > file.$db; \
	mv file.$db $db.one2ones.toga.my_anno.dict_GENES.tsv; \
done
```
