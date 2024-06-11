# calcGDV

Tool for **calc**ulating **g**enotypic **d**iagnostic **v**alues for discrete traits or classes from phenotypic and variant calling data to aid in marker development

## Description

Get predictive value (PV), true positive rate (TPR), false positive rate (FPR) and positive likelihood ratio (LR) of genotypes (GTs) as potential markers for discrete traits. Takes [flapjack](https://github.com/cropgeeks/flapjack) formatted genotype, phenotype and map files.
+ For a given marker site, the PV (a.k.a precision) of a GT for a particular trait is: _PV=n(GT|trait)/n(GT) =TP/(TP+FP)_ where a positive test represents a known GT call - i.e. the proportion of all samples with a particular GT that actually have the trait.
+ The TPR (_=TP/(TP+FN)_) can also be considered the reverse PV - i.e. the proportion of all samples having the trait actually positive for the GT (a.k.a. recall or sensitivity).
+ _FPR=FP/(FP+TN)_, which is the proportion of FPs out of all samples without the trait.
+ False discovery rate (FDR), false negative rate (FNR) and true negative rate (TNR) are reciprocals of PV, TPR and FPR (_FDR=1-PV, FNR=1-TPR, TNR=1-FPR_).
+ _LR=TPR/FPR_ and is used for assessing the value of a positive GT call in usefully changing the odds that the trait exists in a test sample.
+ An ideal marker should have both PV and TPR close to 1, with a LR of 10 or more (assuming a sufficient number and diversity of input samples).

Genotype file must first be transposed with _transpose.awk_ script. Alternate orders of alleles (i.e. A/T or T/A) are treated as distinct. Nucleotide groupings are also tested where applicable, e.g. "hasC" comprises all GTs containing "C" (C, C/T, T/C) for a given site. CalcGDV Does not adjust for population stratification and/or linkage disequilibrium (LD) which the user should factor into their marker/sample selection beforehand. Writes to STDOUT.

## Example Workflow

Here we will use a toy multisample VCF file (variants.vcf.gz) containing ~300 variant calls for 182 samples. A phenotype matrix (phenos.tsv) is also provided containing sample rows and two trait columns labelled "symptoms" and "sex".

### preprocessing (optional)

_variants.vcf.gz_ has already been quality filtered but we will run additional selection using [BCFtools](https://github.com/samtools/bcftools) to inlcude only biallelic SNPs with no more than 10% missing calls.
```
bcftools view -m2 -M2 -v snps -i 'F_MISSING < 0.1' --output-type z -o variants.qc.vcf.gz variants.vcf.gz
```
>new vcf is written to _variants.qc.vcf.gz_<br />


### convert VCF to transposed genotype matrix

We first need to convert the vcf to a flapjack formatted genotype file. In the absence of [flapjack](https://github.com/cropgeeks/flapjack) software, open source tools such as [NGSEPcore](https://github.com/NGSEP/NGSEPcore) can be used. A _map_ file is also generated containing marker position information.
```
java -jar NGSEPcore_5.0.0.jar VCFConverter -flapjack -i variants.qc.vcf.gz -o genos
```
>outputs written to _genos_fj.gen_ and _genos_fj.map_. NGSEPcore automatically generates marker names for variants<br />

_calcGDV_ requires the _.gen_ file to be transposed (markers as rows instead of columns). We can do this using the provided _transpose.awk_ script to generate _genos_fj.tr.gen_
```
./transpose.awk genos_fj.gen > genos_fj.tr.gen
```

### run calcGDV with genotype and phenotype files

Now we can run _calcGDV.py_, trying out different parameters. Use `python calcGDV.py -h` to see a full list of options.<br />
<br />

**i) run with minimum arguments** - this will output values for every genotype in _genos_fj.tr.gen_ using the "symptoms" column of _phenos.tsv_ 
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms > gdv_out.tsv
```
output for a single marker (i.e. variant site) might look something like:
```
...
#markerID=Marker143,GTavail=182,GTmissing=0
#trait=immune,GTavail=91,GTmissing=0
#trait=mild,GTavail=47,GTmissing=0
#trait=susceptible,GTavail=44,GTmissing=0
#<ID>   <GT>    <n>     <trait> <PV>    <TPR>   <FPR>   <LR>
Marker143       G       97      immune  0.701   0.747   0.319   2.345
Marker143       G       97      mild    0.237   0.489   0.548   0.893
Marker143       G       97      susceptible     0.062   0.136   0.659   0.207
Marker143       T       19      mild    0.053   0.021   0.133   0.16
Marker143       T       19      susceptible     0.947   0.409   0.007   56.455
Marker143       G/T     66      immune  0.348   0.253   0.473   0.535
Marker143       G/T     66      mild    0.348   0.489   0.319   1.536
Marker143       G/T     66      susceptible     0.303   0.455   0.333   1.364
Marker143       hasG    163     immune  0.558   1.0     0.791   1.264
Marker143       hasG    163     mild    0.282   0.979   0.867   1.129
Marker143       hasG    163     susceptible     0.16    0.591   0.993   0.595
Marker143       hasT    85      immune  0.271   0.253   0.681   0.371
Marker143       hasT    85      mild    0.282   0.511   0.452   1.13
Marker143       hasT    85      susceptible     0.447   0.864   0.341   2.536
...
```
>in this example, Marker143 comprises the genotypes 'G', 'T' and 'G/T'. Lines with 'hasG' give the predictiveness of either 'G' or 'G/T' genotypes (likewise for 'hasT'). Values for each of the three symptom traits are shown for every genotype. We can see that genotype 'T' has good precision and specificity (high PV and low FPR) for the "susceptible" trait, but has poor sensitivity (low TPR). Thus, with a LR of 56.455, a positive test for 'T' may be useful in a diagnostic context (when symptoms are already present) but not ideal for screening tests (missing true cases).<br />
<br />

**ii) run on selection of markers** - supply a comma separated list of marker IDs to `-s` option
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -s Marker142,Marker143 > gdv_out.tsv
```
>this will output only results for Marker142 and Marker143<br />

alternatively, a text file containing a single column list of marker IDs can be provided to `-s`
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -s markers.txt > gdv_out.tsv
```
<br />

**iii) filtering output with cutoffs** - output only the most significant genotypes by applying thresholds to PV, TPR, FPR and LR
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -s markers.txt -n 10 -lr 10 -tpr 0.6 > gdv_out.tsv
```
>here `-n 10` indicates the min occurences of a genotype required to be reported, `-lr 10` specifies a LR of 10 or greater, `-tpr 0.6` specifies a TPR of 0.6 or greater<br />

this now reports just a single genotype satisfying all cutoffs:
```
#<ID>   <GT>    <n>     <trait> <PV>    <TPR>   <FPR>   <LR>
Marker142       T       32      susceptible     0.938   0.682   0.014   47.045
```
<br />

**iv) disregard traits or missing pheno data** - label for missing or undesired trait can be added the `-t` string to omit samples with that label
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms,mild -s markers.txt -lr 10 -n 10 -tpr 0.6 > gdv_out.tsv
```
>samples with the "mild" symptom are not included in calculating diagnostic values<br />

using the same cutoffs as before, the output has now changed:
```
#<ID>   <GT>    <n>     <trait> <PV>    <TPR>   <FPR>   <LR>
Marker142       A       65      immune  0.985   0.703   0.023   30.945
Marker142       T       30      susceptible     1.0     0.682   0.0     nan
```
>omitting "mild" samples results in zero false positives for 'T', making the LR undefined (treated as infinite)<br />
<br />

**v) run on markers within specific genomic ranges** - as an alternative to `-s` selection, a referece sequence range can be given with option `-r` (requires map file)
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -r genos_fj.map,chr1_A:4022300-4022500 -lr 10 -n 10 -tpr 0.6 > gdv_out.tsv
```
>first argument of `-r` is the flapjack map file _genos_fj.map_ followed by single genomic range in the format _seqid:start-end_<br />

output lines contain additional seqid and position fields:
```
#<ID>   <GT>    <n>     <trait> <PV>    <TPR>   <FPR>   <LR>    <seq>   <pos>
Marker142       T       32      susceptible     0.938   0.682   0.014   47.045  chr1_A  4022348
...
Marker143       T       19      susceptible     0.947   0.409   0.007   56.455  chr1_A  4022497
```
instead of a single range, a BED file containing multiple ranges can be supplied in place of _seqid:start-end_ coordinate
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -r genos_fj.map,regions.bed > gdv_out.tsv
```
<br />

**vi) range and marker selection can be combined** - even if their respective regions don't overlap
```
egrep "Marker142|Marker143" genos_fj.map
    --> Marker142       chr1_A  4022348
        Marker143       chr1_A  4022497

python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -s Marker142 -r genos_fj.map,chr1_A:4022400-4022500 -pv 0.9 > gdv_out.tsv
```
<br />

**vii) sample selection can be applied using label from another trait column** - diagnostic values are derived from only samples of a particular group specified in `-f`
```
python calcGDV.py -g genos_fj.tr.gen -t phenos.tsv,symptoms -lr 50 -n 10 -f sex,male > gdv_out.tsv
```
>here `-f sex,male` tells calcGDV to only consider samples labelled "male" in the "sex" column<br />
<br />
