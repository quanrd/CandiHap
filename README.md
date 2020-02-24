# CandiHaplotypes
Genome-wide association study (GWAS) is widely utilized to identify genes involved in plants, animals and human complex traits. The identified SNP is not necessarily the causal variant, but it is rather in linkage disequilibrium (LD). One key challenge for GWAS data interpretation is to fast identify causal genes and provide profound evidence on how they affect the trait. Researches want to identify candidate causal variants from the most significant SNPs of GWAS in any species and on their local computer, while to complete these tasks are to be time-consuming, laborious and prone to errors and omission. To our knowledge, so far there is no tool available to solve the challenge for GWAS data very quickly. CandiHaplotypes is developed to identify candidate causal SNPs and genes from GWAS by integrating LD analysis, functional SNP annotation, haplotype analyses and traits statistics of haplotypes. CandiHaplotypes provides a very fast preselection form GWAS result to candidate causal variants

# CandiHaplotypes: a tool for fast identifying candidate causal genes from genome-wide association study
![+](CandiHaplotypes.png)

## Description
__`SNP_Gene_Phenotype2haplotypes.pl`__ can be used to draw the Gene structure and haplotypes quickly in `PDF` format based on the VCF, GFF and gene's ID. It will become a useful tool for `drawing gene structure` with the advantages of a simple input data format, easily modified output and very good portability. __`SNP_Gene_Phenotype2haplotypes.pl`__ is open source, so all we ask is that you cite our most recent paper in any publications that use this script:</br>
> 基于重测序vcf数据计算群体的单倍型频率及作图.</br>
> A tool to draw Gene structure and haplotypes based on perl and R language.</br>

## Dependencies
__`perl 5`__ and all __`R`__ versions. </br></br>

## Getting started
Put __`vcf2hmp.pl`__  test.gff, test.vcf, and genome.fa files in a same dir, then run annovar (table_annovar.pl):</br>
```
     gffread  test.gff   -T -o test.gtf
     gtfToGenePred -genePredExt test.gtf  si_refGene.txt
     retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa
     
     table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test --buildver  si --protocol refGene --operation g -remove
     
     perl  vcf2hmp.pl  test.vcf  test.si_multianno.txt
```
</br>

Put __`SNP_Gene_Phenotype2haplotypes.pl`__ and Phenotype.txt files in a same dir, then run:</br>
```
     perl  SNP_Gene_Phenotype2haplotypes.pl  ./Your.hmp  ./Phenotype.txt  ./genome.gff  Your_gene_ID
e.g. perl  SNP_Gene_Phenotype2haplotypes.pl  ./haplotypes.hmp   ./Phenotype.txt  ./test.gff  Si9g49990
```
</br>

By the way, if you want do All gene in LD region of a position, then run:</br>
```
     perl  GWAS2LD_haplotypes.pl   ./genome.gff  ./ann.hmp  ./Phenotype.txt  50kb  Chr:position
e.g. perl  perl  GWAS2LD_haplotypes.pl   ./your.gff  ./haplotypes.hmp  ./Phenotype.txtt  50kb  9:54605172
```
</br>

## Contact information
For any questions please contact xukai_li@sxau.edu.cn or xukai_li@qq.com</br>
