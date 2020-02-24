There are mainly three steps included in the CandiHaplotypes analytical through command lines, and the test data files can freely download at https://github.com/xukaili/SNP_Gene_Phenotype2haplotypes.pl.
1.	To annotate the vcf by ANNOVAR: 
1.1 gffread  test.gff   -T -o test.gtf
1.2 gtfToGenePred -genePredExt test.gtf  si_refGene.txt
1.3 retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa
1.4 table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test --buildver  si --protocol refGene --operation g -remove

2.	To convert the txt result of annovar to hapmap format (0.1 means the minor allele frequency (MAF)): 
perl  vcf2hmp.pl  test.vcf  test.si_multianno.txt  0.1

3.	To run CandiHaplotypes: 
perl  GWAS_LD2haplotypes.pl   ./test.gff  ./haplotypes.hmp   ./Phenotype.txt   50kb  9:54583294
Or to run CandiHaplotypes by one gene: 
perl  CandiHap.pl ./haplotypes.hmp   ./Phenotype.txt  ./test.gff  Si9g49990















sed -n '4119763, 4120870p' all.filt.recode.vcf  > test1.vcf
head -45 all.filt.recode.vcf > tt.head.txt
cat tt.head.txt test1.vcf > test.vcf
rm test1.vcf tt.head.txt

gffread  test.gff   -T -o test.gtf
gtfToGenePred -genePredExt test.gtf  si_refGene.txt
retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa

table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test --buildver  si --protocol refGene --operation g -remove
perl  vcf2hmp.pl  test.vcf  test.si_multianno.txt  0.1
rm test.gtf test.avinput test.si_multianno.vcf
 
perl  GWAS_LD2haplotypes.pl   ./test.gff  ./haplotypes.hmp   ./Phenotype.txt   50kb  9:54583294
perl  CandiHap.pl ./haplotypes.hmp   ./Phenotype.txt  ./test.gff  Si9g49990

