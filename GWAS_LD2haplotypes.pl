#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  ./genome.gff  ./ann.hmp  ./Phenotype.txt   LDkb  Chr:position
e.g. perl  GWAS_LD2haplotypes.pl   /data/sidb/Si.Ann_RNA.gff  /data/lxk/GWAS/2.call.variation/4.annotation/All_ann.hmp  /data/lxk/GWAS/2.call.variation/4.annotation/haplotypes/TG-m0022.mGWAS.txt  50kb  9:54605172
e.g. perl  GWAS_LD2haplotypes.pl   ./test.gff  ./haplotypes.hmp   ./Phenotype.txt  50kb  9:54583294 

USAGE
print $usage if(@ARGV==0);
exit if(@ARGV==0);

$ARGV[2] =~ /.*\/(\S+)\.txt/;
$Phenotype = $1;
print "$Phenotype\n";

$ARGV[3] =~ /(\d+)/;
$kb = $1 * 1000;
print "The LD is $ARGV[3] --> $kb bp\n";

open GFF ,"$ARGV[0]" or die "$!";
while (<GFF>) {
    if ($_ =~ /\tgene\t/) {
        $_ =~ s/\n*|\r*//g;
        @F = split(/\t/,$_);
        $F[8] =~ /ID=(S\w+);?(.*)/;
        $lxk{$1} = "$F[0]\t$F[3]\t$F[4]\t$1\t$2";
    }
}
close GFF;
@F='';

open RES,">>LD_$ARGV[3]_gene_$ARGV[4].txt" or die $!;

@b = split(/:/,$ARGV[4]);
$start = $b[1] - $kb;
$end   = $b[1] + $kb;
for (keys %lxk) {
    @F = split(/\t/,$lxk{$_});
    if (($F[0] eq $b[0]) and ($F[1] > $start) and ($F[2] < $end)) {
       #print "@b:\t($F[0] eq $b[0]) and ($F[1] > $start) and ($F[2] < $end)\n";
        print  RES "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\n";
        system("perl  CandiHap.pl  $ARGV[1]  $ARGV[2]  $ARGV[0]  $F[3]");
    }
}

open OUTR,">Plot_hist-$Phenotype.R" or die $!;
print OUTR <<EOF;
x = read.table("$ARGV[2]")
y = x\$V2
Pheno = y[2:length(y)]
Phenotype =  as.numeric(as.character(Pheno))
pdf(file="Plot_hist-$Phenotype.pdf", width=10, height=5)
par(mfrow=c(1,2))
h = hist(Phenotype, col='gray' ,xlab='Phenotype' , main="Histogram of $Phenotype")
rug(jitter(Phenotype))
xfit<-seq (min(Phenotype) ,max (Phenotype) )
yfit<-dnorm (xfit, mean=mean(Phenotype) , sd=sd(Phenotype))
yfit<-yfit*diff (h\$mids[1:2]) *length (Phenotype)
lines (xfit,yfit, col='blue' , lwd=1)

hist(Phenotype, freq =FALSE, col='gray' ,xlab='Phenotype' , main="Histogram of $Phenotype")
rug(jitter(Phenotype))
lines(density(Phenotype),col='red',lwd=1)
curve(dnorm(x, mean=mean(Phenotype), sd=sd(Phenotype)), add=TRUE, col="darkblue", lwd=1)
dev.off()
EOF
system("Rscript Plot_hist-$Phenotype.R");

system("rm -rf *.R");
