#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  Your.hmp  Phenotype.txt  genome.gff  Your_gene_ID
e.g.
     perl  CandiHap.pl  ./haplotypes.hmp   ./Phenotype.txt  ./test.gff  Si9g49990
     perl  CandiHap.pl  /data/lxk/GWAS/2.call.variation/4.annotation/haplotypes.hmp  /data/lxk/GWAS/mGWAS/TG/TG-mt1201.mGWAS.txt /data/sidb/Si.Ann_RNA.gff  Si9g49940
USAGE
print $usage if(@ARGV==0);
exit if(@ARGV==0);

$ARGV[1] =~ /.*\/(\S+)\.txt/;
$Phenotype = $1;
print "$Phenotype\n";

$mychr = $mystart = $myend =0;
system("grep  '$ARGV[3]' $ARGV[2] > gene_gff-$ARGV[3]-$Phenotype.txt");
open GFF ,"gene_gff-$ARGV[3]-$Phenotype.txt" or die "$!";
while (<GFF>) {
    $_ =~ s/\n*|\r*//g;
    if($_ =~ /\texon\t/){
        @F = split(/\t/,$_);
        $F[8]=~/Parent=(.*?);/;
        $cds_pos{$1}{"$F[3] $F[4] $F[7]"}=1;
        $contig{$1}=$F[0];
        $stream{$1}=$F[6];
    }
    if ($_ =~ /\tgene\t/){
        $_ =~ /;Note=(\d+) /;
        $gene_n = $1;
        @gff = split(/\t/,$_);
        $mychr = $gff[0];
        if ($gff[6] eq '+') {
            $mystart = $gff[3] - 2000;
            $myend   = $gff[4] + 500;
        }
        elsif ($gff[6] eq '-') {
            $mystart = $gff[3] - 500;
            $myend   = $gff[4] + 2000;
        }
    }
}
close GFF;
$gene_n = `grep 'mRNA' gene_gff-$ARGV[3]-$Phenotype.txt |wc -l`;
$gene_n =~ s/\s+//g;

print "Query:  $ARGV[3] $gene_n --> $mychr  $mystart  $myend\n";

open PHENOTYPE ,"$ARGV[1]" or die "$!";
while(<PHENOTYPE>){
    $_ =~ s/\n*|\r*//g;
    @F = split(/\t/,$_);
    $name{$F[0]} = $F[1] if $_ ne 'NaN';
}
close PHENOTYPE;


#open VCF, "$ARGV[0]" or die $!;
# #open(VCF,"gzip -dc $ARGV[0]|") or die ("can not open $ARGV[0]\n");
#while (<VCF>) {
#    $_ =~ s/\r*\n*//g;
#    @names = split(/\t/,$_) if $_ =~ /CHROM/;
#    $n  = "$names[1]\t$names[3]\t$names[4]\t$names[7]\t" . join ("\t",@names[9..$#names]);
#    next if $_ !~ /^$mychr/;
#    @a = split(/\t/,$_);
#    next if ($a[1] < $mystart or $a[1] > $myend);
#    #next if ($a[3] =~ /../ or $a[4] =~ /../);
#    $BP{$a[1]} = $a[3] . '->'. $a[4];
#    for ($i = 9; $i <= $#a; $i++){
#        $seq_alt{$names[$i]} .= "$a[1] $a[4] " if ($a[$i] =~ /1\/1/);
#        if ($a[$i] =~ /\.\/\./) {
#            $a[$i] = 'N/N' ;
#        }
#        elsif ($a[$i] =~ /(.)(.)(.)/ ) {
#            $Fir = $1;
#            $Sen = $3;
#            @alt = split(/,/,$a[4]);
#            unshift( @alt, $a[3]); 
#            $a[$i] = $alt[$Fir].'/'.$alt[$Sen];
#        }
#    }
#    #$a[7] =~ /EFF=(.*)/;
#    $l = "$a[1]\t$a[3]\t$a[4]\t$a[7]\t" . join ("\t",@a[9..$#a]);
#    push @lxk , $l;
#}
#system("grep -E '$ARGV[3]|CHROM' $ARGV[0]$mychr > hmp-$ARGV[3]-$Phenotype.txt");
system("grep -E '$ARGV[3]|CHROM' $ARGV[0] > hmp-$ARGV[3]-$Phenotype.txt");

open HMP, "hmp-$ARGV[3]-$Phenotype.txt" or die $!;
while (<HMP>) {
    $_ =~ s/\r*\n*//g;
    if ($_ =~ /CHROM/) {
        @names = split(/\t/,$_);
        $n  = "$names[1]\t$names[2]\t$names[3]\t$names[4]\t$names[5]\t" . join ("\t",@names[6..$#names]);
    }
    next if $_ !~ /^$mychr\t/;
    @a = split(/\t/,$_);
    next if ($a[1] < $mystart or $a[1] > $myend);
    $BP{$a[1]} = $a[2] . '->'. $a[3];
    #$a[7] =~ /EFF=(.*)/;
    $l = "$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t" . join ("\t",@a[6..$#a]);
    push @lxk , $l;
}
push @lxk , $n;


open OUT,">tmp-$ARGV[3]-$Phenotype.txt" or die $!;
for (@lxk) {
    my @a = split /\t/, $_;
    for (my $i=0;$i<@a;$i++){
        $out{$i}{$a[0]} = $a[$i];
    }
}
foreach my $k (sort {$a<=>$b} keys %out ){
    my $out;
    foreach $d (sort {$a<=>$b} keys %{$out{$k}}){
        $out .= "$out{$k}{$d}\t";
    }
    $out =~ s/\t$//g;
    print OUT "$out\n";
}
close IN;
close OUT;

@a ='';
open HAP,">Hap_result-$ARGV[3]-$Phenotype.txt" or die $!;
open BOX,">Boxplot-$Phenotype-$ARGV[3].txt" or die $!;
open FILE, "tmp-$ARGV[3]-$Phenotype.txt" or die $!;
while (<FILE>) {
    if ($_ =~ /^POS/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $POS = "SNP positions\t". join ("\t",@a[1..$#a]) . "\tNumber of varieties\tVarieties ID\tAverage\tStdev";
        next;
    }
    if ($_ =~ /^REF/  ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $ref = "References allele\t". join ("\t",@a[1..$#a]);
        next;
    }
   if ($_ =~ /^ALT/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $alt = "Alternative allele\t". join ("\t",@a[1..$#a]);
        next;
    }
    if ($_ =~ /^INFO/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $EFF = "SNP annotation\t". join ("\t",@a[1..$#a]);
        next;
    }
    if ($_ =~ /^Allele/ ) {
        $_ =~ s/\r*\n*//g;
        @a = split(/\t/,$_);
        $All = "Allele Frequency\t". join ("\t",@a[1..$#a]);
        next;
    }
    $_ =~ s/\r*\n*//g;
    @a = split(/\t/,$_);
    $n = join ("\t",@a[1..$#a]);
    $lxk{$n} .= "$a[0] : $name{$a[0]}, ";
    $ave{$n} .= "$name{$a[0]} " if $name{$a[0]} =~ /\d/;
    $lxk2{$n} ++;
}
print HAP "$ref\n$alt\n$All\n$EFF\n$POS\n";
$i =1 ;
foreach my $k (sort {$lxk2{$b} <=> $lxk2{$a} } keys %lxk2 ){
    $lxk{$k} =~ s/, $//g;
    $hap = Hap_ . $i;
    @ave2 = split(/ /,$ave{$k});
    if (@ave2 > 1) {
        $ave{$k} = &average(\@ave2);
        $std{$k} = &stdev(\@ave2);
    }
    if (@ave2 > 2) {
        for (@ave2) {
            print BOX "$hap\t$_\n";
        }
        #print BOX "$hap\t@ave2\n";
    }
    print HAP "$hap\t$k\t$lxk2{$k}\t$lxk{$k}\t$ave{$k}\t$std{$k}\n";
    $i ++;
}

##$gene_n
for ($i =1; $i <= $gene_n; $i++) {
  $new_gff = $ARGV[3] . '.' . $i;
  #system("grep  'gene'     gene_gff-$ARGV[3]-$Phenotype.txt >  gene_gff-$new_gff-$Phenotype.txt");
  system("grep  '$new_gff' gene_gff-$ARGV[3]-$Phenotype.txt > gene_gff-$new_gff-$Phenotype.txt");
}

@hash_size = keys %BP;
$hash_size = @hash_size;


sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}


open OUTR,">Plot_gene-$ARGV[3]-$Phenotype.R" or die $!;
print OUTR <<EOF;
mutation_plot<-function(start, stop, text="", drop=-0.15, col="red") {
  rect(start, -0.08, stop, drop, col=col, border=col)
  text(stop+70, drop-0.15, text, cex=0.7, col=col, pos=4, offset=-1)
}
genemodel_plot<-function(model, xaxis=TRUE, drop=0) {
  par(mar=c(1,1,3,1), cex=1)
  colnames(model) =c('chr','db','feature','start','end','tmp1','orientation','tmp3','gene')
  orientation = model\$orientation[1]
  start <- min(c(model\$start,model\$end))
  end <- max(c(model\$start,model\$end))
  tmp_min=min(c(model\$start,model\$end))
  model\$start=model\$start-tmp_min+1
  model\$end=model\$end-tmp_min+1
  tmp_max=max(c(model\$start,model\$end))
  tmp_min=min(c(model\$start,model\$end))
  model<-cbind(as.character(model[,3]), as.numeric(model[,4]), as.numeric(model[,5]) )
  model<-as.data.frame(model)
  colnames(model)<-c("feature", "start", "bpstop")
  model\$start<-as.numeric(as.character(model\$start));model\$bpstop<-as.numeric(as.character(model\$bpstop))
  length<-tmp_max-tmp_min
  if (orientation=="-") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-2010, end+2500), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="CDS") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="mRNA") {
        arrows(x0=model\$newstop[i]+2000,y0=0.1 - drop,x1=model\$newstart[i]-500,y1=0.1 - drop, lwd=1, col="dodgerblue4", length = 0.1)
        #arrows(x0=model\$newstop[i]+.2*length,y0=.1,x1=model\$newstart[i]-.2*length,y1=.1, lwd=1, col="dodgerblue4", length = 0.1)
      }
      if (type=="three_prime_UTR") {
        x<-c(model\$newstop[i], model\$newstop[i], model\$newstart[i]+.02*length, start, model\$newstart[i]+.02*length)
        y<-c(0 - drop,.2 - drop,.2 - drop,.1 - drop,0 - drop)
        polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (orientation=="+") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-2500, end+2010), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="mRNA") {
        arrows(x0=model\$newstart[i]-2000,y0=0.1 - drop,x1=model\$newstop[i]+500,y1=0.1 - drop, lwd=1, col="dodgerblue4", length = 0.1)
        #segments(x0=model\$newstart[i],y0=.1,x1=model\$newstop[i],y1=.1, lwd=1, col="dodgerblue4")
      }
      if (type=="CDS") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="three_prime_UTR") {
      	x<-c(model\$newstart[i], model\$newstart[i], model\$newstop[i]-.02*length, end, model\$newstop[i]-.02*length)
      	y<-c(0 - drop,.2 - drop,.2 - drop,.1 - drop,0 - drop)
      	polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  #par(new = TRUE)
  if (xaxis==T)   {
    Axis(side=3, labels=T, cex.axis=0.7)
  }
}
pdf(file="Plot_gene-$ARGV[3]-$Phenotype.pdf") 
EOF

for ($i = 1; $i <= $gene_n; $i++) {
  $new_gff = $ARGV[3] . '.' . $i;
  #$j = ($i)*0.5;
  print OUTR "$new_gff <- read.table('gene_gff-$new_gff-$Phenotype.txt',stringsAsFactors = F, header = F,comment.char = \"#\",sep = '\\t')\n";
  print OUTR "genemodel_plot(model=$new_gff, xaxis=T, drop=0)\n";
  print OUTR "text($mystart + 1000, 0.5, \"$new_gff\", cex=0.7, col=\"black\", pos=4, offset=-1)\n";
  $w =1;$old=0;
  foreach $k (sort {$a<=>$b} keys %BP ){
    $w++ if ($k < $old + 0.1* ($myend - $mystart));
    $w=1 if ($k > $old + 0.15* ($myend - $mystart));
    $d= $w * (-0.2);
    print OUTR "mutation_plot($k, $k, \"$BP{$k}\", $d, col=\"black\")\n";
    $old = $k;
  }
}

print OUTR "dev.off()\n";


open OUTR,">Plot_hist-$Phenotype.R" or die $!;
print OUTR <<EOF;
x = read.table("$ARGV[1]")
y = x\$V2
Pheno = y[2:length(y)]
Phenoty =  as.numeric(as.character(Pheno))
Phenotype = Phenoty[!is.nan(Phenoty)]
pdf(file="Plot_hist-$Phenotype.pdf", width=10, height=5)
par(mfrow=c(1,2))
h = hist(Phenotype, col='gray' ,xlab='Phenotype' , main="Histogram of $Phenotype")
rug(jitter(Phenotype))
xfit<-seq (min(Phenotype) ,max (Phenotype) )
yfit<-dnorm (xfit, mean=mean(Phenotype) , sd=sd(Phenotype))
yfit<-yfit*diff (h\$mids[1:2]) * length(Phenotype)
lines (xfit,yfit, col='blue' , lwd=1)

hist(Phenotype, freq =FALSE, col='gray' ,xlab='Phenotype' , main="Histogram of $Phenotype")
rug(jitter(Phenotype))
lines(density(Phenotype),col='red',lwd=1)
curve(dnorm(x, mean=mean(Phenotype), sd=sd(Phenotype)), add=TRUE, col="darkblue", lwd=1)
dev.off()


library("ggplot2")
df<-read.delim("Boxplot-$Phenotype-$ARGV[3].txt",sep="\t",header=F)
df\$V1<-factor(df\$V1,levels=unique(df\$V1))
pdf(file="Plot_boxplot-$Phenotype-$ARGV[3].pdf", width=10, height=5)
ggplot(data = df,aes(x = df\$V1, y = df\$V2, group= df\$V1)) +  geom_boxplot(fill=rainbow(length(levels(factor(df\$V1)))), position=position_dodge(1)) +  geom_point(size=0.7,position = position_jitter(0.1)) +labs(title="Boxplot of $ARGV[3] Gene's Haplotypes", x="Haplotypes", y="Value of Phenotype") #+ scale_x_continuous(breaks=seq(0, length(levels(factor(df\$V1))), 1)) 

boxplot(V2~V1,data=df, col=rainbow(length(levels(factor(df\$V1)))), xlab = "Haplotypes", ylab = "Value of Phenotype", main = "Boxplot of $ARGV[3] Gene's Haplotypes")
dev.off()
EOF

system("Rscript Plot_hist-$Phenotype.R");
system("Rscript Plot_gene-$ARGV[3]-$Phenotype.R");
system("rm -rf tmp-$ARGV[3]-$Phenotype.txt");
system("rm -rf gene_gff-$ARGV[3]*");
system("rm -rf *.R");
system("rm -rf hmp-$ARGV[3]-*");
system("rm -rf Boxplot*txt");
