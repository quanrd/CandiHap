#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     gffread  test.gff   -T -o test.gtf
     gtfToGenePred -genePredExt test.gtf  si_refGene.txt
     retrieve_seq_from_fasta.pl --format refGene --seqfile  genome.fa  si_refGene.txt --outfile si_refGeneMrna.fa
     
     The txt result of annovar (table_annovar.pl) is test.si_multianno.txt
     table_annovar.pl  test.vcf  ./  --vcfinput --outfile  test --buildver  si --protocol refGene --operation g -remove

     perl  $0  Your.vcf Your_multianno.txt   maf
e.g.
     perl  vcf2hmp.pl  test.vcf  test.si_multianno.txt  0.1
USAGE

print $usage if(@ARGV==0);
exit if(@ARGV==0);


open OUT,">haplotypes.hmp.txt" or die $!;
open VCF, "$ARGV[0]" or die $!;
while (<VCF>) {
    next if $_ =~ /^##/;
    $_ =~ s/\r*\n*//g;
    if ($_ =~ /CHROM/) {
        @names = split(/\t/,$_);
        print OUT  "$names[0]\t$names[1]\t$names[3]\t$names[4]\t$names[7]\tAllele Frequency\t" . join ("\t",@names[9..$#names]) . "\n";
    }
    last;
}
close VCF;

open ANNOVAR,$ARGV[1] or die $!;
while(<ANNOVAR>){
    #next if $_ =~ /^Chr/;
    $_ =~ s/\r*\n*//g;
    @a = split(/\t/,$_);
    #next if $a[5] =~ /intergenic/;
    $chrpos = "$a[0]:$a[1]";
    $ann = "$a[5] :: $a[6] :: $a[7] :: $a[8] :: $a[9]";
    $ann =~ s/ :: \.//g;
    $h{$chrpos} = $ann;
    for ($i = 22; $i <= $#a; $i++){
        if ($a[$i] =~ /\.\/\./) {
            $a[$i] = 'N/N' ;
        }
        elsif ($a[$i] =~ /(.)(.)(.)/ ) {
            $Fir = $1;
            $Sen = $3;
            @alt = split(/,/,$a[4]);
            unshift( @alt, $a[3]); 
            $a[$i] = $alt[$Fir].'/'.$alt[$Sen];
        }
    }
    $a[20] =~ /;AF=(\d+\.\d+);/;
    $AF = $1;
    next if $AF <= $ARGV[2];
    $l = "$a[0]\t$a[1]\t$a[3]\t$a[4]\t$a[5]\t$AF\t" . join ("\t",@a[22..$#a]);
    print OUT "$l\n";
}
close ANNOVAR;


open OUT,">haplotypes.hmp" or die $!;
@a ='';
open HMP, "haplotypes.hmp.txt" or die $!;
while(<HMP>) {
    print OUT if $_ =~ /#CHROM/;
    $_ =~ s/\n*|\r*//g;
    @a=split /\t/, $_;
    $lxk = $a[0] . ':' . $a[1];
    if( exists $h{$lxk}) {
        $a[4] = $h{$lxk};
        print OUT join ("\t",@a),"\n";
    }
}

system("rm -rf haplotypes.hmp.txt");
