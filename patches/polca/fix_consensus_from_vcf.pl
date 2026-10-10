#!/usr/bin/env perl
# From MaSuRCA 4.1.4 (https://github.com/alekseyzimin/masurca), GNU General Public License v3.0;
# bundled unmodified with AAFTF for aaftf_polca.sh.
#this code fixes errors in the consensus called in vcf file by freebayes
#first we read the sequences into memory
#
my $ref_contigs=$ARGV[0];
my $ctg="",$seq="";
my %rseq =();
open(FILE,$ref_contigs);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $rseq{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $rseq{$ctg}=$seq;
}

#we assume the vcf file is sorted by contig+coordinate; read in all vcf entries for the first contig
$ctg="";
my @fixes=();
my @originals=();
my $offsets=();
while($line=<STDIN>){
next if($line =~ /^\#/);
my @f=split(/\s+/,$line);
next if($f[4] =~ /,/ || not(defined($rseq{$f[0]})));
if(not($f[0] eq $ctg)){
  if($#fixes>=0){
    die("sequence $ctg not found in the input fasta file") if(not(defined($rseq{$ctg})));
    my $oldseq=$rseq{$ctg};
  #proceed with fixing
    for(my $i=$#fixes;$i>-1;$i--){#going in reverse order to avoid shifting sequence due to indels
      #first we check if the sequence at given offset matches the original variant
      my $original_seq=substr($oldseq,$offsets[$i]-1,length($originals[$i]));
      #print "$i $ctg $offsets[$i] $originals[$i] $original_seq $fixes[$i]\n";
      if(($original_seq =~ /a|c|g|t|n|A|C|G|T|N/) && not(uc($originals[$i]) eq uc($original_seq))){
        print STDERR "WARNING! sequence does not match the original $ctg $original_seq $originals[$i] $offsets[$i]\n";
      }else{
        #then substitute
        substr($oldseq,$offsets[$i]-1,length($originals[$i]),$fixes[$i]);
      }
    }
    $rseq{$ctg}=$oldseq;
  }
  @fixes=();
  @originals=();
  @offsets=();
  $ctg=$f[0];
}
my @ff=split(/:/,$f[9]);
#print "$f[9]\n";
if($ff[5]>1 && $ff[5]>=2*$ff[3]){
  push(@fixes,$f[4]);
  push(@originals,$f[3]);
  push(@offsets,$f[1]);
}
}

if($#fixes>=0){
#proceed with fixing
  my $oldseq=$rseq{$ctg};
  for(my $i=$#fixes;$i>-1;$i--){#going in reverse order to avoid shifting sequence due to indels
#first we check if the sequence at given offset matches the original variant
    die("sequence $ctg not found in the input fasta file") if(not(defined($rseq{$ctg})));
    my $original_seq=substr($oldseq,$offsets[$i]-1,length($originals[$i]));
    #print "$offsets[$i] $originals[$i] $original_seq $fixes[$i]\n";
    if(($original_seq =~ /a|c|g|t|n|A|C|G|T|N/) && not(uc($originals[$i]) eq uc($original_seq))){
      print STDERR "WARNING! sequence does not match the original $ctg $original_seq $originals[$i] $offsets[$i]\n";
    }else{
#then substitute
      substr($oldseq,$offsets[$i]-1,length($originals[$i]),$fixes[$i]);
    }
  }
  $rseq{$ctg}=$oldseq;
}

foreach my $c(keys %rseq){
  print ">$c\n$rseq{$c}\n";
}
