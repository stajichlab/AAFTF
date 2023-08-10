#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
export PATH=$MYPATH:$PATH;
set -o pipefail
export NUM_THREADS=4
export MEM=1G
export FIX=1
export BATCH_SIZE=5000000
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage {
echo "Usage:"
echo "polca.sh [arguments]"
echo "-a <assembly contigs or scaffolds>"
echo "-r <'polishing_reads_fastq1 polishing_reads_fastq2'> can use any number of fastq files, polishing reads must be in fastq format!"
echo "-t <number of threads, default:1>"
echo "-n <optional: do not polish, just create vcf file, evaluate the assembly and exit>"
echo "-m <optional: memory per thread to use in samtools sort, set to 2G or more for large genomes>"
echo ""
echo "External dependency: must have bwa available on the PATH"
which bwa
}

#parsing arguments
if [[ $# -eq 0 ]];then
usage
exit 1
fi

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -b|--batch)
            export BATCH_SIZE="$2"
            shift
            ;;
        -t|--threads)
            export NUM_THREADS="$2"
            shift
            ;;
        -n|--nofix)
            export FIX=0
            ;;
        -a|--assembly)
            export ASM="$2"
            shift
            ;;
        -r|--reads)
            READS="$2";
            shift
            ;;
        -m|--memory)
            export MEM="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [[ ! -e $ASM ]];then
echo "Input file $ASM not found or not specified!"
usage
exit 1
fi

which bwa || error_exit "bwa not found on the PATH, please install bwa aligner"
ls $MYPATH/freebayes || error_exit "freebayes not found in MaSuRCA bin, please check your MaSuRCA install"
ls $MYPATH/samtools || error_exit "samtools not found in MaSuRCA bin, please check your MaSuRCA install"

ASMPATH="`dirname \"$ASM\"`"
ASMPATH="`( cd \"$ASMPATH\" && pwd )`"
export BASM=`basename $ASM`
export BWA=`which bwa`
export FREEBAYES=$MYPATH/freebayes
export SAMTOOLS=$MYPATH/samtools
rm -f bwa.err samtools.err

if [ ! -e $BASM.index.success ];then
  log "Creating BWA index for $ASM"
  rm -f $BASM.map.success
  $BWA index $ASM -p $BASM.bwa 2>>bwa.err && touch $BASM.index.success || error_exit "Creating BWA index for $ASM failed"
fi

if [ ! -e $BASM.map.success ];then
  log "Aligning reads to $ASM"
  rm -f $BASM.sort.success
  zcat -f $(echo $READS) | $BWA mem -SP -t $NUM_THREADS $BASM.bwa /dev/stdin 1>$BASM.unSorted.sam 2>>bwa.err && \
  touch $BASM.map.success
  if [ ! -e $BASM.map.success ];then
    error_exit "Aligning reads to $ASM failed"
  fi
fi

if [ ! -e $BASM.sort.success ];then
  log "Sorting and indexing alignment file"
  rm -f $BASM.vc.success
  $SAMTOOLS sort -m $MEM -T $BASM -@ $NUM_THREADS $BASM.unSorted.sam -o $BASM.alignSorted.bam 2>>samtools.err && \
  $SAMTOOLS index $BASM.alignSorted.bam 2>>samtools.err && \
  $SAMTOOLS faidx $ASM  2>>samtools.err && \
  touch  $BASM.sort.success

  if [ ! -e $BASM.sort.success ];then
    error_exit "Sorting and indexing alignment file failed"
  fi
fi

#here we are doing variant calling in parallel, per input contig/scaffold
if [ ! -e $BASM.vc.success ];then
  rm -f  $BASM.report.success $BASM.fix.success
  log "Calling variants in $BASM"
  mkdir -p $BASM.vc
  ufasta sizes -H $ASM | awk 'BEGIN{bs='$BATCH_SIZE';bi=1;btot=0;}{if($2>bs){for(n=0;n<$2;n+=bs){max=n+bs;if(max>$2) max=$2;print bi" "$1" "n" "max;bi++;}btot=0}else{if(btot+$2>bs){bi++;btot=$2}else{btot+=$2}print bi" "$1" 0 "$2}}' > $BASM.batches
  export BATCHES=`awk '{print $1}' $BASM.batches|uniq|wc -l |awk '{print $1}'`
#begin subshell execution
  (
    cd $BASM.vc
    rm -f $BASM.vc.success $BASM.fix.success
    echo "#!/bin/bash" > commands.sh
    echo "if [ ! -e \$1.vc.success ];then" >> commands.sh
    echo "awk '{if(\$1=='\$1') print \$2\" \"\$3\" \"\$4}' ../$BASM.batches > batch.\$1" >> commands.sh
    echo "  $FREEBAYES -t batch.\$1 -m 0 --min-coverage 3 -R 0 -p 1 -F 0.2 -E 0 -b ../$BASM.alignSorted.bam -v \$1.vcf -f $ASMPATH/$BASM && touch \$1.vc.success" >> commands.sh
    echo 'fi' >> commands.sh
    chmod 0755 commands.sh
    seq 1 $BATCHES |xargs -P $NUM_THREADS -I % ./commands.sh %

#checking if jobs finished properly

    for f in $(seq 1 $BATCH);do
      if [ ! -e $f.vc.success ];then
        error_exit "Variant calling failed on batch $f in $BASM.vc"
      fi
    done
    touch $BASM.vc.success
  );

#end subshell execution
  if [ -e ./$BASM.vc/$BASM.vc.success ];then
    seq 1 $BATCHES | xargs -I % ls $BASM.vc/%.vcf |xargs cat  | grep '^##' |grep -v commandline | sort -S 10% |uniq > $BASM.vcf.header.tmp && mv $BASM.vcf.header.tmp $BASM.vcf.header1 && \
    seq 1 $BATCHES | xargs -I % ls $BASM.vc/%.vcf |xargs cat  | grep '^#C' |sort -S 10% |uniq > $BASM.vcf.header.tmp && mv $BASM.vcf.header.tmp $BASM.vcf.header2 && \
    seq 1 $BATCHES | xargs -I % ls $BASM.vc/%.vcf |xargs cat  | grep -v '^#' > $BASM.vcf.body.tmp && mv $BASM.vcf.body.tmp $BASM.vcf.body && \
    cat $BASM.vcf.header1 $BASM.vcf.header2 $BASM.vcf.body > $BASM.vcf.tmp && mv $BASM.vcf.tmp $BASM.vcf && \
    touch $BASM.vc.success && rm -rf $BASM.vc $BASM.vcf.body $BASM.vcf.header{1,2}
  else
    error_exit "Variant calling failed in ./$BASM.vc"
  fi
fi

if [ ! -e $BASM.fix.success ];then
  if [ $FIX -gt 0 ];then #fixing in parallel from $BASM.vcf
  mkdir -p $BASM.fix
  ufasta sizes -H $ASM |sort -S 10% -k2 |awk '{print $1}' > $BASM.names
#begin subshell
  (
    cd $BASM.fix
    CONTIGS=`wc -l ../$BASM.names | awk '{print $1}'`;
    BATCH_SIZE=$(($CONTIGS / $NUM_THREADS+1));
    if [ $BATCH_SIZE -lt 1 ];then
      $BATCH_SIZE=1
    fi
    echo "Processing $BATCH_SIZE scaffold(s) per batch"
    BATCH=1
    INDEX=1
    rm -f $BATCH.names
    for f in $(cat ../$BASM.names);do
      echo $f >> $BATCH.names
      let INDEX=$INDEX+1
      if [ $INDEX -gt $BATCH_SIZE ];then
        INDEX=1
        let BATCH=$BATCH+1
        rm -f $BATCH.names
      fi
    done

    if [ $INDEX -lt 2 ];then
      let BATCH=$BATCH-1
    fi

    echo "#!/bin/bash" > commands.sh
    echo "if [ ! -e \$1.fix.success ];then" >> commands.sh
    echo "    $MYPATH/fix_consensus_from_vcf.pl <($MYPATH/ufasta extract -f \$1.names $ASMPATH/$BASM) < ../$BASM.vcf 1>\$1.fixed.tmp 2>$1.fixed.err && mv \$1.fixed.tmp \$1.fixed && touch \$1.fix.success"  >> commands.sh
    echo "fi" >> commands.sh
    chmod 0755 commands.sh && \
    seq 1 $BATCH |xargs -P $NUM_THREADS -I % ./commands.sh %

    for f in $(seq 1 $BATCH);do
      if [ ! -e $f.fix.success ];then
        error_exit "Fixing errors failed on batch $f in $BASM.fix"
      fi
    done
    touch $BASM.fix.success
  );

#checking if jobs finished properly
  if [ -e ./$BASM.fix/$BASM.fix.success ];then
    cat ./$BASM.fix/*.fixed  | ufasta format > $BASM.masurca.tmp && mv $BASM.masurca.tmp $BASM.PolcaCorrected.fa && touch $BASM.fix.success && rm -rf $BASM.fix
  else
    error_exit "Fixing consensus failed in ./$BASM.fix"
  fi
  fi
fi


if [ ! -e $BASM.report.success ];then
  log "Creating report file"
  NUMSUB=`grep --text -v '^#'  $BASM.vcf  |perl -ane '{if(length($F[3])==1 && length($F[4])==1){ print "$F[9]:1\n";}}' | awk -F ':' 'BEGIN{nerr=0}{if($4==0 && $6>1) nerr+=$NF}END{print nerr}'`
  NUMIND=`grep --text -v '^#' $BASM.vcf  |perl -ane '{if(length($F[3])>1 || length($F[4])>1){$nerr=abs(length($F[3])-length($F[4]));print "$F[9]:$nerr\n";}}' | awk -F ':' 'BEGIN{nerr=0}{if($4==0 && $6>1) nerr+=$NF}END{print nerr}'`
  ASMSIZE=`ufasta n50 -S $ASM | awk '{print $2}'`
  NUMERR=$(($NUMSUB+$NUMIND))
  QUAL=`echo $NUMERR $ASMSIZE | awk '{print 100-$1/$2*100}'`
  QV=`perl -e '{$erate='$NUMERR'/'$ASMSIZE';printf("%.2f\n", -10*log($erate+0.0000000001)/log(10))}'`
  echo "Stats BEFORE polishing:" > $BASM.report
  echo "Substitution Errors Found: $NUMSUB" >> $BASM.report
  echo "Insertion/Deletion Errors Found: $NUMIND" >> $BASM.report
  echo "Assembly Size: $ASMSIZE" >> $BASM.report
  echo "Consensus Quality Before Polishing: $QUAL" >> $BASM.report
  echo "Consensus QV Before Polishing: $QV" >> $BASM.report
  touch $BASM.report.success
fi

cat $BASM.report

if [ $FIX -gt 0 ];then
  log "Success! Final report is in $BASM.report; polished assembly is in $BASM.PolcaCorrected.fa"
else
  log "Success! Final report is in $BASM.report"
fi
