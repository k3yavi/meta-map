PUFFPATH=$1
BASEPATH="/mnt/scratch2/avi/meta-map/"
INDEX=$BASEPATH"pufIndexEuk"
R1=$BASEPATH"reads/A1_1.fastq"
R2=$BASEPATH"reads/A1_2.fastq"
OUT=$PWD/output

echo "Using/Creating output directory for output"
echo $OUT
mkdir -p $OUT

echo "Running PuffMap"
$PUFFPATH align -i $INDEX --mate1 $R1 --mate2 $R2 -p 15 -m -o $OUT/A1.sam

echo "Getting stats"
python $BASEPATH/pipeline/get_stats.py --sam $OUT/A1.sam --fq $R1
