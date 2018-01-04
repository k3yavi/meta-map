PUFFPATH=$1
BASEPATH="/mnt/scratch2/avi/meta-map/kraken/"
INDEX=$BASEPATH"puff/index"
R1=$BASEPATH"reads/Huttenhower_HC1.fasta"
OUT=$PWD/output

echo "Using/Creating output directory for output"
echo $OUT
mkdir -p $OUT

echo "Running PuffMap"
$PUFFPATH align -i $INDEX --mate1 $R1 -p 15 -m -o $OUT/A1.sam

echo "Getting stats"
python $BASEPATH/pipeline/get_mason_stats.py --sam $OUT/A1.sam --fq $R1
