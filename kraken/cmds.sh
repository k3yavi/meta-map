for x in `ls KrakenDB/library/Bacteria`
do
echo -n $x >> taxa.txt
for y in `ls KrakenDB/library/Bacteria/$x`
        do
        echo -n ,`head -1 KrakenDB/library/Bacteria/$x/$y | awk '{print $1}'` >> taxa.txt
        done
echo >> taxa.txt
done
