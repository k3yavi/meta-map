#for x in `ls KrakenDB/library/Bacteria`
#do
#echo -n $x >> taxa.txt
#for y in `ls KrakenDB/library/Bacteria/$x`
#        do
#        echo -n ,`head -1 KrakenDB/library/Bacteria/$x/$y | awk '{print $1}'` >> taxa.txt
#        done
#echo >> taxa.txt
#done

#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_HC1.fasta.gz -p 15 -m -o puff/dmps/HC1_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_HC2.fasta.gz -p 15 -m -o puff/dmps/HC2_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC1.fasta.gz -p 15 -m -o puff/dmps/LC1_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC2.fasta.gz -p 15 -m -o puff/dmps/LC2_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC3.fasta.gz -p 15 -m -o puff/dmps/LC3_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC4.fasta.gz -p 15 -m -o puff/dmps/LC4_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC5.fasta.gz -p 15 -m -o puff/dmps/LC5_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC6.fasta.gz -p 15 -m -o puff/dmps/LC6_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC7.fasta.gz -p 15 -m -o puff/dmps/LC7_unq.dmp  -k 
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC8.fasta.gz -p 15 -m -o puff/dmps/LC8_unq.dmp  -k 


../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_HC1.fasta.gz -p 15 -m -o puff/N_dmps/HC1.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_HC2.fasta.gz -p 15 -m -o puff/N_dmps/HC2.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC1.fasta.gz -p 15 -m -o puff/N_dmps/LC1.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC2.fasta.gz -p 15 -m -o puff/N_dmps/LC2.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC3.fasta.gz -p 15 -m -o puff/N_dmps/LC3.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC4.fasta.gz -p 15 -m -o puff/N_dmps/LC4.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC5.fasta.gz -p 15 -m -o puff/N_dmps/LC5.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC6.fasta.gz -p 15 -m -o puff/N_dmps/LC6.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC7.fasta.gz -p 15 -m -o puff/N_dmps/LC7.dmp  -k --scoreRatio 1.0
../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC8.fasta.gz -p 15 -m -o puff/N_dmps/LC8.dmp  -k --scoreRatio 1.0

#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_HC1.fasta.gz --gzip-compressed > krakOut/HC1.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_HC2.fasta.gz --gzip-compressed > krakOut/HC2.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC1.fasta.gz --gzip-compressed > krakOut/LC1.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC2.fasta.gz --gzip-compressed > krakOut/LC2.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC3.fasta.gz --gzip-compressed > krakOut/LC3.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC4.fasta.gz --gzip-compressed > krakOut/LC4.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC5.fasta.gz --gzip-compressed > krakOut/LC5.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC6.fasta.gz --gzip-compressed > krakOut/LC6.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC7.fasta.gz --gzip-compressed > krakOut/LC7.krk
#kraken --db ./KrakenDB --threads 10 --fasta-input reads/Huttenhower_LC8.fasta.gz --gzip-compressed > krakOut/LC8.krk

#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_HC1.fasta.gz -p 15 -m -o puff/dmps/HC1_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_HC2.fasta.gz -p 15 -m -o puff/dmps/HC2_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC1.fasta.gz -p 15 -m -o puff/dmps/LC1_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC2.fasta.gz -p 15 -m -o puff/dmps/LC2_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC3.fasta.gz -p 15 -m -o puff/dmps/LC3_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC4.fasta.gz -p 15 -m -o puff/dmps/LC4_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC5.fasta.gz -p 15 -m -o puff/dmps/LC5_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC6.fasta.gz -p 15 -m -o puff/dmps/LC6_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC7.fasta.gz -p 15 -m -o puff/dmps/LC7_unq.dmp -k
#../../pufferfish/build/src/pufferfish align -i puff/index --read reads/Huttenhower_LC8.fasta.gz -p 15 -m -o puff/dmps/LC8_unq.dmp -k
