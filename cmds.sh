# command to download the refs (no Eukaryotes yet)
parallel -j 10 "wget https://www.ebi.ac.uk/ena/data/view/{}\&display\=fasta -O {}.fa" ::: `cat ./refs/ref_list.txt`

#run TwoPaCo
mkdir -p gfa
mkdir -p temp
./bin/TwoPaCo/build/graphconstructor/twopaco -k 25 -f 20 -t 15 --tmpdir temp/ --outfile gfa refs/fasta/*

./bin/TwoPaCo/build/graphdump/graphdump -k 25 <(cat ./bin/TwoPaCo/2paco-s.format) -f gfa1 ./gfa/de_bruijn.bin > gfa/k25.gfa
