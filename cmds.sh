# command to download the refs (no Eukaryotes yet)
#parallel -j 10 "wget https://www.ebi.ac.uk/ena/data/view/{}\&display\=fasta -O {}.fa" ::: `cat ./refs/ref_list.txt`

./bin/pufferfish/build/src/fixFasta -k 31 -i refs/combine.fasta -o refs/corr.fasta

#run TwoPaCo
echo "1"
mkdir -p gfa
mkdir -p temp
./bin/TwoPaCo/build/graphconstructor/twopaco -k 31 -f 20 -t 20 --tmpdir temp/ --outfile gfa/de_bruijn_31.bin refs/corr.fasta

echo "2"
./bin/TwoPaCo/build/graphdump/graphdump -k 31 -s refs/corr.fasta -f gfa1 ./gfa/de_bruijn_31.bin > gfa/k31.gfa

echo "3"
./bin/pufferfish/build/src/pufferize -k 31 -g gfa/k31.gfa -f refs/corr.fasta -o gfa/puff_31.gfa > gfa/puffize_log 2>&1

echo "4"
mkdir -p pufIndex
./bin/pufferfish/build/src/pufferfish index -k 31 -o pufIndex -g gfa/puff_31.gfa -r refs/corr.fasta  > gfa/index_log 2>&1
