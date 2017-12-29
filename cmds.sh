# command to download the refs (no Eukaryotes yet)
parallel -j 10 "wget https://www.ebi.ac.uk/ena/data/view/{}\&display\=fasta -O {}.fa" ::: `cat ./refs/ref_list.txt`
