ksize = config["ksize"]
bwa = config["bwa"]
debga = config["debga"]
kallisto = config["kallisto"]
twopaco = config["twopaco"]
pufferize = config["pufferize"]
puffer = config["pufferfish"]
data_path  = config["data_path"]
output_path  = config["output_path"]

de_bruijn = os.path.sep.join([data_path, "de_bruijn.bin"])

#human_txome_ref = config["human_txome_ref"]
#human_genome_ref = config["human_genome_ref"]
bacterial_genome_ref = config["bacterial_genome_ref"]
datasets = [bacterial_genome_ref]#human_txome_ref, human_genome_ref]

#human_txome_read = config["human_txome_read"]
#human_genome_read = config["human_genome_read"]
bacterial_genome_read = config["bacterial_genome_read"]



rule all:
     input:
      expand("{out}/k{k}_n_{outfiles}.bwa_idx", out=output_path, k=ksize, outfiles=datasets),
      expand("{out}/k{k}_n_{outfiles}.kallisto_idx", out=output_path, k=ksize, outfiles=datasets),
      expand("{out}/k{k}_n_{outfiles}.puffer_idx/seq.bin", out=output_path, outfiles=datasets, k=ksize),
      expand("{out}/k{k}_n_{outfiles}.puffer_sparse_idx/seq.bin", out=output_path, outfiles=datasets, k=ksize),
      
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.bwa.lookup.benchmark.txt", out=output_path, ref=human_txome_ref, read=human_txome_read, k=ksize),
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.bwa.lookup.benchmark.txt", out=output_path, ref=human_genome_ref, read=human_genome_read, k=ksize),
      expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.bwa.lookup.benchmark.txt", out=output_path, ref=bacterial_genome_ref, read=bacterial_genome_read, k=ksize),
      
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.kallisto.lookup.benchmark.txt", out=output_path, ref=human_txome_ref, read=human_txome_read, k=ksize),
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.kallisto.lookup.benchmark.txt", out=output_path, ref=human_genome_ref, read=human_genome_read, k=ksize),
      expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.kallisto.lookup.benchmark.txt", out=output_path, ref=bacterial_genome_ref, read=bacterial_genome_read, k=ksize),
      
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.puffer.lookup.benchmark.txt", out=output_path, ref=human_txome_ref, read=human_txome_read, k=ksize),
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.puffer.lookup.benchmark.txt", out=output_path, ref=human_genome_ref, read=human_genome_read, k=ksize),
      expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.puffer.lookup.benchmark.txt", out=output_path, ref=bacterial_genome_ref, read=bacterial_genome_read, k=ksize),
      
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.puffer.sparse.lookup.benchmark.txt", out=output_path, ref=human_txome_ref, read=human_txome_read, k=ksize),
#expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.puffer.sparse.lookup.benchmark.txt", out=output_path, ref=human_genome_ref, read=human_genome_read, k=ksize)
      expand("{out}/benchmarks/k{k}_n_{ref}_vs_{read}.puffer.sparse.lookup.benchmark.txt", out=output_path, ref=bacterial_genome_ref, read=bacterial_genome_read, k=ksize)


rule bwa_lookup:
     input :
           index = os.path.sep.join([output_path, "k{ksize}_n_{ref}.bwa_idx"]),
           reads = os.path.sep.join([data_path, "{reads}.fa"])
     output :
           os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.bwa.lookup.benchmark.txt"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.bwa.lookup.benchmark.txt"])
     message:
          bwa + " fastmap {input.index} {input.reads}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}_vs_{reads}.bwa.lookup.log"])
     shell :      
          bwa + " fastmap {input.index} {input.reads} > {log} 2>&1"


rule kallisto_lookup:
     input :
           index = os.path.sep.join([output_path, "k{ksize}_n_{ref}.kallisto_idx"]),
           reads = os.path.sep.join([data_path, "{reads}.fa"])
     output:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.kallisto.lookup.benchmark.txt"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.kallisto.lookup.benchmark.txt"])
     message:
          kallisto + " lookup -i {input.index} {input.reads}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}_vs_{reads}.kallisto.lookup.log"])
     shell :
          kallisto + " lookup -i {input.index} {input.reads} > {log} 2>&1"


rule puffer_lookup:
     input :
           index = os.path.sep.join([output_path, "k{ksize}_n_{ref}.puffer_idx/seq.bin"]),
           reads = os.path.sep.join([data_path, "{reads}.fa"])
     output:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.puffer.lookup.benchmark.txt"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.puffer.lookup.benchmark.txt"])
     message:
          puffer + " lookup -i {input.index} -r {input.reads}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}_vs_{reads}.puffer.lookup.log"])
     run :
          index_dir = str(input.index).rsplit("/",1)[0]
          shell("{puffer} lookup -i {index_dir} -r {input.reads} > {log} 2>&1")

rule puffer_sparse_lookup:
     input :
           index = os.path.sep.join([output_path, "k{ksize}_n_{ref}.puffer_sparse_idx/seq.bin"]),
           reads = os.path.sep.join([data_path, "{reads}.fa"])
     output:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.puffer.sparse.lookup.benchmark.txt"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}_vs_{reads}.puffer.sparse.lookup.benchmark.txt"])
     message:
          puffer + " lookup -i {input.index} -r {input.reads}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}_vs_{reads}.puffer.sparse.lookup.log"])
     run :
          index_dir = str(input.index).rsplit("/",1)[0]
          shell("{puffer} lookup -i {index_dir} -r {input.reads} > {log} 2>&1")


rule bwa_index:
     input :
           os.path.sep.join([data_path, "{ref}.fa"])
     output :
           os.path.sep.join([output_path, "k{ksize}_n_{ref}.bwa_idx"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}.bwa.index.benchmark.txt"])
     message:
          bwa + " index  -p {output} {input}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}.bwa.index.log"])
     shell :
           "touch {output}; {bwa} index  -p {output} {input} > {log} 2>&1"


rule kallisto_index:
     input :
           os.path.sep.join([data_path, "{ref}.fa"])
     output :
           os.path.sep.join([output_path, "k{ksize}_n_{ref}.kallisto_idx"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}.kallisto.index.benchmark.txt"])
     message:
          kallisto + " index -k {ksize} -i {output} {input}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}.kallisto.index.log"])
     shell :
           kallisto + " index -k {ksize} -i {output} {input} > {log} 2>&1"

rule puffer_twopaco:
     input :
           fastafile = os.path.sep.join([data_path, "{ref}.fa"]),
           de_bruijn = os.path.sep.join([output_path, "de_bruijn.bin"])		   
     output :
           os.path.sep.join([data_path, "k{ksize}_n_{ref}.gfa"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}.puffer.twopaco.benchmark.txt"])
     message:
          "{twopaco}/graphconstructor/twopaco -k {ksize} -t 15 -f 20 {input.fastafile} --outfile {input.de_bruijn} --tmpdir {output_path}/twopacoTmp\n{twopaco}/graphdump/graphdump -k {ksize} -s {input.fastafile} -f gfa1 {input.de_bruijn} > {output}"
     shell :
          "mkdir -p {output_path}/twopacoTmp &&"          
          "{twopaco}/graphconstructor/twopaco -k {ksize} -t 15 -f 20 {input.fastafile} --outfile {input.de_bruijn} --tmpdir {output_path}/twopacoTmp  && "
          "{twopaco}/graphdump/graphdump -k {ksize} -s {input.fastafile} -f gfa1 {input.de_bruijn} > {output}"

rule puffer_pufferize:
     input :
           gfafile = os.path.sep.join([data_path, "k{ksize}_n_{ref}.gfa"]),
           fastafile = os.path.sep.join([data_path, "{ref}.fa"])
     output :
           os.path.sep.join([data_path, "k{ksize}_n_{ref}.pufferized.gfa"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}.puffer.pufferize.benchmark.txt"])
     message:
           "{pufferize} -k {ksize} -g {input.gfafile} -f {input.fastafile} -o {output} > {log} 2>&1"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}.puffer.pufferize.log"])
     shell :
           "{pufferize} -k {ksize} -g {input.gfafile} -f {input.fastafile} -o {output} > {log} 2>&1"

rule puffer_index:
     input :
           os.path.sep.join([data_path, "k{ksize}_n_{ref}.pufferized.gfa"])
     output :
           os.path.sep.join([output_path, "k{ksize}_n_{ref}.puffer_idx","seq.bin"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}.puffer.index.benchmark.txt"])
     message:
          puffer + " index -k {ksize} -o {output} -g {input}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}.puffer.index.log"])
     run :
          output_dir= str(output).rsplit("/",1)[0]
          shell("rm -rf {output}; {puffer} index -k {ksize} -o {output_dir} -g {input} > {log} 2>&1")

rule puffer_index_sparse:
     input :
           os.path.sep.join([data_path, "k{ksize}_n_{ref}.pufferized.gfa"])
     output :
           os.path.sep.join([output_path, "k{ksize}_n_{ref}.puffer_sparse_idx","seq.bin"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_n_{ref}.puffer.index.sparse.benchmark.txt"])
     message:
          puffer + " index -s -k {ksize} -o {output} -g {input}"
     log:
          os.path.sep.join([output_path, "logs/k{ksize}_n_{ref}.puffer.index.sparse.log"])
     run :
          output_dir= str(output).rsplit("/",1)[0]
          shell("rm -rf {output}; {puffer} index -s -k {ksize} -o {output_dir} -g {input} > {log} 2>&1")

rule debga_index:
     input :
           os.path.sep.join([data_path, "{ref}.fa"])
     output :
           os.path.sep.join([output_path, "k{ksize}_{ref}.debga_idx"])
     benchmark:
          os.path.sep.join([output_path, "benchmarks/k{ksize}_{ref}.debga.index.benchmark.txt"])
     message:
          debga + " index -k {ksize} {input} {output}"
     shell :
           debga + " index -k {ksize} {input} {output}"


