# pipelines

Paths:
export root="/rugpfs/fs0/vgl/store/vglshare/edwin"
export VGP_PIPELINE="/rugpfs/fs0/vgl/store/vglshare/tools/VGP-pipeline"
export tools="/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools"

### genomescope
	conda activate VGP
	assumes *_R?_001.fastq.gz in same directory
	sh $VGP_PIPELINE/meryl2/_submit_meryl2_10x.sh 31 genomeId

### N50 QC
	conda activate VGP
	$VGP_PIPELINE/stats/asm_stats.sh <fasta file> <genome size bp> c/p
or for simple output:
	/rugpfs/fs0/vgl/store/ofedrigo/tools/assembly-stats_bin/assembly-stats <fasta file>

### MashMap
	conda activate mash
	sh $VGP_PIPELINE/mashmap/_submit_mashmap_genome_to_genome.sh <genome1 fasta> <genome2 fasta>

### Mash
	conda activate base
	python /rugpfs/fs0/vgl/store/vglshare/edwin/scripts/fetch_PB.py <genomeId>
	conda activate VGP
	ls *.bam > input.fofn
	sh $VGP_PIPELINE/mash/_submit_mash.sh <genomeId> vgl

### BUSCO
	conda activate busco
	sh $VGP_PIPELINE/busco/_submit_busco.sh <genome fasta>

### Minimap2 (raw PB reads)
	conda activate VGP
	ls *.fasta.gz > input.fofn
	sh $VGP_PIPELINE/minimap2/_submit_minimap2.sh <reference fasta>

### Bionano Solve
	conda activate bionano
	ln -s <assembly.fasta> asm.fasta
	ln -s <bionano.cmap> DLE1.cmap
	sh $VGP_PIPELINE/bionano/_submit_hybrid_scaffold_dle1.sh <genomeId>
