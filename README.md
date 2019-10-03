# PIPELINES

Partitions:
hpc
vgl
vgl_bigmem

Paths:
export root="/rugpfs/fs0/vgl/store/vglshare/edwin"
export VGP_PIPELINE="/rugpfs/fs0/vgl/store/vglshare/tools/VGP-pipeline"
export tools="/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools"

# VGP pipeline

## purge dups
	sh $VGP_PIPELINE/purge_dups/_submit_purge_dups.sh <asm> <path to *.subreads.bam files> <partition>

## scaff10x
	conda activate VGP
	sh $VGP_PIPELINE/scaff10x/_submit_scaff10x.sh <asm> <path to 10X fastq files> <partition>
	
## solve
	conda activate bionano
	sh $VGP_PIPELINE/bionano/_submit_hybrid_scaffold_dle1.sh <asm> <cmap> <partition>
	sh $VGP_PIPELINE/bionano/_submit_trimNs.sh <asm> <partition>

## salsa
	conda activate VGP
	sh $VGP_PIPELINE/salsa/_submit_salsa2.sh <asm> <path to Hi-C fastq files> <partition>

## arrow polishing
	sh $VGP_PIPELINE/arrow/_submit_arrow.sh <asm> <path to *.subreads.bam files> <partition> <cpus>

# Data and asm QC	
## genomescope
	conda activate VGP
	assumes *_R?_001.fastq.gz in same directory
	sh $VGP_PIPELINE/meryl2/_submit_meryl2_10x.sh 31 <genomeId> <partition> <cpus>
## N50 QC
	conda activate VGP
	$VGP_PIPELINE/stats/asm_stats.sh <fasta file> <genome size bp> c/p
## MashMap
	conda activate mash
	sh $VGP_PIPELINE/mashmap/_submit_mashmap_genome_to_genome.sh <genome1 fasta> <genome2 fasta>
## Mash
	conda activate base
	python /rugpfs/fs0/vgl/store/vglshare/edwin/scripts/fetch_PB.py <genomeId>
	conda activate mash
	ls *.bam > input.fofn
	sh $VGP_PIPELINE/mash/_submit_mash.sh <genomeId> <partition>
## BUSCO
	conda activate busco
	sh $VGP_PIPELINE/busco/_submit_busco.sh <genome fasta>

## Blast a subset of reads
	conda activate VGP
	sh $VGP_PIPELINE/blast/_submit_blast.sh <bam file> <partition> <cpus> <#reads>
	cat *blast_results_* > <blast output file>
	python $VGP_PIPELINE/blast/parseblast.py <blast output file> <list of taxonomic groups>
	
	###example:
	sh $VGP_PIPELINE/blast/_submit_blast.sh /ru-auth/local/home/smrtanalysis2/store/data_root/r64055_20190930_191719/1_A01/m64055_190930_192559.subreads.bam hpc 24 1000
	cat *blast_results_* > m64055_190930_192559.subreads_blast_results.tb
	python $VGP_PIPELINE/blast/parseblast.py m64055_190930_192559.subreads_blast_results.tb Chondrichthyes,Teleostei,Coelacanthiforme,Tetrapoda,Platyhelminthes,Protostomia,Viridiplantae,Fungi,Bacteria
	
# Misc tools
## Minimap2 (raw PB reads)
	conda activate VGP
	ls *.fasta.gz > input.fofn
	sh $VGP_PIPELINE/minimap2/_submit_minimap2.sh <reference fasta>