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

## falcon and falcon_unzip
	tmux new -s mysession
	conda activate VGP
	bam2fasta -o projectName *.subreads.bam
	conda activate denovo_asm
	ls *.fasta.gz > input.fofn
	fc_run fc_run.cf
	ls *.subreads.bam > input_bam.fofn
	fc_unzip.py fc_unzip.cfg 

## purge dups
	conda activate VGP
	sh $VGP_PIPELINE/purge_dups/_submit_purge_dups.sh <asm> <path to *.fasta files> <ploidy_mode> <rm_OVLP_only> <partition> <cpus>
	# a purge_dups_ccs.sh script is available when working with ccs reads. need to make input.fofn first before using it
	ls /path/to/files/*.fastq.gz > input.fofn
	sh $VGP_PIPELINE/purge_dups/_submit_purge_dups_ccs.sh <asm> <path to *.fasta files> <ploidy_mode> <rm_OVLP_only> <partition> <cpus>

### example:
	sh $VGP_PIPELINE/purge_dups/_submit_purge_dups.sh ../sCarCar2_p1_arrowed.fasta /vggpfs/fs3/vgl/scratch/vglshare/sandbox/ofedrigo/sCarCar2/assembly_vgp/intermediates/purge_dups/fasta/ diploid false vgl 32

## scaff10x
	conda activate VGP
	sh $VGP_PIPELINE/scaff10x/_submit_scaff10x.sh <asm> <partition> <path to 10X fastq files>
	
## solve
	conda activate bionano
	sh $VGP_PIPELINE/bionano/_submit_hybrid_scaffold_dle1.sh <asm> <cmap> <partition>
	sh $VGP_PIPELINE/bionano/_submit_trimNs.sh <asm> <partition>

## salsa
	conda activate VGP
	sed 's/:/_/g' <asm> > <asm_renamed>
	sh $VGP_PIPELINE/salsa/_submit_salsa_2.2.sh <asm_renamed> <path to Hi-C fastq files> <partition> <cpus>

## arrow polishing
	conda activate VGP
	sh $VGP_PIPELINE/arrow/_submit_arrow.sh <asm> <path to *.subreads.bam files> <partition> <cpus>

# Data and asm QC	

## fastqc
	conda activate VGP
	sh $VGP_PIPELINE/qc/_submit_fastqc.sh <directory where fastq.gz files are> <partition>
	
## trimming Illumina reads
	sh $VGP_PIPELINE/qc/_submit_trimming.sh <directory where fastq.gz files are> <partition>

## plot read length distribution
    conda activate VGP
    sh $tools/plots/_submit_readlength.sh <fasta file>

## bam to fastq
	conda activate VGP
	sbatch -e err.log -o out.log <<"EOF"
	#!/bin/bash
	#partition=vgl
	#SBATCH -J bam2fastq
	#SBATCH -n 10
	#SBATCH -t 10:00:00
	for qry in *.bam; do
	bam2fastq $qry -o $(basename "$qry"); 
	done;
	EOF
	
## genomescope
	conda activate VGP
	assumes *_R?_001.fastq.gz in same directory
	sh $VGP_PIPELINE/meryl2/_submit_meryl2_10x.sh 31 <genomeId> <partition> <cpus>

## N50 QC
	conda activate VGP
	$VGP_PIPELINE/stats/asm_stats.sh <fasta file> <genome size bp> c/p

## kat
	conda activate kat
	ls <path to 10X data>/* > input.fofn
	sh $VGP_PIPELINE/kat/_submit_kat_comp.sh <asm fasta file>

## MERQURY	
	conda activate merqury
	export MERQURY="/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools/merqury"
	export PATH=$PATH:/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools/meryl/Linux-amd64/bin	

	### Build meryl DB
	WGS Illumina data
	ls *.fastq.gz > input.fofn # important: no full path in symlink!
	sbatch --partition=vgl $tools/merqury/_submit_build.sh 21 input.fofn <out_prefix> mem=F
	or
	10X data
	ls *_R1_001.fastq.gz > > R1.fofn # important: no full path in symlink!
	ls *_R2_001.fastq.gz > > R2.fofn # important: no full path in symlink!
	sbatch --partition=vgl $tools/merqury/_submit_build_10x.sh 21 R1.fofn R2.fofn <out_prefix> mem=F

	### run merqury
	sbatch --partition=vgl $tools/merqury/_submit_merqury.sh <out_prefix>.meryl pri.fasta [alt.fasta] <output>

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
	sh $VGP_PIPELINE/busco/_submit_busco.sh <genome fasta> <partition> <cpus>

## Blast a subset of reads
	conda activate VGP
	sh $VGP_PIPELINE/blast/_submit_blast.sh <bam file> <partition> <cpus> <#reads>
	if you want to directly blast a fasta file: sh $VGP_PIPELINE/blast/_submit_blast_fasta_.sh <fasta file> <partition> <cpus> <#reads>
	cat *blast_results_* > <blast output file>
	python $VGP_PIPELINE/blast/parseblast.py <blast output file> <list of taxonomic groups>
	
	
### example:
	sh $VGP_PIPELINE/blast/_submit_blast.sh /ru-auth/local/home/smrtanalysis2/store/data_root/r64055_20190930_191719/1_A01/m64055_190930_192559.subreads.bam hpc 24 1000
	cat *blast_results_* > m64055_190930_192559.subreads_blast_results.tb
	python $VGP_PIPELINE/blast/parseblast.py m64055_190930_192559.subreads_blast_results.tb Chondrichthyes,Teleostei,Coelacanthiforme,Tetrapoda,Platyhelminthes,Protostomia,Viridiplantae,Fungi,Bacteria
	
# Misc tools
## Minimap2 (raw PB reads or Iso-seq)
	conda activate VGP
	ls *.fast?.gz > input.fofn
	sh $VGP_PIPELINE/minimap2/_submit_minimap2.sh <reference fasta> <datatype: genome or isoseq> <cpus> <partition>

## RepeatMasker
	conda activate VGP
	sh $VGP_PIPELINE/repeatmasker/_submit_repeatmasker.sh <reference fasta> <repeat modeler fasta file>
	
## Isoseq
