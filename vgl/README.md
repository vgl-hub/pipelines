all code blocks assume i am using my conda unless otherwise specified

VGP conda: `$VGPCONDA=/rugpfs/fs0/vgl/store/vglshare/tools/miniconda3/bin/conda`

# cluster/general commands
sbatch, request one node, all 32 cores on it, in the vgl partition, and name the job xyz: `sbatch -p vgl -n 1 -c 32 -J xyz`

squeue vgl partition with jobID, user, partition, nodelist(reason), start_time, state, and name

`squeue -p vgl -o "%A %u %P %R %S %T %j"` add a `%Z` for job working directory to be nosy. pipe into `tr " " "\t"` for tab delim

check disk usage in human readable format: `du -sh`

check disk free in human readable format: `df -h`

ssh listening command from laptop: `ssh -f -N -L 8080:localhost:8080 rocku`

run galaxy detached: `nohup sh galaxy/run.sh > nohup_outs/20220308.1850.out 2>&1 &`

view top for my jobs only: `top -u labueg` then hit `c` to expand command

jump to line 123 in vim: `vim +123 file.txt` and then delete next four words `4dw` and join next three lines `3J`

fix the "subseq:" parts of bionano AGP: `sed 's/W\t\(.*\)_subseq_\([0-9]*\):\([0-9]*\)\t[0-9]*\t[0-9]*\t\(.\)/W\t\1\t\2\t\3\t\4/g' old.agp > new.agp`

check if two files are the same: `cmp --silent $old $new || echo "files are different"`

create conda env with more recent c++ compilers
````bash
conda install -c conda-forge gxx_linux-64
conda install -c conda-forge cxx-compiler
````

aws commands
````bash
## upload ./test.txt to the fScoJap1 folder
## don't forget the last slash, so it knows to put it INTO the folder, instead of creating a new thing
aws s3 cp ./test.txt s3://genomeark/species/Scomber_japonicus/fScoJap1/
aws s3 cp ./bTheCae1.alt.cur.nopipe.20210525.fasta.gz s3://genomeark/working/NCBI_nopipe_asms/bTheCae1.alt.cur.nopipe.20210525.fasta.gz

## download only the fastq.gz files from these folders
aws s3 cp --recursive --exclude "*" --include "*.fastq.gz" s3://genomeark/species/Taeniopygia_guttata/bTaeGut2/genomic_data/hifi/ . --no-sign-request
aws s3 cp --recursive --exclude "*" --include "*.fq.gz" s3://genomeark/species/Taeniopygia_guttata/bTaeGut2/genomic_data/arima/ . --no-sign-request

## download everything from this folder
aws s3 cp --recursive s3://genomeark/species/Taeniopygia_guttata/bTaeGut2/genomic_data/bionano/ . --no-sign-request

## upload all fasta.gz files from current dir to genomeark
aws s3 cp --recursive --exclude "*" --include "*.fasta.gz" . s3://genomeark/species/Scomber_japonicus/fScoJap1/assembly_vgp_standard_2.0/intermediates/
````

evaluating gap fills from TGS gap closer
````bash
## generate BED file and parseable file from prefix.gap_fill_detail file output from TGS gap closer
sh parseTGSdetailfile.sh prefix.gap_fill_detail

## generate snapshot batch script for igv using 20 random gap fills from the bed file
shuf -n 20 <(grep "GapFill" prefix.gap_fill_detail.bed) > prefix.20random.bed
sh bed2snapshotscript.sh prefix.20random.bed 250 250 ./snapshots/ prefix
````

# program commands

## hifiasm default + HiC mode
````bash
$STORE/programs/hifiasm/hifiasm -o rEryReg1 --h1 ../genomic_data/arima/rEryReg1_Royal_Ground_Snake_R1.fastq.gz --h2 ../genomic_data/arima/rEryReg1_Royal_Ground_Snake_R2.fastq.gz ../rEryReg1.trimmed.cat.fastq.gz 2> rEryReg1.log
````

## hifiasm trio mode
aPseCor1 = mom, aPseCor2 = dad, aPseCor3 = offspring

1) meryl database for mom, dad, child (including hapmer DB for parents)
2) yak for mom, dad
3) hifiasm on child, using mom and dad yak
4) merqury on child assemblies, using mom and dad hapmers

````bash
## MERYL DATABASES FOR ALL
# i don't activate the VGP conda env because that one has a different meryl version than this script needs

# parents:
mkdir mom_meryl
cd mom_meryl
ls $parentdir/genomic_data/mom_*.fastq.gz > input.fofn
sh $VGP_PIPELINE/meryl2/_submit_meryl2_build.sh 21 input.fofn <out_prefix> vgl
cd ..
cd dad_meryl
ls $parentdir/genomic_data/dad_*.fastq.gz > input.fofn
sh $VGP_PIPELINE/meryl2/_submit_meryl2_build.sh 21 input.fofn <out_prefix> vgl

# child:
same as above lol
sh $VGP_PIPELINE/meryl2/_submit_meryl2_build_ccs.sh 21 input.fofn <out_prefix> vgl

## MERYL HAPMER DATABASE FOR PARENTS
sh $tools/merqury/_submit_hapmers.sh <mat.meryl> <pat.meryl> [child.meryl]
# will create mat.hapmer.meryl and pat.hapmer.meryl, this is what you need to evaluate trios

## YAK DATABASES FOR MOM AND DAD
# concatenate the parental reads into one file first
sbatch --partition=vgl --nodes=1 --cpus-per-task=32 --wrap="$STORE/programs/yak/yak count -b37 -t32 -o aPseCor1.yak aPseCor1_cat.fastq.gz"
sbatch --partition=vgl --nodes=1 --cpus-per-task=32 --wrap="$STORE/programs/yak/yak count -b37 -t32 -o aPseCor2.yak aPseCor2_cat.fastq.gz"

## HIFIASM IN TRIO MODE
$STORE/programs/hifiasm/hifiasm -o aPseCor3.20220409 -1 ../../parents/2a_yak_aPseCor1/aPseCor1.yak -2 ../../parents/2b_yak_aPseCor2/aPseCor2 -t 32 <child reads concatenated>.fastq.gz 2> aPseCor3.20220409.log

## MERQURY IN TRIO MODE
# _submit_merqury.sh <read-db.meryl> [<mat.meryl> <pat.meryl>] <asm1.fasta> [asm2.fasta] <out>
sbatch $tools/merqury/_submit_merqury.sh ../../1b_meryl/aPseCor3.meryl \
../../../parents/1c_merylDB_hapmers/aPseCor1.hapmer.meryl/ \
../../../parents/1c_merylDB_hapmers/aPseCor2.hapmer.meryl/ \
../aPseCor3.20220409.dip.hap1.p_ctg.fasta ../aPseCor3.20220409.dip.hap2.p_ctg.fasta \
aPseCor3.20220409.trio
sbatch $tools/merqury/_submit_merqury.sh ../../1b_meryl/aPseCor3.meryl ../../../parents/1c_merylDB_hapmers/aPseCor1.hapmer.meryl/ ../../../parents/1c_merylDB_hapmers/aPseCor2.hapmer.meryl/ ../aPseCor3.20220409.dip.hap1.p_ctg.fasta ../aPseCor3.20220409.dip.hap2.p_ctg.fasta aPseCor3.20220409.trio

## command for coloring nodes in the graph by hapmer content
## COMMAND FROM ANN MCCARTNEY!
sh $MERQURY/trio/hap_blob.sh maternal.k30.hapmer.meryl paternal.k30.hapmer.meryl assembly.fa prefix_output
cat prefix_output.count |grep -v Assembly |awk '{if ($3+$4 < 25) print $2"\t0\t0"; else print $2"\t"$3/($3+$4)*100"\t"$4/($3+$4)*100"\t"$5}'|awk '{if (NR == 1) print "Name,chr,color"; if ($2 > 95) print $1",maternal,khaki"; else if ($3>95) print $1",paternal,darkred"; else print $1",unknown,gray"}' > trio.csv


````

## cutadapt
command taken from galaxy pipeline 12 april 2022
````bash
This is cutadapt 3.5 with Python 3.9.12
Command line parameters: -j=32 -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT --output=out1.fq.gz --error-rate=0.1 --times=1 --overlap=3 --action=trim --revcomp --discard-trimmed m64330e_211110_060034_hifi_reads_fastq_gz.fq.gz
````

## main pipeline QC
MERYL

some of the build scripts in the pipeline dirs don't work.

`sh $tools/merqury/_submit_build.sh 21 input.fofn aPseCor3 T` -> error: invalid partition name (norm) specified

`sh $VGP_PIPELINE/meryl2/_submit_meryl2_ccs.sh 21 input.fofn aPseCor3 vgl` -> error: cannot access *fastq

probably works if your fastq files are ungzipped and in the wd, but it seems the only ccs diff is the genomescope script? 


````bash
ls genomic_data/pacbio/*.fastq.gz > input.fofn
sh $VGP_PIPELINE/meryl2/_submit_meryl2_build.sh 21 input.fofn aPseCor3 vgl
````


MERQURY
````bash
## run Merqury -- needs specific R function
## also using own Meryl (v 1.3), because the shared VGP one is v 1.7
## RENAME THE (ASM).FA TO (ASM).FASTA, THE JAVA SCRIPTS DO NOT LIKE (ASM).FA!!!!
## merqury script will not take paths in the output_prefix
#### cd to the directory you want to put it in, first
$MYCONDA init
conda activate R_3_6
export MERQURY="/rugpfs/fs0/vgl/store/vglshare/tools/VGP-tools/merqury"
sbatch --partition=vgl $tools/merqury/_submit_merqury.sh child.meryl [<mat.hapmer.meryl> <pat.hapmer.meryl>] pri.fasta [alt.fasta] <output>
````

BUSCO
````bash
conda activate busco_test
busco					# gives "no module found" error
which busco				# points to conda busco
conda deactivate		# puts you in base conda
busco					# shouldn't work
which busco				# shouldn't find anything
conda deactivate		# turns conda off entirely
busco					# still shouldn't work
conda activate busco_test
busco					# now it works for some ungodly reason
busco -i <input.fasta> -o <outprefix> -m geno -l $STORE/busco_DBs/vertebrata_odb10 -c 32
````

## misc
nucmer and mummerplot
````bash
## nucmer first
# -l 10000 turns down the noise in the plot, useful when aligning (expected) exact sequences
srun -p vgl -c 32 $STORE/programs/mummer4/bin/nucmer --mum <reference.fasta> <query.fasta> -t 32 -p <prefix> -l 10000

## mummerplot to visualize the nucmer
conda activate gnuplot5
# this one seems to somehow sort the fastas
$STORE/programs/mummer4/bin/mummerplot -large -layout <prefix>.delta --fat -f -t png
# this will have fasta sequences in the order they are in the file
$STORE/programs/mummer4/bin/mummerplot <prefix>.delta -t png

## zooming into specfic coordinates. reads fasta as long linear sequence, so to get the starting
## point of the Z chromosome i summed the lengths of all chromosomes preceding it
# zooms x axis into the Z chromosome for the bTaeGut1.4 reference 
$STORE/programs/mummer4/bin/mummerplot bTaeGut2.hap2.s2.ref.delta -t png -x [937009338:1012405514]
# zooms x axis into the W chromosome for the bTaeGut1.4 reference
$STORE/programs/mummer4/bin/mummerplot bTaeGut2.hap2.s2.ref.delta -t png -x [1012405514:1033251908]
````

blast
````bash
conda activate blast
makeblastdb -in ${W_URL}/db/${ASM%.*} -parse_seqids -dbtype nucl -out ${W_URL}/db/${ASM%.*}.db
blastn -outfmt 6 -query ${W_URL}/reference/${REF} -db ${W_URL}/db/${ASM%.*}.db -out ${W_URL}/${ASM%.*.*}_in.out
````

smudgeplot
````bash
conda activate smudgeplot
cd $ASSEMBLIES/rEryReg1
mkdir smudgeplot && cd smudgeplot
mkdir tmp
ls ../genomic_data/*.fastq.gz > FILES
# kmer 21, 32 threads, 64G of memory, counting kmer coverages between 1 and 10000x
$STORE/programs/KMC/bin/kmc -k21 -t32 -m64 -ci1 -cs10000 @FILES kmcdb tmp
$STORE/programs/KMC/bin/kmc_tools transform kmcdb histogram kmcdb_k21.hist -cx10000
## "The next step is to extract genomic kmers using reasonable coverage thresholds. You can either inspect the kmer spectra and 
## choose the L (lower) and U (upper) coverage thresholds via visual inspection, or you can estimate them using command 
## smudgeplot.py cutoff <kmcdb_k21.hist> <L/U>."
smudgeplot.py cutoff kmcdb_k21.hist L
# Running smudgeplot v0.2.2
# Task: cutoff
# 10
# Done!
smudgeplot.py cutoff kmcdb_k21.hist U
# Running smudgeplot v0.2.2
# Task: cutoff
# 790
# Done!

## alternatively:
L=$(smudgeplot.py cutoff kmcdb_k21.hist L)
U=$(smudgeplot.py cutoff kmcdb_k21.hist U)
echo $L $U # these need to be sane values

## using instructions that use KMC installed from this repo: https://github.com/tbenavi1/KMC
## will have smudge_pairs installed in the bin
## compile it with gxx11 conda env
kmc_tools transform kmcdb -ci"$L" -cx"$U" reduce kmcdb_L"$L"_U"$U"
smudge_pairs kmcdb_L"$L"_U"$U" kmcdb_L"$L"_U"$U"_coverages.tsv kmcdb_L"$L"_U"$U"_pairs.tsv > kmcdb_L"$L"_U"$U"_familysizes.tsv
# now generate the plot
smudgeplot.py plot kmcdb_L"$L"_U"$U"_coverages.tsv
# see plot parameters
smudgeplot.py plot --help
````

using meryl to find contigs that are 1) only in hap2, and 2) are 2-copy
````bash
meryl print difference bTaeGut2_hap2_count bTaeGut2_hap1_count output bTaeGut2_hap2only
meryl equal-to 2 bTaeGut2_hap2_only output bTaeGut2_hap2_only_equalto2
meryl-lookup -existence -sequence bTaeGut2_trio.asm.dip.hap2.p_ctg.fa -mers bTaeGut2_hap2_only_equalto2 > bTaeGut2_hap2_only_equalto2.tsv
````

pretextsnapshot (on my macbook)
````bash
conda activate pretext
PretextSnapshot -m ../bTaeGut2.trio.trim.rebin.s1.hap1.pretext --sequences "=full"
````

exporting a conda environment with URL to install packages in another environment
````bash
conda list -n bionano --explicit > bionano-package-list.expl.txt
conda create -n bionano2 --file bionano-package-list.expl.txt
````

running mitohifi locally
````bash
while read paths; do
	export PATH=$PATH:/lustre/fs5/vgl/store/labueg/programs/mitohifi_deps/$paths/ ;
done < <(ls)
while read paths; do
	export PATH=$PATH:/lustre/fs5/vgl/store/labueg/programs/mitohifi_deps/$paths/bin/ ;
done < <(ls)
sb --wrap='python $STORE/programs/MitoHiFi/mitohifi.py -r "$BTAEGUT2/genomic_data/hifi/m54306U_210519_154448.hifi_reads.fastq.gz $BTAEGUT2/genomic_data/hifi/m54306U_210521_004211.hifi_reads.fastq.gz $BTAEGUT2/genomic_data/hifi/m54306Ue_210629_211205.hifi_reads.fastq.gz $BTAEGUT2/genomic_data/hifi/m54306Ue_210719_083927.hifi_reads.fastq.gz $BTAEGUT2/genomic_data/hifi/m64055e_210624_223222.hifi_reads.fastq.gz" -f DQ453514.1.fasta -g DQ453514.1.gb -t 32 -o 2
````

winnowmap for hg002 (~6 Gb)
```
$STORE/programs/Winnowmap-2.03/bin/winnowmap -I 7G -W ../repetitive_k15.txt -ax map-pb ../../../assembly/assembly.v0.2.fasta ../test.fastq.gz > output.sam ; samtools view -S -b output.sam | samtools sort - -o output.sorted.bam
```

installing verkko
```
curl -L https://github.com/marbl/verkko/releases/download/v1.1/verkko-v1.1.tar.gz --output verkko-v1.1.tar.gz
tar -xzf verkko-v1.1.tar.gz
cd verkko-v1.1/src
conda create -n verkko_src
conda activate verkko_src
conda install gxx
conda install rust
conda install -c bioconda graphaligner
conda install -c bioconda mbg
conda install snakemake
srun -c 32 -p vgl make -j 32 

sb --job-name=verkko --wrap="$STORE/programs/verkko-v1.1/bin/verkko -d asm/ --hifi ../genomic_data/pacbio/demultiplex.bc1001--bc1001.hifi_reads.fastq.gz ../genomic_data/pacbio/m64330e_211008_165718.hifi_reads.fastq.gz --threads 32"
```

getting an agp overlay GFA into bandageNG (17 nov 2022)
```
~/Programs/gfastats/build/bin/gfastats rLiaOli1_hap2_contigs.gfa --discover -a rLiaOli1_hap2_s1.agp -o rLiaOli1_hap2_s1.gfa
~/Programs/gfastats/build/bin/gfastats rLiaOli1_hap2_s1.gfa -a rLiaOli1_hap2_s2.agp -o rLiaOli1_hap2_s2.gfa
# bandageNG doesn't want path lines
grep -v "^P" rLiaOli1_hap2_s2.gfa > rLiaOli1_hap2_s2.NOP.gfa
# add perfect overlaps back from hifiasm, threshold = 1000 bases
~/Programs/gfastats/build/bin/gfastats --discover-perfect-overlaps rLiaOli1_hap2_s2.NOP.gfa -o rLiaOli1_hap2_s2.NOP.perfoverlap.gfa
# can then edit the J lines to be L lines so it will be prettier
### change the identifier at beginning of line; add an M to end of 6th field; and optionally add a tag for line customization, see https://github.com/asl/BandageNG/wiki/Custom-GFA-tags#styling
```

running linux-only docker container on m1 mac
```
docker pull quay.io/biocontainers/fastk:1.0--h3e8787d_1
docker run --platform linux/amd64 -ti -d quay.io/biocontainers/fastk:1.0--h3e8787d_1
# 530cad7d51935ebfdd69af4971066ead64d7795d729358bf043a92ee7727b426
docker image ls
# REPOSITORY                                                                 TAG                                          IMAGE ID       CREATED         SIZE
# quay.io/biocontainers/fastk                                                1.0--h3e8787d_1                              6167f900d4ed   8 days ago      56.5MB
# quay.io/biocontainers/mulled-v2-c9dd01c51b75826c5cd3fe1c076d72edf5c05e93   de0e66755f7dfb7ec893f2c8e64f7cf9db099a96-0   bddaf71cecd3   7 weeks ago     192MB
# quay.io/biocontainers/yahs                                                 1.2a.2--h7132678_0                           885fcd7a9acf   3 months ago    15.2MB
# continuumio/miniconda3                                                     latest                                       3e9fe0c4644a   5 months ago    403MB
# docker/getting-started                                                     latest                                       157095baba98   7 months ago    27.4MB
# bgruening/galaxy-stable                                                    latest                                       e0ca6537abaa   19 months ago   2.07GB
docker container ls
# CONTAINER ID   IMAGE                                         COMMAND                  CREATED         STATUS         PORTS     NAMES
# 530cad7d5193   quay.io/biocontainers/fastk:1.0--h3e8787d_1   "/usr/local/env-exec…"   7 seconds ago   Up 7 seconds             gallant_johnson
[linelle@blaziken ~/Documents/_sandbox/fastk]$ docker exec -it 530cad7d5193 /bin/bash
root@530cad7d5193:/#
```

pulling & running a docker tool and giving it a directory to work in 
```
docker run --rm -ti -v ~/Documents/_sandbox/vcfstats_test/:/data justold/vcfstats:latest vcfstats --vcf /data/masta.ddoc.SNPandINDEL.vcf -o /data/ -f 'COUNT(1) ~ CONTIG' --title test
# -v ~/Documents/_sandbox/vcfstats_test/:/data sets the directory before the colon as the "shortcut" name of /data
```

```
~/Programs/planemo/.venv/bin/planemo serve --biocontainers --galaxy_root ./galaxy/ --galaxy_email user@user.com [tools].xml
```

gfastats transform newline to tab for spreadsheets:
```
gfastats -t mEubGla1.HiC.hap2.20220517.fasta.gz 2661810066 | cut -f 2 | tail -n +2 |  tr "\n" "\t" | pbcopy
```

installing picard on cluster
```
cd $STORE/programs
wget https://download.oracle.com/java/17/latest/jdk-17_linux-x64_bin.tar.gz
tar -xvzf jdk-17_linux-x64_bin.tar.gz
cd BINARIES
ln -s $STORE/programs/jdk-17.0.6/bin/java
ln -s $STORE/programs/jdk-17.0.6/bin/jar
cd ~
which java
cd $STORE/programs/picard/
./gradlew shadowJar
```

promethION / ONT notes: computer is running ubuntu linux x86_64, graphics card is NVIDI GV100GL (QuadroGV100)
1. install some dev packages needed by ONT programs: `sudo apt-get update && sudo apt-get install python3-dev zlib1g-dev libbz2-dev liblzma-dev`
2. install dorado binary from https://github.com/nanoporetech/dorado
3. install duplex_tools in a venv from https://github.com/nanoporetech/duplex-tools
