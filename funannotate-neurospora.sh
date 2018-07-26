#!/bin/bash
#SBATCH --ntasks 32 --mem 32G --time 48:00:00 --nodes 1 --out neurospora_predict.%A.log -J funNeuroPred

# rewrite as a makefile/snakemake?

module load funannotate/git-live
module load python/2.7.12

#shell script to run funannotate on Neurospora crassa
CPUS=32
HOMEDIR=`pwd`
SED=sed
m=$(which gsed 2> /dev/null)
if [ $m ]; then
    SEd=$m
fi

mkdir run_Neurospora38
pushd run_Neurospora38

#download Neurospora reference from FungiDB
if [ ! -f FungiDB-38_NcrassaOR74A_Genome.fasta ]; then
    wget http://fungidb.org/common/downloads/Current_Release/NcrassaOR74A/fasta/data/FungiDB-38_NcrassaOR74A_Genome.fasta
fi
if [ ! -f FungiDB-38_NcrassaOR74A.gff ]; then
    wget http://fungidb.org/common/downloads/Current_Release/NcrassaOR74A/gff/data/FungiDB-38_NcrassaOR74A.gff
fi

#clean up the genome and annotation from FungiDB
#the data contains both a rRNA contig and a mitochondrial contig, neither can be annotated
#with funannotate, so will remove them and then clean up rest of fasta header
if [ ! -f Ncrassa.genome.fa ]; then
    echo -e 'KC683708\nKI440765\nKI440766\nKI440767\nKI440768\nKI440769\nKI440770\nKI440771\nKI440772\nKI440773\nKI440774\nKI440775\nKI440776\nKI440777\n' > remove.list
    sed 's/ | .*$//g' FungiDB-38_NcrassaOR74A_Genome.fasta > FungiDB-38_NcrassaOR74A_Genome_edit.fasta
    python $HOMEDIR/fasta_remove.py FungiDB-38_NcrassaOR74A_Genome_edit.fasta remove.list > Ncrassa.genome.fa
fi

#clean up gff file to remove annotations from mito and rRNA
if [ ! -f Ncrassa.clean.gff3 ]; then
    grep -E -v '^(KC683708|KI440765|KI440766|KI440767|KI440768|KI440769|KI440770|KI440771|KI440772|KI440773|KI440774|KI440775|KI440776|KI440777)' FungiDB-38_NcrassaOR74A.gff > Ncrassa.clean.gff3
fi
echo -e "Now annotate the genome 4 ways 
	\t1) funannotate busco mediated training
	\t2) funannotate RNA-seq data
	\t3) maker base
	\t4) maker with RNA-seq gene models
	-------------------------------------------"

#first run repeatmodeler/repeat masker
if [ ! -f Ncrassa.masked.fa ]; then
    funannotate mask -i Ncrassa.genome.fa --cpus $CPUS -o Ncrassa.masked.fa
fi
REPEATLIB=$(find . -name "repeatmodeler-library*" | sed 's,./,,g')

 
#get UniProt/SwissProt proteins, remove those from Ncrassa as to not bias predictions
if [ ! -f neurospora_crassa.list ]; then
    grep 'OS=Neurospora crassa' $FUNANNOTATE_DB/uniprot_sprot.fasta | $SED 's/>//g' > neurospora_crassa.list
    python $HOMEDIR/fasta_remove.py $FUNANNOTATE_DB/uniprot_sprot.fasta neurospora_crassa.list > uniprot_sprot.no-neurospora.fasta
fi

#map protein evidence to genome so can be re-used for the different funannotate runs (will save some compute time)
if [ ! -f ncrassa_uniprot.alignments.gff3 ]; then
    funannotate util prot2genome -g Ncrassa.masked.fa -p uniprot_sprot.no-neurospora.fasta -o ncrassa_uniprot.alignments.gff3 --cpus $CPUS
fi

#BUSCO only
funannotate predict -i Ncrassa.masked.fa \
					-o busco_only \
					-s "Neurospora crassa" \
					--isolate OR74A_busco \
					--busco_seed_species neurospora_crassa \
					--protein_alignments ncrassa_uniprot.alignments.gff3 \
					--name BUSCO_ \
					--cpus $CPUS

#now run with RNA-seq data. 
#Download data from SRA as 9 singled ended stranded Illumina from different growth conditions

# do we want to do strand-specific paired end too at some point? --js
echo "Now downloading RNA-seq data from NCBI SRA"
for file in SRR957218 SRR5591013 SRR5591014 SRR5591015 SRR5591016 SRR5591017 SRR5591018 SRR5591019 SRR5591020; do
    if [ ! -f $file\_1.fastq.gz ];
	fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' --split-files $file
#	pigz $file\_1.fastq
    fi
done

#data sets from different growth conditions
funannotate train -i Ncrassa.masked.fa \
				  -o rna_seq \
				  -s *.fastq.gz \
				  --stranded R \
				  --jaccard_clip \
				  --species "Neurospora crassa" \
				  --isolate OR74A_rna \
				  --cpus $CPUS

# now run prediction step using the data generated from training
# now train with PASA models
funannotate predict -i Ncrassa.masked.fa \
    -o rna_seq \
    -s "Neurospora crassa" \
    --isolate OR74A_pasa \
    --name PASA_ \
    --rna_bam rna_seq/training/funannotate_train.coordSorted.bam \
    --pasa_gff rna_seq/training/funannotate_train.pasa.gff3 \
    --transcript_evidence rna_seq/training/funannotate_train.trinity-GG.fasta \
    --transcript_alignments rna_seq/training/funannotate_train.transcripts.gff3 \
    --protein_alignments ncrassa_uniprot.alignments.gff3 \
    --cpus $CPUS

#recover isoforms and update gene models using PASA/Kallisto
funannotate update -i rna_seq --cpus $CPUS

#want to compare this to a Maker run
exit

modue load maker

#setup MAKER run
maker -CTL
MAKER_EXE="$HOMEDIR/maker_exe.ctl"
MAKER_BOPTS="$HOMEDIR/maker_bopts.ctl"
MAKER_OPTS="$HOMEDIR/maker_opts.ctl"
sed "s/rmlib= #/rmlib=$REPEATLIB #/g" maker_base_opts.ctl > $MAKER_OPTS
maker -base maker_base $MAKER_OPTS $MAKER_BOPTS $MAKER_EXE
gff3_merge -d maker_base.maker.output/maker_base_master_datastore_index.log -s -g -n > maker_base.gff3
maker_map_ids --prefix MBASE_ --justify 6 --suffix -T --iterate 1 maker_base.gff3 > maker_base.map.ids
map_gff_ids maker_base.map.ids maker_base.gff3

#now use RNA-seq PASA data to run Maker using those data
#use PASA high quality gene models to train SNAP
mkdir snap_train
cd snap_train
maker2zff rna_seq/training/funannotate_train.pasa.gff3
#since PASA GFF doesn't have fasta file appended like Maker, just copy of genome file
cp ../Ncrassa.masked.fa genome.dna
#now train SNAP
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Ncrassa . > ../Ncrassa.snap.hmm
cd $HOMEDIR

sed "s/rmlib= #/rmlib=$REPEATLIB #/g" maker_rna_opts.ctl > $MAKER_OPTS
maker -base maker_rna $MAKER_OPTS $MAKER_BOPTS $MAKER_EXE
gff3_merge -d maker_rna.maker.output/maker_rna_master_datastore_index.log -s -g -n > maker_rna.gff3
maker_map_ids --prefix MBASE_ --justify 6 --suffix -T --iterate 1 maker_rna.gff3 > maker_rna.map.ids
map_gff_ids maker_rna.map.ids maker_rna.gff3


#now run annotation comparison

ln -s $HOMEDIR/rna_seq/update_results/Neurospora_crassa_OR74A_rna.gff3 $HOMEDIR/fun_rna-seq.gff3
ln -s $HOMEDIR/busco_only/predict_results/Neurospora_crassa_OR74A_busco.gff3 $HOMEDIR/fun_busco.gff3

funannotate util compare -f Ncrassa.masked.fa -r Ncrassa.clean.gff3 -q maker_base.gff3 maker_rna.gff3 fun_busco.gff3 fun_rna-seq.gff3 -o compare -c

popd
