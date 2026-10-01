
# DATA SOURCE

## blast-2.15.0+
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz

## Escherichia coli str. K-12 substr. MG1655 reference genome
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

## UniProtKB/Swiss-Prot 2025_02 of UniProtKB, published on Wed Apr 23 2025 (573,230 sequence entries)
    wget https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

--

# SHUFFLE-BASED STATISTICAL SIGNIFICANCE EVALUATION

## Computer-generated amino acid sequences with the average size and composition displayed by UniProtKB/Swiss-Prot sequences
    perl shuffle.pl uniprot_sprot.fasta

    #!/usr/bin/perl

    # SUMMARY
    # Takes a multi-FASTA file and creates "n" artificial sequences
    # with the average size and amino acid composition representing
    # the whole set of sequences in this file ("n" = number of seqs
    # in the multi-FASTA file)

    use strict;
    use List::Util qw (shuffle);

    my $infile = shift;

    unless (-f $infile) {
    # script usage
    print <<EOF;

    Usage: shuffle.pl <Multi-FASTA file>

    EOF
    exit;
    }

    my (%MULTI_FASTA,$id);
    my $c = 0;
    my $number_of_seqs  = 0;
    my $total_seq_size  = 0;
    my @ARTIFICIAL_SEQ  = ();
    my %AA              = ();
    my %AVERAGE_AA_FREQ = ();
    my %NUMBER_OF_AA    = ();

    open (READ, "$infile") || die "cannot open $infile: $!";
    while (<READ>){
    	if (/^(>.+)$/ || /^\n/){
    		$c++;
    		$number_of_seqs++;
    		if ($c > 1) {
    			foreach my $key (keys %MULTI_FASTA) {
    				# Counting aminoacids
    				for (my $pos = 0; $pos < length ($MULTI_FASTA{$key}); $pos ++) {
    					my $aa = substr ($MULTI_FASTA{$key}, $pos, 1);
    					$AA{$aa}++;
    				}
    				$total_seq_size = $total_seq_size + length ($MULTI_FASTA{$key});
    			}
    			%MULTI_FASTA = ();
    		}
    		$id = $1;
    	}
    	elsif (/^([A-Z]+)$/i){
    		$MULTI_FASTA{$id} .= $1;
    	}
    }
    close(READ);

    $number_of_seqs = $number_of_seqs -1;
    print "\nNUMBER OF SEQUENCES: $number_of_seqs\n";

    # Calculating the average sequence length for the whole set of sequences
    my $average_seq_size = $total_seq_size / $number_of_seqs;
    #print $average_seq_size," ",$number_of_seqs,"\n";

    # Calculate the average frequence of each aminoacid and the number of times
    # each aminoacid must be included in the artificial sequence so that it
    # holds the average aminoacid frequence calculated earlier
    foreach my $key (keys %AA) {
    	$AVERAGE_AA_FREQ{$key} = $AA{$key} / $total_seq_size;
    	$NUMBER_OF_AA{$key} = sprintf("%.f", ($AVERAGE_AA_FREQ{$key} * $average_seq_size));
    }
    print "\nNUMBER OF AMINO ACIDS: $total_seq_size\n\n";
    foreach my $key (keys (%AVERAGE_AA_FREQ)) {print "$key -> $AVERAGE_AA_FREQ{$key}\n";}
    print "\nAVERAGE SEQUENCED SIZE: $average_seq_size\n\n";
    foreach my $key (keys (%NUMBER_OF_AA)) {print "$key -> $NUMBER_OF_AA{$key}\n";}

    # Create an artificial sequence with the average sequence size and
    # the average aminoacid frequence calculated before
    foreach my $key (keys %AA) {
    	for (1 .. $NUMBER_OF_AA{$key}) {
    		push (@ARTIFICIAL_SEQ, $key);
    	}
    }
    print "\nARTIFICIAL SEQUENCE\n\n",join ("", @ARTIFICIAL_SEQ),"\n\n";

    # Shuffle the artificial sequence and print it
    open (WRITE, ">$infile.shuffle");
    for (1 .. 10000) { # $number_of_seqs
    	my @RANDOM_SEQ = shuffle (@ARTIFICIAL_SEQ);
    	print WRITE ">shuffled_$_\n" , join("",@RANDOM_SEQ),"\n";
    }
    close (WRITE);

## Pairwise sequence comparison between UniProtKB/Swiss-Prot real sequences against 10,000 Swiss-Prot artificial protein sequences (NULL MODEL)
    ncbi-blast-2.15.0+/bin/makeblastdb -in uniprot_sprot.fasta.shuffle -dbtype prot -out uniprot_sprot_shuffle
    nohup ncbi-blast-2.15.0+/bin/blastp -num_threads 4 -word_size 3 -gapopen 11 -gapextend 1 -matrix BLOSUM62 -threshold 13 -comp_based_stats 2 -seg yes -soft_masking true -lcase_masking -evalue 10 -max_target_seqs 1000000 -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qcovs qcovhsp" -query uniprot_sprot.fasta -db uniprot_sprot_shuffle -out uniprot_sprot_shuffle.blastp.out &

## Parse results searching for the most statistically significant sequence similarity and the highest BLAST score between artificial x real sequences (BIOLOGICAL SIGNIFICANCE TRESHOLD)
    awk '{print $13,$11}' uniprot_sprot_shuffle.blastp.out | LC_ALL=C sort -gr -k 2 | tail -n 1

--

# GENOME-WIDE HOMOLOGY-BASED SEARCHING

    ncbi-blast-2.15.0+/bin/makeblastdb -in GCF_000005845.2_ASM584v2_genomic.fna -dbtype nucl -out ecolik12db
    nohup ncbi-blast-2.15.0+/bin/tblastn -num_threads 4 -word_size 3 -gapopen 11 -gapextend 1 -matrix BLOSUM62 -threshold 13 -comp_based_stats 2 -seg yes -soft_masking true -lcase_masking -evalue 10 -max_target_seqs 1000000 -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sframe sstrand qcovs qcovhsp" -query uniprot_sprot.fasta -db ecolik12db -out genome.tblastn.out &

--

# NEXT STEPS

* For each alignment obtained in the previous step, re-align the pair of sequences with SSEARCH version 36 <https://fasta.bioch.virginia.edu/wrpearson/fasta/CURRENT/fasta36-linux64.tar.gz> to calculate the best rigorous (dynamic programming) alignment between them, retaining statistically and biologically significant pairwise alignments for subsequent analyzes. To estimate the statistical significance of each pairwise alignment in this step, set SSEARCH to shuffle the functional sequence of the pair randomly, then align this artificially created sequence with the similar genomic sequence of the pair, and repeat this process 500 times, calculating the E-value for the original pair of sequences based on an effective (apparent) database size of 10,000 sequences and the distribution of alignment scores obtained by these random sequences: ssearch36 -b 1 -d 1 -k 500 -f -11 -g 1 -m 10 -q -s BL62 -z 11 -Z 10000 (-b: number of best scores to show, -d: number of best alignments to show, -k: specify the number of shuffles for statistical parameter estimation, -f: penalty for opening a gap, -g: penalty for additional residues in a gap, -m: alignment display options, -q: quiet option, -s: specify substitution matrix, -z: specify statistical calculation, -Z: set the apparent database size used for expectation value calculations)

* Merge adjacent and overlapping UniProtKB/Swiss-Prot high scoring pairs (HSPs) in each mapped genomic regions with BEDTools, obtaining the coordinates of all non-overlapping merged segments displaying statistically and biologically significant similarity to UniProtKB/Swiss-Prot real proteins: bedtools merge -c 2,4,4,4,5 – count, mean, min, max, collapse

* Parse the SSEARCH alignment outputs previously obtained for each of these merged regions, calculating the score density produced by each UniProtKB/Swiss-Prot protein hitting one or more DNA segments inside a merged region as the summation of the scores obtained by their high scoring pairs (HSPs), divided by the summation of their sequences' length, obtaining normalized scores

* Select the UniProtKB/Swiss-Prot protein sequence producing the highest score density as the unique parent 'functional' protein, simultaneously reconstructing and annotating these genomic regions

* Next, having the HSPs of each reconstructed/annotated genomic region aligned with their parent 'functional' protein, scan the pairwise alignment using the UniProtKB/Swiss-Prot protein as a reference, searching for pseudogene-like features in the genomic sequence, that is, disablements rendered by non-synonymous substitutions, in-frame insertions/deletions, frameshifts, loss of initiation and/or termination codons, and/or internal termination codons, calculating the total number of mutations as the summation of all disablements

