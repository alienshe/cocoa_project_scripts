This is a collection of scripts which are used to analyse nanopore amplicon sequence data. 



## Summary of scripts

**demultiplex_nested_barcodes.sh**- Demultiplexes reads that were multiplexed using "nested" barcodes. This means that an inner barcode was applied to the reads using tagged PCR primers, and then a second outer barcode was applied using the ONT Native barcoding kit.

**haplotyping_pipeline.sh** - Takes reads and sorts them into haplotypes based on variant SNPs and indels (depending on which output you choose)

**batch_run_haplotyping_pipeline.sh** - Uses an input .tsv file to run the haplotyping pipeline script over a large number of samples, also summarises the outputs.


## demultiplex_nested_barcodes.sh

#### Usage
```bash
./demultiplex_nested_barcodes.sh <barcodes.fa> <reads.fastq.gz> <output_directory> optional: <min_read_length> <max_read_length>
```

#### Inputs
- a single .fastq.gz file containing the raw reads
- .fasta file containing FULL COMBINED BARCODE sequences. Works best if full sequence is >40bp. Should contain the ONT barcode, adapter, and inner barcode (i also like to include the primer)
- Minimum and maximum read lengths. Any reads outside of this range will be filtered out. Defaults are 3000 and 9000 bp.

eg:
```
      >Tc1318
      CCAAACCCAACAACCTAGATAGGCCAGCACCTTAGGCGAAAAGAGATTGCCGGTCGTTGT
      >Tc1320
      CCAAACCCAACAACCTAGATAGGCCAGCACCTGCGAGAATTGACAAGTTGGCCAGTCGTT
```
#### What does it do?
1. Reads are aligned to barcodes.fa file with minimap2
2. Hits from minimap2 output that are from reads that are not between 3 and 9 kb, and hits that have barcodes in the middle, are removed.
3. The best hit for each read is kept
4. Based on this, reads are sorted into separate files, one for each barcode.
5. Porechop is used to trim barcodes/adapters/primers from the reads

     NOTE: for this step to work, any custom barcodes/primers must be added to porechop's adapters.py file

#### Outputs
A .fastq file is output for each barcode found in the reads. The reads trimmed with porechop can be found in the sorted_trimmed_fastq_files folder. The untrimmed reads are in the intermediate_files fodler.
```
├── intermediate_files
├── run.log
└── sorted_trimmed_fastq_files
    ├── trimmed_Tc1318.fastq
    ├── trimmed_Tc1320.fastq
    └── trimmed_unclassified.fastq
```


## haplotyping_pipeline.sh

#### Usage
```bash
./haplotyping_pipeline.sh <home/user/input_dir> <trimmed_reads.fastq> <reference.fasta> </home/user/output_dir> <gene_name>
```
#### Inputs

1. Input directory: This directory must contain any other input files (reference fasta +reads). MUST BE NON RELATIVE PATH FROM HOME DIR
2. Trimmed reads: path relative to input directory. .fastq file containing reads with primers/barcodes already trimmed with something like porechop.
3. Reference fasta. (note: must indexed and .fai file must be present in the same directory)
4. Output directory: MUST BE NON RELATIVE PATH FROM HOME DIR
5. Gene name: This is used to name output files neatly.

#### What does it do?

1. Reads are aligned against the reference fasta with minimap2
2. Clair3 identifies variant SNPs and indels using the alignment and the reference sequences. Outputs phased .vcf file.
3. Haplotyping: (sorting reads into haplotype groups based on variant alleles from .vcf) Done by both devider and whatshap. You can choose which output you want to use. In either case, a tagged .bam file will be output (tag: HP) which indicates which reads belong to which haplotypes. There will also be consensus sequences output for each haplotype.
   - devider: optimised for nanopore data, doesnt take indels into account, doesnt require you to define ploidy (so can work for mixed samples)
   - whatshap: less sophisticated, takes indels into account, requires defined ploidy (default: 2).

#### Outputs
```
.
├── alignments        #contains reads aligned against reference
├── clair3_output     #contains .vcf files identifying variant alleles
├── devider_output    #contains reads sorted into haplotypes and consensus sequences for each haplotype
├── log.txt
└── whatshap_output   #contains reads sorted into haplotypes and consensus sequences for each haplotype
```

#### Example

##### Example input directory structure
```
├── inputs
    ├── trimmed_reads
    │   ├── trimmed_Tc1318.fastq
    │   ├── trimmed_Tc1320.fastq
    │   ├── trimmed_Tc2391.fastq
    │   └── trimmed_unclassified.fastq
    └── reference_files
        ├── Tc1318_ref_chr8.fasta
        ├── Tc1318_ref_chr8.fasta.fai
        ├── Tc1320_ref_chr8.fasta
        ├── Tc1320_ref_chr8.fasta.fai
        ├── Tc2391_ref_chr3.fasta
        └── Tc2391_ref_chr3.fasta.fai
```
##### Example command
```bash
./haplotyping_pipeline.sh home/user/inputs trimmed_reads/trimmed_Tc1318.fastq reference_files/Tc1318_ref_chr8.fasta /home/user/output_dir Tc1318
```

##### Example output directory structure
```
├── output_dir
│   ├── alignments
│   │   ├── Tc1318_alignment.bam
│   │   ├── Tc1318_alignment.sam
│   │   ├── Tc1318_alignment_sorted.bam
│   │   ├── Tc1318_alignment_sorted.bam.bai
│   │   ├── Tc1318_alignment_sorted.bam.tagged.bam
│   │   └── Tc1318_alignment_sorted.bam.tagged.bam.bai
│   ├── clair3_output
│   │   ├── full_alignment.vcf.gz
│   │   ├── full_alignment.vcf.gz.tbi
│   │   ├── log
│   │   ├── merge_output.vcf.gz
│   │   ├── merge_output.vcf.gz.tbi
│   │   ├── phased_merge_output.vcf.gz
│   │   ├── phased_merge_output.vcf.gz.tbi
│   │   ├── pileup.vcf.gz
│   │   ├── pileup.vcf.gz.tbi
│   │   ├── run_clair3.log
│   │   └── tmp
│   ├── devider_output
│   │   ├── hap_info.txt
│   │   ├── ids.txt
│   │   ├── intermediate
│   │   ├── majority_vote_haplotypes.fasta
│   │   └── snp_haplotypes.fasta
│   ├── log.txt
│   └── whatshap_output
│   │   ├── split_alignments
│   │   ├── Tc1318_alignment_haplotagged.bam
│   │   ├── Tc1318_haplotype_0_consensus.fasta
│   │   ├── Tc1318_haplotype_1_consensus.fasta
│   │   └── Tc1318_haplotype_unassigned_consensus.fasta
```

## batch_run_haplotyping_pipeline.sh

#### Usage
```bash
./batch_run_devider_pipeline.sh <file_list.tsv> <home/user/input_dir> <home/user/output_dir> <query_dir> <reference_dir> <minimum_abundance>
```

#### Inputs
1. file list tsv: .tsv file which tells the script the reads and corresponding reference files you want to input.
         Column 1 should be the name of the .fastq file containing your trimmed reads (in the query directory)
         Column 2 should be the name of the reference .fasta file (in the reference directory)
         Column 3 is the name of the sample (used for neatly naming output files)
   Example:
```
trimmed_Tc3425.fastq	Tc3425_ref_chr9.fasta	Tc3425
trimmed_Tc3078.fastq	Tc3078_ref_chr9.fasta	Tc3078
trimmed_Tc1320.fastq	Tc1320_ref_chr8.fasta	Tc1320
trimmed_Tc2391.fastq	Tc2391_ref_chr3.fasta	Tc2391
trimmed_Tc1318.fastq	Tc1318_ref_chr8.fasta	Tc1318
```
2. Input directory: This directory must contain any other input files (reference fastas + reads). MUST BE NON RELATIVE PATH FROM HOME DIR
3. Output directory: MUST BE NON RELATIVE PATH FROM HOME DIR
4. Query directory: directory containing all of your trimmed .fastq files containing your reads. You should have 1 file per sample (path relative to input directory)
5. Reference directory: directory containing all reference .fasta files (path relative to input directory)
6. Minimum abundance: integer from 0-100. passed to devider. Defines the minimum abundance required for each devider haplotype.

##### Example 
###### Command
```bash
~batch_run_haplotyping_pipeline.sh \
      ~/inputs/reference_files/NLRs.tsv \   #tsv which tell the script what files to process
      /home/sylvie/inputs \                 #objective path to input directory
      /home/sylvie/NLR_haplotyping_output \ #objective path to output directory
      trimmed_reads \                       #directory (relative to input dir) containing .fastq files trimmed with porechop
      reference_files \                     #directory (relative to input dir) containing reference .fasta files
      0                                     #minimum abundance for devider haplotypes
```
```
└── inputs
    ├── trimmed_reads
    │   ├── trimmed_Tc1318.fastq
    │   ├── trimmed_Tc1320.fastq
    │   ├── trimmed_Tc2391.fastq
    │   ├── trimmed_Tc3078.fastq
    │   ├── trimmed_Tc3425.fastq
    │   ├── trimmed_Tc3618.fastq
    └── reference_files
        ├── NLRs.tsv
        ├── Tc1318_ref_chr8.fasta
        ├── Tc1318_ref_chr8.fasta.fai
        ├── Tc1320_ref_chr8.fasta
        ├── Tc1320_ref_chr8.fasta.fai
        ├── Tc2391_ref_chr3.fasta
        ├── Tc2391_ref_chr3.fasta.fai
        ├── Tc3078_ref_chr9.fasta
        ├── Tc3078_ref_chr9.fasta.fai
        ├── Tc3425_ref_chr9.fasta
        ├── Tc3425_ref_chr9.fasta.fai
        ├── Tc3618_ref_chr9.fasta
        ├── Tc3618_ref_chr9.fasta.fai

```


#### What does it do?

1. Iterates through file list tsv and runs the haplotyping_pipeline.sh script on the files from each row.
2. Summarises devider output, telling you how many haplotypes were found etc.
3. Collects devider haplotagged alignments from all samples into one directory
   
#### Outputs

```
.
├── collected_alignments
├── haplotype_summary.txt
├── Tc1318
│   ├── clair3_output
│   ├── devider_output
│   ├── log.txt
│   └── whatshap_output
├── Tc1320
│   ├── alignments
│   ├── clair3_output
│   ├── devider_output
│   ├── log.txt
│   └── whatshap_output
├── Tc2391
│   ├── alignments
│   ├── clair3_output
│   ├── devider_output
│   ├── log.txt
│   └── whatshap_output
├── Tc3078
│   ├── alignments
│   ├── clair3_output
│   ├── devider_output
│   ├── log.txt
│   └── whatshap_output
└── Tc3425
    ├── alignments
    ├── clair3_output
    ├── devider_output
    ├── log.txt
    └── whatshap_output

```
