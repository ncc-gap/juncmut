# juncmut

The function for analyzing splicing junction associated variant (generated by STAR or bamTojunction, .SJ.out.tab files )

## Dependency

### Binary programs

[samtools](http://www.htslib.org/), [bedtools](https://github.com/arq5x/bedtools2)

### Python

[pysam](https://github.com/pysam-developers/pysam), 
[junc_utils](https://github.com/friend1ws/junc_utils), 
[annot_utils](https://github.com/friend1ws/annot_utils), 
[edlib](https://github.com/Martinsos/edlib)

## Install

```
pip install git+https://github.com/ncc-gap/juncmut.git
```

## Commands
### get

Requirment of RNA bam
```
usage: juncmut get [-h]
 -input_file input_file
 -output_file output_file
 -reference reference
 -rna_bam rna_bam 
 -genecode_gene_file genecode_gene_file
 -output_bam output_bam
 [-control_file [CONTROL_FILE [CONTROL_FILE ...]]] 
 [-read_num_thres READ_NUM_THRES]
 [-freq_thres FREQ_THRES] 
 [-mut_num_thres MUT_NUM_THRES]
 [-mut_freq_thres MUT_FREQ_THRES]
 [-debug debug]

optional arguments:
  -h, --help            show this help message and exit
  -input_file input_file
                        Input file(such as sample.SJ.out.tab file generated by STAR)
  -output_file output_file
                        Output file
  -reference reference  Path to reference genome
  -rna_bam rna_bam      Path to RNA bam file
  -genecode_gene_file genecode_gene_file
                        Path to genecode gene file
  -output_bam output_bam
                        Output bam
  -control_file [CONTROL_FILE [CONTROL_FILE ...]]
                        Path to control data created by merge_control (default: None)
  -read_num_thres READ_NUM_THRES
                        Splicing junctions with reads >= read_num_thres is saved (default: 3)
  -freq_thres FREQ_THRES
                        Splicing junctions with reads >= freq_thres is saved (default: 0.05)
  -mut_num_thres MUT_NUM_THRES
                        A mutation with mutation alleles >= mut_num_thres is a true candidate (default: 1)
  -mut_freq_thres MUT_FREQ_THRES
                        A mutation with frequency >= mut_freq_thres is a true candidate (default: 0.05)
  -debug debug          True keeps the intermediate files.
```

### annot
```
usage: juncmut annot [-h]
 -input_file input_file
 -output_file output_file
 -reference reference
 -gnomad gnomad
 [-debug debug]

optional arguments:
  -h, --help            show this help message and exit
  -input_file input_file
                        Input file(such as juncmut.txt generated by juncmut)
  -output_file output_file
                        Output file
  -reference reference  Path to reference genome
  -gnomad gnomad        Path to gnomad vcf file
  -debug debug          True keeps the intermediate files.
```

## Tutorial

See this [page](https://github.com/ncc-gap/juncmut/wiki/Tutorial)

