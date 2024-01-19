# FASTAhandler
Python scripts designed for efficiently managing various types of FASTA file formats


## Brief Background
**PoreQC** is a [Nextflow](https://github.com/nextflow-io/nextflow) pipeline for Oxford nanopore reads ([Slow5](https://github.com/hasindu2008/slow5tools), [Pod5](https://github.com/nanoporetech/pod5-file-format) and [Fastq](https://en.wikipedia.org/wiki/FASTQ_format)). Integrating with [Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview), [Dorado](https://github.com/nanoporetech/dorado), [Buttery-eel](https://github.com/Psy-Fer/buttery-eel), [Cutadapt](https://github.com/marcelm/cutadapt), and [Sequali](https://github.com/rhpvorderman/sequali), the automated pipeline can work for basecalling, quality control and removal adapters. We (Hyungtaek Jung and the [National Centre for Indigenous Genomics](https://ncig.anu.edu.au/) at [The Australian National University](https://www.anu.edu.au/), Australia) initially started this project to provide comprehensive data management at the [National Computational Infrastructure](https://nci.org.au/) for biologists. As a command-line interface (CLI) application, we have tested it for ONT long-read data focusing on whole genome shotgun datasets that can be widely used by the greater research community. However, please note that basecalling and visualising a big dataset would require large computational resources on HPC or Cloud. 


## Citation
**Hyungtaek Jung**, Kirat Alreja, Kosar Hooshmand, Hadi Nazem-Bokaee, Hasindu Gamaarachchi, **Hardip Patel**: **PoreQC**: An automated nextflow pipeline for oxford nanopore read basecalling, quality control and adapter removal, [PLoS Comp Biol Submitted](https://www.biorxiv.org/XXXX).


## Contents:
+ STABLE
+ INSTALLATION
+ LICENSE 
+ GETTING STARTED
+ FAQ
+ WIKI PAGE
+ AUTHOR
+ COPYRIGHT


## STABLE (version 0.0.XXX)
Release date: January 2024
**PoreQC** comprises two key features (basecalling and quality control) and four interactive steps with open-source programs (See LICENSE). 


## INSTALLATION
Please download the program from [this link](https://github.com/OZTaekOppa/PoreQC)
!!! Please note, that programs and dependencies can also be installed via Bioconda. For any other issues, we highly encourage users to use the [Issues](https://github.com/OZTaekOppa/PoreQC/issues).

    ~~~
    Create the virtual environment
    Need to be updated
    
    Get source
    Need to be updated
    
    Install packages
    Need to be updated
    
    Run
    Need to be updated
    ~~~

## License

**PoreQC** is provided under the MIT license and is based on other open-source software:

[Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview) for  basecalling and processing raw signal data from nanopore sequencing devices, providing accurate DNA sequence information.

[Dorado](https://github.com/nanoporetech/dorado) for Oxford Nanopore long-read sequencing data, offering enhanced accuracy in detecting structural variants and single nucleotide variants.

[Buttery-eel](https://github.com/Psy-Fer/buttery-eel) for a Slow5 file reader and basecalling wrapper for Guppy and Dorado.

[Cutadapt](https://github.com/marcelm/cutadapt) for removing adapters, primers, and other unwanted sequences from high-throughput sequencing data.

[Sequali](https://github.com/rhpvorderman/sequali) for evaluating the quality of sequencing data through the generation of comprehensive metrics and visualizations.

[Nextflow](https://github.com/nextflow-io/nextflow) for a data-driven computational workflow engine designed to facilitate scalable and reproducible scientific workflows.

[In-house Perl Script] for calculating basic statistics of a FASTQ file in-house.


### Tested Datasets
Reference genome(https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.3/#/st)
Oxford Nanopore reads(https://ngdc.cncb.ac.cn/gsa/browse/CRA004538) and (https://www.sciencedirect.com/science/article/pii/S1672022921001741)

## GETTING STARTED

**PoreQC**, integrated with Nextflow, has two specific features: a basecalling (Slow5) and a result summary and visualisation of quality control (Fastq). The data input/output enables end-to-end file selection. The result summary and visualisation are mainly designed to visualise the outcome for quality control. Please note that all required input files (e.g. Slow5) must be prepared from [Slow5tools](https://github.com/hasindu2008/slow5tools) to have a seamless experience of **PoreQC**. However, users can use Fastq files for quick quality control. 

**FASTAhandler**, mainly written in Python 3.12+ and ??, has two specific features: a data input module and a result visualisation window. The data input module enables end-to-end file selection. The result visualisation window is mainly designed to visulise the outcome selected by a user. Please note that all required input files (e.g. vcf, gff, gtf, bam, and fasta) must be prepared from third party programs before running the iVPSV. And, all input files must be in your local drive on your desktops.


### multi2singleline
- To convert a multi-fasta (multiline) into a single-line fasta.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A multiline fasta file.
	+ Output: A single-line fasta.

Example usage
```
python3 multi2singleline.py --input-seq test_dna.fasta --out test_output_sl.fasta --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. multi2singleline.py:  Call multi2singleline module
	1. python3 multi2singleline.py --help: Check help menu
		+ --input-seq: Indicate an input multi-fasta (multiline) fasta file and its path
		+ --out: Indicate an output single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only only with Gb size)

### renameid
- To rename prefix IDs and headers from a single-line fasta.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A single-line fasta with a new prefix ID name.


Example usage
```
python3 renameid.py --input-seq test_dna.fasta --new-name FunNGS --out output_reN.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, please use multi2singleline module first. 

- Parameter explanation
	1. python 3: Call python 3
	1. renameid.py:  Call renameid module
	1. python3 renameid.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --new-name: Indicate a new prefix ID/header name (accept both integer and strings but no space)
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### idextract
- To extract matched IDs and their corresponding sequences.
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.

Example usage
```
python3 idextract.py --input-seq test_dna.fasta --input-header header_id.txt --out output_extracted.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 
	
- Parameter explanation
	1. python 3: Call python 3
	1. idextract.py:  Call idextract module
	1. python3 idextract.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --input-header: Indicate an input ID and header (without ">") text file and its path
		+ --out: Indicate an output ID matched and extracted single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### idextractlocation
- To extract matched IDs, locations and their corresponding sequences (focused on a single ID)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.

Example usage
```
python3 idextractlocation.py --input-seq test_dna.fasta --header-id test3_3%week --start 2 --end 10 --out output_test.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. idextractlocation.py:  Call idextractlocation module
	1. python3 idextractlocation.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --header-id: Indicate an input ID and header (without ">") name and pattern
		+ --start: Indicate a start position to extract (please use -1 value due to the python index)
		+ --end: Indicate en end position to extract (please use -1 value due to the python index)
		+ --out: Indicate an output ID matched and extracted single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)
	
### idextractlocamulti
- To extract matched IDs, locations and their corresponding sequences (focused on multiple IDs)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.

Example usage
```
python3 idextractlocamulti.py --input-seq test_dna.fasta --input-extract input_extract.txt --out output_extest.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. idextractlocamulti.py:  Call idextractlocamulti module
	1. python3 idextractlocamulti.py --help: Check help menu
	--input-seq: Indicate an input single-line fasta file and its path
	--input-extract: Indicate an input ID and header (without ">") text file (a tap-separated) and its path including start and end positions (please use -1 value due to the python index)
	--out: Indicate an output ID matched and extracted single-line fasta file and its path
	--t: Specify thread numbers (intergr only)
	--mem: Specify memory numbers (integer only with Gb size)

### revcomplement
- To make a reverse complement sequence
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A reverse complement single-line fasta.

Example usage
```
python3 revcomplement.py --input-seq test_dna.fasta --out output_revctest.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. revcomplement.py:  Call revcomplement module
	1. python3 revcomplement.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --out: Indicate a reverse complement converted output single-line fasta file and its path
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### findcountdupl
- To find and count the duplicated IDs and sequences from a multi-fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A text with a duplication number.

Example usage
```
python3 findcountdupl.py --input-seq test_dna2.fasta --out output_files.txt --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. findcountdupl.py:  Call findcountdupl module
	1. python3 findcountdupl.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --out: Indicate an output text file after finding and counting the duplicated IDs and sequences (only for both matched IDs and their corresponding sequences)
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### removedupl
- To remove the duplicated IDs and sequences from a multi-fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A single-line fasta after removing updated IDs and sequences.

Example usage
```
python3 removedupl.py --input-seq test_dna2.fasta --outfasta output_testdupl.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. removedupl.py:  Call removedupl module
	1. python3 removedupl.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --outfasta: Indicate an output fasta file after removing the duplicated IDs and sequences (only for both matched IDs and their corresponding sequences)
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### subsetfasta
- To make a subset of data with a sequence length filter option
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A subseted single-line fasta.

Example usage
```
python3 subsetfasta.py --input-seq test_mRNA1.fasta --filter 50 --out output_subset.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. subsetfasta.py:  Call subsetfasta module
	1. python3 subsetfasta.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --filter: Indicate a length size to filter out
		+ --out: Indicate an output fasta file after filtering the length size (intergr only)
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### extractpattern
- To make a subset of data with find, filter and extract options
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of patterns.
	+ Output: A extracted single-line fasta with IDs and their corresponding sequences.

Example usage
```
python3 extractpattern.py --input-seq test_dna1.fasta --input-pattern seq_pattern.txt --input-length 45 --out output_pattern.fasta --t 1 --mem 2
```
	+ If the input fasta file is not a single-line fasta, the embedded pipeline will automatically convert multiline fasta files into a single-line fasta. 

- Parameter explanation
	1. python 3: Call python 3
	1. extractpattern.py:  Call extractpattern module
	1. python3 extractpattern.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --input-pattern: Indicate an input text file and its path to find and filter a specific sequence pattern (accept reverse complement sequences)
		+ --input-length: Indicate a length size to filter out
		+ --out: Indicate an output fasta file after filtering the length size (intergr only)
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)


### concatenate (여기 Unlimited Input으로 변경. ASMstats와 비슷한 방법으로. 변경 했음. 확인 요망)
- To make a subset of data with find, filter and extract options
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multiple fasta files with the same prefix IDs.
	+ Output: A single-line fasta with concatenated IDs and their corresponding sequences.

Example usage
```
python concatenate.py --input-seq test_dna1.fasta test_dna2.fasta test_dna3.fasta --out output_concat.fasta --t 1 --mem 2
```
	+ To maximise this module, the multiple input fasta files must have the same prefix IDs and headers along with a single-line fasta before concatenating.
	+ If not, please use renameid and multi2singleline modules before using concatenate module.

- Parameter explanation
	1. python 3: Call python 3
	1. concatenate.py:  Call concatenate module
	1. python3 concatenate.py --help: Check help menu
		+ --input-seq: Indicate input single-line fasta files and their path
		+ --out: Indicate a concatenated output fasta file and its path
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)

 


## FAQ

We encourage users to use the [Issues](https://github.com/OZTaekOppa/FASTAhandler/issues).


## WIKI PAGE

Please see GitHub page.


## AUTHOR

**Hyungtaek Jung** and the [**National Centre for Indigenous Genomics**](https://ncig.anu.edu.au/).


## COPYRIGHT

The full **FASTAhandler** is distributed under the MIT license. 
