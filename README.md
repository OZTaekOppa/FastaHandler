# FASTAhandler
Python scripts designed for efficiently managing various types of FASTA file formats


## Brief Background
**FASTAhandler** is designed for analysing and manipulating FASTA data efficiently. With 14 work modules, it simplifies input/output indications, covers diverse aspects of FASTA data analysis, and supports post-processing, filtering, and format conversion. We (Hyungtaek Jung and the [National Centre for Indigenous Genomics](https://ncig.anu.edu.au/) at [The Australian National University](https://www.anu.edu.au/), Australia) initially started this project to provide comprehensive data management at the [National Computational Infrastructure](https://nci.org.au/) for biologists. As a command-line interface (CLI) application, we have tested it for various FASTA file formats focusing on life science datasets so that the greater research community can widely use it. However, please note that analysing and manipulating a big dataset would require large computational resources on HPC or Cloud. 


## Citation
Hyungtaek Jung, Kirat Alreja, Kosar Hooshmand, Hadi Nazem-Bokaee, Hardip Patel: **FASTAhandler**: An easy Python-based tool set for handling FASTA files, [PLoS Comp Biol Submitted](https://www.biorxiv.org/XXXX).


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
**FASTAhandler** is a standalone Python application with 14 modules for manipulating FASTA files via interactive steps with open-source programs (See LICENSE). 


## INSTALLATION
Please download the program from [this link](https://github.com/OZTaekOppa/FASTAhandler)
!!! Please note, that programs and dependencies can also be installed via Bioconda. For any other issues, we highly encourage users to use the [Issues](https://github.com/OZTaekOppa/FASTAhandler/issues).

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

**FASTAhandler** is provided under the MIT license and is based on other open-source software. Please see the manuscript for the full details including Python packages, modules and libraries integrated into **FASTAhandler** and their applications. 


### Tested Datasets
Please see the example dataset folder. 


## GETTING STARTED
**FASTAhandler**, mainly written in Python 3.12+ and ??, has 14 modules. The data input and output via CLI enables end-to-end file selection. Please note that all required input files (e.g. fasta and txt) must be prepared to have a seamless experience of **FASTAhandler**. 


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
