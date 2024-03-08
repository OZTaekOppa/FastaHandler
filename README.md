# FastaHandler
A collection of Python scripts designed for the efficient management of various FASTA file formats.

![FASTAhandler Logo](https://github.com/OZTaekOppa/FastaHandler/blob/main/images/FASTAhandler_Logo.png)


## Brief Background
**FastaHandler**, created by Hyungtaek Jung and the team at the [National Centre for Indigenous Genomics](https://ncig.anu.edu.au/) at [The Australian National University](https://www.anu.edu.au/), is a Python script suite for efficient FASTA file management. It boasts 19 work modules to ease input/output processes, covering various aspects of FASTA data analysis including post-processing and format conversion. Optimized for life science datasets, **FastaHandler** is a CLI application tested across different FASTA formats. Users should note that processing large datasets may require substantial computational resources on Linux, HPC or Cloud platforms.


## Citation
Hyungtaek Jung et al. 2024: **FastaHandler**: An easy Python-based toolset for handling fasta files, [PLoS Comp Biol Submitted](https://www.biorxiv.org/XXXX).


## Contents:
+ [STABLE](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#stable-version-101)
+ [INSTALLATION](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#installation)
+ [LICENSE](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#license) 
+ [GETTING STARTED](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#getting-started)
+ [FAQ](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#faq)
+ [WIKI PAGE](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#wiki-page)
+ [AUTHORS](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#authors)
+ [COPYRIGHT](https://github.com/OZTaekOppa/FASTAhandler/blob/main/README.md#copyright)


## STABLE (version 1.0.1)
- Release date: January 2024
- **FastaHandler** is a standalone Python application equipped with 18 sub-modules for interactive FASTA file manipulation, available as open-source (see [LICENSE](https://github.com/OZTaekOppa/FastaHandler/blob/main/README.md#license)).


## INSTALLATION
- Access the program via [GitHub](https://github.com/OZTaekOppa/FASTAhandler)
- Installation options include Bioconda or Python pip. For support, refer to [Issues on GitHub](https://github.com/OZTaekOppa/FastaHandler/issues).
- **FastaHandler** requires no separate installation process.
- Just clone this repository, and run
```
git clone https://github.com/OZTaekOppa/FASTAhandler/
python3 {path}/fastahandler.py
```

## LICENSE
**FastaHandler** is available under the MIT license and incorporates various open-source software. For detailed information on the integrated Python packages, modules, and libraries, and their specific applications within **FastaHandler**, please refer to the [manuscript](https://www.biorxiv.org/XXXX)


### Tested Datasets
Please refer to the example dataset folder for sample data and usage demonstrations.

## GETTING STARTED
**FastaHandler** is developed primarily in Python 3.9+ and Biopython and features 19 modules. It facilitates data input and output through a Command-Line Interface (CLI), ensuring smooth end-to-end file handling. To optimize the use of **FastaHandler**, users should prepare all necessary input files, such as FASTA and TXT formats, in advance.

![FASTAhandler Workflow](https://github.com/OZTaekOppa/FastaHandler/blob/main/images/FASTAhandler_Workflow.png)

### General Usage
```
FastaHandler: Fasta File Manipulation Toolkit
version 1.0.1

Usage: python3 fastahandler.py <module> <parameters>

Modules:
Multi2Single    		| m2s   	Convert a multi-fasta (multiline) into a single-line fasta.
Gfa2Fasta	    		| g2a   	Convert a gfa (Graphical Fragment Assembly) into a single-line fasta.
RenameId        		| rid   	Rename prefix IDs and headers.
PrefixRename    		| prn   	Rename prefix IDs and headers with a user’s input.
PrefixSelectRename      	| psr   	Rename prefix IDs and headers with a user’s input (Only).
PrefixFindReplaceRename      	| pfr   	Replace and rename prefix IDs and headers with a user’s input (Only).
IdExtract       		| idx   	Extract matched IDs and their corresponding sequences.
IdExtractLocation       	| iel   	Extract matched IDs, locations and their corresponding sequences.
IdExtractLocationMultiple       | iem   	Extract matched IDs, locations and their corresponding sequences (Multiple).
ReverseComplement       	| rcp   	Make a reverse complement sequence.
FindCountDuplication    	| fcd   	Find and count the duplicated IDs and sequences.
RemoveDuplication       	| rvd   	Remove the duplicated IDs and sequences.
SubsetFasta     		| ssf   	Make a subset of data with a sequence length filter.
ExtractPattern  		| xpt   	Make a subset of data with find, filter and extract.
EachFastaStats  		| efs   	Generate each line fasta statistic for a multi-line fasta.
AllFastaStats   		| afs   	Generate a summary of multi-line fasta statistics.
MultipleFastaStats      	| mfs   	Generate a summary of multi-line fasta statistics (Multiple).
ConcatenateFasta        	| ccf   	Make a concatenated fasta file (Multiple).
TranslateSequence       	| tls   	Find the translated sequences as a protein and open reading frames (ORFs).

Use <module> --help for module usage.
```

### multi2single (m2s)
- To convert a multi-fasta (multiline) into a single-line fasta.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A multiline fasta file.
	+ Output: A single-line fasta.
	+ Example file: [mltseq2sl](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/mltseq2sl.fa) in the "example_data" folder. 

Example usage
```
python3 multi2single.py --input-seq test_dna.fasta --out test_output_sl.fasta --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. multi2single.py:  Call multi2single module
	1. python3 multi2singleline.py --help: Check help menu
		+ --input-seq: Indicate an input multi-fasta (multiline) fasta file and its path
		+ --out: Indicate an output single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only only with Gb size)

### gfa2fa (g2a)
- To convert a gfa (Graphical Fragment Assembly) into a single-line fasta.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A gfa file.
	+ Output: A single-line fasta.
	+ Example file: [gfa2fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/gfa2fa.gfa) in the "example_data" folder. 

Example usage
```
python3 gfa2fa.py --input-gfa test_dna.gfa --out test_output_sl.fasta --t 1 --mem 2
```

- Parameter explanation
	1. python 3: Call python 3
	1. gfa2fa.py:  Call gfa2fa module
	1. python3 gfa2fa.py --help: Check help menu
		+ --input-gfa: Indicate an input gfa file and its path
		+ --out: Indicate an output single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only only with Gb size)
  
### renameid (rid)
- To rename prefix IDs and headers from a single-line fasta.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A single-line fasta with a new prefix ID name.
 	+ Example file: [renameid_seq.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/renameid_seq.fa) in the "example_data" folder. 


Example usage
```
python3 renameid.py --input-seq test_dna.fasta --new-name FunNGS --out output_reN.fasta --t 1 --mem 2
```
+ Use the multi2single module first if your input FASTA file isn't in single-line format.

- Parameter explanation
	1. python 3: Call python 3
	1. renameid.py:  Call renameid module
	1. python3 renameid.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --new-name: Indicate a new prefix ID/header name (accept both integer and strings but no space)
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### prfxrename (prn)
- To rename prefix IDs and headers from a single-line fasta with a user's input text file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a tap-separated text file (e.g. old_IDs	new_IDs)
	+ Output: A single-line fasta with a new prefix ID name based on the user's input text file. Unmatched IDs will be also produced with its original IDs.
	+ Example file: [header_id_only.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/header_id_only.txt) in the "example_data" folder.

Example usage
```
python3 prfxrename.py --input-seq test_dna.fasta --input-id new_ids.txt --out output_reN.fasta --t 1 --mem 2
```
+ Although the script can process multiline FASTA files, if you're unsure about your input FASTA file format, it's recommended to first convert it to single-line format using the multi2single module.
+ The script is ideal for renaming specific parts of IDs/headers in FASTA files, such as those from genome assemblies acquired from public databases like NCBI and EBI, based on user input.

- Parameter explanation
	1. python 3: Call python 3
	1. prfxrename.py:  Call prfxrename module
	1. python3 renameid.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --input-id: Indicate a new tap-separated prefix ID/header name file (accept both integer and strings but no space)
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### prfxselrename (psr)
- To rename prefix IDs and headers from a single-line fasta with a user's input text file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a tap-separated text file (e.g. old_IDs	new_IDs)
	+ Output: A single-line fasta with a new prefix ID name based on the user's input text file. Unmatched IDs will be discarded.
 	+ Example file: [header_id_only.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/header_id_only.txt) in the "example_data" folder.

Example usage
```
python3 prfxselrename.py --input-seq test_dna.fasta --input-id new_ids.txt --out output_reN.fasta --t 1 --mem 2
```
+ Although the script can process multiline FASTA files, if you're unsure about your input FASTA file format, it's recommended to first convert it to single-line format using the multi2single module.
+ The script is effective for selectively renaming IDs/headers in FASTA files (with a user's specific input), such as genome assemblies from databases, where you want to exclude and not rename certain elements like scaffolds.

- Parameter explanation
	1. python 3: Call python 3
	1. prfxselrename.py:  Call prfxselrename module
	1. python3 prfxselrename.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --input-id: Indicate a new tap-separated prefix ID/header name file (accept both integer and strings but no space)
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### prfxfindreplace (pfr)
- To rename prefix IDs and headers from a single-line fasta with a user's input text file.
	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a find and replace pattern (e.g. old_IDs	new_IDs)
	+ Output: A single-line fasta with a new prefix ID name based on the user's input text file. Worked as "sed" command.
 	+ Example file: [header_id_only.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/header_id_only.txt) in the "example_data" folder.

Example usage
```
python3 prfxfindreplace.py --input-seq --find-ptrn hifiasm --replace Assembly --out output_reN.fasta --t 1 --mem 2
```
+ Although the script can process multiline FASTA files, if you're unsure about your input FASTA file format, it's recommended to first convert it to single-line format using the multi2single module.
+ The script is effective for selectively renaming IDs/headers in FASTA files (with a user's specific input), such as genome assemblies from databases, where you want to exclude and not rename certain elements like scaffolds.

- Parameter explanation
	1. python 3: Call python 3
	1. prfxfindreplace.py:  Call prfxfindreplace module
	1. python3 prfxfindreplace.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --find-ptrn: Indicate a pattern in prefix ID/header name file (accept both integer and strings but no space)
  		+ --replace: Indicate a new prefix ID/header name for replacement
		+ --out: Indicate an output renamed single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### idextract (idx)
- To extract matched IDs and their corresponding sequences.
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.
 	+ Example file: [header_id.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/header_id.txt) and [idextract.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/idextract.fa) in the "example_data" folder.

Example usage
```
python3 idextract.py --input-seq test_dna.fasta --input-header header_id.txt --out output_extracted.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).
	
- Parameter explanation
	1. python 3: Call python 3
	1. idextract.py:  Call idextract module
	1. python3 idextract.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --input-header: Indicate an input ID and header (without ">") text file and its path
		+ --out: Indicate an output ID matched and extracted single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### idextloct (iel)
- To extract matched IDs, locations and their corresponding sequences (focused on a single ID)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.
	+ Example file: [header_id.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/header_id.txt) and [idextloct.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/idextloct.fa) in the "example_data" folder.

Example usage
```
python3 idextloct.py --input-seq test_dna.fasta --header-id test3_3%week --start 2 --end 10 --out output_test.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. idextloct.py:  Call idextloct module
	1. python3 idextloct.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --header-id: Indicate an input ID and header (without ">") name and pattern
		+ --start: Indicate a start position to extract (please use -1 value due to the python index)
		+ --end: Indicate en end position to extract (please use -1 value due to the python index)
		+ --out: Indicate an output ID matched and extracted single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)
	
### idextloctmlt (iem)
- To extract matched IDs, locations and their corresponding sequences (focused on multiple IDs)
  	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of IDs.
	+ Output: A single-line fasta with extracted IDs and their corresponding sequences.
 	+ Example file: [hdrmulti_id.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/hdrmulti_id.txt), [idextract.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/idextract.fa) and [idextloct.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/idextloct.fa) in the "example_data" folder.

Example usage
```
python3 idextloctmlt.py --input-seq test_dna.fasta --input-extract input_extract.txt --out output_extest.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. idextloctmlt.py:  Call idextloctmlt module
	1. python3 idextloctmlt.py --help: Check help menu
	--input-seq: Indicate an input single-line fasta file and its path
	--input-extract: Indicate an input ID and header (without ">") text file (a tap-separated) and its path including start and end positions (please use -1 value due to the python index)
	--out: Indicate an output ID matched and extracted single-line fasta file and its path
	--t: Specify thread numbers (integer only)
	--mem: Specify memory numbers (integer only with Gb size)

### revcomplt (rcp)
- To make a reverse complement sequence
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A reverse complement single-line fasta.
 	+ Example file: [revscomplt.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/revscomplt.fa) in the "example_data" folder.

Example usage
```
python3 revcomplt.py --input-seq test_dna.fasta --out output_revctest.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. revcomplt.py:  Call revcomplt module
	1. python3 revcomplt.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path
		+ --out: Indicate a reverse complement converted output single-line fasta file and its path
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### findcntdupl (fcd)
- To find and count the duplicated IDs and sequences from a multi-fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A text with a duplication number.
	+ Example file: [dupids.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/dupids.fa) in the "example_data" folder.

Example usage
```
python3 findcntdupl.py --input-seq test_dna2.fasta --out output_files.txt --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. findcntdupl.py:  Call findcntdupl module
	1. python3 findcntdupl.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --out: Indicate an output text file after finding and counting the duplicated IDs and sequences (only for both matched IDs and their corresponding sequences)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### removedupl (rvp)
- To remove the duplicated IDs and sequences from a multi-fasta file
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A single-line fasta after removing updated IDs and sequences.
 	+ Example file: [rmvdup.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/rmvdup.fa) in the "example_data" folder.

Example usage
```
python3 removedupl.py --input-seq test_dna2.fasta --outfasta output_testdupl.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. removedupl.py:  Call removedupl module
	1. python3 removedupl.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --outfasta: Indicate an output fasta file after removing the duplicated IDs and sequences (only for both matched IDs and their corresponding sequences)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### subsetfa (ssf)
- To make a subset of data with a sequence length filter option
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file.
	+ Output: A subseted single-line fasta.
 	+ Example file: [subsetfa.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/subsetfa.fa) in the "example_data" folder.

Example usage
```
python3 subsetfa.py --input-seq test_mRNA1.fasta --filter 50 --out output_subset.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. subsetfa.py:  Call subsetfa module
	1. python3 subsetfa.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --filter: Indicate a length size to filter out
		+ --out: Indicate an output fasta file after filtering the length size (integer only)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### extractptrn (xpt)
- To make a subset of data with find, filter and extract options
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A fasta file and a list of patterns.
	+ Output: A extracted single-line fasta with IDs and their corresponding sequences.
 	+ Example file: [find_ptrn.txt](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/find_ptrn.txt), [find_ptrn_mlt.txt
](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/find_ptrn_mlt.txt), [find_ptrn_mlt_seq.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/find_ptrn_mlt_seq.fa), and [find_ptrn_seq.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/find_ptrn_seq.fa) in the "example_data" folder.

Example usage
```
python3 extractptrn.py --input-seq test_dna1.fasta --input-pattern seq_pattern.txt --input-length 45 --out output_pattern.fasta --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. extractptrn.py:  Call extractptrn module
	1. python3 extractptrn.py --help: Check help menu
		+ --input-seq: Indicate an input single-line fasta file and its path (accept reverse complement sequences)
		+ --input-pattern: Indicate an input text file and its path to find and filter a specific sequence pattern (accept reverse complement sequences)
		+ --input-length: Indicate a length size to filter out
		+ --out: Indicate an output fasta file after filtering the length size (integer only)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### eachfastats (efs)
- To generate each line fasta statistic for a multi-line fasta
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A multi-line fasta file.
	+ Output: A summary of single-line fasta with its corresponding sequence length.
	+ Example file: [all_stat_asm.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/all_stat_asm.fa) and [stat_asm_mlt_unlm.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/stat_asm_mlt_unlm.fa) in the "example_data" folder.

Example usage
```
python3 eachfastats.py --input-seq test_dna.fasta --out output_dna.txt --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).

- Parameter explanation
	1. python 3: Call python 3
	1. eachfastats.py:  Call eachfastats module
	1. python3 eachfastats.py --help: Check help menu
		+ --input-seq: Indicate an input multi-line fasta file and its path
		+ --out: Indicate an output text file with the sequence length
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### allfastats (afs)
- To generate a summary of multi-line fasta statistics. 
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: A multi-line fasta file.
	+ Output: A summary of multi-line fasta with its diverse statistics (e.g. sequence length, GC content, N50, and more).
  	+ Example file: [all_stat_asm.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/all_stat_asm.fa) and [stat_asm_mlt_unlm.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/stat_asm_mlt_unlm.fa) in the "example_data" folder.

Example usage
```
python3 allfastats.py --input-seq test_dna.fasta --out output_dna.txt --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).
+ Both eachfastats and allfastasts modules use the same multi-line FASTA input. Eachfastats focuses on individual FASTA lines and their sequence lengths, while allfastasts provides a summary of assembly statistics, such as for genomes or transcriptomes.

- Parameter explanation
	1. python 3: Call python 3
	1. allfastats.py:  Call allfastats module
	1. python3 allfastats.py --help: Check help menu
		+ --input-seq: Indicate an input multi-line fasta file and its path
		+ --out: Indicate an output text file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### asmstatsunlm (mfs)
- To generate a summary of multi-line fasta statistics for unlimited fasta files. An extended version of the allfastats for multiple fasta files. 
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Unlimited multi-line fasta files.
	+ Output: A summary of multi-line fasta with its diverse statistics (e.g. sequence length, GC content, N50, and more).
  	+ Example file: [all_stat_asm.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/all_stat_asm.fa) and [stat_asm_mlt_unlm.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/stat_asm_mlt_unlm.fa) in the "example_data" folder.

Example usage
```
python3 asmstatsunlm.py --input-seqs test_dna1.fasta test_dna2.fasta test_dna3.fasta test_dna4.fasta test_dna5.fasta --out output_dna.txt --t 1 --mem 2
```
+ If your input FASTA file is in multi-line format, the script will automatically convert it to single-line format for processing (embedded pipeline).
+ Both allfastasts and asmstatsunlm modules take the same multi-line FASTA files. Allfastats generates a summary of assembly statistics from a single FASTA file, while asmstatsunlm does so for multiple FASTA files, useful for genome or transcriptome analysis.

- Parameter explanation
	1. python 3: Call python 3
	1. asmstatsunlm.py:  Call asmstatsunlm module
	1. python3 asmstatsunlm.py --help: Check help menu
		+ --input-seqs: Indicate an input multi-line fasta file and its path (Separated by a space for each new fasta file)
		+ --out: Indicate an output text file with the summary of statistics
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)

### concatenate (ccf)
- To make a concatenated fasta file for unlimited fasta files. 
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multiple fasta files with the same prefix IDs.
	+ Output: A single-line fasta with concatenated IDs and their corresponding sequences.
  	+ Example file: [concat_seq1.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/concat_seq1.fa), [concat_seq2.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/concat_seq2.fa), [concat_seq3.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/concat_seq3.fa), [concat_seq4.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/concat_seq4.fa), and [concat_seq1.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/concat_seq5.fa) in the "example_data" folder.

Example usage
```
python concatenate.py --input-seqs test_dna1.fasta test_dna2.fasta test_dna3.fasta test_dna4.fasta test_dna5.fasta --out output_concat.fasta --t 1 --mem 2
```
+ For optimal use of this module, ensure that all input FASTA files have matching prefix IDs and headers and are formatted as single-line FASTA before concatenation.
+ Before using the concatenate module, ensure your FASTA files are in single-line format by using the renameid and multi2singleline modules, even though the pipeline automatically converts multi-line FASTA files.

- Parameter explanation
	1. python 3: Call python 3
	1. concatenate.py:  Call concatenate module
	1. python3 concatenate.py --help: Check help menu
		+ --input-seqs: Indicate input multiple single-line fasta files and their path (Separated by a space for each new fasta file)
		+ --out: Indicate a concatenated output fasta file and its path
		+ --t: Specify thread numbers (intergr only)
		+ --mem: Specify memory numbers (integer only with Gb size)
 
### translated (tls)
- To find the translated sequences as a protein including the open reading frames (ORFs)
 	+ Requirement: The script of Python/bash requires a Python library.
	+ Input: Multline fasta files.
	+ Output: A single-line fasta with translation and selected ORFs (nucleotide and protein sequences).
	+ Example file: [trnsldna.fa](https://github.com/OZTaekOppa/FastaHandler/blob/main/example_data/trnsldna.fa) in the "example_data" folder.

Example usage
```
python translatedna.py --input-seq test_dna.fasta --out output/folder --t 1 --mem 2
```
+ For optimal use of this module, ensure that all input FASTA files have matching prefix IDs and headers and are formatted as single-line FASTA before concatenation.
+ Before using the translatedna module, ensure your FASTA files are in single-line format by using the renameid and multi2singleline modules, even though the pipeline automatically converts multi-line FASTA files.

- Parameter explanation
	1. python 3: Call python 3
	1. translatedna.py:  Call translatedna module
	1. python3 translatedna.py --help: Check help menu
		+ --input-seq: Indicate input single-line fasta files and their path
		+ --out: Indicate a translated output fasta file and its path (a total of four fasta files for both nucleotide and protein sequences)
		+ --t: Specify thread numbers (integer only)
		+ --mem: Specify memory numbers (integer only with Gb size)


## FAQ

We encourage users to use the [Issues](https://github.com/OZTaekOppa/FastaHandler/issues).


## WIKI PAGE

Please see the GitHub page.


## AUTHORS

**Hyungtaek Jung** and the [**National Centre for Indigenous Genomics**](https://ncig.anu.edu.au/).


## COPYRIGHT

The full **FastaHandler** is distributed under the MIT license. 
