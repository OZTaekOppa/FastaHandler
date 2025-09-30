# A main module for FastaHandler created by Hyungtaek Jung

#!/usr/bin/env python3

import sys
import subprocess

usage = '''FastaHandler: Fasta File Manipulation Toolkit
version 1.0.1

Usage: python3 fastahandler.py <module> <parameters>

Modules:
AllFastaStats\t| all_fa_stats\tGenerate a summary of multi-line fasta statistics.
AssemblyStatsUnlimit\t| asm_stats_unlimit\tGenerate a summary of multi-line fasta statistics (Multiple).
ChrPanSpecNameExtract\t| chr_pansn_extract\tMake sequence partition based on the name of prefix (CL) in the fasta header.
ConcatenateFasta\t| concatenate_fa\tMake a concatenated fasta file (Multiple).
EachFastaStats\t| each_fa_stats\tGenerate each line fasta statistic for a multi-line fasta.
ExtractPattern\t| extract_pattern\tMake a subset of data with find, filter and extract.
FindAnchorTrim\t| find_anchor_trim\tMake a subset of fasta with find, anchor, trim, and extract.
FindCountDuplication\t| find_count_duplicate\tFind and count the duplicated IDs and sequences.
FindMergeFasta\t| find_merge_fa\tFind and merge fasta files if they matched with a certain pattern.
Gfa2Fasta\t| gfa2fa\tConvert a gfa into a single-line fasta.
IdExtractLocation\t| id_extract_location\tExtract matched IDs, locations and their corresponding sequences.
IdExtractLocationMultiple\t| id_extract_multi_location\tExtract matched IDs, locations and their corresponding sequences (Multiple).
IdExtract\t| id_extract\tExtract matched IDs and their corresponding sequences.
Multiple2Each\t| multi2each\tConvert a multi-fasta file to a each (split) fasta.
Multiple2Single\t| multi2single\tConvert a multi-fasta (multiline) file to a single-line fasta.
OverlapSplit\t| overlap_split\tConvert a multi-fasta (multiline) file to a single-line fasta.
PangenomeIdRename\t| pangenome_id_rename\tRename prefix IDs and headers for a Pangenome format.
PrefixPatternReplace\t| prefix_pattern_replace\tReplace and rename prefix IDs and headers with a user’s input (Only).
PrefixRename\t| prefix_rename\tRename prefix IDs and headers with a user’s input.
PrefixSelectRename\t| prefix_select_rename\tRename prefix IDs and headers with a user’s input (Only).
RemoveDuplication\t| remove_duplicate\tRemove the duplicated IDs and sequences.
RenameId\t| rename_id\tRename prefix IDs and headers.
ReverseComplement\t| reverse_complement\tMake a reverse complement sequence.
SizePatternSearch\t| size_pattern_search\tFind unique and similariy sequecnes against its own input sequecne
SubsetFasta\t| subset_fa\tMake a subset of data with a sequence length filter.
TranslateSequence\t| translate_dna\tFind the translated sequences as a protein and open reading frames (ORFs).

Use <module> --help for module usage.'''

module_map = {
    'AssemblyStatsUnlimit': 'asm_stats_unlimit.py',
    'ChrPanSpecNameExtract': 'chr_pansn_extract.py',
    'ConcatenateFasta': 'concatenate_fa.py',
    'EachFastaStats': 'each_fa_stats.py',
    'ExtractPattern': 'extract_pattern.py',
    'FindAnchorTrim': 'find_anchor_trim.py',
    'FindCountDuplication': 'find_count_duplicate.py',
    'FindMergeFasta': 'find_merge_fa.py',
    'Gfa2Fasta': 'gfa2fa.py',
    'IdExtractLocation': 'id_extract_location.py',
    'IdExtractLocationMultiple': 'id_extract_multi_location.py',
    'IdExtract': 'id_extract.py',
    'Multiple2Each': 'multi2each.py',
    'Multiple2Single': 'multi2single.py',
    'OverlapSplit': 'overlap_split.py',
    'PangenomeIdRename': 'pangenome_id_rename.py',
    'PrefixPatternReplace': 'prefix_pattern_replace.py',
    'PrefixRename': 'prefix_rename.py',
    'PrefixSelectRename': 'prefix_select_rename.py',
    'RemoveDuplication': 'remove_duplicate.py',
    'RenameId': 'rename_id.py',
    'ReverseComplement': 'reverse_complement.py',
    'SizePatternSearch': 'size_pattern_search.py',
    'SubsetFasta': 'subset_fa.py',
    'TranslateSequence': 'translate_dna.py'
}


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] in ['--h', '--help']:
        print(usage)
        sys.exit(0)

    module = sys.argv[1]
    if module not in module_map:
        print('Unexpected module. Use --h for help.')
        sys.exit(0)

    script_path = f'{sys.path[0]}/scripts/{module_map[module]}'
    parameters = ' '.join(sys.argv[2:])
    subprocess.run(f'python3 {script_path} {parameters}', shell=True)

