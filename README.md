# PAMPHLET

PAM Prediction HomoLogous Enhancement Toolkit (DEV version)

## Dependencies

Python packages:  `argparse`, `Biopython`, `pyfaidx`, `urllib`, `requests`, `func_timeout`

Softwares: `WebLOGO`(need `GhostScript`), `seqkit` (ADD TO PATH)

NEED INTERNET CONNECTION

## Parameters

`-s`,`--spacer`: Input spacer sequences, required.

`-r`,`--repeat`: Input repeat sequences, need to be used together with `-p` or `-P` and `-r`, not required.

`-p`,`--protein`: Relative Cas protein sequences, need to be used together with `-r`, not required.

`-o`,`--outdir`: Output directory, required. Need to specify a path that does not exist。

`-O`,`--orientation`: CRISPR array orientation, choices = `positive` or `negative`. Default is `positive`

`-u`,`--unique`: Unique mode, if set, only revise unique spacer.

`-P`,`--proteinblast`: Protein blastp output file in outfmt 6. Do not co-submit with `-p`

`-l`,`--proteinlen`: Protein length, need to be used together with `-P`

`-d`,`--spacerdb`: Spacer blast database, choices = `phage` or `prokaryote`. Default is `prokaryote`

`-L`,`--flanklen`: Length of flank sequence, default is `12`

`-R`,`--refmode`: Reference mode, choices = `r` or `nr` or `a`.  `r`  means reference mode, which use taxid to revise homologue protein. `nr` is non-reference mode (use most-redundant taxid as homologue protein reference) and `a` is all-mode (use all significant hit as homologue protein). Default is `r`

`-f`,`--freqmode`: Base frequency calculation mode, choices = `linear` and `sigmoid`, which effect weblogo. Default is `sigmoid`

`-b`,`--blastmode`: Spacer blastn mode. Common mode means use default blastn parameters and strict mode means use specific parameters, which could get more putative protospacers but also could cause false positives. Choices = `relax` and `common`, default is `common`

`--pcovs`: Minimum percent coverage of spacer sequence, default is `0.9`.

`--pident`: Minimum percent identity of spacer sequence, default is `0.9`.

`--rident`: Minimum percent identity of repeat sequence, default is `0.8`.

`--MaxProteinNumber`: Maximum number of protein homologs, default is `20`. If the size is too large, NCBI will forbidden your request.

## INSTALLATION

Install the required Python packages:
```bash
pip install argparse biopython pyfaidx urllib3 requests func_timeout
```

Install WebLOGO and seqkit, and ensure they are added to your PATH.

## USAGE

### BASIC USAGE

```bash
python pamphlet.py -s spacers.fa -o output_directory
```

### ADVANCED USAGE

Specify additional parameters as needed, here is an example with input protein alignment tabular txt file
```bash
python pamphlet.py -s spacers.fasta -r GCTAGTGTAGCTGTCAGTCGATGTCAC -P proteins_alignment.tabular.txt -l 1092 -O positive -u -R a -b relax -d prokaryote --pcovs 0.9 --pident 0.9 --rident 0.9 --MaxProteinNum 20
```

Specify additional parameters as needed, here is an example with input protein sequence fasta file
```bash
python pamphlet.py -s spacers.fasta -r GCTAGTGTAGCTGTCAGTCGATGTCAC -p protein_sequence.fasta -O positive -u -R a -b relax -d prokaryote --pcovs 0.9 --pident 0.9 --rident 0.9 --MaxProteinNum 20
```

## NOTICE

We recommend preparing the protein BLAST tabular text file prior to executing PAMPHLET, as this practice can save computational time. While utilizing PAMPHLET for blastp is also feasible, repeated use may lead to increased time for individual BLAST operations. Therefore, to conserve online BLAST resources, please complete the protein homology search using BLAST in advance.
