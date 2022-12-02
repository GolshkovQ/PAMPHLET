# PAMPHLET

PAM Prediction HomoLogous Enhancement Toolkit

## Dependencies

Python packages:  argparse, Biopython, pyfaidx, urllib, requests, 

Softwares: WebLOGO, seqkit (ADD TO PATH)

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
