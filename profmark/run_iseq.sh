#/bin/bash

iseq pscan2 Pfam-A_small.hmm AE014075.1_nucl.fasta --model 2
iseq gff-filter output.gff --max-e-value 1e-10 > output2.gff
