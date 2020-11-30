# HT-recruit Analyze
Analysis of high-throughput recruitment screens, from FASTQ to element enrichments.

Paper: <https://www.biorxiv.org/content/10.1101/2020.09.09.288324v1> (Tycko, J., ..., Bassik, M.C.#, Bintu, L.#, 2020)

# Activate environment
Use `conda` to install (list here: )? in a python3 environment

You will also need to install `bowtie` for read alignment.

# Align reads to an index of the library
First, create your reference index. The oligo file should have 2 columns: `label, oligo sequence`. It is recommended to trim the cloning adapters off the oligo sequence by using the `-e` and `-s` parameters of `makeIndices.py`. Optionally, if your NGS reads do not cover the full length of the domain, it is recommended to further trim so the index is only the part of the domain that will be covered by the read (and no shorter). Being precise with the index trimming helps decrease ambiguous alignments of a read to multiple domains that may have homology in some positions (e.g. in a tiling screen). If your reads are too short, there is a greater risk of counting reads from clones with problematic mutations (e.g. a deletion outside the sequenced portion such that the domain is truncated).

`makeIndices.py oligoFile.csv shortLibraryName fullName`

Then, align your reads to that index. Use the trim option `-t` to remove the constant region in the beginning of the sequencing read, and the length option `-l` to determine the remaining length of read that will be used for the alignment (e.g. removing poor quality positions from the end of the read).

`makeCounts.py FASTQfile OutputName shortLibraryName -t TrimLength -m NumMismatches -l ReadLength`

It is helpful to optimize the parameters on a small FASTQ (e.g. a subsample of 250,000 reads) with the goal of a high unique alignment percentage and a low ambiguous alignment percentage. The optimal parameters will depend on library composition. Due to errors in oligo synthesis and library sequencing, the optimized alignment percentage is often ~40-60% and the ambiguous percentage is <0.5%.

# Compute OFF:ON ratios and combine bio-replicates
Compare the read counts between your OFF and ON samples for a given bio-replicate. 

`makeRhos.py *`

You can then do this optional step to check the diversity and uniformity of the elements' count distributions. You can list as many count file names as needed.
`plotDist.py OutputName countFile1 countFile2...`

Finally, combine the values from bio-replicates.

`combineRhos.py *`
