# HT-recruit Analyze
Analysis of high-throughput recruitment screens, from FASTQ to element enrichments.

Paper: <https://www.biorxiv.org/content/10.1101/2020.09.09.288324v1> (Tycko, J., ..., Bassik, M.C.#, Bintu, L.#, 2020)

# Activate environment
Use `conda` to install (list here: )? in a python3 environment

You will also need to install `bowtie` for read alignment.

# Align reads to an index of the library
First, create your reference index.

`makeIndices.py oligoFile.csv shortScreenName fullName`

Then, align your reads to that index.

`makeCounts.py FASTQfile OutputName shortScreenName -t TrimLength -m NumMismatches -l ReadLength`

It is helpful to optimize the parameters on a small FASTQ (e.g. a subsample of 250,000 reads) with the goal of a high unique alignment percentage and a low ambiguous alignment percentage. The optimal parameters will depend on library composition. Due to errors in oligo synthesis and library sequencing, the optimized alignment percentage is often ~40-60% and the ambiguous percentage is <0.5%.

# Compute OFF:ON ratios and combine bio-replicates
Compare the read counts between your OFF and ON samples for a given bio-replicate. 

`makeRhos.py *`

You can then do this optional step to check the diversity and uniformity of the elements' count distributions. You can list as many count file names as needed.
`plotDist.py OutputName countFile1 countFile2...`

Finally, combine the values from bio-replicates.

`combineRhos.py *`