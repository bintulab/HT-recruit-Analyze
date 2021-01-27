# HT-recruit Analyze
Analysis of high-throughput recruitment screens, from FASTQ to element enrichments.

Paper: <https://www.biorxiv.org/content/10.1101/2020.09.09.288324v1> (Tycko, J., ..., Bassik, M.C.#, Bintu, L.#, 2020)

# Install dependencies
Install `bowtie` (We use the original Bowtie 1 software for read alignment)

# Clone the github repo
`git clone https://github.com/bintulab/HT-recruit-Analyze`

Navigate into the directory. The pipeline depends on having this directory structure.

`cd HT-recruit-Analyze`

# Prepare Python environment
Use `conda` to install the required packages in a python3 environment:

`conda create -n HTrec3 python=3.7`

`conda install -n HTrec3 numpy pandas seaborn Cython`

# Activate environment
`source activate HTrec3`

`(HTrec3)` should appear on the command line

If this is the first time using the environment, finish installations into your environment with `pip`:
`conda install pip`
`pip install pysam HTSeq`

# Create an index of the library
First, create your reference index for alignment, using the oligo file. The oligo file should have 2 columns: `label, oligo sequence`.

It is recommended to trim the cloning adapters off the oligo sequence by using the `-e` and `-s` parameters of `makeIndices.py`. 

Optionally, if your NGS reads do not cover the full length of the domain, it is recommended to further trim so the index is only inclusive of the sequences that will be covered by the read (and no shorter). Being precise with the index trimming decreases ambiguous alignments of a read to multiple elements that may have homology in some positions (e.g. in a tiling screen). If your reads are too short, there is a greater risk of counting reads from clones with problematic mutations (e.g. a deletion outside the sequenced portion such that the domain is truncated).

`makeIndices.py -o -t oligoFile.csv shortLibraryName fullName`

# Align reads to an index of the library
Then, align your reads to that index.

Use the trim option `-t` to remove the constant region in the beginning of the sequencing read. Note that `-t` will vary depending on the length of stagger sequence in the primer. Use the length option `-l` to determine the remaining length of read that will be used for the alignment (e.g. removing extraneous sequence from the end of the read if it exceeds the length of the elements). Use `-p` for parallel processing on a desired number of processors.

`makeCounts.py FASTQfile OutputName shortLibraryName -t TrimLength -m NumMismatches -l ReadLength -p Processors`

It is helpful to optimize the parameters on a small FASTQ (e.g. a subsample of 250,000 reads) with the goal of a high unique alignment percentage and a low ambiguous alignment percentage. The optimal parameters will depend on library composition. Due to errors in oligo synthesis and library sequencing, the optimized alignment percentage is often ~40-60% and the ambiguous percentage is <0.5%.

If there is an additional FASTQ file for the same sample (e.g. it was re-sequenced at a later date for more depth), you can add it to the counts with `-a PathToFASTQ`.

# Check the diversity and uniformity of the count distributions
`plotDist.py OutputName Data/Sample1_Bound Data/Sample1_Unbound... -l Label1 Label2 ...`

You can list as many count file names as needed.

# Compute OFF:ON ratios 
Compare the read counts between your OFF and ON samples for a given bio-replicate. If, for example, you want to use random sequence controls as the negative control set and shift the ratios such that the median of negative controls is 0, then you can use `-n 'Random'`. In this case, the first word in the `_` delimited label of your random sequence controls, when making the index, must be `'Random'`.

`makeRhos.py Data/Sample1_Bound_counts.csv Data/Sample1_Unbound_counts.csv SampleRep1 -n 'Random' -b 'none'`

# Combine bio-replicates
Finally, combine the values from bio-replicates.

`combineRhos.py Data/SampleRep1_rhos.csv Data/SampleRep2_rhos.csv SampleName`
