# ngs-analysis

Intended for analysis of sequencing reads that span multiple DNA or protein parts. For instance, given a library of protein variants linked to DNA barcodes, it can answer questions like:

- How accurate are the variant sequences, at the DNA or protein level?
- How frequently is the same barcode linked to two different variants?
- Which reads contain parts required for function (e.g., a kozak start sequence, or a fused protein tag)?

This kind of analysis often involves parsing raw sequencing reads for DNA and/or protein sub-sequences (parts), then mapping the parts to a reference of anticipated part combinations. This package offers a simple workflow: 

1. Define how to parse reads into parts using plain text expressions (no code)
2. Test the parser on simulated DNA sequences (e.g., your vector map)
3. Parse a batch of sequencing samples
4. Map the (combination of) parts found in each read to your reference

It’s been tested with Illumina paired-end reads and Oxford Nanopore long reads. Under the hood it uses [NGmerge](https://github.com/jsh58/NGmerge) to merge paired reads and [MMseqs2](https://github.com/soedinglab/MMseqs2) for sequencing mapping. It is moderately performant: 1 million paired-end reads can be mapped to a reference of 100,000 variant-barcode pairs in ~1 minute.

# Installation

```bash
pip install ngs-analysis
```

Tested on Linux and MacOS (Apple Silicon).

