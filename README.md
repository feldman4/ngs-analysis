# ngs-analysis

Convenient analysis of sequencing reads that span multiple DNA or protein parts. For instance, given a library of protein variants linked to DNA barcodes, this tool can answer questions like:

- How accurate are the variant sequences, at the DNA or protein level?
- How frequently is the same barcode linked to two different variants?
- Which reads contain parts required for function (e.g., a kozak start sequence, or a fused protein tag)?

This kind of analysis often involves parsing raw sequencing reads for DNA and/or protein sub-sequences (parts), then mapping the parts to a reference of anticipated part combinations. This package offers a simple workflow: 

1. Define how to parse reads into parts using plain text expressions (no code)
2. Test the parser on simulated DNA sequences (e.g., your vector map)
3. Parse a batch of sequencing samples
4. Map the (combination of) parts found in each read to your reference

It’s been tested with Illumina paired-end reads and Oxford Nanopore long reads. Under the hood it uses [NGmerge](https://github.com/jsh58/NGmerge) to merge paired reads and [MMseqs2](https://github.com/soedinglab/MMseqs2) for sequencing mapping. It is moderately performant: 1 million paired-end reads can be mapped to a reference of 100,000 variant-barcode pairs in ~1 minute.

# Workflow

A cartoon example with two reference sequences, each consisting of a variant linked to a barcode:

<img src="examples/sequences.png" alt="sequences" width="400"/>

Here's the analysis workflow and outputs:

<img src="examples/workflow.png" alt="analysis workflow" width="400"/>

Check out the [example notebook for paired end reads](examples/paired_reads/paired_read_example.ipynb) for details.

### **TL;DR**

Run `ngs-analysis --help` to see available commands.

1. Make an empty directory, add `config.yaml` and `samples.csv` based on the example.
2. Add `reference_dna.csv` with anticipated DNA sequences (including adapters).
3. Run `ngs-analysis setup`. Add `--clean` to start the analysis from scratch.
4. Run `ngs-analysis dna_to_designs`. Check that `designs.csv` is accurate; if not, fix `config.yaml`.
5. 
    - If you have paired-end data, put it in `0_paired_reads/` and run `ngs-analysis merge_read_pairs <sample>`.
    - If you have single-end data (e.g., nanopore), put it in `1_reads/`.
6. Run `ngs-analysis parse_reads <sample>`. Check that `2_parsed/<sample>.parsed.pq` looks alright (with pandas, use `pd.read_parquet`)
7. Run `ngs-analysis map_parsed_reads <sample>`. Results are in `3_mapped/<sample>.mapped.csv`


# Install

```bash
pip install ngs-analysis
```

Make sure that the `mmseqs` and `NGmerge` executables are available (NGmerge only needed for paired reads).

Tested on Linux and MacOS (Apple Silicon).

