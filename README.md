# Flongle Assemble
ONT assembly pipeline with conservative polishing to prevent systematic sequence errors.

## Requirements

Create a mamba/conda environment with:
```bash
mamba create -n flye -c conda-forge -c bioconda \
  flye minimap2 racon seqkit ca-certificates openssl \
  canu miniasm minipolish raven-assembler nanopolish \
  samtools seqtk unicycler quast bandage mummer \
  biopython pandas numpy openjdk=17

mamba activate flye
```

## Usage
```bash
nextflow run ryandward/flongle_assemble --fastq ./ --ref_fasta reference.fa
```
