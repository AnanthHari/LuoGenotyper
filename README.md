# Luo-genotyper


[![DOI](https://zenodo.org/badge/493383966.svg)](https://zenodo.org/badge/latestdoi/493383966)


A Immunoglobulin Heavy Chain Variable Gene genotyping tool for WGS short reads based on Luo et al. (2019). Link to the [paper](www.doi.org/10.26508/lsa.201800221).


## Installation

### Prerequisites

This tool needs the following software packages to be installed and the path to the binaries to be present in `$PATH`:

| Software | Version |
| -------- |:-------:|
| samtools | 1.15 |
| BWA | 0.7 |
| Bowtie2 | 2.4 |
| Minimap2 | 2.24 |
| ClustalW | 2.1 |
| SPAdes | 3.15 |
| IgBLAST | 1.16 |
| GATK | 4.2 |

Once you have installed the above dependencies, just download the latest release binary (see right toolbar) and install with pip (needs a `python3.8` environment):

Then run `pip install <binary.whl>`

### Environment and Dependencies

Installing LuoGenotyper using `pip` will automatically install the following dependencies:

- [logbook](https://logbook.readthedocs.io/en/stable/)
- [dill](https://pypi.org/project/dill/)
- [biopython](https://biopython.org/)
- [pysam](https://pysam.readthedocs.io/en/latest/api.html)

### Installing from source

If the binary fails to install for whatever reason, you can build the tool from source as follows:

```
conda create -n luo-genotyper python=3.8
conda activate luo-genotyper
git clone git@github.com:AnanthHari/LuoGenotyper.git ./LuoGenotyper
cd LuoGenotyper
python -m pip install --upgrade  build
python -m build
pip install dist/<.tar.gz or .whl build>
```


## Running LuoGenotyper

After installation with pip, simply use the command `luo_genotyper`. The only required input is a bam file that consists of the sample's WGS paired-end short reads mapped against the functional IGHV allele database using `minimap2`:
- If only the reads are known, they could be provided using the arguments `--reads1 READS1 --reads2 READS2`, and the alignments are written to the BAM provided as the argument.

The default output directory is the current working directory. The sequencing depth of the original sample (when mapped against a standard human genome reference like GRCh37/GRCh38) needs to be provided as an input.

```
$ luo_genotyper -h
usage: luo_genotyper [-h] [--reads1 READS1] [--reads2 READS2] [--output_dir OUTPUT_DIR] [--sequencing_depth SEQUENCING_DEPTH] [--read_length READ_LENGTH] [--threads THREADS] bam_file

luo_genotyper: IGHV Genotyping using Short Read WGS -- luo et al (2019)

positional arguments:
  bam_file              Input BAM file: Reads mapped against the functional allele database

optional arguments:
  -h, --help            show this help message and exit
  --reads1 READS1       Path to file containing the first mate of the paired-end reads
  --reads2 READS2       Path to file containing the second mate of the paired-end reads
  --output_dir OUTPUT_DIR
                        Path to output directory
  --sequencing_depth SEQUENCING_DEPTH
                        Sequencing depth
  --read_length READ_LENGTH
                        Read length
  --threads THREADS     Max number of threads to use
```
