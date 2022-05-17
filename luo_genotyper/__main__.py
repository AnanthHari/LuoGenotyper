import argparse, os
import statistics
from posixpath import splitext
from .allele_database import ImgtNovelAlleleDatabase as AlleleDatabase
from .common import log, initialize_logger, resource_path, script_path
from .luo_ighv_genotyper import *

initialize_logger()

ALLELE_DATABASE_PATH = resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-no_duplicates.fa')
FUNC_ALLELE_DATABASE_PATH = resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-no_duplicates-functional.fa')
IGBLAST_V_ALLELE_DATABASE_PATH = resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-no_duplicates-trimmed.fa')
D_ALLELE_DATABASE_PATH = resource_path('IGHD_igblast.fasta')
J_ALLELE_DATABASE_PATH = resource_path('IGHJ_igblast.fasta')
AUX_FILE = resource_path('human_gl.aux')
KGP_COVERAGES = resource_path('1kgp-coverages.tsv')

TWO_CNV_ANALYSIS_SCRIPT = script_path('two_copy_gene_analysis.sh')

parser = argparse.ArgumentParser(description='luo_genotyper: IGHV Genotyping using Short Read WGS -- luo et al (2019)')

parser.add_argument(
    'bam_file',
    type=str,
    help='Input BAM file: Reads mapped against the functional allele database'
)

parser.add_argument(
    '--reads1',
    type=str,
    default='',
    help='Path to file containing the first mate of the paired-end reads'
)

parser.add_argument(
    '--reads2',
    type=str,
    default='',
    help='Path to file containing the second mate of the paired-end reads'
)

parser.add_argument(
    '--output_dir',
    type=str,
    default=os.getcwd(),
    help='Path to output directory'
)

parser.add_argument(
    '--sequencing_depth',
    type=float,
    default=0.0,
    help='Sequencing depth'
)

parser.add_argument(
    '--read_length',
    type=int,
    default=150,
    help='Read length'
)

parser.add_argument(
    '--threads',
    type=int,
    default=6,
    help='Max number of threads to use'
)

def main():
    args = parser.parse_args()
    run_luo_genotyper(args.bam_file, args.reads1, args.reads2, args.output_dir,  args.sequencing_depth, args.read_length, args.threads)

def run_luo_genotyper(bam_file: str='',
                      reads1: str='',
                      reads2: str='',
                      output_base_dir: str='',
                      sequencing_depth: int=0,
                      read_length: int=0,
                      threads: int=6):

    
    EXPERIMENT_NAME = os.path.splitext(os.path.basename(bam_file))[0]
    OUT_DIR = os.path.join(output_base_dir, EXPERIMENT_NAME)
    MINIMAP_EXTRACTED_READS_FA = os.path.join(OUT_DIR, 'minimap2_out', f'{EXPERIMENT_NAME}-extracted_reads.fa')
    MINIMAP_EXTRACTED_READS_FQ = os.path.join(OUT_DIR, 'minimap2_out', f'{EXPERIMENT_NAME}-extracted_reads.fq')
    
    IGBLAST_OUT_PREFIX = f'{OUT_DIR}/igblast_out'
    SPADES_OUT_PREFIX = f'{OUT_DIR}/spades_out'
    TWO_CNV_GENES_OUTPUT = f'{OUT_DIR}/spades_out/two_cnv_genes.txt'
    ALL_GENES_OUTPUT = f'{OUT_DIR}/spades_out/all_genes.txt'
    BOWTIE2_OUT_PREFIX = f'{OUT_DIR}/bowtie2_out'
    GATK_OUT_PREFIX = f'{OUT_DIR}/gatk_out'
    HAPCUT2_OUT_PREFIX = f'{OUT_DIR}/hapcut2_out'
    
    allele_db = AlleleDatabase(**{'db_fasta_path': resource_path('IMGT_IGHV_reference_allele_db-updated+oscar_novel-aligned.fasta'),
                                  'consensus_path': resource_path('V-QUEST-reference-allele-db+no-period-references+consensus.clustalw.consensus-seq.fasta'),
                                  'ignored_alleles_path': resource_path('ignored_alleles.txt')})

    if [reads1, reads2] == ['',''] and not os.path.exists(bam_file):
        log.error('Please provide a bam file of the input reads mapped against the IGHV functional allele db or provide paths to fastqs of both mates of the paired-end reads.')
        os._exit(1)

    # 1. minimap2 returns alignment in SAM format. Convert it to bam and extract reads in both FASTA and FASTQ formats
    if os.path.exists(bam_file):
        log.info(f'Input bam file is {bam_file}. Not running minimap2...')
    else:
        log.info('---------- Mapping input reads against functional allele db using minimap2 ----------')
        map_paired_reads_allele_db(threads, FUNC_ALLELE_DATABASE_PATH, reads1, reads2, bam_file)
    log.info('---------- Extracting reads from minimap2-produced BAM ----------')
    extract_reads(bam_file, MINIMAP_EXTRACTED_READS_FA, MINIMAP_EXTRACTED_READS_FQ)

    # 2. Instantiate the genotyper
    log.info('---------- Instantiating LuoGenotyper instance ----------')
    genotyper = LuoGenotyper(MINIMAP_EXTRACTED_READS_FQ, allele_db, IGBLAST_OUT_PREFIX, SPADES_OUT_PREFIX, HAPCUT2_OUT_PREFIX, ALL_GENES_OUTPUT, TWO_CNV_GENES_OUTPUT)

    # 3. Run igBLAST and extract V segment mappings
    log.info('---------- Running IgBLAST ----------')
    run_igblast(IGBLAST_V_ALLELE_DATABASE_PATH, D_ALLELE_DATABASE_PATH, J_ALLELE_DATABASE_PATH, AUX_FILE, MINIMAP_EXTRACTED_READS_FA, IGBLAST_OUT_PREFIX)

    # 4. Process igBLAST output and output reads assigned to each gene
    log.info('---------- Processing IgBLAST output ----------')
    genotyper.extract_igblast_mappings()

    # 5. Run SPAdes assembly
    log.info('---------- Running SPAdes assembly ----------')
    genotyper.run_spades(threads=threads)

    # 6. Get copy number estimate
    log.info('---------- Getting CN estimate ----------')
    genotyper.get_cnv_estimate(sequencing_depth=sequencing_depth, read_length=read_length, seq_depth_adjustment=0.9)
    y = genotyper.spades_copy_number_estimates
    log.info(f'Summary statistics of CNV esimates:\n Mean = {statistics.mean(y)}, Median = {statistics.median(y)}, STDEV = {statistics.stdev(y)}, Min = {min(y)}, Max = {max(y)}')
    genotyper.get_point_cnv_estimate()

    # 7. Proceed further for genes with two copies
    log.info(f' List of SPAdes-assembled genes and their esitmated copy numbers: {sorted(list(zip(genotyper.spades_assembled_genes, genotyper.point_cnv_estimates)), key=lambda x: x[1])}')
    genotyper.get_two_copy_genes()
    
    log.info('---------- Analyzing two copy genes ----------')
    # cmd = f'cat {TWO_CNV_GENES_OUTPUT} | xargs -I % {TWO_CNV_ANALYSIS_SCRIPT} % {EXPERIMENT_NAME} {IGBLAST_OUT_PREFIX} {SPADES_OUT_PREFIX} {BOWTIE2_OUT_PREFIX} {GATK_OUT_PREFIX} {HAPCUT2_OUT_PREFIX}'
    # print(cmd)
    # ! {cmd}

    with open(TWO_CNV_GENES_OUTPUT) as f:
        two_cnv_genes = f.read().strip().split('\n')
    for g in two_cnv_genes:
        if not os.path.exists(os.path.join(HAPCUT2_OUT_PREFIX, g, f'{g}-haplotype.txt')):
            log.info(f'---------- Running two copy analysis for gene: {g} ----------')
            cmd = f'{TWO_CNV_ANALYSIS_SCRIPT} {g} {EXPERIMENT_NAME} {IGBLAST_OUT_PREFIX} {SPADES_OUT_PREFIX} {BOWTIE2_OUT_PREFIX} {GATK_OUT_PREFIX} {HAPCUT2_OUT_PREFIX}'
            sp.run(shlex.split(cmd))
        else:
            log.info(f'---------- Not running HAPCUT2 on gene {g} ----------')

if __name__ == '__main__':
	main() 
