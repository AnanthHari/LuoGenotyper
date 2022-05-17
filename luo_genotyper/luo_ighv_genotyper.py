import os, shlex, itertools, random, pysam
from pathlib import Path
import subprocess as sp
from .common import log, fasta_from_seq, SeqRecord
from .msa_wrappers import ClustalWWrapper
from Bio import Seq, SeqIO

# Bio - SeqRecord max id length
MAX_ID_LENGTH = 30
# 2-copy genes
TWO_COPY_GENES = ['IGHV6-1', 'IGHV1-18', 'IGHV3-20', 'IGHV1-24', 'IGHV2-26', 'IGHV1-45', 'IGHV5-51', 'IGHV1-58', 'IGHV3-72', 'IGHV3-73', 'IGHV3-74']

AMBIGUOUS_GENE_SET = ['IGHV4-4', 'IGHV4-59', 'IGHV4-61', 'IGHV4-38-2']

# def extract_reads_minimap(allele_db_fasta, reads_fasta, extracted_reads_fasta):
#     ''' Map reads_fasta against allele_db_fasta and return a file with mapped reads '''
#     with open(extracted_reads_fasta, 'w') as f:
#         p1 = sp.Popen(['minimap', '-w1', '-f1e-9',  f'{allele_db_fasta}', f'{reads_fasta}'], stdout=sp.PIPE)
#         p2 = sp.Popen(['cut', '-f', '1'], stdin=p1.stdout, stdout=f)
#         p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
#         p2.communicate()

def map_paired_reads_allele_db(threads, func_allele_db, fq1, fq2, bam_file):
    cmd1 = shlex.split(f'minimap2 -ax sr -w1 -f1e-9 -t {threads} --sam-hit-only {func_allele_db} {fq1} {fq2}')
    cmd2 = shlex.split(f'samtools sort -b -')
    bai_file = bam_file.replace('.bam', '.bai')
    with open(bam_file, 'w') as f:
        p1 = sp.Popen(cmd1, stdout=sp.PIPE)
        p2 = sp.Popen(cmd2, stdin=p1.stdout, stdout=f)
        p1.stdout.close()
        p2.communicate()
    pysam.index(bam_file, bai_file)

def extract_reads(bam_file, reads_fa, reads_fq):
    if os.path.exists(reads_fa):
        log.info(f'extract_reads: {reads_fa} exists.. Not extracting reads')
    else:
        pysam.fasta('-o', reads_fa, bam_file)
    if os.path.exists(reads_fq):
        log.info(f'extract_reads: {reads_fq} exists.. Not extracting reads')
    else:
        pysam.fastq('-o', reads_fq, bam_file)

def run_igblast(v_allele_db_path, d_allele_db_path, j_allele_db_path, aux_file, query_seq_fasta, igblast_output_prefix):
    # "igblastn -germline_db_V {IGBLAST_V_ALLELE_DATABASE_PATH} -germline_db_D {D_ALLELE_DATABASE_PATH} -germline_db_J {J_ALLELE_DATABASE_PATH} -auxiliary_data {AUX_FILE} -organism human -query {MINIMAP_EXTRACTED_READS_FA} -outfmt 7 | grep '^V' > {IGBLAST_TRIMMED_OUTPUT}"
    out_file = os.path.join(igblast_output_prefix, 'igblast-trimmed.tsv')
    if not os.path.exists(out_file):
        os.makedirs(igblast_output_prefix, exist_ok=True)
        with open(out_file, 'w') as f:
            p1 = sp.Popen(['igblastn', '-germline_db_V', v_allele_db_path, '-germline_db_D', d_allele_db_path, '-germline_db_J', j_allele_db_path, '-auxiliary_data', aux_file, '-organism', 'human', '-query', query_seq_fasta, '-outfmt', '7'], stdout=sp.PIPE)
            p2 = sp.Popen(['grep', r"^V"], stdin=p1.stdout, stdout=f)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            p2.communicate()
    else:
        log.warn('IgBLAST output is already present. Not running...')

# --------------- read filtering --------------- 
class Allele:
    def __init__(self, allele):
        self.allele = allele
        self.gene = self.allele.split('*')[0]

class igblast_mapping:
    def __init__(self, igblast_line):
        ighv_segment, self.query_id, self.subject_id, self.identity, self.alignment_length, self.mismatches, self.gap_opens, self.gaps, self.q_start, self.q_end, self.s_start, self.s_end, self.evalue, self.bit_score = igblast_line.split()
        assert ighv_segment == 'V', f'IGHV segment is not V. It is {ighv_segment}'
        self.query_id = self.query_id.lstrip('reversed|')
        self.ref = Allele(self.subject_id)
        self.identity = float(self.identity)
        self.s_start, self.s_end = int(self.s_start), int(self.s_end)


def toss_fair_coin():
    ''' Returns True/False '''
    return random.random() < 0.5

def get_cov_proportion(cov_gene, mappings_per_gene):
    ''' coverage per gene, mappings per gene '''
    total_gene_cov, len_gene_cov = sum(cov_gene), len(cov_gene)
    len_of_cov_outside_mapping = [len(cov_gene) - (m.s_end - m.s_start + 1) for m in mappings_per_gene]
    cov_per_allele = [(total_gene_cov - sum(cov_gene[m.s_start-1:m.s_end]))/(len_gene_cov - m.s_end + m.s_start - 1) for m in mappings_per_gene]
    return sum(cov_per_allele)/len(cov_per_allele)

def get_kmer_cov(contigs_fasta):
    ''' Outputs average k-mer coverage of the assembled contig '''
    with open(contigs_fasta, 'r') as f:
        contigs = [l.strip().split('_') for l in f.readlines() if l.startswith('>')]
    total_len, total_cov = 0, 0.0
    for i in contigs:
        contig_len, contig_cov = int(i[-3]), float(i[-1])
        total_cov += contig_len * contig_cov
        total_len += contig_len
    return total_cov/total_len

# def fasta_from_seq(id_seq_tuple):
#     ''' Builds a fasta sequence out of (id, seq) '''
#     return f'>{id_seq_tuple[0]}\n{id_seq_tuple[1]}'


class LuoGenotyper:
    def __init__(self, input_reads_fastq, allele_db, igblast_output_prefix, spades_output_prefix, hapcut2_output_prefix, all_gene_out, two_cnv_gene_out, kmer_size=21, check_cache=True):
        self.allele_db = allele_db
        self.check_cache = check_cache
        # igblast
        self.input_reads_fastq = input_reads_fastq
        self.igblast_output_prefix = igblast_output_prefix
        # spades
        self.kmer_size = kmer_size
        self.spades_output_prefix = spades_output_prefix
        self.all_gene_out = all_gene_out
        self.two_cnv_gene_out = two_cnv_gene_out
        # # bowtie2
        # self.bowtie2_output_prefix = bowtie2_output_prefix
        # # gatk
        # self.gatk_output_prefix = gatk_output_prefix
        # two_cnv_alleles
        self.hapcut2_output_prefix = hapcut2_output_prefix

        self.functional_gene_list = set(a.gene for a in self.allele_db.alleles if a.is_functional)
        self.gene_read_pairings = {g: [] for g in self.functional_gene_list}
        self.ambiguous_read_mappings = {}

    # --------------- read filtering --------------- 
    def process_mappings(self, read, mappings):
        mappings.sort(key=lambda x: x.identity, reverse=True)
        max_identity = mappings[0].identity
        top_hits = [i.ref for i in mappings if i.identity == max_identity]
        
        UNIQ_FUNC_GENE_SET = ['IGHV6-1', 'IGHV3-9', 'IGHV5-10-1', 'IGHV1-18', 'IGHV3-20', 'IGHV1-24', 'IGHV2-26', 'IGHV3-43', 'IGHV3-43D', 'IGHV1-45', 'IGHV3-49', 'IGHV5-51', 'IGHV1-58', 'IGHV1-69', 'IGHV1-69D', 'IGHV1-69-2', 'IGHV3-72', 'IGHV3-73']
        UNIQ_TOP_HIT_SET = ['IGHV1-2', 'IGHV1-3', 'IGHV7-4-1', 'IGHV2-5', 'IGHV3-7', 'IGHV1-8', 'IGHV3-11', 'IGHV3-13', 'IGHV3-15', 'IGHV3-21', 'IGHV3-23', 'IGHV3-23D', 'IGHV4-30-2', 'IGHV4-30-4', 'IGHV3-30', 'IGHV3-30-3', 'IGHV3-30-5', 'IGHV3-33', 'IGHV4-34', 'IGHV1-46', 'IGHV3-48', 'IGHV3-53', 'IGHV3-64', 'IGHV3-64D', 'IGHV3-66', 'IGHV2-70', 'IGHV2-70D', 'IGHV3-NL1']
        NON_UNIQ_TOP_HIT_SET = ['IGHV4-31', 'IGHV4-39', 'IGHV4-28']
        TOSS_COIN_GENE_SET = [('IGHV3-74', 'IGHV1-OR/16-13*01')]
        SAMPLE_COV_GENE_SET = ['IGHV4-4', 'IGHV4-59', 'IGHV4-61', 'IGHV4-38-2']
        all_genes = UNIQ_FUNC_GENE_SET + UNIQ_TOP_HIT_SET + NON_UNIQ_TOP_HIT_SET + SAMPLE_COV_GENE_SET + ['IGHV3-74']
        GENES_TO_BE_IGNORED = ['IGHV3-35', 'IGHV3-62']
    
        gene_set = set(i.gene for i in top_hits)
        state = 0
        if len(gene_set) == 1:
            top_hit_gene = gene_set.pop()
            
            # discard read if all top hits belong to the same pesudogene
            if top_hit_gene not in self.functional_gene_list:
                state = 1
                log.debug(f'Discarding read {read}: Maps to a pseudogene')
            # else save the mapping
            elif top_hit_gene in UNIQ_FUNC_GENE_SET+UNIQ_TOP_HIT_SET:
                state = 2
                self.gene_read_pairings[top_hit_gene].append(read)
            else:
                state = 3
                assert top_hit_gene in all_genes+GENES_TO_BE_IGNORED, f'Gene unknown to Luo: {top_hit_gene}'
                # log.warn(f'WARNING: Uniquely mapped: {top_hit_gene}')
                self.gene_read_pairings[top_hit_gene].append(read)
        else:
            # non-unique top hits
            for i in NON_UNIQ_TOP_HIT_SET:
                if i in gene_set:
                    state = 4
                    self.gene_read_pairings[i].append(read)
            if state != 4:
                # |gene_set| > 1, read is unassigned
                for func_gene, pseudogene in TOSS_COIN_GENE_SET:
                    if func_gene in gene_set and pseudogene in gene_set:
                        state = 5
                        if toss_fair_coin():
                            self.gene_read_pairings[i].append(read)
                # SAMPLE_COV_GENE_SET
                for gene in AMBIGUOUS_GENE_SET:
                    if gene in gene_set:
                        top_gene_mappings = [m for m in mappings if m.ref.gene == gene and m.identity == max_identity]
                        try:
                            self.ambiguous_read_mappings[read][gene] = top_gene_mappings
                        except KeyError:
                            self.ambiguous_read_mappings[read] = {gene: top_gene_mappings}

    def ambiguous_gene_proportional_assignment(self):
        ''' sample and assign '''
        cov_genes = {g: [0] * max(len(i) for i in self.allele_db.alleles if i.gene == g) for g in AMBIGUOUS_GENE_SET}
    
        for read_dict in self.ambiguous_read_mappings.values():
            for gene, mappings in read_dict.items():
                m = max(mappings, key=lambda x: x.s_end - x.s_start) # get the mapping with the longest alignment
                for i in range(m.s_start-1, m.s_end):
                    cov_genes[gene][i] += 1
                #cov = [(1 if (i >= m.s_start-1 and i < m.s_end) else 0) for i in range(len(cov_genes[gene]))]
                #cov_genes[gene] += cov
    
        # remove genes which have 0 coverage
        remove_gene_ids = [k for k, v in cov_genes.items() if not sum(v)]
        for k in remove_gene_ids:
            _ = cov_genes.pop(k)
    
        for i, (read, mappings) in enumerate(self.ambiguous_read_mappings.items()):
            # if not i%100:
            #     log.warn(f'Processing reads {i}/{len(self.ambiguous_read_mappings)}')
            # process each read
            cov_proportions = {g: get_cov_proportion(cov_genes[g], mappings[g]) for g in mappings.keys()}
            genes, weights = zip(*cov_proportions.items())
            gene_choice = random.choices(list(genes), weights=list(weights), k=1)[0]
            self.gene_read_pairings[gene_choice].append(read)

    def output_reads_per_gene(self):
        reads = {record.id: record for record in SeqIO.parse(self.input_reads_fastq, "fastq")}
        for gene, list_of_reads in self.gene_read_pairings.items():
            sequences = [reads[r] for r in list_of_reads]
            if self.check_cache:
                if not os.path.exists(f'{self.igblast_output_prefix}/{gene}.fq'):
                    SeqIO.write(sequences, f'{self.igblast_output_prefix}/{gene}.fq', 'fastq')
                else:
                    log.warn('Reads for at least one gene are already written. Not running...')
                    break
            else:
                SeqIO.write(sequences, f'{self.igblast_output_prefix}/{gene}.fq', 'fastq')

    def extract_igblast_mappings(self):  # igblast-trimmed.tsv, input_reads_fastq, allele_db, output_prefix):
        ''' igblast-trimmed.tsv should contain only lines starting with V, others are ignored '''
        rev = 'reversed|'
        old_read_name = ''
        with open(os.path.join(self.igblast_output_prefix, 'igblast-trimmed.tsv'), 'r') as f:
            for l in f:
                # create an igblast_mapping
                mapping = igblast_mapping(l)
    
                read_name = mapping.query_id
                if read_name != old_read_name:
                    if old_read_name:
                        self.process_mappings(old_read_name, list_of_mappings)
                    list_of_mappings = []
                list_of_mappings.append(mapping)
                old_read_name = read_name
            # the last pair
            self.process_mappings(old_read_name, list_of_mappings)
    
        self.ambiguous_gene_proportional_assignment()
    
        # remove genes which have 0 reads assigned
        unassigned_genes = [g for g, v in self.gene_read_pairings.items() if not v]
        _ = [self.gene_read_pairings.pop(g) for g in unassigned_genes]
    
        # output reads per gene
        self.output_reads_per_gene()
    
        self.igblast_assigned_genes = sorted(list(self.gene_read_pairings.keys()))

    # --------------- contig assembly --------------- 
    def run_spades(self, threads=1):
        ''' Outputs assmebled contigs from the filtered reads '''
        # spades.py –k 21 –-careful –s IGHV6-1.fastq –o contigs/IGHV6-1
        if not os.path.exists(self.spades_output_prefix):
            os.makedirs(self.spades_output_prefix, exist_ok=True)
        genes = self.igblast_assigned_genes
        output_genes = [1] * len(genes)
        for idx, gene in enumerate(genes):
            out_dir = f'{self.spades_output_prefix}/{gene}'
            # run spades if the directory does not exist
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
                completed_process = sp.run(['spades.py', '-k', f'{self.kmer_size}', '-t', f'{threads}', '--careful', '-s', f'{self.igblast_output_prefix}/{gene}.fq', '-o', f'{out_dir}'], capture_output=True)
                try:
                    completed_process.check_returncode()
                except sp.CalledProcessError as e:
                    output_genes[idx] = 0
                    log.warn(f'SPAdes assembly failed for {gene}: {completed_process.stderr}')
            else:
                # populate output_genes if already run
                if not os.path.exists(f'{out_dir}/contigs.fasta'):
                    output_genes[idx] = 0
                    log.warn(f'SPAdes assembly failed (earlier) for {gene}')
        # return all genes for which assembly has been successful
        self.spades_assembled_genes = [g for idx, g in enumerate(genes) if output_genes[idx]]
    
    def get_cnv_estimate(self, sequencing_depth, read_length, seq_depth_adjustment=0.9):
        ''' copy number from k-mer coverage '''
        kmer_cov = [get_kmer_cov(f'{self.spades_output_prefix}/{g}/contigs.fasta') for g in self.spades_assembled_genes]
        per_base_cov_factor = read_length/(read_length - self.kmer_size + 1)
        self.spades_copy_number_estimates = [c * per_base_cov_factor * 8/(3 * sequencing_depth * seq_depth_adjustment) for c in kmer_cov]
    
    def get_point_cnv_estimate(self):
        ''' Integer point estimates for copy numbers'''
        # TODO: special cases: ["IGHV7-4-1", "IGHV3-7", "IGHV3-49", "IGHV1-69-2"]
        self.point_cnv_estimates = [round(i) for i in self.spades_copy_number_estimates]
        if not os.path.exists(self.all_gene_out):
            with open(self.all_gene_out, 'w') as f:
                print(*[f'{i} {j}' for i,j in zip(self.spades_assembled_genes, self.point_cnv_estimates)], sep='\n', file=f)

    def get_two_copy_genes(self):
        self.two_copy_genes = [g for idx, g in enumerate(self.spades_assembled_genes) if self.point_cnv_estimates[idx] == 2]
        if not os.path.exists(self.two_cnv_gene_out):
            with open(self.two_cnv_gene_out, 'w') as f:
                print(*self.two_copy_genes, sep='\n', file=f)

    # # --------------- Haplotype phasing and allele/SNV calls for 2 copy segments --------------- 
    # def apply_variants(self):
    #     ''' Apply phased variants to the contig '''
    #     for gene in self.two_copy_genes:
    #         all_contigs = {record.id: record for record in SeqIO.parse(f'{self.spades_output_prefix}/{gene}/contigs.fasta', "fasta")}
    #         max_contig_id, max_contig_len = max([(i.id, len(i.seq)) for i in list(all_contigs.values())], key=lambda x: x[1])
    #         contig_id, contig_seq = all_contigs[max_contig_id].id, str(all_contigs[max_contig_id].seq)
    #         homozygous, no_variant_copy = False, False
    #         copy = ['', '']

    #         # check if variants are present
    #         variants = list(pysam.VariantFile(f'{self.gatk_output_prefix}/{gene}/{gene}-vs-ref.vcf').fetch())
    #         variants = [i for i in variants if i.contig == contig_id]
    #         if not variants:
    #             # homozygous: same copy as the contig for both copies
    #             copy = [fasta_from_seq(contig_id, contig_seq)] * 2
    #         else:
    #             # get variants from HapCUT2
    #             with open(f'{self.hapcut2_output_prefix}/{gene}/{gene}-haplotype.txt', 'r') as hap_f:
    #                 haplotype = hap_f.readlines()
    #                 # homozygous if the same variant is present in both copies
    #                 no_variant_copy = not haplotype
    #                 # if no phase file is present, apply all variants to one copy only
    #                 if no_variant_copy:
    #                     copy[0] = fasta_from_seq(contig_id, contig_seq)
    #                     # apply all variants to the second copy
    #                     contig_str_list = list(contig_seq)
    #                     for var in variants:
    #                         # check if reference position matches
    #                         ref, allele = var.alleles
    #                         assert contig_str_list[var.pos - 1] == ref, f'Reference present in contig for gene {gene} is {contig_seq[var.pos - 1]} while VCF claims ref is {ref}'
    #                         # apply variant
    #                         contig_str_list[var.pos - 1] = allele
    #                     # assemble the string again
    #                     copy[1] = fasta_from_seq(contig_id, ''.join(contig_str_list))
    #                 else:
    #                     # heterozygous
    #                     pass
    #         # apply variants and write to a new file

    #         # TODO: check if length of the longest contig is at least as long as the IGHV allele
    #         # assert 

class Variant:
    def __init__(self,  line):
        self.unphased = True
        self.pos, self.ref, self.allele, self.genotype = int(line[1])-1, line[3], line[4], line[-1][:3]
        if '|' in self.genotype:
            self.unphased = False
            self.genotype = [int(x) for x in self.genotype.split('|')]
            assert (len(self.ref), len(self.allele)) == (1, 1), f'Variant is phased but ref, allele == {self.ref, self.allele}'
        else:
            self.genotype = [int(x) for x in self.genotype.split('/')]

class PhasedVariants:
    def __init__(self, variants):
        ''' Input: list of variants outputs: variants on either copy, unphased '''
        _variants = [Variant(i) for i in variants]
        self.phased, self.unphased = [], []
        for i in _variants:
            if i.unphased:
                self.unphased.append(i)
            else:
                self.phased.append(i)
        self.on_copy = [[idx for idx, v in enumerate(self.phased) if v.genotype[copy]] for copy in range(2)]

class AlignStruct:
    def __init__(self, id_str):
        self.start, self.end, self.NM, self.allele_similarity = -1, -1, 0, 0
        self.id = id_str

    def __repr__(self):
        return f'id: {self.id} - start: {self.start} end: {self.end} NM: {self.NM} similarity: {self.allele_similarity}'

    def compute_ed(self):
        # self.allele_similarity = self.end - self.start + 1 - self.NM
        self.edit_distance = self.NM

def align_seq_code(contig, allele):
    ''' 1: match with char
        2: mismatch with char
        3: match with -
        4: only contig has -
        5: only allele has -
    '''
    if (contig, allele) == ('-', '-'):
        return 3
    if contig == '-':
        return 4
    if allele == '-':
        return 5
    if contig == allele:
        return 1
    return 2

class VariantAnalyzer():
    ''' Apply phased variants to the contig '''
    def __init__(self, gene, fasta_prefix, phased_vcf_prefix, allele_db, msa_wrapper=ClustalWWrapper):
        self.gene = gene
        # allele_db
        self.all_alleles = [i for i in allele_db.alleles if i.gene == gene]
        PHASED_VCF_FILE = f'{phased_vcf_prefix}/{gene}/{gene}-haplotype.txt.phased.VCF'
        # contig
        CONTIGS_FILE = f'{fasta_prefix}/{gene}/contigs.fasta'
        all_contigs = {record.id: record for record in SeqIO.parse(CONTIGS_FILE, "fasta")}
        max_contig_id, max_contig_len = max([(i.id, len(i.seq)) for i in list(all_contigs.values())], key=lambda x: x[1])
        self.contig_id, self.contig_seq = all_contigs[max_contig_id].id, str(all_contigs[max_contig_id].seq)
        self.copies = [list(self.contig_seq) for _ in range(2)]
        self.aligned_copies = [''] * 2
        self.closest_alleles = [''] * 2
        self.closest_alleles_list = [[]] * 2

        # get_variants
        self.get_variants(PHASED_VCF_FILE)

        # apply phased variants
        self.apply_phased_variants()

        # run msa to get the closet allele
        self.msa_wrapper = msa_wrapper()
        self.get_msa_output()
        self.get_closest_allele()

    def get_variants(self, vcf_file):
        # manually read the vcf file
        vcf = [i.split() for i in Path(vcf_file).read_text().split('\n') if i.strip() and not i.startswith('#')]
        vcf = [i for i in vcf if i[0] == self.contig_id]
        self.variants = PhasedVariants(vcf)

    def apply_phased_variants(self):
        for c in range(2):
            for p_idx in self.variants.on_copy[c]:
                variant = self.variants.phased[p_idx]
                # check if ref is the same
                assert self.copies[c][variant.pos] == variant.ref, f'Reference present in contig for gene {gene} is {self.copies[c]} while VCF claims ref is {variant.ref}'
                # apply variant
                self.copies[c][variant.pos] = variant.allele
        self.copy_records = [SeqRecord(self.contig_id, ''.join(self.copies[c])) for c in range(2)]

    def get_msa_output(self):
        self.msa_output = [''] * 2
        # align the sequences
        for c in range(2):
            # fasta = fasta_from_seq(*zip(*(self.all_alleles + [(self.contig_id, ''.join(self.copies[c]))])))
            seqs = self.all_alleles + [self.copy_records[c]]
            aligned_seqs = self.msa_wrapper.align(seqs)
            self.msa_output[c] = list(aligned_seqs)
            self.aligned_copies[c] = self.msa_output[c].pop([i.id for i in self.msa_output[c]].index(self.contig_id[:MAX_ID_LENGTH]))
        # TODO unphased
        if self.variants.unphased:
            pass

    def get_closest_allele(self):
        for c in range(2):
            genotype_copy = self.aligned_copies[c]
            align_struct = {i.id: AlignStruct(i.id) for i in self.msa_output[c]}
            idx_last_contig_char = idx_last_allele_char = -1
            for allele in self.msa_output[c]:
                _align_struct = align_struct[allele.id]
                for i in range(len(genotype_copy.seq)):
                    code = align_seq_code(genotype_copy.seq[i], allele.seq[i])
                    if _align_struct.start == -1:
                        if code in [1, 2, 4]:
                            _align_struct.start = i
                    if _align_struct.start >= 0:
                        if code == 1:
                            # match with char
                            _align_struct.end = i
                            idx_last_contig_char = idx_last_allele_char = i
                        elif code == 2:
                            # mismatch with char
                            _align_struct.end = i
                            _align_struct.NM += 1
                            idx_last_contig_char = idx_last_allele_char = i
                        elif code == 4:
                            # only contig has -
                            _align_struct.end = i
                            _align_struct.NM += 1
                            idx_last_allele_char = i
                        elif code == 5:
                            # only allele has -
                            _align_struct.NM += 1
                            idx_last_contig_char = i
                    if i == len(genotype_copy.seq) - 1:
                        if code == 5:
                            _align_struct.NM -= (idx_last_contig_char - idx_last_allele_char)
                        elif code == 3:
                            if idx_last_contig_char > idx_last_allele_char:
                                _align_struct.NM -= (idx_last_contig_char - idx_last_allele_char)
                # compute allele_similarity
                _align_struct.compute_ed()
            # get closest_alleles
            # alleles = sorted(list(align_struct.values()), key=lambda x: x.allele_similarity, reverse=True)
            alleles = sorted(list(align_struct.values()), key=lambda x: x.edit_distance)
            # self.closest_alleles[c] = alleles[0]
            self.closest_alleles[c] = [i for i in alleles if i.edit_distance == alleles[0].edit_distance]
            self.closest_alleles_list[c] = alleles

class TwoCNVPhaser():
    ''' Get phased variants and output alleles for each gene '''
    def __init__(self, allele_output, two_cnv_gene_out, fasta_prefix, phased_vcf_prefix, allele_db):
        self.two_cnv_genes = [i for i in Path(two_cnv_gene_out).read_text().split() if i in TWO_COPY_GENES]
        self.gene_variants = {i: VariantAnalyzer(i, fasta_prefix, phased_vcf_prefix, allele_db) for i in self.two_cnv_genes}
        self.output_alleles(allele_output)

    def output_alleles(self, allele_output):
        write_val = [f'{gene}\t{",".join([i.id for i in val.closest_alleles[0]])}\t{val.closest_alleles[0][0].edit_distance}\t{",".join([i.id for i in val.closest_alleles[0]])}\t{val.closest_alleles[1][0].edit_distance}' for gene, val in self.gene_variants.items()]
        #write_val = [f'{gene}\t{val.closest_alleles[0].id}\t{val.closest_alleles[0].edit_distance}\t{val.closest_alleles[1].id}\t{val.closest_alleles[1].edit_distance}' for gene, val in self.gene_variants.items()]
        write_val = 'Gene\tAllele_copy1\tEdit_distance1\tAllele_cop2\tEdit_distance2\n' + '\n'.join(write_val)
        Path(allele_output).write_text(write_val)

