#!/usr/bin/env python
# 786

import os, sys
import pkg_resources
import logbook, logbook.more
from Bio.Seq import Seq

DATA_DIR = 'luo_genotyper.data'
SCRIPT_DIR = 'luo_genotyper'

def resource_path(key, data_dir_path=DATA_DIR):
    """
    Args:
        key (str): resource to be extracted from the ``data'' directory.
        data_dir_path (str): Relative path from the package root directory.
    Returns:
        str: Full path of the resource.
    """
    return pkg_resources.resource_filename(data_dir_path, key)

def script_path(key, script_dir_path=SCRIPT_DIR):
    """
    Args:
        key (str): resource to be extracted from the ``script'' directory.
        script_dir_path (str): Relative path from the package root directory.
    Returns:
        str: Full path of the resource.
    """
    return pkg_resources.resource_filename(script_dir_path, key)

def header(string):
    return '\n\n' + '-'*len(string) + string + '-'*len(string) + '\n\n'

def colorize(text, color='green'):
    return logbook._termcolors.colorize(color, text)

log = logbook.Logger('Cypiripi')
# log.level = logbook.DEBUG

def initialize_logger(debug_log_path='luo_genotyper-debug'):
    LOG_FORMAT = '{record.message}'
    if debug_log_path:
        debug_log_path = debug_log_path+'.log'
        if os.path.exists(debug_log_path):
            os.remove(debug_log_path)
        
        handler = logbook.NestedSetup([logbook.NullHandler(),
                                    logbook.FileHandler(debug_log_path, level='DEBUG', format_string=LOG_FORMAT),
                                    logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='INFO', bubble=True)])
    else:
        handler = logbook.NestedSetup([logbook.NullHandler(),
                                    logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='INFO', bubble=True)])

    handler.push_application()



def log_msg(msg):
    for m in msg:
        f = m[0]
        f(m[1])

mean = lambda x: round(float(sum(x))/len(x), 3)
median = lambda x: sorted(x)[len(x)/2]
def mean_median(x):
    return '{} mean, {} median'.format(mean(x), median(x))


def fasta_from_seq(name, seq):
    ## Input: Name, seq can be str or iterable yielding str
    result = []
    if not isinstance(name, str):
        try:
            for n, s in zip(name, seq):
                result.append('>{}\n{}'.format(n, s))
        except TypeError:
            log.error('Please provide a iterable or string')
            raise TypeError
    else:
        result = ['>{}\n{}'.format(name, seq)]
    return '\n'.join(result)


def get_columns(data):
    col_width = max(len(str(word)) for row in data for word in row) + 2  # padding
    result = []
    for row in data:
        result.append("".join(str(word).ljust(col_width) for word in row))
    return '\n'.join(result)

def print_columns(data):
    print(get_columns(data))


class Read():
    def __init__(self, iden, seq):
        self.id = iden
        self.seq = seq
        self.mappings = []
        self.mappings_dict = {}

    def __getitem__(self, key):
        if isinstance(key, slice):
            indices = range(*key.indices(len(self.seq)))
            return ''.join(self.seq[i] for i in indices)
        return self.seq[key]
    
    def __str__(self):
        return str(self.seq)
    
    def __len__(self):
        return len(self.seq)

    @property
    def reverse_complement(self):
        return str(Seq(self.seq).reverse_complement())
    
    @property
    def is_discarded(self):
        return True if self.discard_var.X > 1 else False
    
    def add_mapping(self, mapping):
        self.mappings.append(mapping)
        self.mappings_dict[self.reference_id_parser(mapping.reference_name)] = mapping
        try:
            self.primary_mapping = self.get_primary_mapping(self.mappings)
        except ValueError:
            pass

    def get_distance(self, allele):
        '''Gets NM tag edit distance from mapping to provided allele'''
        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        return self.mappings_dict[allele_id].get_tag('NM')

    def get_start(self, allele):
        '''Gets reference_start from mapping to provided allele'''

        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        return self.mappings_dict[allele_id].reference_start
    
    def get_end(self, allele):
        '''Gets reference_end from mapping to provided allele'''
        if isinstance(allele, str):
            allele_id = allele
        else:
            allele_id = allele.id
        return self.mappings_dict[allele_id].reference_end
    
    def covers_position(self, pos, allele):
        '''Returns true if mapping to allele covers pos'''
        return True if (self.get_start(allele) <= pos and pos <= self.get_end(allele)) else False

    @staticmethod
    def reference_id_parser(raw_id):
        if '|' in raw_id:
            return raw_id.split('|')[1]
        else:
            return raw_id

    @staticmethod
    def get_primary_mapping(mapping_list):
        '''Takes a list of pysam.AlignedSegment objects, returns the one that is primary mapping
        Raises ValueError if it does not have,or has multiple, primary mappings'''
        primary_mapping = [x for x in mapping_list if not x.is_secondary and not x.is_supplementary]
        if len(primary_mapping) > 1:
            s = "\n".join([str(x)+f'\t Reference = {x.reference_name}' for x in primary_mapping])
            log.error(f'More than 1 primary mapping:\n{s}')
        if len(mapping_list) and not primary_mapping:
            raise ValueError('Provided mappings have no primary mapping'.format(str([str(x) for x in mapping_list])))
        return None if not primary_mapping else primary_mapping[0]


class SeqRecord(Read):
    '''Implemented for backwards compatibility'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    

class NamedFunction():
    def __init__(self, name, f):
        self.f = f
        self.name = name

    def __call__(self, *args, **kwargs):
        return self.f(*args, **kwargs)

    def __str__(self):
        return self.name


class PrintRedirect:
    def __init__(self, output_path=None):
        if not output_path: 
            self.output_path = os.devnull
        else:
            self.output_path = output_path
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(self.output_path, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
