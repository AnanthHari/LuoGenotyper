import os, abc
from os.path import exists
from .common import log, fasta_from_seq
from Bio import AlignIO

class MSAWrapper():
    """
    Attr:
        params (str): 			saved bwa parameters to be used with run
        src (str): 				path / alias to bwa executable
        output_path (str):		mapping save location for caching / not using temp files
    """
    __metaclass__ = abc.ABCMeta

    params = ''			## command arguments
    src = 'clustalw' 	## path / command to access executabke
    output_path =  None
    default_format = None

    def __init__(self, src=None, params=None, output_path=None):
        self.params = params
        if src: self.src = src
        if output_path: self.output_path = output_path

    # def __str__(self):
    #     return(self.src + self.params)

    @staticmethod
    def create_temp_file(write_data=None):
        import tempfile
        result = tempfile.NamedTemporaryFile(delete=True, mode="w")
        if write_data:
            result.seek(0)
            result.write(write_data)
            result.truncate()
            result.flush()
        return result
    
    def align(self, sequences, src=None, params="-type=dna", parser=None, output_path=None, output_format=None, *args, **kwargs):
        """Runs src using os.system. Depends on tool-specific implementation of self.build_command()

        Args:
            src (str)			path to bwa executable. self.src if None
            params (str):		string of bwa parameters. self.params if None
            parser (func(x)):	parser func for bwa stdout result. bwaWrapper.paf_parser if None
            output_path (str):	cache path to save mapping result to
        
        Note:
            read-like requires 'id' and 'seq' attributes

        Returns:
            output: 			result of parser
        """

        
        ## Check type(query), make temp file and write query seqs as needed
        sequences_file = None
        if isinstance(sequences, str):
            if not exists(sequences):
                log.error('Provided query path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
            sequences_path = sequences
        else:
            sequences_file = self.create_temp_file(write_data=fasta_from_seq(*zip(*[(x.id, x.seq) for x in sequences])))
            sequences_path = sequences_file.name

        if not src:
            src = self.src
        if not params:
            params = self.params if self.params else ''
        if not output_path:
            output_path = self.output_path
        if not parser:
            parser = self.biopython_parser	

        # make output file if needed
        if not output_path:
            output_file = self.create_temp_file()
            output_path = output_file.name

        if not output_format:
            output_format = self.default_format

        # run commmand
        command = self.build_command(src, params, sequences_path, output_path, output_format, *args, **kwargs)

        log.debug('Running {}:\n{}'.format(self.src, command))
        print('Running {}:\n{}'.format(self.src, command))
        os.system(command)

        if sequences_file:
            sequences_file.close()

        result = parser(output_path, output_format)
        if output_file:
            output_file.close()
        return result
        
    @abc.abstractmethod
    def build_command(self, src, params, sequences_path, output_path, output_format, *args, **kwargs):
        raise NotImplementedError

    @staticmethod
    def biopython_parser(output_path, format):
        file = AlignIO.read(output_path, format.lower())
        return file

class ClustalWWrapper(MSAWrapper):
    default_format = 'CLUSTAL'
    def build_command(self, src, params, sequences_path, output_path, output_format, *args, **kwargs):
        sequences_path =f"-INFILE={sequences_path}"
        if output_format:
            output_format = f"-OUTPUT={output_format}"
        else:
            output_format = f"-OUTPUT={self.default_format}"
        output_path = f"-OUTFILE={output_path}"
        return ' '.join([x for x in (src, '-align', params, sequences_path, output_path, output_format) if x])


from collections import Counter
def make_consensus(seqs):
    '''Takes a list of strings or objects with .str attribute representing MSA aligned sequnences
    Returns a string of the majority consensus
    Raises TypeError if not string or does not have .seq attribute
    Raises ValueError if seqs are not the same length'''
    if not isinstance(seqs[0], str):
        seqs = [x.seq for x in seqs]

    if any([len(x) != len(seqs[0]) for x in seqs]):
        raise ValueError('Sequences are not the same length!')
    # list(print(zip([list(x) for x in seqs])))
    return ''.join([Counter(p).most_common(1)[0][0] for p in zip(*[list(x) for x in seqs])])

        
