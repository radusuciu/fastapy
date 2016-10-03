"""Simple methods for reading and writing fasta files."""

import os

LINE_LENGTH = 80
ENTRY_DELIM_CHAR = '>'
HEADER_SEP_CHAR = '|'

VALID_NUCLEIC_ACIDS = 'ACGTNUKSYMWRBDHV-'
VALID_AMINO_ACIDS = 'APBQCRDSETFUGVHWIYKZLXMN*-'

NUCLEIC_DELCHARS = str.maketrans({ord(c): None for c in VALID_NUCLEIC_ACIDS})
AA_DELCHARS = str.maketrans({ord(c): None for c in VALID_AMINO_ACIDS})


def read_file(filename):
    """Read a multi-sequence fasta file into a list of Sequence instances."""
    fasta_file = FastaFile(filename)
    return fasta_file.read()


def write_file(filename, sequences):
    """Write a list of Sequence instances to disk in the fasta format."""
    fasta_file = FastaFile(filename)
    return fasta_file.write(sequences)


class FastaFile():
    """Simple class for interacting with multi-sequence fasta files."""

    def __init__(self, filename):
        """Initialize with filename."""
        self.filename = filename

    def read(self):
        """Get Sequences from fasta file."""
        with open(self.filename, 'r') as f:
            data = f.read().split(ENTRY_DELIM_CHAR)
            raw = filter(None, data)

        return [Sequence.validate_fasta(ENTRY_DELIM_CHAR + f) for f in raw]

    def write(self, sequences):
        """Write a Sequence or list of Sequences in fasta format."""
        # take single instances as arguments but conver to list for convenience
        if isinstance(sequences, Sequence):
            sequences = [sequences]

        with open(self.filename, 'w') as f:
            f.writelines(seq.fasta + '\n' for seq in sequences)

    def __repr__(self):
        return 'FastaFile(filename={})'.format({
            self.filename
        })


class Sequence():
    """Container for sequence information."""

    def __init__(self, header, sequence):
        """Initialize new Sequence with header and raw sequence info."""
        self.header = header
        self.sequence = sequence
        self._fasta = None

    @property
    def fasta(self):
        """Get fasta representation of sequence."""
        if self._fasta is None:
            self._fasta = self._to_fasta()
        return self._fasta

    @fasta.setter
    def fasta(self, fasta_sequence):
        """Get Sequence instance from fasta sequence."""
        sequence = self.validate(fasta_sequence)

        if not sequence:
            raise InvalidFastaSequenceException

        self.header = sequence.header
        self.sequence = sequence.sequence

    def _to_fasta(self):
        """Output sequence to fasta format."""
        header = ENTRY_DELIM_CHAR + self.header
        sequence = self.split_text(self.sequence)

        return '{}\n{}'.format(
            header,
            sequence
        )

    @staticmethod
    def split_text(string='', join_char=os.linesep, line_length=LINE_LENGTH):
        """Utility method for splitting string every n characters."""
        return join_char.join(string[i:i + line_length] for i in range(0, len(string), line_length))

    @staticmethod
    def validate_fasta(entry):
        """Check whether sequence is in valid fasta format."""
        lines = entry.splitlines()

        # a valid fasta file must have a header, followed by a sequence
        if len(lines) < 2:
            return False

        # the first line must begin with a header
        if not lines[0].startswith(ENTRY_DELIM_CHAR):
            return False

        header = lines[0][1:]
        sequence = ''.join(lines[1:])

        def validate_sequence(sequence, valid_chars):
            """Check whether sequence is valid nucleotide or protein sequence."""
            return len(sequence.upper().translate(valid_chars)) == 0

        # if sequence contains exclusively nucleic or exclusively amino acids
        # it is valid
        if validate_sequence(sequence, NUCLEIC_DELCHARS) or validate_sequence(sequence, AA_DELCHARS):
            return Sequence(header, sequence)
        else:
            return False

    def __repr__(self):
        return 'Sequence(header={}, sequence={})'.format(
            self.header,
            self.sequence
        )


class InvalidFastaSequenceException(Exception):
    """Sequence is not a valid fasta sequence."""
