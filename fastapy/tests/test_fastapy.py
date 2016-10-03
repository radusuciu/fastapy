"""Tests for fastapy reading and writing of files."""

import fastapy.fastapy as fastapy
import unittest
import os

VALID_NUCLEIC_SEQUENCE = fastapy.VALID_NUCLEIC_ACIDS
VALID_PROTEIN_SEQUENCE = fastapy.VALID_AMINO_ACIDS

_dir_path = os.path.dirname(os.path.realpath(__file__))


class TestCase(unittest.TestCase):
    """Implement fastapy test suite using unittest."""

    def test_read_file(self):
        """Test reading of multi-sequence fasta file."""
        test_filename = os.path.join(_dir_path, 'test.fasta')
        sequences = fastapy.read_file(test_filename)
        self.assertEqual(len(sequences), 3)
        self.assertEqual(sequences[0].sequence, 'ACD')

    def test_write_file(self):
        """Test writing of multi-sequence fasta file."""
        temp_filename = os.path.join(_dir_path, 'temp.fasta')
        test_sequence = 'ACD'
        header = 'test'
        sequences = [fastapy.Sequence(header, test_sequence)] * 2
        fastapy.write_file(temp_filename, sequences)

        with open(temp_filename, 'r', newline=os.linesep) as f:
            raw = f.read()

        content_should_be = os.linesep.join([self._build_fasta('>' + header, test_sequence)] * 2) + os.linesep
        self.assertEqual(raw, content_should_be)
        self.addCleanup(lambda: os.remove(temp_filename))

    def test_validate_fasta(self):
        """Test the fasta validation method with valid and invalid sequences."""
        valid_header = '> I am valid!'

        # J is not a valid nucleic acid or amino acid and should thus throw an exception
        invalid_dna_fasta = self._build_fasta(valid_header, VALID_NUCLEIC_SEQUENCE + 'J')
        invalid_protein_fasta = self._build_fasta(valid_header, VALID_PROTEIN_SEQUENCE + 'J')
        self.assertEqual(fastapy.Sequence.validate_fasta(invalid_dna_fasta), False)
        self.assertEqual(fastapy.Sequence.validate_fasta(invalid_protein_fasta), False)

        valid_dna_fasta = self._build_fasta(valid_header, VALID_NUCLEIC_SEQUENCE)
        valid_protein_fasta = self._build_fasta(valid_header, VALID_PROTEIN_SEQUENCE)
        self.assertNotEqual(fastapy.Sequence.validate_fasta(valid_dna_fasta), False)
        self.assertNotEqual(fastapy.Sequence.validate_fasta(valid_protein_fasta), False)

    def test_invalid_headers(self):
        """Test that sequence does not validate when an invalid header is given."""
        invalid_header = 'invalid'
        invalid_fasta = self._build_fasta(invalid_header, fastapy.VALID_NUCLEIC_ACIDS)
        self.assertEqual(fastapy.Sequence.validate_fasta(invalid_fasta), False)

    def _build_fasta(self, header, sequence):
        return '{header}{linesep}{sequence}'.format(
            header=header,
            linesep=os.linesep,
            sequence=sequence
        )
