"""
Sequence of events:

    Read config file
    Read stems
    Read matures
    For each stem:
        if has corresponding mature:
            blast stem to genome
            for each match:
                blast mature to match
                if results not in panda dataframe:
                    add results to panda dataframe
    write dataframe to file

"""
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import config


class MirFinder():

    def __init__(self):
        self.genome_database = config.genome_db
        self.stem_file = config.stemloop_query_file
        self.mature_file = config.mature_query_file
        self.species_abreviation = config.species_abreviation
        self.e_value_threshold = config.e_value_threshold

    def read_reference_stems(self, stem_file=None):
        if stem_file is None:
            stem_file = self.stem_file

        reference_stems = []
        for record in SeqIO.parse(stem_file, 'fasta'):
            reference_stems.append(record)
        return reference_stems

    def read_reference_matures(self, mature_file=None):
        if mature_file is None:
            mature_file = self.mature_file

        reference_matures = []
        for record in SeqIO.parse(mature_file, 'fasta'):
            reference_matures.append(record)
        return reference_matures

    def get_mir_name(self, fasta_entry):
        return fasta_entry.id

    def get_corr_reference_matures(self, stem_entry, reference_matures):
        corr_reference_matures = []
        stem_name = self.get_mir_name(stem_entry)
        mature_stem_names = [stem_name + '-3p', stem_name + '-5p']
        for ref_mature in reference_matures:
            ref_mature_name = self.get_mir_name(ref_mature)
            if ref_mature_name.lower() in mature_stem_names:
                corr_reference_matures.append(ref_mature)
        return corr_reference_matures

    def blast_stem_vs_genome(self, stem_entry, genome_database=None,
                             stem_file=None, blast_output_file=None):
        potential_stems = []

        if genome_database is None:
            genome_database = self.genome_database

        if stem_file is None:
            stem_file = './tmp/stem-' + stem_entry.id + '.fasta'

        if blast_output_file is None:
            blast_output_file = '/tmp/blast-stem-gen-' + stem_entry.id + '.xml'

        SeqIO.write(stem_entry, stem_file, 'fasta')
        blastn_cline = NcbiblastnCommandline(cmd='blastn', query=stem_file,
                                             db=self.genome_database,
                                             out=blast_output_file,
                                             outfmt=5, word_size=16)

        blastn_cline()

        results_handle = open(blast_output_file)
        blast_records = NCBIXML.parse(results_handle)

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.e_value_threshold:
                        potential_stem_id = self.species_abreviation + \
                                            stem_entry.id[3:] + \
                                            '_potential_stem_' + \
                                            alignment.hit_id + '_' + \
                                            str(hsp.query_start) + \
                                            '_' + str(hsp.query_end)

                        record = SeqRecord(Seq(hsp.sbjct.replace('-', '')),
                                           id=potential_stem_id)

                        potential_stems.append(record)

        return potential_stems

    def blast_mature_vs_potential_stem(self, mature_record,
                                       potential_stem_record,
                                       mature_record_file=None,
                                       pot_stem_record_file=None,
                                       multi_blast_output_file=None):
        potential_matures = []

        if mature_record_file is None:
            mature_record_file = './tmp/mature-' + mature_record.id + '.fasta'

        if pot_stem_record_file is None:
            pot_stem_record_file = '/tmp/pot-stem-' + potential_stem_record.id + '.fasta'

        if multi_blast_output_file is None:
            multi_blast_output_file = '/tmp/multi-blast-' + \
                                      mature_record.id + \
                                      '-' + potential_stem_record.id + '.xml'

        SeqIO.write(mature_record, mature_record_file, 'fasta')
        SeqIO.write(potential_stem_record, pot_stem_record_file, 'fasta')

        multi_blast = NcbiblastnCommandline(cmd='blastn',
                                            query=pot_stem_record_file,
                                            subject=mature_record_file,
                                            out=multi_blast_output_file,
                                            outfmt=5, word_size=16)

        multi_blast()
        results_handle = open(multi_blast_output_file)
        multi_blast_records = NCBIXML.parse(results_handle)
        for multi_blast_record in multi_blast_records:
            for alignment in multi_blast_record.alignments:
                for hsp in alignment.hsps:
                    potential_mature_id = self.species_abreviation + \
                                          mature_record.id[3:]

                    record = SeqRecord(Seq(hsp.sbjct.replace('-', '')),
                                       id=potential_mature_id)
                    potential_matures.append(record)
        return potential_matures

    def add_to_dataframe(self, potential_matures, corr_ref_mature):
        # TODO: Check if potential_mature is already in df, if not add to df
        pass

    def write_potential_matures(self, potential_matures_df):
        # TODO: write potential_matures_df to file
        pass


if __name__ == '__main__':
    mf = MirFinder()
    reference_stems = mf.read_reference_stems()
    reference_matures = mf.read_reference_matures()
    for ref_stem in reference_stems:
        corr_reference_matures = mf.get_corr_reference_matures(ref_stem, reference_matures)
        if len(corr_reference_matures) > 0:
            potential_stems = mf.blast_stem_vs_genome(ref_stem)
            for pot_stem in potential_stems:
                for corr_ref_mature in corr_reference_matures:
                    potential_matures = mf.blast_mature_vs_potential_stem(corr_ref_mature, pot_stem)
                    for i in potential_matures:
                        print(i.id)
                        print(i.seq)
                    # mf.add_to_datafram(potential_matures, corr_ref_mature)
    # mf.write_potential_matures(potential_matures_df)
