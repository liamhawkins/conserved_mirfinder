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
        stem_name = self.get_mir_name(self, stem_entry)
        mature_stem_names = [stem_name + '-3p', stem_name + '-5p']
        for ref_mature in reference_matures:
            ref_mature_name = self.get_mir_name(self, ref_mature)
            if ref_mature_name in mature_stem_names:
                corr_reference_matures.append(ref_mature)
        return corr_reference_matures

    def blast_stem_vs_genome(self, stem_entry, genome_database, output_file):
        '''
        TODO: Write stem_entry to file, blast stem_entry_file against
            genome_database, read in output_file to potential_stems list
            and return potential_stems containing fasta entries
        '''
        # return potential_stems
        pass

    def blast_mature_vs_potential_stem(self, mature_entry, potential_stem, output_file):
        '''
        TODO: Write mature_entry and potential_stem to file, blast
            mature_entry against potential_stem, read in output_file
            to potential_matures list and return potential_matures
        '''
        # return potential_matures 
        pass

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
            potential_stems = mf.blast_stem_vs_genome(ref_stem, genome_database, output_file)
            for pot_stem in potential_stems:
                for corr_ref_mature in corr_reference_matures:
                    potential_matures = mf.blast_mature_vs_potential_stem(corr_ref_mature, pot_stem, output_file)
                    mf.add_to_datafram(potential_matures, corr_ref_mature)
    mf.write_potential_matures(potential_matures_df)
