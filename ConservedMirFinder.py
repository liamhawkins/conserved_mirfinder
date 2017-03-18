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
                add results to list

"""
from Bio import SeqIO

class MirFinder():

    def __init__(self):
        pass

    def read_reference_stems(self, stem_file):
        reference_stems = []
        for record in SeqIO.parse(stem_file, 'fasta'):
            reference_stems.append(record)
        return reference_stems

    def read_reference_matures(self, mature_file):
        reference_matures = []
        for record in SeqIO.parse(mature_file, 'fasta'):
            reference_matures.append(record)
        return reference_matures

    def get_mir_name(self, fasta_entry):
        # TODO: Pass in a fasta entry and return the name of the miRNA
        # return mir_name
        pass

    def get_corr_reference_matures(self, stem_entry, reference_matures):
        corr_reference_matures = []
        stem_name = self.get_mir_name(self, stem_entry)
        for ref_mature in reference_matures:
            if stem_name == self.get_mir_name(self, ref_mature):
                corr_reference_matures.append(ref_mature)
        return corr_reference_matures

    def blast_stem_vs_genome(self, stem_entry, genome_database, output_file):
        # TODO: Blast stem entry against genome and return potential_stems
        # return potential_stems
        pass

    def blast_mature_vs_potential_stem(self, mature_entry, potential_stem, output_file):
        # TODO: Blast mature entry against potential_stems and return potential_matures
        # return potential_matures 
        pass

    def write_potential_matures(self, potential_matures, corr_ref_mature):
        # TODO: write potential matures with corresponding reference mature to file
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
                    mf.write_potential_matures(potential_matures, corr_ref_mature)
