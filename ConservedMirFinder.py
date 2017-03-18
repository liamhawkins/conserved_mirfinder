"""
Sequence of events:

    Read config file
    For each stemloop:
        if has corresponding mature:
            blast stemloop to genome
            for each match:
                blast mature to match
                add results to list

"""

class MirFinder():

    def __init__(self):
        pass

    def read_stemloops(self, stemloop_file):
        pass
    
    def read_mature(self, mature_file):
        pass

    def blast_stemloops_vs_genome(self, stemloop, genome_database, output_file):
        pass

    def blast_mature_vs_potential_stemloops(self, stemloop, mature, output_file):
        pass


