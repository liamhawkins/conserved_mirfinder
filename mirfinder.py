from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

QUERY_NAME = 'mmu-miR-15a'
STEMLOOP_OUT_FILE = './output/stemloop-blast/stemloop-blast-' + QUERY_NAME + '.xml'
STEMLOOP_QUERY_FILE = './queries/stemloop/' + QUERY_NAME + '.fasta'
MATURE_QUERY_FILE = './queries/matures/' + QUERY_NAME + '.fasta'
DB_FILES = './db/xen/xen'
SPECIES_ABV = 'xen'
E_VALUE_THRESH = 0.04


def blast_stemloops_vs_genome(query, db, out):
    blastn_cline = NcbiblastnCommandline(cmd="blastn", out=out, outfmt=5, query=query, db=db, word_size=16)
    blastn_cline()


def blast_mature_vs_potential_stemloops(query, subject, out):
    blastn_cline = NcbiblastnCommandline(cmd='blastn', out=out, outfmt=5, query=query, subject=subject, word_size=19)
    blastn_cline()

def get_mature_sequences():
    # TODO pass in reference stemloop fasta entry, return fasta entry(ies) for reference mature sequences
    pass


blast_stemloops_vs_genome(STEMLOOP_QUERY_FILE, DB_FILES, STEMLOOP_OUT_FILE)
result_handle = open(STEMLOOP_OUT_FILE)
blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    seq_list = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                seq_list.append(SeqRecord(Seq(hsp.sbjct.replace('-', '')), id=alignment.title))

    pstemloops_output_file = './output/potential-stemloops/' + blast_record.query + '-potential-stemloops.fasta'
    SeqIO.write(seq_list, pstemloops_output_file, 'fasta')
    blast_mature_vs_potential_stemloops(pstemloops_output_file, MATURE_QUERY_FILE, STEMLOOP_OUT_FILE)

    mature_handle = open(STEMLOOP_OUT_FILE)
    mature_records = NCBIXML.parse(mature_handle)

    for mature_record in mature_records:
        for mature_alignment in mature_record.alignments:
            for mature_hsp in mature_alignment.hsps:
                query_mir_name = mature_alignment.hit_def
                query_mir_name = query_mir_name.split(' ')
                print('Original miRNA:', mature_alignment.hit_def)
                print('Possible conserved miRNA:', SPECIES_ABV + query_mir_name[0][3:])
                print(mature_hsp.sbjct)
