from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

QUERY_NAME = 'mmu-miR-24-1'
STEMLOOP_OUT_FILE = './output/stemloop-blast/stemloop-blast-' + QUERY_NAME + '.xml'
STEMLOOP_QUERY_FILE = './queries/' + QUERY_NAME + '.fasta'
MATURE_QUERY_FILE = './queries/mature/' + QUERY_NAME + '.fasta'
DB_FILES = './db/xen/xen'
E_VALUE_THRESH = 0.04


def blast_stemloops_vs_genome(query, db, out):
    blastn_cline = NcbiblastnCommandline(cmd="blastn", out=out, outfmt=5, query=query, db=db, word_size=16)
    blastn_cline()


def blast_mature_vs_potential_stemloops(query, subject, out):
    blastn_cline = NcbiblastnCommandline(cmd='blastn', out=out, outfmt=5, query=query, subject=subject, word_size=10)
    blastn_cline()

blast_stemloops_vs_genome(STEMLOOP_QUERY_FILE, DB_FILES, STEMLOOP_OUT_FILE)

result_handle = open(STEMLOOP_OUT_FILE)
blast_records = NCBIXML.parse(result_handle)
# blast_records = list(blast_records)

for blast_record in blast_records:
    seq_list = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                # print('****Alignment***')
                # print('Sequence:', alignment.title)
                # print('Length:', alignment.length)
                # print('e value:', hsp.expect)
                # print(hsp.query)
                # print(hsp.match)
                # print(hsp.sbjct)
                # print('xen hp for', blast_record.query[0:-10] + ':', hsp.sbjct.replace('-', ''))
                seq_list.append(SeqRecord(Seq(hsp.sbjct.replace('-','')), id=alignment.title))
                phairpins_output_file = './output/potential-hairpins/' + blast_record.query + '-potential-hairpins.fasta'
                SeqIO.write(seq_list, phairpins_output_file, 'fasta')
                blast_mature_vs_potential_stemloops(phairpins_output_file, MATURE_QUERY_FILE, 'mature.xml')

                mature_handle = open('mature.xml')
                mature_records = NCBIXML.parse(mature_handle)

                for mature_record in mature_records:
                    for mature_alignment in mature_record.alignments:
                        for mature_hsp in mature_alignment.hsps:
                            if mature_hsp.align_length > (mature_alignment.length -3 ):
                                print(mature_alignment.title)
                                print(mature_hsp.query)
                                print(mature_hsp.match)
                                print(mature_hsp.sbjct)
