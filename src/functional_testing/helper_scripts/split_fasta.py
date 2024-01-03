from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class FastaModifier:
    
    @staticmethod
    def read_fasta(file_path):
        with open(file_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                header = ">" + record.id
                sequence = str(record.seq)
        return sequence, header

    @staticmethod
    def create_fasta(output_file,sequences, headers):
        records = []
        for header, sequence in zip(headers, sequences):
            seq_record = SeqRecord(Seq(sequence), id=header[1:], description="")
            records.append(seq_record)
        with open(output_file, "w") as fasta_file:
            SeqIO.write(records, fasta_file, "fasta")
            
    @staticmethod
    def calc_positions(json_data_true):
        seq_variations_true = json_data_true['seq_regions']
        seq_variations_true = seq_variations_true[next(iter(seq_variations_true))]
        start_pos = seq_variations_true['ref_start_pos']
        end_pos = seq_variations_true['ref_end_pos']
        half_pos = int((int(end_pos) - int(start_pos))/2)
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        return start_pos, end_pos, half_pos


    def modify_seqs(self, json_data_true, fasta_file):
        start_pos, end_pos, half_pos = self.calc_positions(json_data_true)
        sequence, header = self.read_fasta(fasta_file)
        first_contig = sequence[:half_pos]
        second_contig = sequence[half_pos:end_pos]
        header2 = header + "_split2"
        both_seqs = [first_contig, second_contig]
        both_headers = [header, header2]
        return both_seqs, both_headers


