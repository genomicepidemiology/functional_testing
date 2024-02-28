from .markdown_builder import MarkdownBuilder
import pandas as pd
import os

class MakeTests:
    def __init__(self, output_json_file = None, species = None, software_command = None):
        """
        :param output_json_file: The path to the json file which was created as base true file for comparing this with the new outputs when the software is changed.
        :param species: The species for which the test is created.
        :param software_command: Is the new command for running the software.
        """
        self.Markdown = MarkdownBuilder()
        #self.add_initial_command(species, software_command)
        self.output_json_file = output_json_file
        self.species = species
        self.software_command = software_command
        self.switch = None
        
    def add_switcher(self, software_command):
        if self.switch != None:
            initial_wd = os.getcwd()
            self.Markdown.add_codeline(f"initial_wd = {initial_wd}")
            code = f"os.chdir({self.switch})"
            self.Markdown.add_codeline(code)
            code = f"subprocess.run({software_command})"
            self.Markdown.add_codeline(code)
            completed_process = f"CompletedProcess(args={software_command}, returncode=0)"
            self.Markdown.add_line(completed_process)
            code = f"os.chdir(initial_wd)"
        else:
            self.Markdown.add_codeline(f"subprocess.run({software_command})")
            completed_process = f"CompletedProcess(args={software_command}, returncode=0)"
            self.Markdown.add_line(completed_process)
    
    @staticmethod
    def change_fasta(software_command, filepath_fasta):
        assert "-ifa" in software_command, "Software command must contain -ifa flag"
        ifa_index = software_command.index("-ifa")
        software_command[ifa_index + 1] = filepath_fasta
        return software_command
        
        
    def add_imports(self):
        self.Markdown.add_header("Imports",)
        self.Markdown.add_codeblock()
        self.Markdown.add_codeline("import subprocess")
        self.Markdown.add_codeline("import json")
        self.Markdown.add_codeline("import os")
        self.Markdown.add_codeline("import pandas as pd")          
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()  
        self.add_fasta_modifier()
    
    def add_initial_command(self, output_json_base):
        self.Markdown.add_header(f"Software Command for {self.species}", level = "##")
        self.Markdown.add_codeblock()
        self.add_switcher(self.software_command)
        code_file = f"json_file = '{output_json_base}'"
        self.Markdown.add_codeline(code_file)
        json_str = "with open(json_file, 'r') as f:"
        self.Markdown.add_codeline(json_str)
        self.Markdown.add_code_intend(code = "json_data_true = json.load(f)")
        code_json_new = f"json_file_new = '{self.output_json_file}'"
        self.Markdown.add_codeline(code_json_new)
        json_str_new = "with open(json_file_new, 'r') as f:"
        self.Markdown.add_codeline(json_str_new)
        self.Markdown.add_code_intend(code = "json_data_new = json.load(f)")
        self.Markdown.add_codeline("if not os.path.isdir('test_temp'):")
        self.Markdown.add_code_intend(code = "os.mkdir('test_temp')")
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        

    def test_base_dic(self):
        self.Markdown.add_header(f"Test keys and base structure of json dic", level = "###")
        self.Markdown.add_codeblock()
        code = f"assert json_file != json_file_new, 'Json files are the same'"
        self.Markdown.add_codeline(code)
        code_keys_true = "keys_true = list(json_data_true.keys())"
        self.Markdown.add_codeline(code_keys_true)
        code_keys_new = "keys_new = list(json_data_new.keys())"
        self.Markdown.add_codeline(code_keys_new)
        test_keys = "assert keys_true == keys_new, 'Keys are not equal'"
        self.Markdown.add_codeline(test_keys)
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        
        
    def test_seq_region(self):
        self.Markdown.add_header(f"Test seq_regions for {self.species}", level = "###")
        self.Markdown.add_codeblock()
        codeline = "seq_regions = json_data_true['seq_regions']"
        self.Markdown.add_codeline(codeline)
        code_table = "seq_regions_table_true = pd.DataFrame(seq_regions)"
        self.Markdown.add_codeline(code_table)
        code_new_region = "seq_region_new = json_data_new['seq_regions']"
        self.Markdown.add_codeline(code_new_region)
        code_table_new = "seq_regions_table_new = pd.DataFrame(seq_region_new)"
        self.Markdown.add_codeline(code_table_new)
        test_phenotypes = "assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'"
        self.Markdown.add_codeline(test_phenotypes)
        test_alignment_length = f"assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for {self.species} is not equal'"
        self.Markdown.add_codeline(test_alignment_length)
        # query start pos does not exist in some - not resolved yet (comment out)
       # test_query_start = f"assert seq_regions_table_true.loc['query_start_pos'].values.tolist() == seq_regions_table_new.loc['query_start_pos'].values.tolist(), f'Query start for {self.species} is not equal'"
       # self.Markdown.add_codeline(test_query_start)
       # test_query_end = f"assert seq_regions_table_true.loc['query_end_pos'].values.tolist() == seq_regions_table_new.loc['query_end_pos'].values.tolist(), f'Query end for {self.species} is not equal'"
       # self.Markdown.add_codeline(test_query_end)
        test_ref_start = f"assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for {self.species} is not equal'"
        self.Markdown.add_codeline(test_ref_start)
        test_ref_end = f"assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for {self.species} is not equal'"
        self.Markdown.add_codeline(test_ref_end)
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        

    def test_mutations(self, mut_type):
        if mut_type == "-c":
            mut_type = "chromosomal"
        elif mut_type == "-acq":
            mut_type = "acquired"
        elif mut_type == "-u":
            mut_type = "unknown"
        self.Markdown.add_header(f"Test mutations for {self.species} and mut_type: {mut_type}")
        self.Markdown.add_codeblock()
        code_table_true = "seq_variations_true = json_data_true['seq_variations']"
        self.Markdown.add_codeline(code_table_true)
        code_pandas_true = "seq_variations_table_true = pd.DataFrame(seq_variations_true)"
        self.Markdown.add_codeline(code_pandas_true)
        code_table_new = "seq_variations_new = json_data_new['seq_variations']"
        self.Markdown.add_codeline(code_table_new)
        code_pandas_new = "seq_variations_table_new = pd.DataFrame(seq_variations_new)"
        self.Markdown.add_codeline(code_pandas_new)
        self.Markdown.add_codeline("assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'")
        
        self.Markdown.add_codeline("if seq_variations_table_true.shape[0] > 0:")
        
        test_key = "assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'" # compares the first element of key (before first ;) of both dataframes"
        self.Markdown.add_code_intend(test_key)
        test_mut_pos = f"assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_code_intend(test_mut_pos)
        test_mut_end = f"assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_code_intend(test_mut_end)
        test_mut_codon = f"assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_code_intend(test_mut_codon)
        test_ref_codon = f"assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_code_intend(test_ref_codon)
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        
    def add_contig_test(self, filepath_fasta):
        self.Markdown.add_header(f"Test if resistant gene is same after splitting in two contigs."  , level = "###")
        self.Markdown.add_codeblock()
        code_read_fasta = f"with open('{filepath_fasta}', 'r') as f:"
        self.Markdown.add_codeline(code_read_fasta)
        self.Markdown.add_code_intend(code = "fasta = f.read()")
        init_class = "ModFasta = FastaModifier()"
        self.Markdown.add_codeline(init_class)
        self.Markdown.add_codeline(f"both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '{filepath_fasta}')")
        self.Markdown.add_codeline(f"ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)") # requires import of ModFasta
        # software command must be changed
        software_command = self.change_fasta(self.software_command, "test_temp/new_fasta.fa")
        self.add_switcher(software_command)
        code_json_new = f"json_file_new = '{self.output_json_file}'"
        self.Markdown.add_codeline(code_json_new)
        json_str_new = "with open(json_file_new, 'r') as f:"
        self.Markdown.add_codeline(json_str_new)
        self.Markdown.add_code_intend(code = "json_data_new = json.load(f)")
        codeline = "seq_variations = json_data_true['seq_variations']"
        self.Markdown.add_codeline(codeline)
        code_table = "seq_variations_table_true = pd.DataFrame(seq_variations)"
        self.Markdown.add_codeline(code_table)
        code_test_regions = "if seq_variations_table_true.shape[0] > 0:"
        self.Markdown.add_codeline(code_test_regions)
        code_new = "seq_variations_new = json_data_new['seq_variations']"
        self.Markdown.add_code_intend(code_new)
        code_table_new = "seq_variations_table_new = pd.DataFrame(seq_variations_new)"
        self.Markdown.add_code_intend(code_table_new)
        test_alignment_length = f"assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for {self.species} is not equal'"
        self.Markdown.add_code_intend(test_alignment_length)
        test_ref_start = f"assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for {self.species} is not equal'"
        self.Markdown.add_code_intend(test_ref_start)
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()        

    def add_custom_test(self, ):
        self.Markdown.add_header(f"Custom test", level = "###")
        self.Markdown.add_codeblock()
        codeline = "seq_regions = json_data_true['seq_regions']"
        self.Markdown.add_codeline(codeline)
        code_table = "seq_regions_table_true = pd.DataFrame(seq_regions)"
        self.Markdown.add_codeline(code_table)
        code_test_regions = "if seq_regions_table_true.shape[0] > 0:"
        self.Markdown.add_codeline(code_test_regions)
        code_new = "seq_region_new = json_data_new['seq_regions']"
        self.Markdown.add_code_intend(code_new)
        code_table_new = "seq_regions_table_new = pd.DataFrame(seq_region_new)"
        self.Markdown.add_code_intend(code_table_new)
        test_phenotypes = "assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'"
        self.Markdown.add_code_intend(test_phenotypes)
        test_alignment_length = f"assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for {self.species} is not equal'"
        self.Markdown.add_code_intend(test_alignment_length)
        test_ref_start = f"assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for {self.species} is not equal'"
        self.Markdown.add_code_intend(test_ref_start)
        test_ref_end = f"assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for {self.species} is not equal'"
        self.Markdown.add_code_intend(test_ref_end)
        self.Markdown.add_codeline("else:")
        self.Markdown.add_code_intend(code = "pass")
        code_test_variations = "seq_variations_true = json_data_true['seq_variations']"
        self.Markdown.add_codeline(code_test_variations)
        code_pandas_true = "seq_variations_table_true = pd.DataFrame(seq_variations_true)"
        self.Markdown.add_codeline(code_pandas_true)
        code_len_variations = "if seq_variations_table_true.shape[0] > 0:"
        self.Markdown.add_codeline(code_len_variations)
        code_new_variations = "seq_variations_new = json_data_new['seq_variations']"
        self.Markdown.add_code_intend(code_new_variations)
        code_pandas_new = "seq_variations_table_new = pd.DataFrame(seq_variations_new)"
        self.Markdown.add_code_intend(code_pandas_new)
        ref_end = "assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), 'Ref end position is not equal'"
        self.Markdown.add_code_intend(ref_end)
        ref_start = "assert seq_variations_table_true.loc['ref_start_pos'].values.tolist() == seq_variations_table_new.loc['ref_start_pos'].values.tolist(), 'Ref start position is not equal'"
        self.Markdown.add_code_intend(ref_start)
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
    
        
    def write_markdown(self,filename, filedir = None):
        if filename != None and filedir != None:
            filename = os.path.join(filedir, filename) 
        elif filename == None and filedir != None:
            filename = os.path.join(filedir, "functional_tests.md")
        elif filename != None and filedir == None:
            filename = os.path.join(os.getcwd(), filename)
        else:
            filename = os.path.join(os.getcwd(), "functional_tests.md")
        with open(filename, 'w') as file:
            file.write(self.Markdown.content)
            
            
    def add_fasta_modifier(self,):
        self.Markdown.add_header("Fasta Modifier Import", level = "###")
        self.Markdown.add_codeblock()
        self.Markdown.add_codeline("from Bio import SeqIO")
        self.Markdown.add_codeline("from Bio.SeqRecord import SeqRecord")
        self.Markdown.add_codeline("from Bio.Seq import Seq")
        self.Markdown.add_codeline("class FastaModifier:")
        self.Markdown.add_code_intend("@staticmethod")
        self.Markdown.add_code_intend("def read_fasta(file_path):")
        self.Markdown.add_code_intend(code = "with open(file_path, 'r') as fasta_file:", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "for record in SeqIO.parse(fasta_file, 'fasta'):", intendation="\t\t\t\t")
        self.Markdown.add_code_intend(code = "header = '>' + record.id", intendation="\t\t\t\t\t")
        self.Markdown.add_code_intend(code = "sequence = str(record.seq)", intendation="\t\t\t\t\t")
        self.Markdown.add_code_intend(code = "return sequence, header", intendation="\t\t\t")
        self.Markdown.add_code_intend("@staticmethod")
        self.Markdown.add_code_intend("def create_fasta(output_file,sequences, headers):")
        self.Markdown.add_code_intend(code = "records = []", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "for header, sequence in zip(headers, sequences):", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "seq_record = SeqRecord(Seq(sequence), id=header[1:], description='')", intendation="\t\t\t\t")
        self.Markdown.add_code_intend(code = "records.append(seq_record)", intendation="\t\t\t\t")
        self.Markdown.add_code_intend(code = "with open(output_file + '.fa', 'w') as fasta_file:", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "SeqIO.write(records, fasta_file, 'fasta')", intendation="\t\t\t\t")
        self.Markdown.add_code_intend("@staticmethod")
        self.Markdown.add_code_intend("def calc_positions(json_data_true):")
        self.Markdown.add_code_intend(code = "seq_variations_true = json_data_true['seq_regions']", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "seq_variations_true = seq_variations_true[next(iter(seq_variations_true))]", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "start_pos = seq_variations_true['ref_start_pos']", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "end_pos = seq_variations_true['ref_end_pos']", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "half_pos = int((int(end_pos) - int(start_pos))/2)", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "start_pos = int(start_pos)", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "end_pos = int(end_pos)", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "return start_pos, end_pos, half_pos", intendation="\t\t\t")
        self.Markdown.add_code_intend("def modify_seqs(self, json_data_true, fasta_file):")
        self.Markdown.add_code_intend(code = "start_pos, end_pos, half_pos = self.calc_positions(json_data_true)", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "sequence, header = self.read_fasta(fasta_file)", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "first_contig = sequence[:half_pos]", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "second_contig = sequence[half_pos:end_pos]", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "header2 = header + '_split2'", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "both_seqs = [first_contig, second_contig]", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "both_headers = [header, header2]", intendation="\t\t\t")
        self.Markdown.add_code_intend(code = "return both_seqs, both_headers", intendation="\t\t\t")
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        
    def remove_file_temp(self):
        index_json = self.software_command.index("-j")
        json_path = self.software_command[index_json + 1]
        self.Markdown.add_codeblock()
        self.Markdown.add_codeline(f"os.remove('{json_path}')")
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        