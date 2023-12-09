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
        
        
    def add_imports(self):
        self.Markdown.add_header("Imports",)
        self.Markdown.add_codeblock()
        self.Markdown.add_codeline("import subprocess")
        self.Markdown.add_codeline("import json")
        self.Markdown.add_codeline("import pandas as pd")  
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()  
    
    def add_initial_command(self, output_json_base):
        self.Markdown.add_header(f"Software Command for {self.species}", level = "##")
        self.Markdown.add_codeblock()
        code = f"subprocess.run({self.software_command})"
        self.Markdown.add_codeline(code)
        code_file = f"json_file = {output_json_base}"
        self.Markdown.add_codeline(code_file)
        json_str = "with open(json_file, 'r') as f:"
        self.Markdown.add_codeline(json_str)
        self.Markdown.add_code_intend(code = "json_data_true = json.load(f)")
        code_json_new = f"json_file_new = {self.output_json_file}"
        self.Markdown.add_codeline(code_json_new)
        json_str_new = "with open(json_file_new, 'r') as f:"
        self.Markdown.add_codeline(json_str_new)
        self.Markdown.add_code_intend(code = "json_data_new = json.load(f)")
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
        self.Markdown.add_header(f"Test seq_region for {self.species}", level = "###")
        self.Markdown.add_codeblock()
        codeline = "seq_region = json_data_true['seq_region']"
        self.Markdown.add_codeline(codeline)
        code_table = "seq_regions_table_true = pd.DataFrame(seq_region)"
        self.Markdown.add_codeline(code_table)
        code_new_region = "seq_region_new = json_data_new['seq_region']"
        self.Markdown.add_codeline(code_new_region)
        code_table_new = "seq_regions_table_new = pd.DataFrame(seq_region_new)"
        self.Markdown.add_codeline(code_table_new)
        test_phenotypes = "assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new['phenotypes'].values.tolist(), 'Phenotypes are not equal'"
        self.Markdown.add_codeline(test_phenotypes)
        test_alignment_length = "assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new['alignment_length'].values.tolist(), 'Alignment length for {self.species} is not equal'"
        self.Markdown.add_codeline(test_alignment_length)
        test_query_start = "assert seq_regions_table_true.loc['query_start_pos'].values.tolist() == seq_regions_table_new['query_start_pos'].values.tolist(), 'Query start for {self.species} is not equal'"
        self.Markdown.add_codeline(test_query_start)
        test_query_end = "assert seq_regions_table_true.loc['query_end_pos'].values.tolist() == seq_regions_table_new['query_end_pos'].values.tolist(), 'Query end for {self.species} is not equal'"
        self.Markdown.add_codeline(test_query_end)
        test_ref_start = "assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new['ref_start_pos'].values.tolist(), 'Ref start for {self.species} is not equal'"
        self.Markdown.add_codeline(test_ref_start)
        test_ref_end = "assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new['ref_end_pos'].values.tolist(), 'Ref end for {self.species} is not equal'"
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
        test_key = "assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'" # compares the first element of key (before first ;) of both dataframes"
        self.Markdown.add_codeline(test_key)
        test_mut_pos = f"assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_codeline(test_mut_pos)
        test_mut_end = f"assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_codeline(test_mut_end)
        test_mut_codon = f"assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_codeline(test_mut_codon)
        test_ref_codon = f"assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for {mut_type} and {self.species} is not equal'"
        self.Markdown.add_codeline(test_ref_codon)
        self.Markdown.add_empty_line()
        self.Markdown.close_codeblock()
        
    def add_contig_test(self, filepath_fasta):
        self.Markdown.add_header(f"Test if resistant gene is same after splitting in two contigs."  , level = "###")
        self.Markdown.add_codeblock()
        code_read_fasta = f"with open('{filepath_fasta}', 'r') as f:"
        self.Markdown.add_codeline(code_read_fasta)
        self.Markdown.add_code_intend(code = "fasta = f.read()")
        code_variations = "seq_variations_true = json_data_true['seq_variations']"
        self.Markdown.add_codeline(code_variations)
        init_class = "ModFasta = FastaModifier()"
        self.Markdown.add_codeline(init_class)
        self.Markdown.add_codeline(f"both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, {filepath_fasta})")
        self.Markdown.add_codeline(f"ModFasta.create_fasta('new_fasta.fa', both_seqs, both_headers, )")
        code = f"subprocess.run({self.software_command})"
        self.Markdown.add_codeline(code)
        code_json_new = f"json_file_new = {self.output_json_file}"
        self.Markdown.add_codeline(code_json_new)
        json_str_new = "with open(json_file_new, 'r') as f:"
        self.Markdown.add_codeline(json_str_new)
        self.Markdown.add_code_intend(code = "json_data_new = json.load(f)")
        code_regions = "seq_region_new = json_data_new['seq_region']"
        self.Markdown.add_codeline(code_regions)
        code_pandas_true = "seq_variations_table_true = pd.DataFrame(seq_variations_true)"
        self.Markdown.add_codeline(code_pandas_true)
        code_pandas_new = "seq_variations_table_new = pd.DataFrame(seq_variations_new)"
        self.Markdown.add_codeline(code_pandas_new)
        query_end = "assert seq_variations_table_true.loc['query_end_pos'].values.tolist() == seq_variations_table_new.loc['query_end_pos'].values.tolist(), 'Query end position is not equal'"
        self.Markdown.add_codeline(query_end)
        query_start = "assert seq_variations_table_true.loc['query_start_pos'].values.tolist() == seq_variations_table_new.loc['query_start_pos'].values.tolist(), 'Query start position is not equal'"
        self.Markdown.add_codeline(query_start)
        ref_end = "assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), 'Ref end position is not equal'"
        self.Markdown.add_codeline(ref_end)
        ref_start = "assert seq_variations_table_true.loc['ref_start_pos'].values.tolist() == seq_variations_table_new.loc['ref_start_pos'].values.tolist(), 'Ref start position is not equal'"
        self.Markdown.add_codeline(ref_start)
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