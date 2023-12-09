import os

class DirController:
    def __init__(self):
        self.cwd = os.getcwd()
        if not os.path.isdir(os.path.join(self.cwd, 'functional_testing_data')):
            os.mkdir(os.path.join(self.cwd, 'functional_testing_data'))
        self.data_storage = os.path.join(self.cwd, 'functional_testing_data')
        self.create_test_temp()
        
    def create_dir_database(self, database):
        os.mkdir(os.path.join(self.data_storage, database))
    
    def create_org_dir(self, species, org_filename):
        """
        species: name of species
        org_filename: filename of organism
        
        """
        species_dir = os.path.join(self.data_storage, species, org_filename)
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)
        else:
            pass
        return species_dir
    
    def create_test_temp(self):
        temp_dir = os.path.join(self.cwd, 'test_temp')
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
        return temp_dir
    
    
    