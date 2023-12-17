import os

class DirController:
    def __init__(self):
        """
        Initializes the DirController instance.
        Sets the current working directory, creates a 'functional_testing_data' directory if it does not exist,
        and sets the 'data_storage' path. Calls create_test_temp method.
        """
        
        self.cwd = os.getcwd()
        if not os.path.isdir(os.path.join(self.cwd, 'functional_testing_data')):
            os.mkdir(os.path.join(self.cwd, 'functional_testing_data'))
        self.data_storage = os.path.join(self.cwd, 'functional_testing_data')
        self.create_test_temp()
        
    def create_dir_database(self, database):
        """
        Creates a directory for a specific database within the 'functional_testing_data' directory.
        
        Parameters:
            database (str): Name of the database.
        """
        os.mkdir(os.path.join(self.data_storage, database))
    
    def create_org_dir(self, species, org_filename):
        """
        Creates a directory for a specific species and organism within the 'functional_testing_data' directory.
        
        Parameters:
            species (str): Name of the species.
            org_filename (str): Filename of the organism.
        
        Returns:
            str: Path to the created directory.
        """
        species_dir = os.path.join(self.data_storage, species, org_filename)
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)
        else:
            pass
        return species_dir
    
    def create_test_temp(self):
        """
        Creates a 'test_temp' directory in the current working directory if it does not exist.
        
        Returns:
            str: Path to the created directory.
        """
        temp_dir = os.path.join(self.cwd, 'test_temp')
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
        return temp_dir
    
    
    