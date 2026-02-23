import json
import linecache
import logging
import os
import subprocess
import traceback

import MDAnalysis as mda
import numpy as np
import pandas as pd

cwd = os.getcwd() + '/'

logging.basicConfig(level=logging.INFO ,format='%(asctime)s -  %(name)s - %(levelname)s - %(message)s')

dens_st_cond = ["reg","sdk","comp"]


    
class CompoundData:
    
    
    def __init__(self, 
                 info_file : str,
                 opt_type : str = 'reg'):
        
        """ The CompoundData class encapsulates the logic for reading and updating compound data from a file.
            It provides a method to read the data and stores the data as instance variables.
            The class also uses a logging object for logging messages and raises appropriate exceptions when errors occur.

            Args:
                
                info_file (str): The file path of the data file.
                opt_type (str): The optimization type. Can be one of 'reg', 'sdk', or 'it'.

            Attributes:
                
                data_file (str): The file path of the data file.
                opt_type (str): The optimization type.
                logging (Logger): The logging object for logging messages.
                data (dict): A dictionary to store the compound data lists.

            Example:
                
                import logging

                logging.basicConfig(level=logging.INFO)

                # Create an instance of the CompoundData class
                compound_data = CompoundData('data.csv', opt_type='reg')

                # Read and update the compound data from the file
                compound_data.read_compound_data()

                # Access the data lists as instance variables of the compound_data object
                print(compound_data.data['MOLECULE'])
                print(compound_data.data['DENSITY'])"""
        
        
        self.data_file = str(info_file)
        self.opt_type = opt_type
        self.data = {
            'MOLECULE': [],
            'CODE': [],
            'UID': [],
            'TEMPERATURE': [],
            'GROUP': []
        }

        if self.opt_type in dens_st_cond:
            
            self.data['DENSITY'] = []
            self.data['SURFACETENSION'] = []
            
            if self.opt_type == 'reg':
            
                self.data['EV'] = []
        
            elif self.opt_type == 'comp':

                self.data['COMP'] = []

        elif self.opt_type == 'it':
            
            self.data['INTERFACIALTENSION'] = []

    def read_compound_data(self):
        
        """
        Read and update the compound data from the file.

        Reads the data from the file specified during initialization and updates the compound data lists.
        The data is stored in the instance variables of the class.

        Raises:
            
            FileNotFoundError: If the data file specified is not found.
            Exception: If any other error occurs during data reading or updating.
        """
        
        try:
            
            logging.info('Reading data from file: {}'.format(self.data_file))
            df = pd.read_csv(self.data_file, sep='\s+')
            logging.info('Updating property data ....')

            columns_map = {
                'reg': ['Density', 'Surfacetension', 'ev'],
                'comp':['Density','Surfacetension','comp'],
                'sdk': ['Density', 'Surfacetension'],
                'it': ['Interfacialtension']
            }
            
            columns = columns_map.get(self.opt_type, [])

            for index, row in df.iterrows():
                
                uid = row['UID']
                
                if uid not in self.data['UID']:
                
                    self.data['UID'].append(uid)
                    
                    for column, value in row.items():
                        if column != 'UID':
                            self.data[column.upper()].append(value)
                    
                    logging.info('Including data of compound ' + str(uid))

            logging.info('Update done ...')

        except FileNotFoundError as e:
            
            logging.error('Error reading data file: {}'.format(e))
            traceback.print_exc()

        except Exception as e:
            
            logging.error('Error updating property data: {}'.format(e))
            traceback.print_exc()

            
   

## The below class requires an update or a bug clearance in the near future.

class TabulatedPotential:
   
    """
    The TabulatedPotential class creates a table of potential energy values for a given set of exponential parameters.
    The table is stored in a pandas DataFrame object and can be saved to a file.
    
    Args:
        
        update_dir (str): The directory path for saving the updated tables. Default is "TABULATED_POTENTIAL".
        lamb_1 (float): The exponential parameter for h(r) in the potential equation. Default is 8.9901.
        lamb_2 (float): The exponential parameter for g(r) in the potential equation. Default is 6.9347.

    Attributes:
        
        update_dir (str): The directory path for saving the updated tables.
        lamb_1 (float): The exponential parameter for h(r) in the potential equation.
        lamb_2 (float): The exponential parameter for g(r) in the potential equation.
        logging (Logger): The logging object for logging messages.
        

    Example:
        
        tp = TabulatedPotential(lamb_1=12.0, lamb_2 =  6.0)
        tp.update(table='specific', group='ALK_ALK')
    """
    
    def __init__(self, 
                 update_dir : str = "TABULATED_POTENTIAL",
                 lamb_1 : float = 9.0, 
                 lamb_2 : float = 6.0):
        
        self.update_dir = update_dir
        self.lamb_1 = lamb_1
        self.lamb_2 = lamb_2
        self.rnm = np.linspace(0.002, 2.5, num=1250)
        self.r = self.rnm * (10 ** (-9))
        
        self.tabledf = pd.DataFrame(np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]),
                                     columns=['range', 'f(r)', '-f\'(r)', 'g(r)', '-g\'(r)', 'h(r)', '-h\'(r)'])
        
        self.Potentialdf = pd.DataFrame({'range': self.rnm, 
                                         'f(r)': 0.0, 
                                         '-f\'(r)': 0.0, 
                                         'g(r)': -1 / self.rnm ** self.lamb_2,
                                         '-g\'(r)': -self.lamb_2 / (self.rnm ** (self.lamb_2 + 1)),
                                         'h(r)': 1 / self.rnm ** self.lamb_1,
                                         '-h\'(r)': self.lamb_1 / (self.rnm ** (self.lamb_1 + 1))})
        
        self.Potentialtotal = pd.concat([self.tabledf, self.Potentialdf])
        self.roundtodecimals = self.Potentialtotal.round(decimals=7)

    def update(self, 
               group : str = None,
               table : str = 'default'):

        """
        Update and save the tabulated potentials.

        Args:
            
            group (str): The group name for the table. Required when table is set to 'specific'.
            table (str): The type of table to update and save. Can be 'default' or 'specific'. Default is 'default'.

        Raises:
            
            FileNotFoundError: If the specified file path is not found.
            Exception: If any other error occurs during table updating or saving.
        """
        
        try:
            logging.info(f"Updating tabulated potentials with LJ_{self.lamb_1}_{self.lamb_2} for {group} ....")

            if table == 'specific':
                file = f"{self.update_dir}/table_{group}.xvg"
            
            elif table == 'default':
                file = f"{self.update_dir}/table.xvg"
            
            else:
                logging.error(f"Table type {table} not recognized.")
                traceback.print_exc()
                return

            self.roundtodecimals.to_csv(file, sep=' ', header=None, index=False, float_format='%.10E')

            logging.info("Update done.")
            
        except FileNotFoundError as e:
            
            logging.error('Error reading data file: {}'.format(e))
            traceback.print_exc()

        except Exception as e:
            
            logging.error('Error updating table: {}'.format(e))
            traceback.print_exc()

mdp_lists=['minim.mdp','eq.mdp','eq1.mdp','eq2.mdp','md.mdp','eqnvt.mdp','eqnvt1.mdp','mdnvt.mdp']
mdp_lists_it=['IT_equil.mdp','IT_prod.mdp']




#Update mdp files and move them into respective directories.
def update_mdp_with_code(input_data : dict,
                         opt_type : str = 'reg',
                         reason : str = 'homo'):

    """
    Updates MDP (molecular dynamics simulation input files for GROMACS 2019.6 or earlier) with specific codes and groups based on the opt_type and reason.

    Args:
        
        input_data (dict): A dictionary containing input data, including "CODE" and "GROUP" lists.
        opt_type (str): The type of simulation. Default is 'reg'.
        reason (str): The reason for the simulation. Default is 'homo'.

    Raises:
        
        Exception: If any error occurs during the MDP file updating process.

    Notes:
        
        The function relies on the availability of the following variables and lists:
        - mdp_lists: A list of MDP files to be copied and updated for regular simulation (opt_type = 'reg' or 'sdk').
        - mdp_lists_it: A list of MDP files to be copied and updated for iterative simulation (opt_type = 'it').
        - cwd: The current working directory before entering the code-specific directories.
        - The function uses the `sed` command to perform the search and replace operation.

    Example:
        
        input_data = {"CODE": ["code1", "code2"], "GROUP": ["group1", "group2"]}
        update_mdp_with_code(input_data, opt_type='reg', reason='homo')
    """
    
    try:
        
        logging.info('Updating mdp files with code ...')
        
        code = input_data["CODE"]
        group = input_data["GROUP"]
        mdp_files = mdp_lists if opt_type in dens_st_cond else mdp_lists_it

        for j, c in enumerate(code):
            
            logging.info(f'Including data of code {c} and transferring files to the respective directory.')
            os.chdir(str(c))

            for l in mdp_files:
                
                os.system(f'cp ../MDP_FILES/{l} .')
                
            if reason == 'homo':
                    
                os.system(f'sed -i "s/INT_GROUPS/{group[j]}/g" *.mdp')
                os.system(f'sed -i "s/INT_TABLES/{group[j]} {group[j]}/g" *.mdp')
                
            elif reason == 'hetero':
                    
                groups = group[j].split('_')
                os.system(f'sed -i "s/INT_GROUPS/{groups[0]} {groups[1]}/g" *.mdp')
                os.system(f'sed -i "s/INT_TABLES/{groups[0]} {groups[0]} {groups[1]} {groups[1]}/g" *.mdp')


            os.chdir(cwd)

        logging.info('Update done ...')
    
    except Exception as e:
        
        logging.error(str(e))
        traceback.print_exc()
        raise

        
        
        
def update_temperature_in_mdp_files(temp_K : float = 298.15):
    
    """
    Updates the temperature value in MDP (molecular dynamics simulation input) files with the specified temperature in Kelvin.

    Args:
        temp_K (float): The temperature value to replace in the MDP files. Default is 298.15 K.

    Notes:
        - The function assumes that there are MDP files present in the current working directory.
        - The function uses the `sed` command to perform the search and replace operation.
        - The string 'TEMP' is used as a placeholder in the MDP files, which will be replaced by the specified temperature.

    Example:
        update_temperature_in_mdp_files(temp_K=310.0)
    """
    
    logging.info('Replacing temperature in mdp files with %s K', temp_K)
    
    os.system(f"sed -i s/TEMP/{temp_K}/g  *.mdp")
        

        

        
        

#Update the forcefield dispersive and repulsive parameters such that 

def update_parameters_of_FF(parameter_groups : dict,
        **kwargs : dict):
    
    """
    Updates specific parameters in a forcefield file based on user-provided input.

    Args:
        
        parameter_groups (dict): A list of tuples representing parameter groups to update. Each tuple consists of three elements:
            
            - prefix (str): The prefix used in the parameter name (e.g., 'disp' or 'rep').
            - suffix (str): The suffix used in the parameter name (e.g., 'w' or 'cm').
            - keys (dict): A dictionary containing the keys used to retrieve corresponding values from the `kwargs` argument. The keys include:
                - 'eps' (str): The key to retrieve the value of the epsilon (eps) parameter for the group.
                - 'sig' (str): The key to retrieve the value of the sigma (sig) parameter for the group.
                - 'lamb_1' (str or float): The key or float value to retrieve the value of the lamb_1 parameter for the group.
                - 'lamb_2' (str or float): The key or float value to retrieve the value of the lamb_2 parameter for the group.

        **kwargs (dict): Additional keyword arguments to specify the values for the parameters and lambdas. The keys in `kwargs` should match the keys specified in the `parameter_groups` argument for eps, sig, lamb_1, and lamb_2.

    Notes:
        
        - The function assumes that the forcefield file 'FORCEFIELD_ITP/Parameters.itp' exists in the current working directory.
        - The function performs parameter updates based on the provided `parameter_groups` and `kwargs`.
        - The function calculates and replaces the 'disp' and 'rep' parameter values in the forcefield file using the specified epsilon (eps), sigma (sig), and lambda (lamb_1, lamb_2) values.
        - The updated forcefield file is saved as 'FORCEFIELD_ITP/Parameters_updated.itp'.

    Example:
        
        parameter_groups = [
            ('prefix1', 'suffix1', {'eps': 'eps_1', 'sig': 'sig_1', 'lamb_1': 'lamb1_1', 'lamb_2': 'lamb2_1'}),
            ('prefix2', 'suffix2', {'eps': 'eps_2', 'sig': 'sig_2', 'lamb_1': 'lamb1_2', 'lamb_2': 'lamb2_2'})
        ]
        kwargs = {
            'eps_1': 1.5,
            'sig_1': 3.0,
            'lamb1_1': 1.8,
            'lamb2_1': 1.2,
            'eps_2': 2.0,
            'sig_2': 2.5,
            'lamb1_2': 2.5,
            'lamb2_2': 3.0
        }
        update_parameters_of_FF(parameter_groups, **kwargs)
    """
    
    os.chdir(cwd)

    logging.info('Updating parameters in the forcefield files ....')

    with open('FORCEFIELD_ITP/Parameters.itp', 'r') as f:
        
        fdata = f.read()

    for group in parameter_groups:
        
        prefix, suffix, keys = group
        eps = kwargs.get(keys['eps'], None)
        sig = kwargs.get(keys['sig'], None)
        
        if eps is not None:
        
            sig_nm =  sig/10
            
            if isinstance(keys['lamb_1'],float):
                
                lamb_1 = keys['lamb_1']
                lamb_2 = keys['lamb_2']
            
            else:
                pass
            
            prefactor = (lamb_1 / (lamb_1 - lamb_2)) * (lamb_1 / lamb_2) ** (lamb_2 / (lamb_1 - lamb_2))
            disp = round(prefactor * eps * (sig_nm) ** lamb_2, 6)
            rep = round(prefactor * eps * (sig_nm) ** lamb_1, 6)
            fdata = fdata.replace(f'disp_{prefix}_{suffix}_', f'{disp}')
            fdata = fdata.replace(f'rep_{prefix}_{suffix}_', f'{rep}')

    with open('FORCEFIELD_ITP/Parameters_updated.itp', 'w') as f:
        
        f.write(fdata)

    logging.info('Update completed of forcefield file completed ...')
    
    




# #Creating directory for simulations of successive iterations.
class FileTransfer:

    """
    The FileTransfer class provides methods for creating directories, transferring files, and executing protocols.

    Args:
        
        dirpath (str): The path of the directory to be created.
        opt_type (str): The type of operation to perform. Defaults to 'reg'.
        reason (str): The reason for the operation. Defaults to 'homo'.

    Attributes:
        
        dirpath (str): The path of the directory to be created.
        opt_type (str): The type of operation to perform.
        reason (str): The reason for the operation.
        logging (Logger): The logging object for logging messages.

    Methods:
        
        _mkdir_and_cd(): Creates a new directory, copies a protocol file based on the opt_type, and changes the current working directory.
        transferfiles(temp_K=298.15): Executes the protocol file, updates the temperature, and changes the current working directory.

    Example:
        
        transfer = FileTransfer(dirpath='/path/to/directory', opt_type='reg', reason='homo')
        transfer.transferfiles(temp_K=300.0)
    """
    
    def __init__(self,
                 dirpath : str,
                 opt_type : str = 'reg',
                 reason : str = 'homo'):
    
        self.dirpath=dirpath
        self.opt_type = opt_type
        self.reason=reason
    
    def _mkdir_and_cd(self):
        
        """Creates a new directory, copies a protocol file, and changes the current working directory."""
        
        logging.info(f'Creating directory in {cwd} ....')
        logging.info(self.dirpath)
        os.system('mkdir ' + str(self.dirpath))

        if self.opt_type in ['reg','sdk','comp']:
            
            os.system('cp PROTOCOLS/file_transfer_protocol.file ' + str(self.dirpath) + '/')

        elif self.opt_type == 'it':
            
            os.system('cp PROTOCOLS/file_transfer_protocol_it.file ' + str(self.dirpath) + '/')

        logging.info('Transfer protocol added to the new directory.')
        os.chdir(cwd)
        logging.info('Task completed ...')

    def transferfiles(self, 
                      temp_K : float = 298.15):
        
        """Executes the protocol file, updates the temperature, and changes the current working directory.

        Args:
            
            temp_K (float): The temperature value in Kelvin. Default temperature is 298.15."""
         
        try:
            
            self._mkdir_and_cd()
        
            if self.opt_type in dens_st_cond:
                
                os.chdir(self.dirpath)
                os.system('./file_transfer_protocol.file')

            elif self.opt_type == 'it':
            
                os.chdir(self.dirpath)
                os.system('./file_transfer_protocol_it.file')

            update_temperature_in_mdp_files(temp_K)
        
            os.chdir(cwd)
            
        except FileNotFoundError as e:
            
            logging.error('Error reading data file: {}'.format(e))
            traceback.print_exc()

        except Exception as e:
            
            logging.error('Error updating table: {}'.format(e))
            traceback.print_exc()
            

    

    
# #Simulating density, surface tension, interfacial tension 
class Simulation:

    """
    The Simulation class provides methods for performing density, surface tension, and interfacial tension simulations.

    Args:
        
        dirpath (str): The directory path.
        opt_type (str): The optional type of simulation. Defaults to 'reg'.
        reason (str): The reason for the simulation. Defaults to 'homo'.

    Attributes:
        
        dirpath (str): The directory path.
        opt_type (str): The optional type of simulation.
        reason (str): The reason for the simulation.
        logging (Logger): The logging object for logging messages.

    Methods:
        
        sim_density(): Performs a density simulation.
        _gas_phase_nvt_box_simulation(): Selects a single atom for gas phase simulation.
        _xyz_npt_nvt(): Returns the dimensions of the simulation box.
        sim_st(): Performs a surface tension simulation.
        sim_it(): Performs an interfacial tension simulation.

    Example:
        
        sim = Simulation(dirpath='/path/to/directory', opt_type='reg', reason='homo')
        sim.sim_density()
        sim.sim_st()
    """
    
    
    def __init__(self,
                 dirpath : str,
                 opt_type : str = 'reg',
                 reason : str = 'homo'):
        
        self.dirpath = dirpath
        self.opt_type = opt_type
        self.reason = reason
        self.logging = logging.getLogger(__name__)

    def sim_density(self):
        
        """
        Performs a density simulation by calling other methods based on the optional type specified.
        """
        
        logging.info('Starting density simulation')
        
        if self.opt_type in dens_st_cond:
            
            os.chdir(str(self.dirpath))
            logging.info('Starting NPT simulation')
            os.system('./nptlog.file')
            logging.info('NPT simulation completed')
            subprocess.run(['gmx', 'trjconv', '-f', 'md.gro', '-s', 'md.tpr', '-pbc', 'whole', '-o', 'md.gro'], capture_output=True, input="System\n0".encode())

        
        elif self.opt_type != 'it':
            
            logging.error('Invalid opt_type')
            traceback.print_exc()
        
        elif self.opt_type == 'it':
            
            logging.error("Optimization procedure chosen for interfacial tension simulation please choose opt_type:reg or sdk")
            traceback.print_exc()

        if self.opt_type == 'reg':
            
            self._gas_phase_nvt_box_simulation()
            os.system('./nvtgas.file')
        
        os.chdir(cwd)
        
        logging.info('Density simulation completed')

    def _gas_phase_nvt_box_simulation(self):
        
        """
        Selects a single atom for gas phase simulation of the given molecule.
        """
        
        logging.info('Selecting single atom for gas phase simulation of the given molecule')
        initialize = mda.Universe('md.gro')
        selection = initialize.select_atoms('resid 750')
        selection.write('gas_phase.gro')
        os.system('gmx editconf -f gas_phase.gro -c -o gas_phase.gro')

    def _xyz_npt_nvt(self):
        
        """
        Returns the dimensions of the simulation box.
        """
        
        if self.opt_type in dens_st_cond:
            line = subprocess.check_output(['tail', '-1', 'md.gro']).split()
            x, y, z= float(line[0]), float(line[1]), 3*float(line[2])
            z_half = z / 2
            x_c, y_c, z_c = x / 2, y / 2, z / 2
            return x, y, z, x_c, y_c, z_c, z_half
        
        elif self.opt_type == 'it':
        
            line = subprocess.check_output(['tail', '-1', 'md_cube_current.gro']).split()
            z_l = float(line[2])
        
            return z_l
        else:
        
            logging.error('Invalid opt_type')
            traceback.print_exc()

    def sim_st(self):
        
        """
        Performs a surface tension simulation by calling other methods based on the optional type specified.
        """
        
        logging.info('Starting surface tension simulation')
        
        if self.opt_type in dens_st_cond:
        
            os.chdir(str(self.dirpath))
            logging.info('Preparing system for NVT simulation')
            subprocess.run(['gmx', 'trjconv', '-s','md.tpr','-f', 'md.gro', '-pbc','whole','-o', 'md.gro'],
                           capture_output=True, input="System\n0".encode(), check=True)
            os.system('gmx editconf -f md.gro -o start_nvt.gro -box ' + str(
                self._xyz_npt_nvt()[0]) + ' ' + str(self._xyz_npt_nvt()[1]) + ' ' + str(
                self._xyz_npt_nvt()[2]) + ' -center ' + str(
                self._xyz_npt_nvt()[3]) + ' ' + str(self._xyz_npt_nvt()[4]) + ' ' + str(
                self._xyz_npt_nvt()[5]))
            logging.info('Starting NVT simulation')
            os.system('./nvtlog.file')
            logging.info('NVT simulation completed')
            os.chdir(cwd)
        
        else:
        
            logging.error('Invalid opt_type')
            traceback.print_exc()
        
        logging.info('Surface tension simulation completed')
        
    def sim_it(self):
        
        """
        Performs an interfacial tension simulation if the optional type is 'it'.
        """
    
        if self.opt_type=="it":

            logging.info('Starting interfacial tension simulation')    
            os.chdir(str(self.dirpath))
            logging.info('Preparing system for IT simulation ....')   	
            os.system('./snptlog.file')    		
            logging.info('IT simulation completed ...')
            os.chdir(cwd)
                      
        elif self.opt_type in dens_st_cond:
                      
            logging.error("This optimization procedure for bulk MD of one kind of molecules please choose opt_type: it")
            traceback.print_exc()
    
        else:
        
            logging.error('Invalid opt_type')
            traceback.print_exc()
        
        logging.info('Interfacial tension simulation completed')


    
    
    

class SimulationData:

    """
    The SimulationData class provides methods for calculating various properties of a molecular simulation.

    Args:
        
        dirpath (str): The directory path of the simulation.

    Attributes:
        
        dirpath (str): The directory path of the simulation.

    Methods:
        
        _run_gmx_energy(edr_file, output_file, input_str): Runs the Gromacs "gmx energy" command with the given arguments and input string.
        _get_line_number(search_string, filename): Searches for a given string in a given file and returns the line number where it was found.
        calculate_mean_st(): Calculates the average surface tension of the simulation.
        calculate_mean_it(): Calculates the average interfacial tension of the simulation.
        calculate_mean_density(): Calculates the average density of the simulation.
        calculate_mean_enthalpy_of_vaporization(code, temp_K): Calculates the average enthalpy of vaporization of the simulation.

    Example:
        
        sim_data = SimulationData(dirpath='/path/to/simulation')
        mean_st = sim_data.calculate_mean_st()
        mean_it = sim_data.calculate_mean_it()
        mean_density = sim_data.calculate_mean_density()
        enthalpy_of_vaporization = sim_data.calculate_mean_enthalpy_of_vaporization(code='XYZ', temp_K=300.0)
    """
    
    
    def __init__(self, 
            dirpath : str):
        
        self.dirpath = dirpath

    
    def _run_gmx_energy(self, 
            edr_file : str,
            output_file : str,
            input_str : str,
            b_time : float = 5):
        
        """
        Runs the Gromacs "gmx energy" command with the given arguments and input string.

        Args:
            
            edr_file (str): The input .edr file.
            output_file (str): The output file to store the results.
            input_str (str): The input string for gmx energy command.

        Raises:
            
            subprocess.CalledProcessError: If the command fails.
        """
        
        try:
            
            subprocess.run(['gmx', 'energy', '-f', edr_file, '-b',str(b_time*1000),'-o', output_file, '-xvg', 'none'],
                           capture_output=True, input=input_str.encode(), check=True)
        
        except subprocess.CalledProcessError as e:
            
            logging.error(f"Error running gmx energy: {e}")
            raise

    def _get_line_number(self, 
            search_string : str,
            filename : str):

        """
        Searches for a given string in a given file and returns the line number where it was found.

        Args:
            
            search_string (str): The string to search for.
            filename (str): The name of the file to search in.

        Returns:
            
            int: The line number where the search string was found.

        Raises:
            
            ValueError: If the string is not found in the file.
        """
        
        try:
            
            with open(filename) as f:
                
                for i, line in enumerate(f, 1):
                    
                    if search_string in line:
                        
                        return i
                
                raise ValueError(f"{search_string} not found in {filename}")
        
        except Exception as e:
            
            logging.error(f"Error getting line number: {e}")
            raise

    
    
    
    def calculate_mean_st(self):
        
        """
        Calculates the average surface tension of the simulation.

        Returns:
            
            float: The mean surface tension.

        Raises:
            
            Exception: If there is an error calculating the mean surface tension.
        """
        
        try:
            
            os.chdir(str(self.dirpath))
            logging.info('Calculating the average surface tension ....')
            
            self._run_gmx_energy('mdnvt.edr', 
                    'surfacetension_mdnvt.xvg', 
                    "#Surf*SurfTen\n0")
            
            data = np.loadtxt('surfacetension_mdnvt.xvg')
            mean_st = np.mean(data[:, 1]) / 20
            logging.info('Calculation and storage of average surface tension concluded ...')
            os.chdir(cwd)
            
            return mean_st
        
        except Exception as e:
            
            logging.error(f"Error calculating mean surface tension: {e}")
            raise

    

    def calculate_mean_it(self):
        
        """
        Calculates the average interfacial tension of the simulation.

        Returns:
            
            float: The mean interfacial tension.

        Raises:
            
            Exception: If there is an error calculating the mean interfacial tension.
        """
        
        try:
            
            os.chdir(str(self.dirpath))
            logging.info('Calculating the average interfacial tension ....')
            
            self._run_gmx_energy('md_cube_current.edr', 
                    'interfacialtension_md.xvg', 
                    "#Surf*SurfTen\n0")
            
            data = np.loadtxt('interfacialtension_md.xvg')
            mean_it = np.mean(data[:, 1]) / 20
            logging.info('Calculation and storage of average interfacial tension concluded ...')
            os.chdir(cwd)
            
            return mean_it
        
        except Exception as e:
            
            logging.error(f"Error calculating mean interfacial tension: {e}")
            
            raise

    def calculate_mean_density(self):
        
        """
        Calculates the average density of the simulation.

        Returns:
            
            float: The mean density.

        Raises:
            
            Exception: If there is an error calculating the mean density.
        """
        
        try:
            
            os.chdir(str(self.dirpath))
            
            self._run_gmx_energy('md.edr', 
                    'average_density.xvg', 
                    "Density\n0")
            
            logging.info('Calculating the average density ....')
            data = np.loadtxt('average_density.xvg')
            mean_density = np.mean(data[:, 1])
            logging.info('Calculation and storage of average density concluded ...')
            os.chdir(cwd)
            
            return mean_density
        
        except Exception as e:
            
            logging.error(f"Error calculating mean density: {e}")
            raise

    


    def calculate_mean_enthalpy_of_vaporization(self, 
            code : str,
            temp_K : float):
        
        """
        Calculates the average enthalpy of vaporization of the simulation.

        Args:
            
            code (str): The code to search for in the topol.top file.
            temp_K (float): The temperature in Kelvin.

        Returns:
            
            float: The mean enthalpy of vaporization.

        Raises:
            
            Exception: If there is an error calculating the mean enthalpy of vaporization.
        """
        
        try:
        
            os.chdir(str(self.dirpath))
            
            self._run_gmx_energy('md.edr', 
                    'potential.xvg', 
                    "Potential\n0")
            
            self._run_gmx_energy('mdgas.edr', 
                    'potential_gas.xvg', 
                    "Potential\n0")
            
            logging.info('Calculating the average enthalpy of vaporization ....')
            data_pot = np.loadtxt("potential.xvg")
            data_pot_gas = np.loadtxt("potential_gas.xvg")
            potential = data_pot[:, 1]
            potential_gas = data_pot_gas[:, 1]
            mean_potential = np.mean(potential)
            mean_potential_gas = np.mean(potential_gas)
            line = linecache.getline("topol.top", self._get_line_number(f"{code} ", "topol.top"))
            number_of_molecules = float(line.split()[1])
            linecache.clearcache()
            R_kJ = 8.314 / 1000
            T = float(temp_K)
            H = mean_potential_gas - mean_potential / number_of_molecules + R_kJ * T
            logging.info('Calculation and storage of average enthalpy of vaporization concluded ...')
            os.chdir(cwd)
        
            return H

        except Exception as e:
        
            logging.error(f"Error calculating mean enthalpy of vaporization: {e}")
            
            raise


    def calculate_mean_compressibility(self,
                                       code : str):

        """
                Calculates the average enthalpy of vaporization of the simulation.

                Returns:

                    float: The mean isothermal compressibility.

                Raises:

                    Exception: If there is an error calculating the mean isothermal compressibility.
        """

        try:

            os.chdir(str(self.dirpath))
            logging.info('Calculating the average compressibility....')

            line = linecache.getline('topol.top', self._get_line_number(f"{code} ", 'topol.top'))
            k = line.split()

            linecache.clearcache()

            with open("compress.txt", "w") as output:

                subprocess.run(
                    ['gmx', 'energy', '-f', 'md.edr', '-fluct_props', '-b','5000' , '-nmol', str(k[1]), '-driftcorr'],
                    input='Volume\nTemperature\n0'.encode(), stdout=output)
                data = linecache.getline('compress.txt', 19)
                data_string = str(data).split()
                mean_val = float(data_string[4])
                comp = mean_val * 1e9
                linecache.clearcache()

            logging.info('Calculation and storage of average isothermal compressibility concluded ...')

            os.chdir(cwd)
            print(comp)

            return comp

        except Exception as e:

            logging.error(f"Error calculating mean isothermal compressibility: {e}")
            
            raise


#Functions to write the output data files.

def data_file_writer_iterations(iteration_output_file : str,
                     cost : float,
                     **kwargs : dict):

    """Writes data from iterations to a file in JSON format.

        This function takes the provided cost value and any additional key-value pairs
        from kwargs and writes them to a file specified by iteration_output_file in
        JSON format. Each iteration's data is appended to the file.

        Args:

            iteration_output_file (str): The file path or name where the iteration data
                                         will be written.

            cost (float): The cost value associated with the current iteration.

            **kwargs (dict): Additional key-value pairs representing data associated
                             with the current iteration. The keys are strings representing the data
                             labels, and the values are the corresponding data values.

        Raises:
            IOError: If there is an error in opening or writing to the output file.

        Example:
            # Write iteration data to 'output.json' file
            data_file_writer_iterations('output.json', 3.14, **kwargs)

        In the above example, the function writes the following data to the 'output.json'
        file:

        kwargs={
            "eps1": 2.5,
            "sig1": 4.3,
            "eps2": 1.3,
            "sig2": 4.4
         }

        Note:

            The function assumes that the file specified by iteration_output_file already
            exists. If the file does not exist, an IOError will be raised.
    """

    values = [x for x in kwargs.values()] + [cost]
    data = {key: value for key, value in zip(kwargs.keys(), values)}
    data["cost"] = cost

    file = cwd + iteration_output_file
    with open(file, 'a') as f:
        json.dump(data, f)
        f.write(',\n')




def percentage_error(sim_value,
                     expt_value):
    
    return (abs(sim_value-expt_value)/expt_value)*100




def relu(x):

    return max(0,x)





#Cost function for optimization.


def cost_function(input_data : dict,
                  density : list = [],
                  surfacetension : list = [],
                  ev : list = [],
                  comp : list = [],
                  interfacialtension : list = [],
                  kind : str = 'hierarchial'):

    """
    Evaluates the performance of a model that predicts the physical properties of chemical compounds.

    Args:

        input_data (dict): The input data containing the predicted values of the properties.
        density (list): The predicted values of density.
        surfacetension (list): The predicted values of surface tension.
        ev (list): The predicted values of evaporation enthalpy.
        comp (list): The predicted values of compressibility.
        interfacialtension (list): The predicted values of interfacial tension.
        kind (str, optional): The type of cost function to be used. Defaults to 'hierarchical'.

    Returns:
        
        float: The objective value representing the cost of the model's predictions.

    Raises:
        
        None

    Docstring Explanation:
    
    The cost function evaluates the predictions made by a model for the physical properties of chemical compounds. 
    It takes the predicted values of density, surface tension, evaporation enthalpy, isothermal compressibility and interfacial tension,
    and computes an objective value that represents the cost of these predictions.

    The cost function can use two types of cost functions: 'hierarchical' and 'scaled_absolute'. 
    The 'hierarchical' cost function penalizes predictions based on the percentage error between the predicted 
    and experimental values. It assigns a reward value of 1000 if the error falls within a certain tolerance percentage, 
    and decreases the reward linearly for errors outside the tolerance. The 'scaled_absolute' cost function scales 
    the absolute error by the mean experimental value of the corresponding property.

    The objective value is calculated by comparing the predicted values with the experimental values and 
    applying the selected cost function. The objective value represents the cost of the model's predictions and 
    is used for parameter optimization.
    """
    
   
    tolerance_percentages = {'density': 1.5, 
                             'surfacetension': 10, 
                             'ev': 5, 
                             'comp': 15,
                             'interfacialtension': 5}
    
    if kind=='hierarchial':
        reward = 1000
        
        def obj_val(sim : float,
                exp : float ,
                tol : float):
            return (1 - relu((percentage_error(sim,exp) - tol) / 100)) * reward
    
    elif kind=='scaled_absolute':
        reward = 0.0
        
        def obj_val(sim : float,
                exp : float,
                scl : float):
            return scl*np.absolute(sim-exp)
        
        
    objective_value = 0
    
    sim_data = {'density': density, 
                'surfacetension': surfacetension, 
                'ev': ev, 
                'comp': comp,
                'interfacialtension': interfacialtension}
    
    expt_data = {'density': [], 
                 'surfacetension': [], 
                 'ev': [], 
                 'comp':[],
                 'interfacialtension': []}
        
    props = [k for k in sim_data.keys() if len(sim_data[k]) != 0]
        
    for i in props:
            
        expt_data[i] = input_data.get(i.upper())
        
        for j in range(len(expt_data[i])):

            if percentage_error(sim_data[i][j], expt_data[i][j]) <= tolerance_percentages[i]:

                if kind=="hierarchial":

                    objective_value += reward + 100*np.sin((tolerance_percentages[i] - percentage_error(sim_data[i][j],expt_data[i][j]))/tolerance_percentages[i])

                elif kind=="scaled_absolute":

                    objective_value += reward
                
            else:
            
                if kind=="hierarchial":

                    objective_value += obj_val(sim_data[i][j],expt_data[i][j],tolerance_percentages[i])
                    
                elif kind=="scaled_absolute":    

                    scale = 1/np.mean(expt_data[i])
                    objective_value += obj_val(sim_data[i][j],expt_data[i][j],scale)
                    
    return objective_value

    




def path(dirpath : str,
         quantity : str):
    
    
    """
    Generates the file path for a given quantity in a specified directory.

    Args:
        
        dirpath (str): The directory path where the file is located.
        quantity (str): The quantity for which the file path is required.

    Returns:
        
        str: The full file path.

    Raises:
        
        ValueError: If the quantity argument is not found in the file_paths dictionary.
    """
    
    
    file_paths = {
        'density': 'average_density.xvg',
        'st': 'surfacetension_mdnvt.xvg',
        'ev': 'average_density.xvg',
        'comp': 'average_density.xvg',
        'it': 'interfacialtension_md.xvg'
    }
    
    file_name = file_paths.get(quantity)
    
    if file_name is None:
        
        raise ValueError('Invalid quantity: {}'.format(quantity))
    
    return os.path.join(os.getcwd(), dirpath, file_name)





fcost = []

def run_simulations_over_properties(iniguess,
                                    input_data : dict,
                                    cost_kind : str = "hierarchial",
                                    opt_type : str = 'reg',
                                    reason : str = 'homo'):

    """
    Runs simulations over properties of compounds based on the provided arguments.

    Args:
        
        iniguess (list): List of initial guesses.
        input_data (dict): Dictionary containing input data.
        cost_kind (str, optional): Type of cost calculation. Defaults to "hierarchial".
        opt_type (str, optional): Type of optimization. Defaults to 'reg'.
        reason (str, optional): Reason for simulation. Defaults to 'homo'.

    Returns:
        
        float: Total cost of the simulations.

    Raises:
        
        FileNotFoundError: If there is an error reading the data file.
        Exception: If there is an error updating the property data.

    Docstring Explanation:
    The run_simulations_over_properties function performs simulations over properties of compounds based on the provided arguments.

    The function takes several arguments:
    
    - iniguess: A list of initial guesses.
    - input_data: A dictionary containing input data.
    - cost_kind: A string representing the type of cost calculation. It defaults to "hierarchial".
    - opt_type: A string representing the type of optimization. It defaults to 'reg'.
    - reason: A string representing the reason for simulation. It defaults to 'homo'.

    The function defines an internal execution function that performs the main simulation steps.
    It initializes empty lists for storing simulation results.
    It updates the parameters of the FF (Force Field) using the update_parameters_of_FF function.
    It iterates over each compound in the compound data dictionary.
    Inside the iteration, it sets up the directory path, creates instances of Simulation and SimulationData classes,
    and evaluates the mean values of density, surface tension, and evaporation enthalpy based on the availability of simulation files.
    The calculated values are appended to the corresponding lists.
    The total cost is calculated using the cost_function with appropriate arguments based on the opt_type and cost_kind.
    Finally, the execution function returns the total cost.

    The main function attempts to execute the execution function and store the result in cost_tot.
    If an error occurs during execution, it changes the current working directory back to the original directory and retries the execution.
    If an error still occurs, it does nothing.
    The cost_tot value is appended to the fcost list and then returned based on the cost_kind.

    The function also includes exception handling for FileNotFoundError and general Exception.
    If a FileNotFoundError occurs during execution, an error message is printed, and traceback is printed.
    If any other Exception occurs, an error message is printed, and traceback is printed.
    """    
    
                  
        
    def execution():
                
        sim_density_over_temp_mol, sim_st_over_temp_mol, sim_ev_over_temp_mol, sim_comp_over_temp_mol, sim_it_over_temp_mol = [], [], [], [], [] 
        
        update_mdp_with_code(input_data,
                             opt_type=opt_type,
                             reason=reason)
       
        for j in range(len(input_data["CODE"])):

            code, temp_K=input_data["CODE"][j],input_data["TEMPERATURE"][j]
            parameters = '_'.join(str(x) for x in iniguess)
            dirpath = f"{code}/Parameter_{parameters}_{temp_K}"
        
            
            file_generator=FileTransfer(dirpath=dirpath,
                                        opt_type=opt_type,
                                        reason=reason)
            
            file_generator.transferfiles(temp_K=temp_K)
            sims = Simulation(dirpath=dirpath,
                              opt_type=opt_type,
                              reason=reason)
            
            means = SimulationData(dirpath=dirpath)
            

            if opt_type in dens_st_cond:

                if os.path.exists(path(dirpath,"density"))==True and os.stat(path(dirpath,"density")).st_size != 0:
                
                    mean_density_compress = means.calculate_mean_density()
                
                elif (os.path.exists(path(dirpath,"density"))==True and os.stat(path(dirpath,"density")).st_size == 0) or os.path.exists(path(dirpath,"density"))==False:
                
                    sims.sim_density()
                    mean_density_compress = means.calculate_mean_density()
                    
                    
                if os.path.exists(path(dirpath,"st"))==True and os.stat(path(dirpath,"st")).st_size != 0:
                
                    mean_surfacetension = means.calculate_mean_st()
                                    
                elif (os.path.exists(path(dirpath,"st"))==True and os.stat(path(dirpath,"st")).st_size == 0) or os.path.exists(path(dirpath,"st"))==False:
                
                    sims.sim_st()
                    mean_surfacetension = means.calculate_mean_st()
                    
                                    
                if opt_type=='reg':
                    
                    if os.path.exists(path(dirpath,"ev"))==True:
                        
                        mean_ev = means.calculate_mean_enthalpy_of_vaporization(code=code,temp_K=temp_K)
                        sim_ev_over_temp_mol.append(mean_ev)
                       
                    elif os.path.exists(path(dirpath,"ev")) == False:
                        
                        sims.sim_density()
                        mean_ev = means.calculate_mean_enthalpy_of_vaporization(code=code,temp_K=temp_K)
                        sim_ev_over_temp_mol.append(mean_ev)

                elif opt_type=='comp':

                    if os.path.exists(path(dirpath,"comp"))==True:

                        mean_comp = means.calculate_mean_compressibility(code="TO")
                        sim_comp_over_temp_mol.append(mean_comp)

                    elif os.path.exists(path(dirpath,"comp")) == False:

                        sims.sim_density()
                        mean_comp = means.calculate_mean_compressibility(code="TO")
                        sim_comp_over_temp_mol.append(mean_comp)

                        
                       
                sim_density_over_temp_mol.append(mean_density_compress)
                sim_st_over_temp_mol.append(mean_surfacetension)
                
            elif opt_type=="it":
                
                if os.path.exists(path(dirpath,"it"))==True and os.stat(path(dirpath,"it")).st_size != 0:
                
                    mean_it = means.calculate_mean_it()
                    
                elif (os.path.exists(path(dirpath,"it"))==True and os.stat(path(dirpath,"it")).st_size == 0) or os.path.exists(path(dirpath,"it"))==False:
                
                    sims.sim_it()
                    mean_it = means.calculate_mean_it()                                    
                    
                sim_it_over_temp_mol.append(mean_it)
                
             
        if opt_type=="reg":
                
            cost_tot = cost_function(input_data=input_data,density=sim_density_over_temp_mol,surfacetension=sim_st_over_temp_mol,ev=sim_ev_over_temp_mol,kind=cost_kind)
        
        elif opt_type=="sdk":
            
            cost_tot = cost_function(input_data=input_data,density=sim_density_over_temp_mol,surfacetension=sim_st_over_temp_mol,kind=cost_kind)

        elif opt_type=="comp":

            cost_tot = cost_function(input_data=input_data,density=sim_density_over_temp_mol,surfacetension=sim_st_over_temp_mol,comp=sim_comp_over_temp_mol,kind=cost_kind)
            
        elif opt_type=="it":
            
            cost_tot = cost_function(input_data=input_data,interfacialtension=sim_it_over_temp_mol,kind=cost_kind)
            
        return cost_tot
    
    try:
    
        try:
        
            
            cost_tot=execution()
        
   
        
        except:
        
            os.chdir(cwd)
        
            try:
            
                cost_tot=execution()
            
            except:
        
                print('Error occured')

                if cost_kind == "hierarchial":
        
                    cost_tot = 500
        
                elif cost_kind == "scaled_absolute":
        
                    cost_tot = 12.0
        
                os.chdir(cwd)
    
        fcost.append(cost_tot)

        if cost_kind == "hierarchial":
        
            return cost_tot
    
        elif cost_kind == "scaled_absolute":
        
            return -cost_tot

    except FileNotFoundError as e:
            
        logging.error('Error reading data file: {}'.format(e))
        traceback.print_exc()

    except Exception as e:
            
        logging.error('Error updating property data: {}'.format(e))
        traceback.print_exc()
