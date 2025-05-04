#!/usr/bin/env python3
import os,sys,re
import traceback

class Table(object):

    def __init__(self, arrays=None, header=None):
        """
        Create a new Table instance
        """
        self.__numCols = 0
        self.__numRows = 0
        self.__content = {}
        self.__header = []
        
        if arrays is not None and len(arrays) >= 1:
            
            # Check for homogeneous array length
            length = len(arrays[0])
            for array in arrays[1:]:
                if len(array) != length:
                    raise TableError('Heterogeneous array length')
            
            if header is not None:
                if len(header) != len(arrays):
                    raise TableError('Header length different from array size')
                
            else:
                # Build header
                header = ['Col'+str(i) for i in range(len(arrays))]
            
            self.__numCols = len(arrays)
            self.__numRows = length
            self.__content = {}
            self.__header = header
            
            for headerName,array in zip(header,arrays):
                self.__content[headerName] = array

    def clone(self):
        """
        Returns a copy of this table
        """
        newTable = Table()
        newTable.__content = self.__content.copy()
        newTable.__header = self.__header[:]
        newTable.__numCols = self.__numCols
        newTable.__numRows = self.__numRows
        return newTable
                
    @staticmethod
    def read(tableFile):
        """
        Read a Table from a given file
        """
        try:
            f = open(tableFile)
            tab = f.readlines()
            f.close()
        except IOError as e:
            raise TableError(str(e))
        
        if tab[0].split()[0] == '#>T':    # for icm table format
            tab[0] = tab[1].replace('#>',' ').replace('-',' ')
            tab[1] = "-----------"
        
        arrays = []
        header = tab[0].split()
        for i in range(len(header)): arrays.append([])
        
        # in case header has two lines (first for columns, second just a separator)
        if len(tab[1].split()) == len(tab[0].split()):
            i_first = 1
        else:     
            i_first = 2

        for i in range(i_first,len(tab[i_first:])+i_first):
            line = tab[i].split()
            for k in range(len(header)):
                try:
                    arrays[k].append(int(line[k]))
                except:
                    try:
                        arrays[k].append(float(line[k]))
                    except:
                        try:
                            arrays[k].append(str(line[k]))
                        except:
                            arrays[k].append(' ')

        return Table(arrays,header)


    def get_column(self, colName):
        """
        Get a column identified by colName
        """
        if colName not in self.__header:
            raise TableError('colName is not a valid column name')

        return self.__content[colName]

    def __getitem__(self, key):
        return self.get_column(key)

class PyDockParameters:
    to_rename = {'HID': 'HIS'}
    renameResidues = True
    pdb_extension = '.pdb'
    ligand_extension = '_lig'
    receptor_extension = '_rec'
    reference_extension = '_ref'
    energy_extension = '.ene'
    rot_extension = '.rot'

def write_atom_line(atom, output, is_hetatm, to_rename=None):
    """
    Writes a PDB file format line to output.
    """
    atom_type = atom.atom_type
    if len(atom_type) < 4 :
        atom_type = ' ' + atom_type
    if is_hetatm:
        atom_class = 'HETATM'
    else:
        atom_class = 'ATOM  '
    try:
        residue_type = to_rename[atom.residue_type]
    except:
        residue_type = atom.residue_type
    line = "%6s%5d %-4s%-1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (atom_class, atom.atom_number,
                                                                           atom_type, atom.atom_alternative,
                                                                           residue_type, atom.chain_id,
                                                                           atom.residue_number, atom.residue_ext,
                                                                           atom.x, atom.y, atom.z,
                                                                           atom.occupancy,
                                                                           atom.b_factor)
    output.write(line)

def write_pdb_to_file(macroMol, fileOut, includeHydrogens=False, force_hetatm=False,
                      to_rename=None, rename_residues=PyDockParameters.renameResidues):
    """
    Writes a MacroMolecule structure to a file in PDB format.
    """
    output = open(fileOut, "w")
    for mol in macroMol.molecules:
        try:
            if mol.is_cofactor_ion and not force_hetatm:
                is_hetatm = True
            else:
                is_hetatm = False
        except:
            is_hetatm = False
        for res in mol.get_residues():
            for atom in res.atoms:
                if not atom.is_hydrogen() or includeHydrogens:
                    if rename_residues:
                        write_atom_line(atom, output, is_hetatm, to_rename)
                    else:
                        write_atom_line(atom, output, is_hetatm)
    output.close()
        
class TableError(Exception):
    """Table Exception class"""
    def __init__(self, value):
        self.parameter = value
       
    def __str__(self):
        return self.parameter

class Residue(object):
    def __init__(self, atoms = None):
        self.no_H_atoms = []
        if atoms:
            self.atoms = atoms
            for atom in atoms:
                if not atom.is_hydrogen():
                    self.no_H_atoms.append(atom)
        else:
            self.atoms = []

    def add_atom(self, atom):
        self.atoms.append(atom)
        if not atom.is_hydrogen():
            self.no_H_atoms.append(atom)

    def clone(self):
        raise NotImplementedError()

class AminoAcid(Residue):
    def __init__(self, atoms=None):
        Residue.__init__(self, atoms)

    def clone(self):
        """Deep copy"""
        new_atoms = []
        for atom in self.atoms:
            new_atoms.append(atom.clone())
        return AminoAcid(new_atoms)

def read_rot_file(rot_file_name):
    """
    Reads a rotations table file
    """
    lines = open(rot_file_name).readlines()
    rot_table = []
    for line_count, line in enumerate(lines):
        try:
            l = line.split()
            r11 = float(l[0])
            r12 = float(l[1])
            r13 = float(l[2])
            r21 = float(l[3])
            r22 = float(l[4])
            r23 = float(l[5])
            r31 = float(l[6])
            r32 = float(l[7])
            r33 = float(l[8])
            t1 = float(l[9])
            t2 = float(l[10])
            t3 = float(l[11])
            i_conf = int(l[12])
            rot_table.append((r11, r12, r13, r21, r22, r23, r31, r32, r33, t1, t2, t3, i_conf))
        except:
            log.warning('Skipping line %d in file %s: Empty or wrong line' % (line_count, rot_file_name))
    return rot_table

def is_int(a):
    """
    Checks if a given variable is an integer
    """
    try:
        int (a)
        return True
    except:
        return False

class Atom:
    """
    Atom class
    """
    
    backbone_atoms = ['CA', 'C', 'N', 'O']
    
    def __init__(self, atom_number=None, atom_type=None, atom_alternative=None, 
                 residue_type=None, chain_id=None, residue_number=None, 
                 residue_ext=None, x=None, y=None, z=None, occupancy=None, 
                 b_factor=None, element=None, mass=0., amber_type=None, charge=0.,
                 radius=0.,epsilon_vdw=None):
        self.atom_number = atom_number 
        self.atom_type = atom_type
        self.atom_alternative = atom_alternative
        self.residue_type = residue_type
        self.chain_id = chain_id
        self.residue_number = residue_number
        self.residue_ext = residue_ext
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.b_factor = b_factor
        self.element = element
        self.mass = mass
        self.amber_type = amber_type
        self.charge = charge
        self.radius = radius
        self.epsilon_vdw = epsilon_vdw

    def is_hydrogen(self):
        """
        Checks if this atom is of hydrogen type
        """
        if self.atom_type[0] == 'H':
            return True
        if len(self.atom_type) > 1 and is_int(self.atom_type[0]) and self.atom_type[1] == 'H':
            return True
        return False

    def clone(self):
        """
        Creates a copy of the current atom
        """
        return Atom(self.atom_number,
                    self.atom_type,
                    self.atom_alternative,
                    self.residue_type,
                    self.chain_id,
                    self.residue_number,
                    self.residue_ext,
                    self.x,
                    self.y,
                    self.z,
                    self.occupancy,
                    self.b_factor,
                    self.element,
                    self.mass,
                    self.amber_type,
                    self.charge,
                    self.radius)
    

class Nucleotide(Residue):
    """Represents standard a nucleotide and ribonucleic acid in
    two code letter"""
    __two_letter = {'ADE': 'DA',
                    'GUA': 'DG',
                    'CYT': 'DC',
                    'THY': 'DT',
                    'URA': 'RU',
                    'A': 'RA',
                    'G': 'RG',
                    'C': 'RC',
                    'U': 'RU'}

    # Molecular weight
    mw = {'DA': 347.0,
          'DC': 323.0,
          'DT': 322.0,
          'DG': 363.0}

    # Complementary nucleotides
    complement = {'DA': 'DT',
                  'DC': 'DG',
                  'DG': 'DC',
                  'DT': 'DA'}

    atoms_to_rename = {"OP1": "O1P",
                       "OP2": "O2P",
                       "C1*": "C1'",
                       "C2*": "C2'",
                       "C3*": "C3'",
                       "C4*": "C4'",
                       "C5*": "C5'",
                       "O1*": "O1'",
                       "O2*": "O2'",
                       "O3*": "O3'",
                       "O4*": "O4'",
                       "O5*": "O5'"}
    # DNA and RNA residues to rename ref: J. Chem. Theory Comput. 2007, 3, 4, 1464-1475
    residues_to_rename = {'DA5': 'DA',
                          'DG5': 'DG',
                          'DC5': 'DC',
                          'DT5': 'DT',
                          'DA3': 'DA',
                          'DG3': 'DG',
                          'DC3': 'DC',
                          'DT3': 'DT',
                          'M2A': 'DMA',
                          'A2M': 'MRA',
                          'MIA': 'SPA',
                          '12A': 'STA',
                          'OMC': 'MRC',
                          'G7M': '7MG',
                          'UBG': 'BUG',
                          'M2G': 'DMG',
                          'OMG': 'MRG',
                          'QUO': 'QUG',
                          'YG': 'WBG',
                          'SUR': '2SU',
                          'S4U': '4SU',
                          'H2U': 'DHU',
                          '2MU': 'MMU',}

    @staticmethod
    def get_nucleic_acids_two_letters():
        """Gets the list of all supported nucleic acids."""
        return Nucleotide.__two_letter.values()

    @staticmethod
    def get_nucleic_acids_three_letters():
        """Gets the list of all supported nucleic acids."""
        return Nucleotide.__two_letter.keys()

    def clone(self):
        """Deep copy"""
        new_atoms = []
        for atom in self.atoms:
            new_atoms.append(atom.clone())
        return Nucleotide(new_atoms)

class PyDockError(Exception):
    """
    PyDock exception base class
    """
    def __init__(self, cause):
        self.cause = cause

    def __str__(self):
        representation = "[%s] %s" % (self.__class__.__name__, self.cause)
        return representation

class PDBParsingError(PyDockError):
    """
    Custom PDB Parsing Exception
    """
    pass

class MacroMolecule:

    def __init__(self, molecules=None):
        if molecules is None:
            self.molecules = []
        else:
            self.molecules = molecules

    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    def add_molecules(self, molecules):
        self.molecules.extend(molecules)

    def remove_molecule(self, molecule):
        self.molecules.remove(molecule)

    def get_residues(self):
        residue_list = []
        for mol in self.molecules:
            residue_list.extend(mol.get_residues())
        return residue_list

    def get_atoms(self, consider_hydrogens=True):
        atom_list = []
        for molecule in self.molecules:
            atom_list.extend(molecule.get_atoms(consider_hydrogens))
        return atom_list

    def clone(self):
        new_molecules = []
        for mol in self.molecules:
            new_molecules.append(mol.clone())
        return MacroMolecule(new_molecules)

    def renumber_atoms(self):
        for atom_id, atom in enumerate(self.get_atoms()):
            atom.atom_number = atom_id+1

    def rotation_and_translation(self, rotation_translation):
        """
        Applies first a rotation and then a translation to this macromolecule.
        """
        for atom in self.get_atoms():
            new_x = (rotation_translation[0] * atom.x +
                         rotation_translation[1] * atom.y +
                         rotation_translation[2] * atom.z +
                         rotation_translation[9])
            new_y = (rotation_translation[3] * atom.x +
                         rotation_translation[4] * atom.y +
                         rotation_translation[5] * atom.z +
                         rotation_translation[10])
            new_z = (rotation_translation[6] * atom.x +
                         rotation_translation[7] * atom.y +
                         rotation_translation[8] * atom.z +
                         rotation_translation[11])
            atom.x = new_x
            atom.y = new_y
            atom.z = new_z

    def move_to_origin(self):
        """
        Moves this macromolecule to the coordinates origin given by the center of coordinates
        """
        center_x, center_y, center_z = self.get_coordinates_center()
        self.translate ([-center_x, -center_y, -center_z])

    def rotate(self, rotation_matrix):
        """
        Applies a rotation to this macromolecule given by the rotation matrix
        """
        for atom in self.get_atoms():
            new_x = (rotation_matrix[0] * atom.x +
                         rotation_matrix[1] * atom.y +
                         rotation_matrix[2] * atom.z)
            new_y = (rotation_matrix[3] * atom.x +
                         rotation_matrix[4] * atom.y +
                         rotation_matrix[5] * atom.z)
            new_z = (rotation_matrix[6] * atom.x +
                         rotation_matrix[7] * atom.y +
                         rotation_matrix[8] * atom.z)
            atom.x = new_x
            atom.y = new_y
            atom.z = new_z

    def translate(self, translation):
        """
        Applies an spatial translation to this macromolecule given by the
        translation vector
        """
        for atom in self.get_atoms():
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]

    def get_coordinates_center(self, consider_hydrogens=True):
        """
        Calculates the center of coordinates of this macromolecule.
        """
        totalX = 0.0
        totalY = 0.0
        totalZ = 0.0
        atoms = self.get_atoms(consider_hydrogens)
        for atom in atoms:
            totalX += atom.x
            totalY += atom.y
            totalZ += atom.z
        numAtoms = float(len(atoms))
        return totalX / numAtoms, totalY / numAtoms, totalZ / numAtoms

    def __str__(self):
        output = ""
        for molecule in self.molecules:
            output += "%s" % (str(molecule))
        return output

class Cofactor(Residue):
    pass

class Molecule:
    def __init__(self, n_chainId, n_residues=None, is_cofactor_ion=False, contains_DNA=False, contains_AAs=False):
        self.__chainId = n_chainId
        if n_residues is None:
            self.residues = []
        else:
            self.residues = n_residues
        self.is_cofactor_ion = is_cofactor_ion
        self.contains_DNA = contains_DNA
        self.contains_AAs = contains_AAs

    def add_residue(self, residue):
        self.residues.append(residue)

    def get_chain_id(self):
        return self.__chainId

    def get_residues(self):
        return self.residues

    def get_atoms(self, consider_hydrogens=True):
        atom_list = []
        for residue in self.residues:
            if consider_hydrogens:
                atom_list.extend(residue.atoms)
            else:
                atom_list.extend(residue.no_H_atoms)
        return atom_list

    def clone(self):
        new_residues = []
        for residue in self.residues:
            new_residues.append(residue.clone())
        return Molecule(self.get_chain_id(), new_residues, self.is_cofactor_ion, self.contains_DNA, self.contains_AAs)


    def __str__(self):
        output = ""
        residues = self.get_residues()
        for residue in residues:
            output += "%s" % (str(residue))
        return output

def cstrip(string):
    """
    Remove unwanted symbols from string.
    """
    return string.strip(' \t\n\r')

def read_atom_line(line, nanometers=False):
    """
    Parses a PDB file line starting with 'ATOM'.
    """
    ELEMENT = re.compile('[A-Za-z ][A-Za-z]')
    elem = cstrip(line[76:78])
    if not ELEMENT.match(elem):
        elem = cstrip(line[12:14])

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except:
        raise PDBParsingError("Error in coordinates")

    # Convert coordinates from A to nm if required
    if nanometers:
        x /= 10.0
        y /= 10.0
        z /= 10.0

    try:
        atom_number = int(line[6:11])
    except ValueError:
        raise PDBParsingError("Error in atom number")

    try:
        residue_number = int(line[22:26])
    except ValueError:
        raise PDBParsingError("Error in residue number")

    try:
        occupancy = float(line[54:60])
    except:
        occupancy = 1.0

    try:
        b_factor = float(line[60:66])
    except:
        b_factor = 0.0

    return Atom (atom_number,        # atomNumber
                 cstrip(line[12:16]),# atomName
                 cstrip(line[16]),   # atAlt
                 cstrip(line[17:21]),# resName
                 cstrip(line[21]),   # chainId --> if no chain ID is specified --> ''
                 residue_number,     # resNumber
                 line[26],           # resExt
                 x,                  # coordX
                 y,                  # coordY
                 z,                  # coordZ
                 occupancy,          # occ
                 b_factor,           # bFactor
                 elem)               # element

def read_pdb_from_file(file_in, nanometers=False, read_hetatm=False):
    """
    Reads and parses a given file_in PDB file.
    """
    lines = open(file_in).readlines()

    macroMol = MacroMolecule()
    last_chain = "dummy"
    last_residue = 9999999999
    last_residue_ext = "!" #dummy
    newMol = 0
    numline = 0
    numModels = 0

    for line in lines:
        numline += 1
        try:
            if line[0:6] == "MODEL ":
                numModels += 1
            if (line[0:4] == "ATOM") and numModels <= 1:
                line = line[:-1]
                a = read_atom_line (line, nanometers)
                if (last_chain != a.chain_id):
                    mol = Molecule(a.chain_id)
                    macroMol.add_molecule(mol)
                    newMol = 1
                if ((last_residue != a.residue_number) or
                    (last_residue_ext != a.residue_ext) or (newMol)):
                    if a.residue_type in Nucleotide.get_nucleic_acids_three_letters() or \
                                    a.residue_type in Nucleotide.get_nucleic_acids_two_letters() or \
                                    a.residue_type in Nucleotide.residues_to_rename.keys():
                        res = Nucleotide()
                        mol.add_residue(res)
                        mol.contains_DNA = True
                        newMol = 0
                    else:
                        res = AminoAcid()
                        mol.add_residue(res)
                        mol.contains_AAs = True
                        newMol = 0
                res.add_atom(a)
                last_residue = a.residue_number
                last_chain = a.chain_id
                last_residue_ext = a.residue_ext
            if read_hetatm and (line[0:6] == "HETATM") and numModels <= 1:
                line = line[:-1]
                a = read_atom_line(line, nanometers)
                if ((last_chain != a.chain_id) or (last_residue != a.residue_number) or (last_residue_ext != a.residue_ext)):
                    mol = Molecule(a.chain_id, is_cofactor_ion=True)
                    macroMol.add_molecule(mol)
                    res = Cofactor()
                    mol.add_residue(res)
                res.add_atom(a)
                last_residue = a.residue_number
                last_chain = a.chain_id
                last_residue_ext = a.residue_ext

        except Exception as e:
            raise PDBParsingError('Problem reading %s file at line %s: %s' %
                                  (str(file_in), str(numline), str(e)))
    return macroMol

class PyDockParameters():

    to_rename = {'HID': 'HIS'}
    renameResidues = True
    pdb_extension = '.pdb'
    ligand_extension = '_lig'
    receptor_extension = '_rec'
    reference_extension = '_ref'
    energy_extension = '.ene'
    rot_extension = '.rot'

def exit_with_error(msg):
    """
    Prints an error message and exits with failure
    """
    log.error(msg)
    log.debug(traceback.format_exc())
    raise SystemExit("%s terminated with error" % _app_tag)

class Logger(object):
    
    (ERROR, WARNING, INFO, DEBUG) = (0,1,2,4)
    
    def __init__(self, tag, file_name='', level=INFO):
        """
        Creates a Logger object.
        If file_name is set, output will be printed to a file called file_name.
        Default level is set to 'info', debug messages will be ignored.
        """
        if file_name:
            self._output = open(file_name, 'a')
        self._level = level
        self._tag = tag
            
    def set_level(self, level):
        """
        Sets logging level
        """
        if level >= Logger.ERROR and level <= Logger.DEBUG:
            self._level = level
    
    def debug(self, message=''):
        """
        Prints to selected output a debug message
        """
        if self._level >= self.DEBUG:
            self._write(message, "DEBUG")
    
    def info(self, message=''):
        """
        Prints to selected output an informative message
        """
        if self._level >= self.INFO:
            self._write(message, "INFO")
    
    def warning(self, message=''):
        """
        Prints to selected output a warning message
        """
        if self._level >= self.WARNING:
            self._write(message, "WARNING")
    
    def error(self, message=''):
        """
        Prints to selected output an error message
        """
        if self._level >= self.ERROR:
            self._write(message, "ERROR")
    
    def _write(self, message, level):
        """
        Outputs log info to the selected channel
        """
        out_message = "[%s] %s: %s" % (self._tag, level, message)
        try:
            self._output.write(out_message+os.linesep)
            self._output.flush()
        except:
            print (out_message)
            
    def __del__(self):
        """
        Frees file descriptor
        """
        try:
            self._output.close()
        except:
            pass

class LoggingManager(object):
    """
    Logger facility
    """
    _loggers = {}
    
    @staticmethod
    def get_logger(tag, file_name='', level=Logger.INFO):
        """
        Gets logger object identified by tag if exists or creates a new one
        """
        try:
            return LoggingManager._loggers[tag]
        except KeyError:
            try:
                level = int(os.environ['PYDOCK_LOGGING'])
            except KeyError:
                pass
            log = Logger(tag,file_name,level)
            LoggingManager._loggers[tag] = log  
            return log

def create_from_energy_table(dock_name, receptor_pdb, ligand_pdb, rot_file,
                             output_path, i_from, i_to, energy_table=None, prefix="",
                             includeHydrogens=False):
    """
    This will make pdb files of the docking solutions numbered from i_from to i_to.
    - i_from is the number of the first conformation to be saved as a PDB file
    - i_to is the number of the last conformation to be saved as a PDB file

    The conformation number can be the one in <dock_name>.rot file (if energy_table=='none')
    or the ranking of the docking solution as sorted in the <energy_table> file.

    """
    log = LoggingManager.get_logger('MakePDB')
    # read receptor and ligand files
    log.info("Reading pdb files...")
    recPDB = read_pdb_from_file(receptor_pdb)
    log.info("Receptor read")
    ligPDB = read_pdb_from_file(ligand_pdb)
    log.info("Ligand read")

    currentComplex = recPDB
    ligPDB.move_to_origin()

    rotTable = read_rot_file(rot_file)

    confToTableIndex = {}
    for i in range(len(rotTable)):
        confNum = rotTable[i][12]
        confToTableIndex[confNum] = i

    if energy_table != None:
        try:
            eneTab = Table.read(energy_table)
            # Configuration number as sorted by energy_table
            confList = [int(i) for i in eneTab["Conf"]] [i_from-1: i_to]
        except:
            raise PyDockError('Problem found reading %s file. Is it a valid energy table file?' % energy_table)
    else:
        confList = range(i_from, i_to+1)

    for i in confList:
        try:
            i_conf = confToTableIndex[i]
        except:
            log.warning("Upper index greater than number of conformations")
            break

        currentLig = ligPDB.clone()
        currentLig.rotation_and_translation(rotTable[i_conf])
        for mol in currentLig.molecules:
            currentComplex.add_molecule(mol)

        pdb_file_name = "%s%s_%s%s" % (prefix, dock_name, i, PyDockParameters.pdb_extension)
        write_pdb_to_file(currentComplex, os.path.normpath(output_path + os.sep + pdb_file_name),
                          includeHydrogens, to_rename=PyDockParameters.to_rename)

        for mol in currentLig.molecules:
            currentComplex.remove_molecule(mol)

        log.info("Conformation %s done " % i)

def search_ene_file(dock_name):
    rootdir = "."
    regex = re.compile(dock_name + '.*.ene')
    print (regex)
    for root, dirs, files in os.walk(rootdir):
        for file in files:
            if regex.match(file):
                print (file, rootdir)
                return(file)


def run_makepdb(dock_name, starting_rank_index, ending_rank_index, energy_table_file=None, prefix=""):
    """
    Creates desired PDB files from rotations file.

    Keyword arguments:
    starting_rank_index -- The first complex considered.
    ending_rank_index -- The last complex considered (included).
    energy_table_file -- Usually a .ene file generated in dockser step.
    prefix -- prefix used in the name of the PDB files generated.

    """
    # File names
    receptor_pdb = dock_name + PyDockParameters.receptor_extension + PyDockParameters.pdb_extension
    ligand_pdb = dock_name + PyDockParameters.ligand_extension + PyDockParameters.pdb_extension
    if not energy_table_file:
        energy_table_file = dock_name + PyDockParameters.energy_extension
    rot_file = dock_name + PyDockParameters.rot_extension

    create_from_energy_table(dock_name, receptor_pdb, ligand_pdb,
                                      rot_file, ".",
                                      starting_rank_index, ending_rank_index,
                                      energy_table_file, prefix=prefix)


if __name__ == "__main__":

    os.environ['PYDOCK_LOGGING'] = str(Logger.INFO)
    _app_tag = 'MakePDB'
    log = LoggingManager.get_logger(_app_tag, level=int(os.environ['PYDOCK_LOGGING']))
    log.info("Running pyDockMakePDB...")
    if (len(sys.argv) > 1) and (len(sys.argv) < 6) :
        dock_name = sys.argv[1]
        try:
            i_startingRank = int(sys.argv[2])
        except (IndexError, ValueError) as e:
            exit_with_error('No valid value for starting_rank')
        try:
            i_endingRank = int(sys.argv[3])
        except (IndexError, ValueError) as e:
            exit_with_error('No valid value for ending_rank')
        try:
            eneTab = sys.argv[4]
        except:
            eneTab = search_ene_file(dock_name)
        try:
            run_makepdb(dock_name, i_startingRank, i_endingRank, eneTab)
        except:
            exit_with_error("Check if file name: "+eneTab+ " is a correct ene file")
    else:
        exit_with_error("usage: %s dock_name starting_rank ending_rank [eneTab] (is not given we use the ene table in the folder)\n"
                        "[MakePDB] INFO: example: %s project902 1 10000 project902.ene" % (sys.argv[0],sys.argv[0]))
