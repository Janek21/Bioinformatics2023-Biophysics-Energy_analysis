class ResiduesDataLib(): # Creating the class

    def __init__(self, fname): # initialitation of the class with closed variables
        self.residue_data = {}

        try:
            fh = open(fname, "r") # Oppening the file for reading

        except OSError: # Searching for errors and if found fiving msg error and quitting
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)

        for line in fh: # Reading the file line per line

            if line[0] == '#': # Checking for the header and continuing
                continue

            data = line.split() # splitting the data
            r = Residue(data) # calling the class Residue to get the residue

            self.residue_data[r.id] = r # Assigning the residue to a local var of the class

        self.nres = len(self.residue_data) # Getting the lenght of the variable and assigining it to a local variable

    def get_params(self, resid, atid): # Defining function to get the parameters
        atom_id = resid + ':' + atid # getting the atom id

        if atom_id in self.residue_data: # Checking if the atom id is inside the residue data and if so returning it
            return self.residue_data[atom_id]
        
        else: # If not in the data getting an error
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue(): # Defining the classs
    
    def __init__(self,data): # Initializing the class
        self.id     = data[0]+':'+data[1] # getting id
        self.at_type = data[2] # 
        self.charge  = float(data[3]) # gettig charge
        
class AtomType(): # Defining class

    def __init__(self, data): # Initializing the class
        # Getting some info and storing it in local variables
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name): # initializing the class
        self.at_types = {}

        try: # oppening the file for reading
            fh = open(file_name, "r")

        except OSError: # checking for errors and quitting
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)

        for line in fh: # Reading the file

            if line[0] == '#': # checking for the headder
                continue

            data = line.split() # splitting the line
            self.at_types[data[0]] = AtomType(data) # getting the atom type by calling the class AtomType

        self.ntypes = len(self.at_types) # assigning the length of the variable to a local var
        fh.close() # clossing