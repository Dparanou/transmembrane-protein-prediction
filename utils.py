class PDB(object):
    """Load a pdb or cif file."""
    def __init__(self, filename, ignore_waters=False, ignore_other_HETATM=False, MultiModels=False):
        self.natoms, self.water_coords, self.HETATM_coords, ciff = 0, [], [], False
        if filename.split('.')[-1] == 'cif' or filename.split('.')[-2] == 'cif':
            ciff = True
        #Open file
        if filename.split('.')[-1] == 'gz':
            with gzip.open(filename, 'rt', encoding='utf-8') as file:
                f = file.readlines()
        else:
            with open(filename, 'r') as file:
                f = file.readlines()
        lineslen = len(f)
        for i in range(lineslen):
            if f[i][:4] == 'ATOM':
                break
        for j in range(i-1, lineslen):
            line = f[j]
            sline = line.split()
            if ciff:
                if not ignore_waters and sline[0] == 'HETATM' and sline[2] == 'O' and ((sline[5]=='HOH') or (sline[5]=='TIP')):
                    #Keep coordinates of water Oxygens.
                    self.water_coords.append([float(sline[10]), float(sline[11]), float(sline[12])])
                    continue
                #Return structure only of Model 1 if there are multiple Models
                if not MultiModels and sline[-1] != '1':
                    break
            else:
                if line[:6] == 'HETATM':
                    if not ignore_waters and line[13] == 'O' and ((line[17:20]=='HOH') or (line[17:20]=='TIP')):
                        self.water_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        continue
                    if not ignore_other_HETATM:
                        self.HETATM_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        continue
                #Return structure only of Model 1 if there are multiple Models
                if sline[0] == 'MODEL' and not ciff and not MultiModels and sline[1] != '1':
                    break
            if line[:4] != 'ATOM':
                # skip other lines
                continue
            self.natoms += 1
        self.SSE = np.zeros((self.natoms), dtype=np.dtype((str,1)))
        self.resolution = 0.
        self.SSEraw = []
        self.atomnum = np.zeros((self.natoms),dtype=int)
        self.atomname = np.zeros((self.natoms),dtype=np.dtype((str,3)))
        self.atomalt = np.zeros((self.natoms),dtype=np.dtype((str,1)))
        self.resname = np.zeros((self.natoms),dtype=np.dtype((str,3)))
        self.resnum = np.zeros((self.natoms),dtype=int)
        self.resalt = np.zeros((self.natoms),dtype=np.dtype((str,1)))
        self.chain = np.zeros((self.natoms),dtype=np.dtype((str,2)))
        self.coords = np.zeros((self.natoms, 3))
        self.water_coords = np.array(self.water_coords)
        self.HETATM_coords = np.array(self.HETATM_coords)
        self.occupancy = np.zeros((self.natoms))
        self.b = np.zeros((self.natoms))
        self.atomtype = np.zeros((self.natoms),dtype=np.dtype((str,2)))
        b_strands = False
        atom = 0
        if ciff:
            for j in range(lineslen):
                line = f[j]
                sline = line.split()
                #Protein atoms info.
                if line[:4] == 'ATOM' and self.natoms != atom:
                    self.atomnum[atom] = int(float(sline[1]))
                    self.atomname[atom] = sline[3]
                    self.atomalt[atom] = line[22]
                    self.resname[atom] = sline[5]
                    atomtype = sline[2]
                    if len(atomtype) == 2:
                        atomtype = atomtype[0].upper() + atomtype[1].lower()
                    self.atomtype[atom] = atomtype
                    self.resnum[atom] = int(float(sline[-5]))
                    self.chain[atom] = sline[6]
                    self.coords[atom, :] = float(sline[10]), float(sline[11]), float(sline[12])
                    self.occupancy[atom] = float(sline[13])
                    self.b[atom] = float(sline[14])
                    atom += 1
                    continue
                #Here you can add more self objects from Header.
                #Resolution info in cif format files.
                if line[:25] == '_reflns.d_resolution_high' or line[:33] == '_em_3d_reconstruction.resolution ':
                    self.resolution = sline[1]
                    continue
                #Helix info in cif format files.
                if line[:6] == 'HELX_P':
                    self.SSEraw.append([range(int(sline[-7]),int(sline[-4])+1), sline[-8], 'H'])
                    continue
                #Beta sheet info in cif format files.
                if line[:35] == '_struct_sheet_range.end_auth_seq_id':
                    b_strands = True
                    continue
                if b_strands:
                    if line[:1] == '#':
                        b_strands=False
                        continue
                    self.SSEraw.append([range(int(sline[-4]),int(sline[-1])+1), sline[-2], 'S'])
                    continue
                #Crystal info in cif format files.
                if line[:5] == '_cell':
                    if sline[0] == '_cell.length_a':
                        self.cella = float(sline[1])
                    if sline[0] == '_cell.length_b':
                        self.cellb = float(sline[1])
                    if sline[0] == '_cell.length_c':
                        self.cellc = float(sline[1])
                    if sline[0] == '_cell.angle_alpha':
                        self.cellalpha = float(sline[1])
                    if sline[0] == '_cell.angle_beta':
                        self.cellbeta = float(sline[1])
                    if sline[0] == '_cell.angle_gamma':
                        self.cellgamma = float(sline[1])
                    continue
        else:
            for j in range(lineslen):
                line = f[j]
                sline = line.split()
                #Protein atoms info.
                if line[:4] == 'ATOM' and self.natoms != atom:
                    self.atomnum[atom] = int(float(sline[1]))
                    self.atomname[atom] = line[12:16].strip()
                    self.atomalt[atom] = line[16]
                    self.resname[atom] = line[17:21].strip()
                    atomtype = sline[-1]
                    if len(atomtype) == 2:
                        atomtype = atomtype[0].upper() + atomtype[1].lower()
                    self.atomtype[atom] = atomtype
                    self.resnum[atom] = int(float(line[22:26]))
                    self.resalt[atom] = line[26]
                    self.chain[atom] = line[21]
                    self.coords[atom, :] = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    self.occupancy[atom] = float(line[56:60])
                    self.b[atom] = float(line[60:66])
                    #self.charge[atom] = line[78:80].strip('\n')
                    #self.nelectrons[atom] = electrons.get(self.atomtype[atom].upper(),6)
                    atom += 1
                    continue
                #Here you can add more self objects from Header.
                #Resolution info in pdb format files.
                if line[:22] == 'REMARK   2 RESOLUTION.':
                    self.resolution = sline[3]
                    continue
                #Helix info in pdb format files.
                if line[:5] == 'HELIX':
                    self.SSEraw.append([range(int(line[20:25]),int(line[32:37])+1), line[19], 'H'])
                    continue
                #Beta sheet info in pdb format files.
                if line[:5] == 'SHEET':
                    self.SSEraw.append([range(int(line[22:26]),int(line[33:37])+1), line[21], 'S'])
                    continue
                #Crystal info in pdb format files.
                if line[:6] == 'CRYST1':
                    self.cella, self.cellb, self.cellc = float(sline[1]), float(sline[2]), float(sline[3])
                    self.cellalpha, self.cellbeta, self.cellgamma = float(sline[4]), float(sline[5]), float(sline[6])
                    continue
        #Return structure when every atom and water atom have been searched.
        try:
            self.resolution = float(self.resolution)
        except:
            pass
        self.SSEraw = np.array(self.SSEraw, dtype=object)
        return

    def remove_waters(self):
        idx = np.where((self.resname=="HOH") | (self.resname=="TIP"))
        self.remove_atoms_from_object(idx)

    def remove_by_atomtype(self, atomtype):
        idx = np.where((self.atomtype==atomtype))
        self.remove_atoms_from_object(idx)

    def remove_by_atomname(self, atomname):
        idx = np.where((self.atomname==atomname))
        self.remove_atoms_from_object(idx)

    def remove_by_atomnum(self, atomnum):
        idx = np.where((self.atomnum==atomnum))
        self.remove_atoms_from_object(idx)

    def remove_by_resname(self, resname):
        idx = np.where((self.resname==resname))
        self.remove_atoms_from_object(idx)

    def remove_by_resnum(self, resnum):
        idx = np.where((self.resnum==resnum))
        self.remove_atoms_from_object(idx)

    def remove_by_chain(self, chain):
        idx = np.where((self.chain==chain))
        self.remove_atoms_from_object(idx)

    def remove_atoms_from_object(self, idx):
        mask = np.ones(self.natoms, dtype=bool)
        mask[idx] = False
        self.atomnum = self.atomnum[mask]
        self.atomname = self.atomname[mask]
        self.atomalt = self.atomalt[mask]
        self.resalt = self.resalt[mask]
        self.resname = self.resname[mask]
        self.resnum = self.resnum[mask]
        self.chain = self.chain[mask]
        self.coords = self.coords[mask]
        self.occupancy = self.occupancy[mask]
        self.b = self.b[mask]
        self.atomtype = self.atomtype[mask]
        self.SSE = self.SSE[mask]
        #self.charge = self.charge[mask]
        #self.nelectrons = self.nelectrons[mask]
        self.natoms = len(self.atomnum)

    def rearrange_resalt(self):
        #Keep indexes where you have added residues
        ind = np.where(self.resalt!=' ')[0]
        if not len(ind):
            return
        #Find chains of these added residues
        resalt_ch = self.chain[ind]
        #Find unique chains
        diff_chains = np.unique(resalt_ch)
        #For each unique chain id
        for ch in diff_chains:
            #Keep indexes of added residues only for your chain
            ind_resalt_ch = ind[np.where(resalt_ch==ch)]
            #Keep last index of your chain's residues
            last_ch_ind = np.where(self.chain==ch)[0][-1]+1
            #For each atom in chain in added residue
            for atom in ind_resalt_ch:
                #Find where added residue starts
                if self.resalt[atom-1]!=self.resalt[atom]:
                    #Add to all consquent residue +1 number
                    self.resnum[atom:last_ch_ind] += 1
                #If added residues are before residue number, change the numbering of the original residue +1
                if self.resalt[atom+1] == ' ' and self.resnum[atom+1] == self.resnum[atom]:
                    self.resnum[atom+1:last_ch_ind] += 1
        #Remove alternative residues indexes
        self.resalt[ind] == ' '

SSEdict = {'L':0., 'H':1., 'S':2.}
inverseSSEdict = {0.:'Loops', 1.:'Helices', 2.:'Sheets'}
def GetSSE(pdb, specific_chain = None):
    for i in range(len(pdb.resnum)):
        if specific_chain and pdb.chain[i] not in specific_chain:
            continue
        for line in pdb.SSEraw:
            if pdb.chain[i] == line[1] and pdb.resnum[i] in line[0]:
                pdb.SSE[i] = line[2]
                break
        if pdb.SSE[i] == '':
            pdb.SSE[i] = 'L'


# read file 
pdb_f = PDB("./data/" + file)

# get secondary structure elements
GetSSE(pdb_f)
pdb_f.SSE
