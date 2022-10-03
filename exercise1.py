"""
exercise1.py

Author: Jackson Sheppard
Last Edit: 10/3/22
"""
import numpy as np

def ReadPDB(PdbFile):
    """
    Function to a read a protein PDB file.
    
    Parameters:
    -----------
    PdbFile - str - the string filename of a PDB file
    
    Returns:
    --------
    Pos, ResNames - tuple
        Pos - numpy.array - An (N, 3) dimensional NumPy array containing the
            position of the alpha carbon in each of the N amino acid residues
        ResNames - list - A length N list containing the three-latter code for
            the type of each of the corresponding amino acid residues.
    """
    Pos = []
    ResNames = []
    ATOM_IDX = 2
    RES_IDX = 3
    X_IDX = 5
    Y_IDX = 6
    Z_IDX = 7
    with open(PdbFile, "r") as f:
        # Read header
        f.readline()
        # Start reading data
        line = f.readline()
        while line != "TER\n":
            line_split = line.split()
            if line_split[ATOM_IDX] == "CA":  # add alpha carbon residues only
                ResNames.append(line_split[RES_IDX])
                Pos.append([float(line_split[X_IDX]),
                            float(line_split[Y_IDX]),
                            float(line_split[Z_IDX])])
            line = f.readline()
    Pos = np.asarray(Pos)
    return Pos, ResNames


def ResHydrophobic(ResNames):
    """
    Determines which residue codes in `ResNames` correspond to hydrophic or
    hydrophillic residues

    Parameters:
    -----------
    ResNames - list - length-N list of three letter residue codes

    Returns:
    --------
    IsPhobic - numpy.array - length-N boolean array of `True` or `False`
        values indicated whether each residue code in `ResNames` is a
        hydrophobic (`True`) or hydrophillic (`False`) residue
    """
    ResNames = np.asarray(ResNames)
    HYDROPHOBIC_RES = np.asarray(["ALA", "CYS", "PHE", "ILE", "LEU", "MET",
                                  "PRO", "VAL", "TRP"])
    return np.isin(ResNames, HYDROPHOBIC_RES)


def RediusOfGyration(Pos):
    """
    Computes the radius of gyration for a set of atomic positions. The radius
    of gyration for a collection of N coordinates is given by:

    R_g^2 = (1/N) * \sum_{i=1}^N |r_i - <r>|^2, <r> = (1/N) * \sum_{i=1}^N r_i

    where r_i is the vector defininig the position of atom i and <r> is the
    average position.

    Parameters:
    -----------
    Pos - numpy.array - dimension (N, 3) array of atomic positions

    Returns:
    --------
    Rg - float - the corresponding radius of gyration for the positions
        defined in `Pos`
    """
    pass


# Main function
# Make a list of all pdb files in the working path (using `glob`)
# For each pdb file:
## 1. Read the file
## 2. Compute Rg for all amino acid residues
## 3. Compute Rg,phobic for only the hydrophobic amino acid residues
## 4. Compute the ration Rg,phobic/Rg
## 5. Print the filename, the number of residues, and items 2-4, in a single
#     line.

# Plot the following:
## Rg,phobic and Rg vs. chain length (number of residues)
## Rg,phobic/Rg vs. chain length (number of residues)

if __name__ == "__main__":
    pos, res_names = ReadPDB("proteins/T0639-D1.pdb")
    print(pos[:10])
    print(res_names[:10])
    print(ResHydrophobic(res_names[:10]))