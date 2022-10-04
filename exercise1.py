"""
exercise1.py

Author: Jackson Sheppard
Last Edit: 10/3/22
"""
import glob
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
    ATOM_MIN_COL = 13
    ATOM_MAX_COL = 16
    RES_MIN_COL = 18
    RES_MAX_COL = 20
    X_MIN_COL = 31
    X_MAX_COL = 38
    Y_MIN_COL = 39
    Y_MAX_COL = 46
    Z_MIN_COL = 47
    Z_MAX_COL = 54
    with open(PdbFile, "r") as f:
        # Read header
        f.readline()
        # Start reading data
        line = f.readline()
        while line != "TER\n":
            # add alpha carbon residues only
            if "CA" in line[ATOM_MIN_COL-1: ATOM_MAX_COL]:
                ResNames.append(line[RES_MIN_COL-1: RES_MAX_COL])
                Pos.append([float(line[X_MIN_COL-1: X_MAX_COL]),
                            float(line[Y_MIN_COL-1: Y_MAX_COL]),
                            float(line[Z_MIN_COL-1: Z_MAX_COL])])
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


def RadiusOfGyration(Pos):
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
    # Compute average atomic position
    r_avg = np.mean(Pos, axis=0)
    # Compute Rg^2, the radius of gyration squared
    Rg_sq = np.mean(np.linalg.norm(Pos - r_avg, axis=1)**2, axis=0)
    Rg = np.sqrt(Rg_sq)
    return Rg


# Main function
# Make a list of all pdb files in the working path (using `glob`)
# For each pdb file:
## 1. Read the file
## 2. Compute Rg for all amino acid residues
## 3. Compute Rg,phobic for only the hydrophobic amino acid residues
## 4. Compute the ratio Rg,phobic/Rg
## 5. Print the filename, the number of residues, and items 2-4, in a single
#     line.

# Plot the following:
## Rg,phobic and Rg vs. chain length (number of residues)
## Rg,phobic/Rg vs. chain length (number of residues)

if __name__ == "__main__":
    # Get all pdb files in working directory
    pdb_files = glob.glob("./**/*.pdb")
    pos, res_names = ReadPDB("proteins/T0639-D1.pdb")
    # Iterate through PDB files
    for fname in pdb_files:
        # read file
        pos, res_names = ReadPDB(fname)
        # Compute Rg for all amino acids
        Rg_all = RadiusOfGyration(pos)
        phobic_mask = ResHydrophobic(res_names)
        Rg_phobic = RadiusOfGyration(pos[phobic_mask])
        Rg_ratio = Rg_phobic/Rg_all
        print(fname, len(res_names), Rg_all, Rg_phobic, Rg_ratio)