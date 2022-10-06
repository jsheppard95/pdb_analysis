"""
exercise1.py

Author: Jackson Sheppard
Last Edit: 10/5/22
"""
import glob
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np


def ReadPDB(PdbFile):
    """
    Reads a protein PDB file and stores alpha carbon resiude names and their
    atomic positions.

    Parameters:
    -----------
    PdbFile - str - the string filename of a PDB file

    Returns:
    --------
    Pos, ResNames - tuple
        Pos - numpy.array - An (N, 3) dimensional NumPy array containing the
            position of the alpha carbon in each of the N amino acid residues
        ResNames - list - A length N list containing the three-latter code for
            the type of each of the corresponding amino acid residues
    """
    Pos = []
    ResNames = []
    ATOM_MIN_COL = 13  # Atom: N, CA, CB, ...
    ATOM_MAX_COL = 16
    RES_MIN_COL = 18  # Residue: SER, HIS, TYR, ...
    RES_MAX_COL = 20
    X_MIN_COL = 31  # x coord in Angstroms
    X_MAX_COL = 38
    Y_MIN_COL = 39  # y coord in Angstroms
    Y_MAX_COL = 46
    Z_MIN_COL = 47  # z coord in Angstroms
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
    hydrophillic residues.

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

    R_g^2 = (1/N) * sum_{i=1}^N |r_i - <r>|^2, <r> = (1/N) * sum_{i=1}^N r_i

    where r_i is the vector defininig the position of atom i and <r> is the
    average position.

    Parameters:
    -----------
    Pos - numpy.array - dimension (N, 3) array of atomic positions

    Returns:
    --------
    Rg - float - the corresponding radius of gyration for the positions
        defined in `Pos` with the same units as that of the coordinates in
        `Pos`
    """
    # Average atomic position: r_avg = np.mean(Pos, axis=0)
    # Rg^2 = the radius of gyration squared
    # Rg^2 = np.mean(np.linalg.norm(Pos - r_avg, axis=1)**2, axis=0)
    # Rg = np.sqrt(Rg^2)
    return np.sqrt(
        np.mean(np.linalg.norm(Pos - np.mean(Pos, axis=0), axis=1)**2, axis=0)
    )


def GetContacts(Pos):
    """
    Finds amino-acid contact pairs. A contact between two amino acid resiudes
    can be declared when the distance between their alpha carbons is less than
    a certain cutoff, here 9 Angstroms.
    Parameters:
    -----------
    Pos - numpy.array - dimension (N, 3) array of atomic positions

    Returns:
    --------
    Contacts - list - list of (i, j) tuples giving all residue-residue
        contacts in Pos, (i, j) indeces correspond to those in `Pos`, returns
        distinct pairs only, i.e returns (i1, j1) but not (j1, i1)
    """
    CONTACT_CUTOFF = 9.0  # contact cutoff distance in Angstroms
    pos_combos = list(combinations(Pos, r=2))
    pos_idxs = list(combinations(np.arange(len(Pos)), r=2))
    Contacts = []
    for i in range(len(pos_combos)):
        pair_dist = np.linalg.norm(pos_combos[i][0] - pos_combos[i][1])
        if pair_dist < CONTACT_CUTOFF:
            Contacts.append(pos_idxs[i])
    return Contacts


def nExtremaContacts(matrix, n, min_max, AMINO_ACIDS):
    """
    Finds the `n` smallest (`min_max`=`np.min`) or largest
    (`min_max`=`np.max`) values in `matrix`
    """
    extrem_contacts = np.empty(n, dtype=object)
    extrem_energies = np.zeros(n)
    matrix_search = np.copy(matrix)
    for i in range(n):
        next_extrem = min_max(matrix_search)
        next_idx = np.where(matrix == next_extrem)
        # Check if we got an off-diagonal element, then take 1st element only,
        # 2nd -> symmetric interaction
        if len(next_idx[0]) == 1:
            idx = [next_idx[0][0], next_idx[0][0]]
        else:
            idx = [next_idx[0][0], next_idx[0][1]]
        interact_str = AMINO_ACIDS[idx[0]] + "-" + AMINO_ACIDS[idx[1]]
        energy = matrix[idx[0], idx[1]]
        extrem_contacts[i] = interact_str
        extrem_energies[i] = energy
        # Update this energy in the search matrix so we can find the next
        # extrema
        # Do this by setting to 0, since lowest energies negative, highest
        # positive
        matrix_search[idx[0], idx[1]] = 0
        matrix_search[idx[1], idx[0]] = 0
    return extrem_contacts, extrem_energies


# Main function
# Make a list of all pdb files in the working path (using `glob`)
# For each pdb file:
# 1. Read the file
# 2. Compute Rg for all amino acid residues
# 3. Compute Rg,phobic for only the hydrophobic amino acid residues
# 4. Compute the ratio Rg,phobic/Rg
# 5. Print the filename, the number of residues, and items 2-4, in a single
#    line.

# Plot the following:
# - Rg,phobic and Rg vs. chain length (number of residues)
# - Rg,phobic/Rg vs. chain length (number of residues)

if __name__ == "__main__":
    AMINO_ACIDS = np.asarray(["ALA", "CYS", "PHE", "ILE", "LEU", "MET", "PRO",
                              "VAL", "TRP", "ASP", "GLU", "GLY", "HIS", "LYS",
                              "ASN", "GLN", "ARG", "SER", "THR", "TYR"])
    N_AA = len(AMINO_ACIDS)
    # N_PAIRS = N_AA!/(N_AA-2)!2! + N interactions
    #         = N_AA*(N_AA-1)/2 + N
    N_PAIRS = N_AA*(N_AA - 1)/2 + N_AA
    kB = 1.987204259e-3  # Boltzmann Constant in kcal/(mol K)
    T = 300  # Temperature in K, kB*T ~ 0.6 kcal/mol
    # Get all pdb files in working directory
    pdb_files = glob.glob("./**/*.pdb")

    # Read PDB files
    n_prot = len(pdb_files)
    Rg_alls = np.zeros(n_prot)
    Rg_phobics = np.zeros(n_prot)
    Rg_ratios = np.zeros(n_prot)
    chain_lengths = np.zeros(n_prot)
    # arrays to hold resiude total counts and contact counts used to compute
    # fractions for statistical interaction potentials
    # indeces -> those in `AMINO_ACIDS`
    aa_counts = np.zeros(N_AA)
    contact_counts = np.zeros((N_AA, N_AA))
    for i in range(len(pdb_files)):
        # read file
        pos, res_names = ReadPDB(pdb_files[i])
        # Compute Rg for all amino acids
        Rg_all = RadiusOfGyration(pos)
        # Get hydrophobic resiudes boolean away
        phobic_mask = ResHydrophobic(res_names)
        # Compute Rg for hydropghobic resiudes only
        Rg_phobic = RadiusOfGyration(pos[phobic_mask])
        # Compute ration Rg,phobic/Rg,all
        Rg_ratio = Rg_phobic/Rg_all
        # Get chain length
        n_res = len(res_names)
        # Output results and save data
        print(pdb_files[i], n_res, Rg_all, Rg_phobic, Rg_ratio)
        Rg_alls[i] = Rg_all
        Rg_phobics[i] = Rg_phobic
        Rg_ratios[i] = Rg_ratio
        chain_lengths[i] = n_res

        # Tally amino acid counts
        for res in res_names:
            res_idx = np.where(AMINO_ACIDS == res)
            aa_counts[res_idx] += 1

        # Get contacts and tally counts
        contacts = GetContacts(pos)
        for contact in contacts:
            # Get the residue names for the indeces return by GetContacts
            contact_res = [res_names[contact[0]], res_names[contact[1]]]
            # Get indeces for these resiudes defined by `AMINO_ACIDS`
            contact_indeces = np.where(np.in1d(AMINO_ACIDS, contact_res))[0]
            if contact_res[0] == contact_res[1]:
                # handle a contact between two of the same residue types
                contact_indeces = np.concatenate(
                    (contact_indeces, contact_indeces)
                )
            # Tally contact counts
            contact_counts[contact_indeces[0], contact_indeces[1]] += 1

    # Get AA fractions f_k over entire data set
    # indeces -> those in AMINO_ACIDS
    f_k = aa_counts / np.sum(aa_counts)

    # Get contact fractions
    c_kl = contact_counts / np.sum(contact_counts)
    # Stack Overflow magic to fill in zeros with their symmetric contact
    c_kl = np.maximum(c_kl, c_kl.transpose())

    # Compute interaction potentials
    # 210 distinct elements:
    # N = 20, computing symmetric elements of the 400 element matrix here
    # u_kl = -k_B*T*log(c_kl/(f_k*f_l))
    # Build matrix F_kl = |f_k><f_l|, then can do element-wise division in log
    F_kl = np.outer(f_k, f_k)
    # Calculate potential energy
    u_kl = -kB*T*np.log(c_kl/F_kl)
    # Normalize s.t mean potential energy is 0
    # Compute mean considering only the 210 distinct interactions, the 210
    # upper elements, and subract from each element in u_kl
    u_avg = np.sum(np.triu(u_kl)) / N_PAIRS
    u_kl_norm = u_kl - u_avg

    # Report 5 most and 5 least favorable interaction potentials
    # Most favorable -> Most negative -> 5 smallest of u_kl_norm
    N_MIN = 5
    min_contacts, min_energies = nExtremaContacts(u_kl_norm, N_MIN, np.min,
                                                  AMINO_ACIDS)
    print("")
    print("Lowest Interaction Energies:")
    print("RES-RES : Energy (kcal/mol)")
    print("---------------------------")
    for i in range(N_MIN):
        print(min_contacts[i] + " : " + str(min_energies[i]))
    print("")

    N_MAX = 5
    max_contacts, max_energies = nExtremaContacts(u_kl_norm, N_MAX, np.max,
                                                  AMINO_ACIDS)
    print("")
    print("Highest Interaction Energies:")
    print("RES-RES : Energy (kcal/mol)")
    print("---------------------------")
    for i in range(N_MAX):
        print(max_contacts[i] + " : " + str(max_energies[i]))
    print("")

    # Plot data
    # Radius of Gyration (All Resiudes and Hydrophobic Resiudes Only)
    f1, ax1 = plt.subplots()
    ax1.plot(chain_lengths, Rg_alls, ".", label=r"All Residues, $R_g$")
    ax1.plot(chain_lengths, Rg_phobics, ".",
             label=r"Hydrophobic Residues Only, $R_{g,phobic}$")
    ax1.legend()
    ax1.set_xlabel("Chain Length (Number of Residues)")
    ax1.set_ylabel(r"Radius of Gyration, $\AA$")
    f1.show()

    # Rg,phobic/Rg,all0
    f2, ax2 = plt.subplots()
    ax2.plot(chain_lengths, Rg_ratios, ".")
    ax2.set_xlabel("Chain Length (Number of Residues)")
    ax2.set_ylabel(r"Ratio Hydrophobic to All, $R_{g,phobic}/R_g$")
    f2.show()

    # potential energy heat map
    f3, ax3 = plt.subplots(figsize=(8.0, 6.0))
    im = ax3.imshow(u_kl_norm, cmap="bone_r")
    ax3.set_xticks(range(N_AA))
    ax3.set_yticks(range(N_AA))
    ax3.set_xticklabels(AMINO_ACIDS, rotation=90)
    ax3.set_yticklabels(AMINO_ACIDS)
    ax3.set_xlabel("Residue")
    ax3.set_ylabel("Residue")
    ax3.set_title("Amino Acid Contact Statistical Interaction Potential\n"
                  r"$u_{kl}=-k_BT\ln{\frac{c_{kl}}{f_kf_l}}$,  "
                  f"T={T} K")
    f3.colorbar(im, label="Potential Energy (kcal/mol)")
    f3.show()

    plt.show()
