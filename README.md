# pdb_analysis

Jackson Sheppard\
CH E 210D, Exercise 1\
10/05/22\

Here we present an analysis of 309 protein sequences, computing structural
quantities such as the radius of gyration and statistical interaction
potentials. The input sequences are stored in
[The Protein Data Bank](www.pdb.org) `*.pdb` format and reside in the
`proteins/` directory at the root of this repository. Each file consists of
a protein X-ray crystal structure and thus includes atomic `X`, `Y`, `Z`
coordinates of each amino acid (residue) comprising the protein structure. In
this analysis, we read the sequences of our 309 proteins and compute the
radius of gyration, $R_g$, defined as follows:

$$
R_g^2 = \frac{1}{N}\sum_{i=1}^N |\vec{r}_i - \langle \vec{r} \rangle|^2, \quad \langle \vec{r} \rangle = \frac{1}{N}\sum_{i=1}^N \vec{r}_i
$$

 and hydrophobic residues
only, $R_{g,phobic}$. We then plot both of these quantities along with their
ratio against the total number of residues in the structure. Finally, we
compute and visualize the "statistical" interation potential for the 20 amino
acid types by considering amino acid contacts present in this data set.

![Rg_length](output/Rg_length.png)

![Rg_ratio_length](output/Rg_ratio_length.png)

![contact_potential](output/contact_potential.png)