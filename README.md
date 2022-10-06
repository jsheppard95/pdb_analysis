# pdb_analysis

Jackson Sheppard\
CH E 210D, Exercise 1\
10/05/22

Here we present an analysis of 309 protein sequences, computing structural
quantities such as the radius of gyration and statistical interaction
potentials. The input sequences are stored in
[The Protein Data Bank](www.pdb.org) `*.pdb` format and reside in the
`proteins/` directory at the root of this repository. Each file consists of
a protein X-ray crystal structure and thus includes atomic `X`, `Y`, `Z`
coordinates of each amino acid (residue) comprising the protein structure. In
this analysis, we read the sequences of our 309 proteins and compute the
radius of gyration considering all residues along with that for hydrophobic
residues only. We then plot both of these quantities along with their ratio
against the total number of residues in the structure. Finally, we compute and
visualize the "statistical" interation potential for the 20 amino acid types
by considering amino acid contacts present in this data set.

## Installation and Usage
Clone [this repository](https://github.com/jsheppard95/pdb_analysis) and
navigate to its root. Install dependencies from the `enironment.yml` file
using `conda`:

```
conda env create --name envname --file=environment.yml
```

If instead files are downloaded individually, ensure the `*.pdb` files are in
the relative path of working directory and that necessary dependencies are
installed. Run the code to generate plots as follows:

```
$ python exercise1.py
```

![Rg_length](output/Rg_length.png)

![Rg_ratio_length](output/Rg_ratio_length.png)

![contact_potential](output/contact_potential.png)