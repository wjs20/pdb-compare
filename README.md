
# About

A tool for registering groups of pdb protein structures into a common reference frame so that they can be compared visually and by automated tools.

It is often useful to jointly analyse multiple crystal structures containing the same of related proteins in a shared reference frame. For example, when analyzing antibody-antigen interactions it can be useful to register all copies of the ligand onto one another, so that different antibody binding modes can be clearly visualized.

This can be difficult however, as crystal structures produced in different studies/by different authors often do not follow standard chain naming schemes. For example the same chain may be referred to as chain A in one file and chain B in another, while both referring to the same underlying protein. In addition, the various copies of a protein structure may each be missing different sets of residues due to modelling issues or ommission of domains/parts of domains. This means that simply identifying the set of unique sequences in your set of structure files will not reliably reveal the latent structural entities. Finally, structures are typically not registered in a common reference frame (corresponding residues do not occupy the same physical coordinates) making comparison visually and programatically difficult.

`pdb-compare` solves these problems by 1) identifying the latent structural entities present in a group of pdb files and 2) aligning them into a common reference frame. These functionalities are exposed through the `cluster` and `align` subcommands respectively.

# Usage

Select the `cluster` subcommand and pass it a plain text file containing the paths to your pdb files, along with clustering configuration options.

Below we select the chains belonging to cluster 0. This will always be the largest cluster, as the clusters are aranged in reverse order with respect to cluster size.
```
$ uv run pdb-compare cluster -i data/group_35.txt -d clusters -t 0.8 -m local | egrep '[[:space:]]0$'
8wrz.pdb|A	0
8ja3.pdb|B	0
8tin.pdb|A	0
8til.pdb|A	0
8vj9.pdb|A	0
8goo.pdb|A	0
8goc.pdb|A	0
8j8r.pdb|A	0
8go9.pdb|F	0
7df9.pdb|A	0
7dfc.pdb|A	0
8gp3.pdb|A	0
8i0q.pdb|A	0
8jaf.pdb|A	0
7r0c.pdb|C	0
7dfb.pdb|A	0
8go8.pdb|A	0
8i0n.pdb|A	0
8zyu.pdb|B	0
7r0j.pdb|C	0
9bt8.pdb|B	0
8j8v.pdb|F	0
7dfa.pdb|A	0
8j97.pdb|A	0
9cx9.pdb|A	0
8wu1.pdb|C	0
8as3.pdb|A	0
6ni2.pdb|B	0
9ii2.pdb|A	0
9ii3.pdb|A	0
8j8z.pdb|A	0
6tko.pdb|B	0
8i10.pdb|A	0
8as2.pdb|A	0
8i0z.pdb|A	0
9e82.pdb|A	0
4jqi.pdb|A	0
8tio.pdb|A	0
6pwc.pdb|A	0
9cx3.pdb|B	0
6u1n.pdb|C	0
8j9k.pdb|C	0
8tii.pdb|A	0
```

We can create a mapping file to be used as input the the alignment subcommand by removing the cluster column name. This file should have the format <filename>:<chain_id>

```
$ uv run pdb-compare cluster -i data/group_35.txt -d clusters -t 0.8 -m local | egrep '[[:space:]]0$' | cut -f1 > mapping-file
```
