CAG: Identigying protein complexes in Protein Protein Interaction network by Core -attachment approach incorporating Gene Expression Profile.

Author: Seketoulie Keretsu

How to use CAG:

step1: compile Protein.java
step2: compile Cag.java
step3: Run   java Cag <PPI_network_data>  <Gene_expression_data>  <reference_data> <density_threshold> <similarity_threshold>
[ex.    java Cag collins2007.txt gene.txt  sgd.txt  0.7  0.1]


[The input files should be kept in the working directory]

parameter values:

density_threshold: A core can be considered only if it has a density greater than or equal to the density_threshold. The default value is 0.6
similarity_threshold:  A protein can be added to form a core only if it has a similarity value with the seed protein greater than the similarity_threshold value. The default value is 0.1.

PPI data:  The PPI network data should contain weighted interactions where the interactions are given by  Protein [space] Protein [space] weight (eg. proein1 Protein2  1.0 
Gene expression data: contains expression values of genes with time course.
sgd : a collection of complexes with each complexes containg protein names seperated by tap or space.  Each complex is given in a seperate line.




