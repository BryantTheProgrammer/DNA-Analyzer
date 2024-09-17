# DNA Analyzer
This project involves writing a program to perform a common analysis on homologous DNA strands, where the objective is to output the original DNA strands, the profile (number of A, C, G, and T bases for each position), and the consensus strand. The consensus strand represents the most likely common ancestor based on the input DNA strands.

## Key Features
**DNA Input Parsing:** Reads homologous DNA strands of the same length from a data file (dnastrands.txt).

**Profile Calculation:** Counts occurrences of each base (A, C, G, T) in each position of the DNA strands.

**Consensus Strand Generation:** Creates a consensus strand based on the most frequent base at each position.

**Error Handling:** Detects unequal strand lengths and file reading issues, terminating the program with appropriate error messages.

**Modular Design:** Utilizes object-oriented principles, reusing existing DNAStrand and DNAStrandSet code, and making use of a flexible vector of vectors for profile computation.
