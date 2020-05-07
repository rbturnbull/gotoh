# gotoh_counts
Aligns two sequences and returns the number of characters which match, mismatch, open gaps or extend gaps.

Perform a global sequence alignment (Needleman-Wunsch) with affine gap penalties (Gotoh) on seq and and 2. It returns a tuple with the number of characters that match, mismatch, open gaps and extend gaps.
The scores/penalties are given as arguments and the defaults correspond to basic Levenshtein distance.

Based almost entirely on Brent Pedersen’s nwalign: https://bitbucket.org/brentp/biostuff/

