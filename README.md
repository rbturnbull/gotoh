# gotoh

Performs a global sequence alignment (Needleman-Wunsch) with affine gap penalties (Gotoh) on seq and and 2. It returns a tuple with the number of characters that match, mismatch, open gaps and extend gaps.
The scores/penalties are given as arguments and the defaults correspond to basic Levenshtein distance.

It can also produce a multiple sequence alignment.

Based on Brent Pedersen’s nwalign: https://bitbucket.org/brentp/biostuff/

For more information, see chapter 7 of Robert Turnbull's thesis 'The Textual History of Codex Sinaiticus Arabicus and its Family'.
