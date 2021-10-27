# gotoh

![pipline](https://github.com/rbturnbull/gotoh/actions/workflows/pipeline.yml/badge.svg)

Performs a global sequence alignment (Needleman-Wunsch) with affine gap penalties (Gotoh). It returns a tuple with the number of characters that match, mismatch, open gaps and extend gaps. The scores/penalties are given as arguments and the defaults correspond to the basic Levenshtein distance.

It can also produce a multiple sequence alignment.

Based on Brent Pedersen’s nwalign which was originally posted at https://bitbucket.org/brentp/biostuff/ but now has disappeared.

For more information, see chapter 7 of Robert Turnbull's thesis 'The Textual History of Codex Sinaiticus Arabicus and its Family'.

Documentation and code clean up to come.