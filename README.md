# gotoh

![pipline](https://github.com/rbturnbull/gotoh/actions/workflows/pipeline.yml/badge.svg)
[<img src="https://img.shields.io/badge/code%20style-black-000000.svg">](<https://github.com/psf/black>)

Performs a global sequence alignment (Needleman-Wunsch) with affine gap penalties (Gotoh). It returns a tuple with the number of characters that match, mismatch, open gaps and extend gaps. The scores/penalties are given as arguments and the defaults correspond to the basic Levenshtein distance.

It can also produce a multiple sequence alignment.

Based on Brent Pedersen’s nwalign: https://bitbucket.org/brentp/biostuff/

For more information, see chapter 7 of Robert Turnbull's thesis 'The Textual History of Codex Sinaiticus Arabicus and its Family'.

Documentation and code clean up to come.