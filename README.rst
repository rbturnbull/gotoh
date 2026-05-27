gotoh
=====

.. start-badges

|pypi badge| |testing badge|

.. |pypi badge| image:: https://img.shields.io/pypi/v/gotoh?color=blue
    :target: https://pypi.org/project/gotoh/

.. |testing badge| image:: https://github.com/rbturnbull/gotoh/actions/workflows/pipeline.yml/badge.svg
    :target: https://github.com/rbturnbull/gotoh/actions

.. end-badges

To install, run:

.. code-block:: bash

    pip install gotoh


This package performs a global sequence alignment (Needleman-Wunsch) with affine gap
penalties (Gotoh). It returns a tuple with the number of characters that match,
mismatch, open gaps and extend gaps. The scores/penalties are given as
arguments and the defaults correspond to the basic Levenshtein distance.

It can also produce a multiple sequence alignment.

Documentation and code clean up to come.

Credits
==================================

.. start-credits

`Robert Turnbull <https://robturnbull.com>`_

Based on Brent Pedersen's nwalign which was originally posted at
https://bitbucket.org/brentp/biostuff/ but now has disappeared.

For more information, see Robert Turnbull's book `Codex Sinaiticus Arabicus and Its Family <https://brill.com/display/title/69300>`_.

.. end-credits
