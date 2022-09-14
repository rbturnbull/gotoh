# -*- coding: utf-8 -*-
from setuptools import setup
import numpy as np

packages = \
['gotoh']

package_data = \
{'': ['*']}

install_requires = \
['numpy>=1.23.3,<2.0.0']

setup_kwargs = {
    'name': 'gotoh',
    'version': '0.1.2',
    'description': 'Sequence Alignment with different penalties for opening gaps and extending them.',
    'long_description': "# gotoh\n\n![pipline](https://github.com/rbturnbull/gotoh/actions/workflows/pipeline.yml/badge.svg)\n\nPerforms a global sequence alignment (Needleman-Wunsch) with affine gap penalties (Gotoh). It returns a tuple with the number of characters that match, mismatch, open gaps and extend gaps. The scores/penalties are given as arguments and the defaults correspond to the basic Levenshtein distance.\n\nIt can also produce a multiple sequence alignment.\n\nBased on Brent Pedersen’s nwalign which was originally posted at https://bitbucket.org/brentp/biostuff/ but now has disappeared.\n\nFor more information, see chapter 7 of Robert Turnbull's thesis 'The Textual History of Codex Sinaiticus Arabicus and its Family'.\n\nDocumentation and code clean up to come.",
    'author': 'Robert Turnbull',
    'author_email': 'robert.turnbull@unimelb.edu.au',
    'maintainer': 'None',
    'maintainer_email': 'None',
    'url': 'None',
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'python_requires': '>=3.8,<4.0',
    'include_dirs': [np.get_include()]
}
from build import *
build(setup_kwargs)

setup(**setup_kwargs)
