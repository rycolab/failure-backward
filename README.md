# Algorithms for Acyclic Weighted Finite-State Automata with Failure Arcs
Code accompanying the EMNLP 2022 publication "Algorithms for Acyclic Weighted Finite-State Automata with Failure Arcs".

The code implements Algorithms 1, 2, 3, 4, 5, and 6 from the paper.
Algorithm 2 is implemented as `expand_phi_arcs` in the `Transformer` class, while the other algorithms are implemented in the `Pathsum` class.

### Running the Code
The code only requires `pytest` to run. To unit-test the code, simply run 
```
pytest test_failure_arcs_backward.py
```
from the `src` directory.

## Cite
```
@inproceedings{svete-etal-2022-algorithms,
    title = "Algorithms for Acyclic Weighted Finite-State Automata with Failure Arcs",
    author = "Svete, Anej  and
      Dayan, Benjamin  and
      Cotterell, Ryan  and
      Vieira, Tim  and
      Eisner, Jason",
    booktitle = "Proceedings of the 2022 Conference on Empirical Methods in Natural Language Processing",
    month = dec,
    year = "2022",
    address = "Abu Dhabi, United Arab Emirates",
    publisher = "Association for Computational Linguistics",
    url = "https://aclanthology.org/2022.emnlp-main.567",
    pages = "8289--8305",
}
```
