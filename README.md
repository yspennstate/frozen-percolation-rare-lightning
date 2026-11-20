# Frozen Percolation with Rare Lightning up to the Critical Time

Yitzchak Shmalo

The paper and the code for its numerical section.

In frozen percolation on a finite graph, edges open at independent uniform times, vertices are struck by lightning at rate λ, and a strike freezes the open cluster of the struck vertex. On graphs of maximum degree at most three, the probability that a given vertex is frozen by the critical time 1/2 is at most (3/4)λ log(1/λ) + O(λ), uniformly in the graph. On the three-regular tree an exact recursion for the freezing time of the root gives the expansion of this probability, and the constant 3/4 cannot be lowered.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.

- `code/residual_check.py`: solves the differential equation of Section 3 at 40 digits (mpmath) for λ from 1e-3 to 1e-12 and prints the quantities in the table of Section 6.
- `code/residual_check_output.json`: its output.

Requires Python 3 and mpmath. Run `python code/residual_check.py` (about a minute).
