# Computes convex hull (CH) of projection of a bounded polytope

Vertex enumeration using rational arithmetic to obtain the convex hull of the projection for a given polytope. We applied the CH method for computing the production envelopes for more than 2 dimensions in metabolic models. This implementation uses rational arithmetic (Fractions variable type in Pyton) to provide exact solutions to the optimization (Linear Programming). Therefore, for genome-scale metabolic models, projections onto less than 5 dimensions is recommended (due to increase in running time). 
 
Based on the algorithm for quantifier elimination in 
[Lassez, Catherine. "Quantifier elimination for conjunctions of linear constraints via a convex hull algorithm." Symbolic and Numerical Computation for Artificial Intelligence (1992): 103-122.]

Author: Sarah N. Galleguillos and Norbet Auer

Austrian Centre of Industrial Biotechnology, Vienna, Austria.

## Requirements
1. Python (tested on python 3.6)
2. Python-QSopt_ex 

~~~bash 
pip install python-qsoptex==0.5
libqsopt-ex-dev
~~~

3. Sympy Python library

## Dependencies

libqsopt-ex-dev/bionic,now 2.5.10.3-1build1

apt install libqsopt-ex-dev



