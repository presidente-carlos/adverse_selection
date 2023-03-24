Adverse Selection in Adaptive Settings (Gonzalez, 2023)
=======================================================

Abstract
---------


Additional Notes to Code
------------------------

_algorithm.R_

1. Parameters Generator Function
* The user may want to adapt parameters according to the specific setting they are working in. Default values of \gamma = 0.02 and B = 10 have been inputed

2. Algorithm 1
* The reader may consider modifiying `exp3_monop` to return additional information which is stored throughout the algorithmic process
* The reader may want to change the number of periods in which an instance of $p_i$ is stored (default = 50)
* The reader might want to set alternative parameters (default include parameters which ensure upper bounds described in Gonzalez, 2023). Particular values might achieve lower regret for particular instances (DGP) of the problem

3. Data preparation and plotting
* 4 workers are imputed for speeding the process through parallel computing. The reader might want to change this attribute according to the possibilities of its machine

4. Evaluation
* The following DGP are evaluated in the paper. All these dgp are generated in _dgp.R_: Uniform linear, Uniform non-linear, Uniform-Uniform, Bern(0.7) linear, Beta(2,1) linear, Beta(1/2, 1/2) linear, Beta(2,1) non-linear
* The reader may consider alternative DGP


_dgp.R_

1. The reader might want to implement additional DGP
* Bounds provided in Gonzalez (2023) hold for any adversarial distribution $F_{U,V}$, thus the reader is encouraged to feed any (bounded) distribution into `exp3_monop()`. However, for the empirical analysis in terms of regret or probability sampling, the user must also provide closed form solutions to \mathbb{E}_{F_{U,V}} S_i(x) and \max_x \mathbb{E}_{F_{U,V}} S_i(x).
