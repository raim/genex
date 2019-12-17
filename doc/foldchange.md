# Estimating Protein Expression Rates

For plasmid-borne expression of a (fluorescent reporter) protein
in a growing culture the intracellular concentration of
the protein (proteins/cell) can be summarized as:

$$
\frac{\text{d}P}{\text{d}t} = n k - (\mu+\delta)P\;,
$$

with protein degradation rate $\delta$, plasmid copy number $n$,
culture growth rate $\mu$, and overall expression rate $k$.
To quantify $P(t)$ from platereader measurements we require to
know the conversion factors fluorescence/protein $f_P$ and and OD/cell
$f_O$. From the measured fluorescence per OD values $F(t)=\frac{f(t)}{OD(t)}$,
we can caculate

$$
P(t) = F(t) \frac{f_0}{f_P}\;.
$$

Without knowing these factors, we can only comparatively analyze
expression between two different experiments, where they cancel out.
If we find a phase of steady state expression, where fluorescence/OD
is stable over a (short) time, the calculations are easy. At
$\frac{\text{d}P_1}{\text{d}t}=\frac{\text{d}P_2}{\text{d}t}=0$ we get:

$$
\frac{n_1 k_1}{n_2 k_2} = \frac{(\mu_1 + \delta) F_1}{(\mu_2 + \delta) F_2}\;.
$$

If we know the protein degradation rate $\delta$ and further assume
constant plasmid copy number between the two conditions, we only need
to establish (a) growth rates $\mu$ and fluorescence values $F(t)$ at
a given experiment time $t$ which sufficiently fulfills our
assumptions (steady state!).

$F(t)$ is directly measured, and to estimate the growth rate
we can use the OD measurements, eg. with the `easylinear`
fitting procedure of R package `growthrates` or with
`dpseg` for more complex growth curves.
Notably, to avoid measurement noise at low cell densities, we
can use the growth model instead of the real OD data for
the normalization of the measured fluorescence by OD. Additionally,
it may be advisable to smooth the fluorescence measurement over
time-points.


#### Dynamic State

If we are NOT in a steady state, we can still get an idea
of expression strength. For that, we need to get both
$P$ and $\frac{\text{d}P}{\text{d}t}$ (the slope), and can
compare

$$
n k = \frac{\text{d}P}{\text{d}t} + (\mu+\delta)P\;.
$$




