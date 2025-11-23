# Similarity and Distance Measures in proxyC

This vignette explains how **proxyC** compute the similarity and
distance measures.

## Notation

$$\begin{matrix}
{\overset{\rightarrow}{x} = \left\lbrack x_{i},x_{i + 1},\ldots,x_{n} \right\rbrack} \\
{\overset{\rightarrow}{y} = \left\lbrack y_{i},y_{i + 1},\ldots,y_{n} \right\rbrack}
\end{matrix}$$

The length of the vector $n = {||}\overset{\rightarrow}{x}{||}$, while
$\left| \overset{\rightarrow}{x} \right|$ is the absolute values of the
elements.

Operations on vectors are element-wise:

$$\begin{matrix}
{\overset{\rightarrow}{z} = \overset{\rightarrow}{x}\overset{\rightarrow}{y}} \\
{n = {||}\overset{\rightarrow}{x}{||} = {||}\overset{\rightarrow}{y}{||} = {||}\overset{\rightarrow}{z}{||}}
\end{matrix}$$

Summation of the elements of vectors is written using sigma without
specifying the range:

$$\sum\overset{\rightarrow}{x} = \sum\limits_{i = 1}^{n}x_{i}$$

When the elements of the vector is compared with a value in a pair of
square brackets, the summation is counting the number of elements that
equal (or unequal) to the value:

$$\sum\left\lbrack \overset{\rightarrow}{x} = 1 \right\rbrack = \sum\limits_{i = 1}^{n}\left\lbrack x_{i} = 1 \right\rbrack$$

## Similarity Measures

Similarity measures are available in
[`proxyC::simil()`](https://koheiw.github.io/proxyC/reference/simil.md).

#### Cosine similarity (“cosine”)

$$simil = \frac{\sum{\overset{\rightarrow}{x}\overset{\rightarrow}{y}}}{\sqrt{\sum{\overset{\rightarrow}{x}}^{2}}\sqrt{\sum{\overset{\rightarrow}{y}}^{2}}}$$

#### Pearson correlation coefficient (“correlation”)

$$simil = \frac{Cov\left( \overset{\rightarrow}{x},\overset{\rightarrow}{y} \right)}{Var\left( \overset{\rightarrow}{x} \right)Var\left( \overset{\rightarrow}{y} \right)}$$

#### Jaccard similarity (“jaccard” and “ejaccard”)

The values of $x$ and $y$ are Boolean for “jaccard”.

$$\begin{matrix}
{e = \sum{\overset{\rightarrow}{x}\overset{\rightarrow}{y}}} \\
{w = \text{user-provided weight}} \\
{simil = \frac{e}{\sum{\overset{\rightarrow}{x}}^{w} + \sum{\overset{\rightarrow}{y}}^{w} - e}}
\end{matrix}$$

#### Fuzzy Jaccard similarity (“fjaccard”)

The values must be $0 \leq x \leq 1.0$ and $0 \leq y \leq 1.0$.

$$simil = \frac{\sum{min\left( \overset{\rightarrow}{x},\overset{\rightarrow}{y} \right)}}{\sum{max\left( \overset{\rightarrow}{x},\overset{\rightarrow}{y} \right)}}$$

#### Dice similarity (“dice” and “edice”)

The values of $x$ and $y$ are Boolean for “dice”.

$$\begin{matrix}
{e = \sum{\overset{\rightarrow}{x}\overset{\rightarrow}{y}}} \\
{w = \text{user-provided weight}} \\
{simil = \frac{2e}{\sum{\overset{\rightarrow}{x}}^{w} + \sum{\overset{\rightarrow}{y}}^{w}}}
\end{matrix}$$

#### Hamann similarity (“hamann”)

$$\begin{matrix}
{e = \sum{\overset{\rightarrow}{x}\overset{\rightarrow}{y}}} \\
{n = {||}\overset{\rightarrow}{x}{||} = {||}\overset{\rightarrow}{y}{||}} \\
{u = n - e} \\
{simil = \frac{e - u}{e + u}}
\end{matrix}$$

#### Faith similarity (“faith”)

$$\begin{matrix}
{t = \sum{\left\lbrack \overset{\rightarrow}{x} = 1 \right\rbrack\left\lbrack \overset{\rightarrow}{y} = 1 \right\rbrack}} \\
{f = \sum{\left\lbrack \overset{\rightarrow}{x} = 0 \right\rbrack\left\lbrack \overset{\rightarrow}{y} = 0 \right\rbrack}} \\
{n = {||}\overset{\rightarrow}{x}{||} = {||}\overset{\rightarrow}{y}{||}} \\
{simil = \frac{t + 0.5f}{n}}
\end{matrix}$$

#### Simple matching (“matching”)

$$simil = \sum\left\lbrack \overset{\rightarrow}{x} = \overset{\rightarrow}{y} \right\rbrack$$

## Distance Measures

Similarity measures are available in
[`proxyC::dist()`](https://koheiw.github.io/proxyC/reference/simil.md).
Smoothing of the vectors can be performed when `method` is “chisquared”,
“kullback”, “jefferys” or “jensen”: the value of `smooth` will be added
to each element of $\overset{\rightarrow}{x}$ and
$\overset{\rightarrow}{y}$.

#### Manhattan distance (“manhattan”)

$$dist = \sum\left| \overset{\rightarrow}{x} - \overset{\rightarrow}{y} \right|$$

#### Canberra distance (“canberra”)

$$dist = \frac{\left| \overset{\rightarrow}{x} - \overset{\rightarrow}{y} \right|}{\left| \overset{\rightarrow}{x} \right| + \left| \overset{\rightarrow}{y} \right|}$$

#### Euclidian (“euclidian”)

$$dist = \sum\sqrt{{\overset{\rightarrow}{x}}^{2} + {\overset{\rightarrow}{y}}^{2}}$$

#### Minkowski distance (“minkowski”)

$$\begin{matrix}
{p = \text{user-provided parameter}} \\
{dist = \left( \sum\left| \overset{\rightarrow}{x} - \overset{\rightarrow}{y} \right|^{p} \right)^{\frac{1}{p}}}
\end{matrix}$$

#### Hamming distance (“hamming”)

$$dist = \sum\left\lbrack \overset{\rightarrow}{x} \neq \overset{\rightarrow}{y} \right\rbrack$$

#### The largest difference between values (“maximum”)

$$dist = \max{\overset{\rightarrow}{x} - \overset{\rightarrow}{y}}$$

#### Chi-squared divergence (“chisquared”)

$$\begin{matrix}
{O_{ij} = {\text{augmented matrix from}\mspace{6mu}}\overset{\rightarrow}{x}{\mspace{6mu}\text{and}\mspace{6mu}}\overset{\rightarrow}{y}} \\
{E_{ij} = {\text{matrix of expected count for}\mspace{6mu}}O_{ij}} \\
{dist = \sum\frac{\left( O_{ij} - E_{ij} \right)^{2}}{E_{ij}}}
\end{matrix}$$

#### Kullback–Leibler divergence (“kullback”)

$$\begin{matrix}
{\overset{\rightarrow}{p} = \frac{\overset{\rightarrow}{x}}{\sum\overset{\rightarrow}{x}}} \\
{\overset{\rightarrow}{q} = \frac{\overset{\rightarrow}{y}}{\sum\overset{\rightarrow}{y}}} \\
{dist = \sum{\overset{\rightarrow}{q}\log_{2}\frac{\overset{\rightarrow}{q}}{\overset{\rightarrow}{p}}}}
\end{matrix}$$

#### Jeffreys divergence (“jeffreys”)

$$\begin{matrix}
{\overset{\rightarrow}{p} = \frac{\overset{\rightarrow}{x}}{\sum\overset{\rightarrow}{x}}} \\
{\overset{\rightarrow}{q} = \frac{\overset{\rightarrow}{y}}{\sum\overset{\rightarrow}{y}}} \\
{dist = \sum{\overset{\rightarrow}{q}\log_{2}\frac{\overset{\rightarrow}{q}}{\overset{\rightarrow}{p}}} + \sum{\overset{\rightarrow}{p}\log_{2}\frac{\overset{\rightarrow}{p}}{\overset{\rightarrow}{q}}}}
\end{matrix}$$

#### Jensen-Shannon divergence (“jensen”)

$$\begin{matrix}
{\overset{\rightarrow}{p} = \frac{\overset{\rightarrow}{x}}{\sum\overset{\rightarrow}{x}}} \\
{\overset{\rightarrow}{q} = \frac{\overset{\rightarrow}{y}}{\sum\overset{\rightarrow}{y}}} \\
{\overset{\rightarrow}{m} = \frac{1}{2}\left( \overset{\rightarrow}{p} + \overset{\rightarrow}{q} \right)} \\
{dist = \frac{1}{2}\sum{\overset{\rightarrow}{q}\log_{2}\frac{\overset{\rightarrow}{q}}{\overset{\rightarrow}{m}}} + \frac{1}{2}\sum{\overset{\rightarrow}{p}\log_{2}\frac{\overset{\rightarrow}{p}}{\overset{\rightarrow}{m}}}}
\end{matrix}$$

## References

- Choi, S., Cha, S., & Tappert, C. C. (2010). A survey of binary
  similarity and distance measures. Journal of Systemics, *Cybernetics
  and Informatics*, 8(1), 43–48.
- Nielsen, F. (2019). On the Jensen–Shannon Symmetrization of Distances
  Relying on Abstract Means. *Entropy*, 21(5), 485.
  <https://doi.org/10.3390/e21050485>
- Jain, G., Mahara, T., & Tripathi, K. N. (2020). A Survey of Similarity
  Measures for Collaborative Filtering-Based Recommender System. In M.
  Pant, T. K. Sharma, O. P. Verma, R. Singla, & A. Sikander (Eds.),
  *Soft Computing: Theories and Applications* (pp. 343–352). Springer.
  <https://doi.org/10.1007/978-981-15-0751-9_32>
- Miyamoto, S. (1990). Hierarchical Cluster Analysis and Fuzzy Sets.
  In S. Miyamoto (Ed.), Fuzzy Sets in Information Retrieval and Cluster
  Analysis (pp. 125–188). Springer Netherlands.
  <https://doi.org/10.1007/978-94-015-7887-5_6>
