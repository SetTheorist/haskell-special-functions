\section{Spence}

Spence's integral for $z\geq0$ is
\[ S(z) = -\int_1^z \frac{\ln t}{t-1}\,dt = -\int_0^{z-1}\frac{ln(1+u)}{z}\,dz \]
and we extend the function via analytic continuation.
Spence's function $S(z)$ is related to the dilogarithm function via $S(z) = Li_2(1-z)$.

\subsection{Preamble}
\begin{code}
module Spence (
    sf_spence,
) where
import Exp
import Util
\end{code}

A useful constant $\verb|pi2_6| = \frac{\pi^2}{6}$
\begin{code}
pi2_6 :: (Value v) => v
pi2_6 = pi^2/6
\end{code}

\subsection{\tt sf\_spence z}
Compute Spence's integral $\verb|sf_spence z| = S(z)$.
We use a variety of transformations to to allow efficient computation with a series.
\begin{eqnarray*}
Li_2(z) + Li_2(\frac{z}{z-1}) &=& -\frac12(\ln(1-z))^2              \qquad{z\in\mathbb{C}\setminus[1,\infty)} \\
Li_2(z) + Li_2(\frac{1}{z}) &=& -\frac{\pi^2}{6}-\frac12(\ln(-z))^2 \qquad{z\in\mathbb{C}\setminus[0,\infty)} \\
Li_2(z) + Li_2(1-z) &=& \frac{\pi^2}{6}-\ln(z)\ln(1-z)              \qquad{0<z<1}
\end{eqnarray*}
(TODO: this code has not be solidly retested after conversion, especially verify complex.)
\begin{code}
sf_spence :: (Value v) => v -> v
sf_spence z
  | is_nan z     = z
  | (re z)<0     = 0/0
  | z == 0       = pi2_6
  | (rabs z)<0.5 = (series z) + (pi2_6 - (sf_log z)*(sf_log (1-z)))
  | (rabs z)<1.0 = -(series (1-z))
  | (rabs z)<2.5 = (series ((z-1)/z)) - (sf_log z)^2/2
  | otherwise    = (series (1/(1-z))) - pi2_6 - (sf_log (z-1))^2/2
\end{code}

\subsection*{*\tt series z}
The series expansion used for Spence's integral:
\[ \verb|series z| = -\sum_{k=1}^{\infty}\frac{z^k}{k^2} \]
\begin{code}
series z = 
  let zk = iterate (*z) z
      terms = zipWith (\ t k -> -t/(#)k^2) zk [1..]
  in ksum terms
\end{code}

