\chapter{Spence}

Spence's integral for $z\geq0$ is
\[ \O{S}(z) = -\int_1^z \frac{\ln t}{t-1}\,dt = -\int_0^{z-1}\frac{ln(1+u)}{z}\,dz \]
and we extend the function via analytic continuation.
Spence's function $\O{S}(z)$ is related to the dilogarithm function via $S(z) = \O{Li}_2(1-z)$.

\section{Preamble}
\begin{titled-frame}{\color{blue}\tt module Spence}
\begin{code}
module Spence (sf_spence) where
import Exp
import Util
\end{code}
\end{titled-frame}

A useful constant $\verb|pi2_6| = \frac{\pi^2}{6}$
\begin{code}
pi2_6 :: (Value v) => v
pi2_6 = pi^2/6
\end{code}

\section{\tt sf\_spence z}
Compute Spence's integral $\verb|sf_spence z| = \O{S}(z)$.
We use a variety of transformations to to allow efficient computation with a series.
\begin{eqnarray*}
\O{Li}_2(z) + \O{Li}_2(\frac{z}{z-1}) &=& -\frac12(\ln(1-z))^2              \qquad{z\in\CC\setminus[1,\infty)} \\
\O{Li}_2(z) + \O{Li}_2(\frac{1}{z}) &=& -\frac{\pi^2}{6}-\frac12(\ln(-z))^2 \qquad{z\in\CC\setminus[0,\infty)} \\
\O{Li}_2(z) + \O{Li}_2(1-z) &=& \frac{\pi^2}{6}-\ln(z)\ln(1-z)              \qquad{0<z<1}
\end{eqnarray*}
(TODO: this code has not be solidly retested after conversion, especially verify complex.)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_spence z} = \O{Li}_2(z)$}
\begin{code}
sf_spence :: (Value v) => v -> v
sf_spence z
  | is_nan z     = z
  | (re z)<0     = 0/0
  | z == 0       = pi2_6
  | (rabs z)<0.5 = (spence__ser z) + (pi2_6 - (sf_log z)*(sf_log (1-z)))
  | (rabs z)<1.0 = -(spence__ser (1-z))
  | (rabs z)<2.5 = (spence__ser ((z-1)/z)) - (sf_log z)^2/2
  | otherwise    = (spence__ser (1/(1-z))) - pi2_6 - (sf_log (z-1))^2/2
\end{code}
\end{titled-frame}

\subsubsection{\tt spence\_\_ser z}
The series expansion used for Spence's integral:
\[ \verb|spence__ser z| = -\sum_{k=1}^{\infty}\frac{z^k}{k^2} \]
\begin{code}
spence__ser z = 
  let zk = iterate (*z) z
      terms = zipWith (\ t k -> -t/(#)k^2) zk [1..]
  in ksum terms
\end{code}

