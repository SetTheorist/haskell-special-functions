\section{Error function}

\subsection{Preamble}

\begin{code}
{-# Language BangPatterns #-}
module Erf (
    sf_erf,
    sf_erfc,
    erfc_asymp_pos,
    erfc_asymp_pos',
) where
import Exp
import Util
\end{code}

\subsection{Error function}

\subsubsection{\tt sf\_erf z}
The error function $\verb|sf_erf z| = \erf z$ where
\[ \erf z = \frac{2}{\sqrt\pi} \int_{-\infty}^{z} e^{-x^2}\,dx \]
For $\Re z<-1$, we transform via $\erf z=-\erf(-z)$ and for $|z|<1$ we use
the power-series expansion, otherwise we use $\erf z=1-\erfc z$.
(TODO: this implementation is not perfect, but workable for now.)
\begin{code}
sf_erf :: (Value v) => v -> v
sf_erf z 
  | (re z)<(-1) = -sf_erf(-z)
  | (rabs z)<1  = erf_series z
  | otherwise   = 1 - sf_erfc z
\end{code}

\subsubsection{\tt sf\_erfc z}
The complementary error-function $\verb|sf_erfc z| = \erfc z$ where
\[ \erfc z = 1-\erf z = \frac{2}{\sqrt\pi} \int_{z}^{\infty} e^{-x^2}\,dx \]
For $\Re z<-1$ we transform via $\erfc z=2-\erf(-z)$ and if $|z|<1$ then we
use $\erfc z=1-\erf z$.  Finally, if $|z|<10$ we use a continued-fraction expansion
and an asymptotic expansion otherwise.
(TODO: there are a few issues with this implementation:
For pure imaginary values and for extremely large values it seems to hang.)
\begin{code}
-- infinite loop when (re z)==0
sf_erfc :: (Value v) => v -> v
sf_erfc z 
  | (re z)<(-1) = 2-(sf_erfc (-z))
  | (rabs z)<1  = 1-(sf_erf z)
  | (rabs z)<10 = erfc_cf_pos1 z
  | otherwise   = erfc_asymp_pos z -- TODO: hangs for very large input
\end{code}

\subsubsection*{\tt erf\_series z}
The series expansion for $\erf z$:
\[ \erf z = \frac{2}{\sqrt\pi}\sum_{n=0}^\infty\frac{(-)^nz^{2n+1}}{n!(2n+1)} \]
There is an alternative expansion
$\erf z=\frac{2}{\sqrt\pi}e^{-z^2}\sum_{n=0}^\infty\frac{2^nz^{2n+1}}{1\cdot3\cdots(2n+1)}$, but
we don't use it here. (TODO: why not?)
\begin{code}
erf_series z =
  let z2 = z^2
      rts = ixiter 1 z $ \n t -> (-t)*z2/(#)n
      terms = zipWith (\ n t ->t/(#)(2*n+1)) [0..] rts
  in (2/sf_sqrt pi)  * ksum terms
\end{code}

\subsubsection*{*\tt sf\_erf z}
This asymptotic expansion for $\erfc z$ is valid as $z\to+\infty$:
\[ \erfc z \sim \frac{e^{-z^2}}{\sqrt\pi}\sum_{n=0}^\infty(-)^n\frac{(1/2)_m}{z^{2m+1}} \]
where the Pochhammer symbol $(1/2)_m$ is given by:
\[ \left(\frac12\right)_m = \frac{1\cdot3\cdot5\cdots(2m-1)}{2^m} = \frac{(2m-1)!}{2^m} \]
TODO: correct the asymptotic term checking (not smallest but pre-smallest term).
\begin{code}
erfc_asymp_pos z =
  let z2 = z^2
      iz2 = 1/2/z2
      terms = ixiter 1 (1/z) $ \n t -> (-t*iz2)*(#)(2*n-1)
      tterms = tk terms
  in (sf_exp (-z2))/(sqrt pi) * ksum tterms
  where tk (a:b:cs) = if (rabs a)<(rabs b) then [a] else a:(tk$b:cs)
\end{code}

\subsubsection*{*\tt erfc\_cf\_pos1 z}
A continued-fraction expansion for $\erfc z$:
\[ \sqrt\pi e^{z^2} \erfc z = \frac{z}{z^2+{}} \frac{1/2}{1+{}} \frac{1}{z^2+{}} \frac{3/2}{1+{}} \cdots\]
\begin{code}
erfc_cf_pos1 z = 
  let z2 = z^2
      as = z:(map fromDouble [1/2,1..])
      bs = 0:cycle [z2,1]
      cf = steeds as bs
  in sf_exp(-z2) / (sqrt pi) * cf
\end{code}

\subsubsection*{*\tt erfc\_cf\_pos1 z}
This is an alternative continued-fraction expansion.
\[ \sqrt\pi e^{z^2} \erfc z = \frac{2z}{2z^2+1-{}} \frac{1\cdot2}{2z^2+5-{}} \frac{3\cdot4}{2z^2+9-{}} \cdots \]
Unused for now.
\begin{code}
erfc_cf_pos2 z = 
  let z2 = z^2
      as = (2*z):(map (\n->(#)$ -(2*n+1)*(2*n+2)) [0..])
      bs = 0:(map (\n->2*z2+(#)4*n+1) [0..])
      cf = steeds as bs
  in sf_exp(-z2) / (sqrt pi) * cf
\end{code}

