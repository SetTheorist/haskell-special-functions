\chapter{Error function}

\section{Preamble}

\begin{titled-frame}{\color{blue}\tt module Erf}
\begin{code}
{-# Language BangPatterns #-}
-- {-# Language BlockArguments #-}
{-# Language ScopedTypeVariables #-}
module Erf (sf_erf, sf_erfc) where
import Exp
import Util
\end{code}
\end{titled-frame}

\section{Error function}

The error function is defined via
\[ \erf(z) = \frac{2}{\sqrt\pi} \int_{0}^{z} e^{-x^2}\,dx \marginnote{$\erf(z)$}\]
and the complementary error function via 
\[ \erfc(z) = \frac{2}{\sqrt\pi} \int_{z}^{\infty} e^{-x^2}\,dx \marginnote{$\erfc(z)$}\]
Thus we have the relation $\erf(z) + \erfc(z) = 1$.

\subsection{\tt sf\_erf z}
The error function $\verb|sf_erf z| = \erf z$ where
\[ \erf(z) = \frac{2}{\sqrt\pi} \int_{-\infty}^{z} e^{-x^2}\,dx \]
For $\Re z<-1$, we transform via $\erf(z)=-\erf(-z)$ and for $|z|<1$ we use
the power-series expansion, otherwise we use $\erf z=1-\erfc z$.
(TODO: this implementation is not perfect, but workable for now.)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_erf z} = \erf(z)$}
\begin{code}
sf_erf :: (Value v) => v -> v
sf_erf z 
  | (re z)<(-1) = -sf_erf(-z)
  | (rabs z)<1  = erf__series z
  | otherwise   = 1 - sf_erfc z
\end{code}
\end{titled-frame}

\subsection{\tt sf\_erfc z}
The complementary error-function $\verb|sf_erfc z| = \erfc z$ where
\[ \erfc z = 1-\erf z = \frac{2}{\sqrt\pi} \int_{z}^{\infty} e^{-x^2}\,dx \]
For $\Re z<-1$ we transform via $\erfc z=2-\erf(-z)$ and if $|z|<1$ then we
use $\erfc z=1-\erf z$.  Finally, if $|z|<10$ we use a continued-fraction expansion
and an asymptotic expansion otherwise.
(TODO: there are a few issues with this implementation:
For pure imaginary values and for extremely large values it seems to hang.)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_erfc z} = \erfc(z)$}
\begin{code}
sf_erfc :: (Value v) => v -> v
sf_erfc z 
  | (re z)<(-1) = 2-(sf_erfc (-z))
  | (rabs z)<1  = 1-(sf_erf z)
  | (rabs z)<10 = erfc_cf_pos1 z
  | otherwise   = erfc_asymp_pos z -- TODO: hangs for very large input
\end{code}
\end{titled-frame}

\subsubsection{\tt erf\_series z}
The series expansion for $\erf z$:
\[ \erf z = \frac{2}{\sqrt\pi}\sum_{n=0}^\infty\frac{(-)^nz^{2n+1}}{n!(2n+1)} \]
There is an alternative expansion
$\erf z=\frac{2}{\sqrt\pi}e^{-z^2}\sum_{n=0}^\infty\frac{2^nz^{2n+1}}{1\cdot3\cdots(2n+1)}$, but
we don't use it here. (TODO: why not?)
\begin{code}
erf__series z =
  let z2 = z^2
      rts = ixiter 1 z $ \n t -> (-t)*z2/(#)n
      terms = zipWith (\n t -> t/(#)(2*n+1)) [0..] rts
  in (2/sf_sqrt pi) * (ksum terms)
\end{code}

\subsubsection{\tt sf\_erf z}
This asymptotic expansion for $\erfc z$ is valid as $z\to+\infty$:
\[ \erfc z \sim \frac{e^{-z^2}}{\sqrt\pi}\sum_{n=0}^\infty(-)^n\frac{(1/2)_m}{z^{2m+1}} \]
where the Pochhammer symbol $(1/2)_m$ is given by:
\[ \left(\frac12\right)_m = \frac{1\cdot3\cdot5\cdots(2m-1)}{2^m} = \frac{(2m)!}{m!2^{2m}} \]
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

\subsubsection{\tt erfc\_cf\_pos1 z}
A continued-fraction expansion for $\erfc z$:
\[ \sqrt\pi e^{z^2} \erfc z = \frac{z}{z^2+{}} \frac{1/2}{1+{}} \frac{1}{z^2+{}} \frac{3/2}{1+{}} \cdots\]
\begin{code}
erfc_cf_pos1 z = 
  let z2 = z^2
      as = z:(map fromDouble [1/2,1..])
      bs = 0:cycle [z2,1]
      cf = sf_cf_steeds as bs
  in sf_exp(-z2) / (sqrt pi) * cf
\end{code}

\subsubsection{\tt erfc\_cf\_pos1 z}
This is an alternative continued-fraction expansion.
\[ \sqrt\pi e^{z^2} \erfc z = \frac{2z}{2z^2+1-{}} \frac{1\cdot2}{2z^2+5-{}} \frac{3\cdot4}{2z^2+9-{}} \cdots \]
Unused for now.
\begin{code}
erfc_cf_pos2 z = 
  let z2 = z^2
      as = (2*z):(map (\n->(#)$ -(2*n+1)*(2*n+2)) [0..])
      bs = 0:(map (\n->2*z2+(#)4*n+1) [0..])
      cf = sf_cf_steeds as bs
  in sf_exp(-z2) / (sqrt pi) * cf
\end{code}

\subsection{Dawson's function}

Dawson's function (or Dawson's integral) is given by
\[ D(z) = e^{-z^2}\int_0^z e^{t^2}\,dt = -\frac{\ii\sqrt\pi}{2} e^{-x^2}\erf(\ii x) \]

\subsubsection{\tt sf\_dawson z}
Compute Dawson's integral $D(z) = e^(-z^2) \int_0^z e^(t^2) dt$ for real z.
(Correct only for reals!)

\begin{code}
sf_dawson :: forall v.(Value v) => v -> v
sf_dawson z
  -- | (rabs z) < 0.5 = (toComplex$sf_exp(-z^2))*(sf_erf((toComplex z)*(0:+1)))*(sf_sqrt(pi)/2/(0:+1))
  | (im z) /= 0    = dawson__seres z
  | (rabs z) < 5   = dawson__contfrac z
  | otherwise      = dawson__contfrac2 z

dawson__seres :: (Value v) => v -> v
dawson__seres z =
  let tterms = ixiter 1 z $ \n t -> t*z^2/(#)n
      terms = zipWith (\n t->t/((#)(2*n+1))) [0..] tterms
      smm = ksum terms
  in (sf_exp(-z^2)) * smm

faddeeva__asymp :: (Value v) => v -> v
faddeeva__asymp z =
  let z' = 1/z
      terms = ixiter 1 z' $ \n t -> t*z'^2*((#)(2*n+1))/2
      smm = ksum terms
  in smm

dawson__contfrac :: (Value v) => v -> v
dawson__contfrac z = undefined

dawson__contfrac2 :: (Value v) => v -> v
dawson__contfrac2 z = undefined

\end{code}
