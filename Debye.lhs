\chapter{Debye functions}

\section{Preamble}
\begin{titled-frame}{\color{blue}\tt module Debye}
\begin{code}
{-# Language BangPatterns #-}
module Debye where
import Data.List(zipWith3)
import Exp
import Numbers
import Util
import Zeta
\end{code}
\end{titled-frame}

\section{Debye functions $D_n(z)$}

The Debye functions $D_n(z)$ for ($n=1,2,\dots$) are defined via
\[ D_n(z) = \int_0^x \frac{t^n}{e^t-1}\,dt \]

\subsection{\tt sf\_debye n z}
We implement this by series for small $|z|\leq2$,

For $\Re z<0$ we use the reflection formula
\[ D_n(-z) = (-)^n D_n(z) + (-)^n \frac{z^{n+1}}{n+1} \]

(TODO: Could try direct quadrature of the defining integral.)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_debye n z} = D_n(z)$}
\begin{code}
sf_debye :: (Value v) => Int -> v -> v
sf_debye n z
  | (re z)<0       = (-1)^n*((sf_debye n (-z)) + (-z)^(n+1)/(#)(n+1))
  | (rabs z)<=2    = debye__ser n z
  | (rabs z)<=(#)n = debye__co2 n z
  | otherwise      = ((#)$factorial n)*(sf_zeta((#)$n+1)) - (debye__coint n z)
\end{code}
\end{titled-frame}

\subsubsection{\tt debye\_\_ser n z}
A series representation for the Debye functions for $|z|<2\pi$ and $n\geq1$
is given by
\[  D_{n}(z)
  = z^n \sum_{k=0}^\infty \frac{B_{k}}{k!} \frac{z^{k}}{k+n}
  = z^n \left( \frac{1}{n} - \frac{z}{2(n+1)} + \sum_{k=1}^\infty \frac{B_{2k}}{(2k)!} \frac{z^{2k}}{2k+n} \right)
  \]
\begin{titled-frame}{$\text{\color{blue}\tt debye\_\_ser n z}$}
\begin{code}
debye__ser !n !z = 
  let !z2    = z^2
      !z2ns  = iterate (*z2) 1
      tk !k  = (fromRational$sf_bernoulli_b_scaled!!k)/(#)(k+n)
      !terms = zipWith (*) (map tk [0,2..]) z2ns
      !smm   = ksum terms
  in (z^n)*(-z/(#)(2*(n+1)) + smm)
\end{code}
\end{titled-frame}

\subsubsection{\tt debye\_\_coint n z}
We can use the series expansion for the complementary integral
(recalling that $\int_0^\infty t^n(e^t-1)\,dt = n!\zeta(n+1)$.)
\[ \int_z^\infty \frac{t^n}{e^t-1}\,dt
    = \sum_{k=1}^\infty e^{-k z}
        \left( \frac{z^n}{k} + \frac{n x^{n-1}}{k^2} + \frac{n(n-1)x^{n-2}}{k^3} + \cdots + \frac{n!}{k^{n+1}} \right)
    \]
\begin{titled-frame}{$\text{\color{blue}\tt debye\_\_coint n z}$}
\begin{code}
debye__coint !n !z =
  let !terms = map (\k -> (sf_exp(-z*(#)k))*(trm k)) [1..]
  in ksum terms
  where
    trm !k = 
      let !terms = take (n+1) $ ixiter 1 (z^n/(#)k) $ \j t -> t*((#)$n+1-j)/(z*(#)k)
      in ksum terms
\end{code}
\end{titled-frame}

\subsubsection{\tt debye\_\_co2 n z}
This is an alternative formulation of the expression
in terms of the complementary integral
which may offer improved numerical stability
(with bracketed terms computed individually with high precision.)
(TODO: write out in terms of functions used...)
\[ D_n(z)
    = n! \zeta(n+1) - \int_z^\infty \frac{t^n}{e^t-1}\,dt
    = n!\left( \left[\zeta(n+1)-1\right]
        + e^{-z}\left[\sum_{j=n+1}^\infty\frac{z^j}{j!}\right]
        - \sum_{k=2}^\infty \frac{e^{-kz}}{k^{n+1}} \left[ \frac{z^j}{j!} \right]
        \right) \]
\begin{titled-frame}{$\text{\color{blue}\tt debye\_\_co2 n z}$}
\begin{code}
debye__co2 !n !z = 
  let !zet = sf_zeta_m1 $ (#)(n+1)
      !eee = (sf_exp(-z)) * (sf_exp_men (n+1) z)
      !terms = map (\k-> -(sf_exp(-z*(#)k))/(((#)k)^(n+1))*(sf_expn n (z*(#)k))) [2..]
      !smm = ksum (zet:eee:terms)
  in smm * ((#)$factorial n)
\end{code}
\end{titled-frame}

\section{Scaled Debye functions}

The scaled Debye functions are
\[ \widetilde{D}_n(z) = \frac{n}{x^n} D_n(z) \]

\subsection{\tt sf\_debye\_scaled n z}

\begin{titled-frame}{$\text{\color{blue}\tt sf\_debye\_scaled n z} = \widetilde{D}_n(z)$}
\begin{code}
sf_debye_scaled :: (Value v) => Int -> v -> v
sf_debye_scaled !n !z
  | (re z) < 0   = -z*(#)n/(#)(n+1) + (sf_debye_scaled n (-z))
  | (rabs z) < 2 = debye_sc__ser n z
  | otherwise    = (#)n*(sf_zeta((#)n+1))*((#)$factorial n)/z^n - (debye_sc__coint n z)
\end{code}
\end{titled-frame}

\subsubsection{\tt debye\_sc\_\_ser n z}
\begin{titled-frame}{$\text{\color{blue}\tt debye\_sc\_\_ser n z}$}
\begin{code}
debye_sc__ser :: (Value v) => Int -> v -> v
debye_sc__ser !n !z = 
  let !n' = (#)n
      !z2 = z^2
      !zps = iterate (*z2) z2
      !bns = map fromRational $ tail . evel $ sf_bernoulli_b_scaled
      !tms = map ((#).(+ n)) [2,4..]
      !tterms = zipWith3 (\ zp bn tm -> zp*bn/tm) zps bns tms
      !terms = (1/n'):(-z/(2*(n'+1))):tterms
  in n' * (ksum terms)
  where evel !(a:b:ts) = a:(evel ts)
\end{code}
\end{titled-frame}

\subsubsection{\tt debye\_sc\_\_coint n z}
\begin{titled-frame}{$\text{\color{blue}\tt debye\_sc\_\_coint n z}$}
\begin{code}
debye_sc__coint :: (Value v) => Int -> v -> v
debye_sc__coint !n !z =
  let !tms = map trm [1..]
      !ees = map (\k->sf_exp $ -z*(#)k) [1..]
      !terms = zipWith (*) tms ees
  in ksum terms
  where
    trm !k =
      let !terms = take (n+1) $ ixiter 1 ((#)n/(#)k) $ \j t -> t*((#)$n+1-j)/(z*(#)k)
      in ksum terms
\end{code}
\end{titled-frame}

