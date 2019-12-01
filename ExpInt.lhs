\chapter{Exponential Integral}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Preamble}
\begin{titled-frame}{\color{blue}\tt module ExpInt}
\begin{code}
{-# Language BangPatterns #-}
module ExpInt(sf_expint_ei, sf_expint_en) where
import Exp
import Gamma
import Util
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Exponential integral $\Ei$}
The exponential integral $\Ei z$ is defined for $x<0$ by
\[ \Ei(z) = -\int_{-x}^\infty \frac{e^{-t}}{t}\,dt \]
It can be defined 

\subsection{\tt sf\_expint\_ei z}
We give only an implementation for $\Re z\geq0$.
We use a series expansion for $|z|<40$ and an asymptotic expansion otherwise.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_expint\_ei z} = \Ei(z)$\marginnote{\tt sf\_expint\_ei}}
\begin{code}
sf_expint_ei :: (Value v) => v -> v
sf_expint_ei !z
  | (re z) < 0.0  = nan
  | z == 0.0      = neg_infty
  | (rabs z) < 40 = expint_ei__series z
  | otherwise     = expint_ei__asymp z
\end{code}
\end{titled-frame}

\subsubsection{\tt expint\_ei\_\_series}
The series expansion is given (for $x>0$)
\[ \Ei(x) = \gamma + \ln x + \sum_{n=1}^\infty \frac{x^n}{n! n} \marginnote{$\Ei(x)$} \]
We evaluate the addition of the two terms with the sum slightly differently
when $\Re z<1/2$ to reduce floating-point cancellation error slightly.
\begin{titled-frame}{\color{blue}\tt expint\_ei\_\_series z}
\begin{code}
expint_ei__series :: (Value v) => v -> v
expint_ei__series !z =
  let !tterms = ixiter 2 z $ \n t -> t*z/(#)n
      !terms = zipWith (\ t n ->t/(#)n) tterms [1..]
      !res = ksum terms
  in if (re z)<0.5
     then sf_log(z * sf_exp(euler_gamma + res))
     else res + sf_log(z) + euler_gamma
\end{code}
\end{titled-frame}

\subsubsection{\tt expint\_ei\_\_asymp}
The asymptotic expansion as $x\to+\infty$ is
\[ \Ei(x) \sim \frac{e^x}{x}\sum_{n=0}^\infty \frac{n!}{x^n} \]
\begin{titled-frame}{\color{blue}\tt expint\_ei\_\_asymp z}
\begin{code}
expint_ei__asymp :: (Value v) => v -> v
expint_ei__asymp !z =
  let !terms = tk $ ixiter 1 1.0 $ \n t -> t/z*(#)n
      !res = ksum terms
  in res * (sf_exp z) / z
  where tk (a:b:cs) = if (rabs a)<(rabs b) then [a] else a:(tk$b:cs)
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Exponential integral $E_n$}
The exponential integrals $E_n(z)$ are defined for $n=0,1,\dots$ and $\Re z>0$ via
\[ E_n(z) = z^{n-1}\int_{z}^\infty \frac{e^{-t}}{t^n}\,dt \marginnote{$E_n(z)$} \]
They satisfy the following relations:
\begin{eqnarray*}
E_0(z) &=& \frac{e^{-z}}{z} \\
E_{n+1}(z) &=& \int_z^\infty E_{n}(t)\,dt \\
\end{eqnarray*}
And they can be expressed in terms of incomplete gamma functions:
\[ E_n(z) = z^{n-1}\Gamma(1-n,z) \]
(which also gives a generalization for non-integer $n$).

\subsection{\tt sf\_expint\_en n z}
We evaluate the exponential integrals $E_n(z)$ by handling the special cases
$n=0,1$ directly, otherwise use a series expansion for $|z|\leq1$
and a continued fraction expansion otherwise.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_expint\_en n z} = E_n(z)$}
\begin{code}
sf_expint_en :: (Value v) => Int -> v -> v
sf_expint_en !n !z | (re z)<0 = nan -- TODO: confirm this
                   | z == 0   = (1/(#)(n-1)) -- TODO: confirm this
sf_expint_en !0 !z = sf_exp(-z) / z
sf_expint_en !1 !z = expint_en__1 z
sf_expint_en !n !z | (rabs z) <= 1.0 = expint_en__series n z
                   | otherwise       = expint_en__contfrac n z
\end{code}
\end{titled-frame}

\subsubsection{\tt expint\_en\_\_1}
We use this series expansion for $E_1(z)$:
\[ E_1(z) = -\gamma - \ln z + \sum_{k=1}^\infty(-)^k\frac{z^k}{k! k}\]
(Note that this will not be good for large positive values of $z$ due to cancellation.)
\begin{titled-frame}{\color{blue}\tt expint\_en\_\_1 z}
\begin{code}
expint_en__1 :: (Value v) => v -> v
expint_en__1 z =
  let !r0 = -euler_gamma - (sf_log z)
      !tterms = ixiter 2 z $ \k t -> -t*z/(#)k
      !terms = zipWith (\t k -> t/(#)k) tterms [1..]
  in ksum (r0:terms)
\end{code}
\end{titled-frame}

\subsubsection{\tt expint\_en\_\_series}
The series expansion for the exponential integral
\[ E_n(z) = \frac{(-z)^{n-1}}{(n-1)!}(-\ln(z) + \psi(n)) - \sum_{m=0,m\neq n}^\infty \frac{(-x)^m}{(m-(n-1))m!} \]
for $n\geq2$, $z\leq1$
\begin{titled-frame}{\color{blue}\tt expint\_en\_\_series n z}
\begin{code}
-- assume n>=2, z<=1
expint_en__series :: (Value v) => Int -> v -> v
expint_en__series n z =
  let !n' = (#)n
      !res = (-(sf_log z) + (sf_digamma n')) * (-z)^(n-1)/(#)(factorial$n-1) + 1/(n'-1)
      !terms' = ixiter 2 (-z) (\m t -> -t*z/(#)m)
      !terms = map (\(m,t)->(-t)/(#)(m-(n-1))) $ filter ((/=(n-1)) . fst) $ zip [1..] terms'
  in ksum (res:terms)
\end{code}
\end{titled-frame}

\subsubsection{\tt expint\_en\_\_contfrac}
The continued fraction expansion for the exponential integral, valid for $z>1$, $n\geq2$.
(TODO: verify for which complex values is this valid?)
\[ e^{-x} \left( \frac{1}{x+n-{}}\ \frac{1\cdot n}{x+(n+2)-{}}\ \frac{2\cdot(n+1)}{x+(n+4)- \cdots} \right) \]
\begin{titled-frame}{\color{blue}\tt expint\_en\_\_contfrac n z}
\begin{code}
expint_en__contfrac :: (Value v) => Int -> v -> v
expint_en__contfrac !n !z =
  let !n' = (#)n
      !an = 1:[-(1+k) * (n'+k) | k'<-[0..], let k=(#)k']
      !bn = 0:[z + n' + 2*k    | k'<-[0..], let k=(#)k']
  in (sf_exp(-z))*(sf_cf_lentz an bn)
\end{code}
\end{titled-frame}
