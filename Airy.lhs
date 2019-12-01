\section{Airy}

The Airy functions $\Ai$ and $\Bi$,
are standard solutions of the ode $y''-zy=0$.

\subsection{Preamble}
\begin{titled-frame}{\color{blue}\tt module Airy}
\begin{code}
{-# Language BangPatterns #-}
module Airy (sf_airy_ai, sf_airy_bi) where
import Exp
import Gamma
import Trig
import Util
\end{code}
\end{titled-frame}

\subsection{Ai}

The solution $\Ai(z)$ of the Airy ODE is given by
\[ \Ai(z) = \frac1\pi \int_0^\infty \cos(\frac{t^3}{3} + xt)\,dt \]
it can be given in terms of Bessel functions, where $\zeta=(2/3)z^{3/2}$
\[ \Ai(z) = \frac{\sqrt{z/3}}{\pi} K_{\pm1/3}(\zeta)
    = \frac{\sqrt{z}}{3} \left( I_{-1/3}(\zeta) - I_{1/3}(\zeta) \right) \]
or
\[ \Ai(-z) = \frac{\sqrt{z}}{3} \left( J_{1/3}(\zeta) - J_{-1/3}(\zeta) \right) \]

\subsubsection{\tt sf\_airy\_ai z}
For now, we use an asymptotic expansion for large values and a series for smaller values.
This gives reasonable results for small-enough or large-enough values, but it has
low accuracy for intermediate values, ({\it e.g.\@} $z=5$).
(TODO: Seems quite bad for complex values)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_airy\_ai z} = \Ai(z)$}
\begin{code}
sf_airy_ai :: (Value v) => v -> v
sf_airy_ai !z
  | (rabs z)>=9 && (re z)>=0 = airy_ai__asympt_pos z
  | (rabs z)>=9 && (re z)< 0 = airy_ai__asympt_neg z
  | otherwise                = airy_ai__series z
\end{code}
\end{titled-frame}

Initial conditions
\begin{eqnarray*}
\Ai(0)  &=& \frac{1}{3^{2/3} \Gamma(2/3)} \\
\Ai'(0) &=& \frac{-1}{3^{1/3} \Gamma(1/3)}
\end{eqnarray*}
\begin{code}
ai0 :: (Value v) => v
ai0  = ( 3**(-2/3)) * (sf_invgamma(2/3))
ai'0 :: (Value v) => v
ai'0 = (-3**(-1/3)) * (sf_invgamma(1/3))
\end{code}

The series expansion for $\Ai(z)$ is given by 
\begin{eqnarray*}
\Ai(z) &=& \Ai(0)f(z) + \Ai'(0)g(z) \\
f(z) &=& \sum_{n=0}^\infty\frac{(3n-2)!!!}{(3n)!}z^{3n}
  = 1 + \frac{1}{3!}z^3 + \frac{1\cdot4}{6!}z^6 + \frac{1\cdot4\cdot7}{9!}z^9 + \cdots \\
g(z) &=& \sum_{n=0}^\infty\frac{(3n-1)!!!}{(3n+1)!}z^{3n+1}
  = z + \frac{2}{4!}z^4 + \frac{2\cdot5}{7!}z^7 + \frac{2\cdot5\cdot8}{10!}z^{10} + \cdots
\end{eqnarray*}
where $n!!!=\max(n,1)$ for $n\leq2$, otherwise $n!!!=n\cdot(n-3)!!!$.
\begin{code}
airy_ai__series !z =
    let !z3 = z^3
        !aiterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        !ai'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in ai0 * (ksum aiterms) + ai'0 * (ksum ai'terms)
\end{code}

The asymptotic expansion for $\Ai(z)$ when $z\to\infty$ with $|\ph z|\leq\pi-\delta$ is given by
\[ \Ai(z) \sim \frac{e^{-\zeta}}{2\sqrt{\pi}z^{1/4}}\sum_{k=0}^\infty (-)^k \frac{u_k}{\zeta^k} \]
where $\zeta=(2/3)z^{3/2}$ and where (with $u_0=1$)
\[ u_{k} = \frac{(2k+1)(2k+3)\cdots(6k-1)}{216^k k!} = \frac{(6k-5)(6k-3)(6k-1)}{(2k-1)216 k} u_{k-1} \]
\begin{code}
airy_ai__asympt_pos :: (Value v) => v -> v
airy_ai__asympt_pos z =
  let !zeta = z**(3/2)*2/3
      !uk = ixiter 1 1 $ \ k u -> let k'=(#)k in u*(6*k'-5)*(6*k'-3)*(6*k'-1)/(2*k'-1)/216/k'
      !zn = iterate (/(-zeta)) 1
      !tterms = zipWith (*) uk zn
      !terms = tk tterms
  in (sf_exp(-zeta))/(2*(sf_sqrt pi)*z**(1/4)) * (ksum terms)
  where tk (a:b:c:ts) = if (rabs b)<(rabs c) then [a] else a:(tk$b:c:ts)
\end{code}

We also have the asymptotic expansion
\[ \Ai(-z) \sim \frac{1}{\sqrt\pi z^{1/4}} \left(
    \cos(\zeta-\frac\pi4) \sum_{k=0}^\infty (-)^k\frac{u_{2k}}{\zeta^{2k}}
   + \sin(\zeta-\frac\pi4) \sum_{k=0}^\infty (-)^k\frac{u_{2k=1}}{\zeta^{2k+1}}
  \right) \]
\begin{code}
airy_ai__asympt_neg :: (Value v) => v -> v
airy_ai__asympt_neg z' =
  let !z = -z'
      !zeta = z**(3/2)*2/3
      !zp4 = zeta - pi/4
      !uk = ixiter 1 1 $ \ k u -> let k'=(#)k in u*(6*k'-5)*(6*k'-3)*(6*k'-1)/(2*k'-1)/216/k'
      !uke = evel uk
      !uko = evel (tail uk)
      !eterms = tk $ zipWith (*) uke (iterate (/(-zeta^2)) 1)
      !oterms = tk $ zipWith (*) uko (iterate (/(-zeta^2)) (1/zeta))
  in ((sf_cos zp4)*(ksum eterms) + (sf_sin zp4)*(ksum oterms))/((sf_sqrt pi)*z**(1/4))
  where tk (a:b:c:ts) = if (rabs b)<(rabs c) then [a] else a:(tk$b:c:ts)
        evel (a:b:cs) = a:(evel cs)
\end{code}

\subsection{Bi}

\subsubsection{\tt sf\_airy\_bi z}
For now, we use an asymptotic expansion for large values and a series for smaller values.
This gives reasonable results for small-enough or large-enough values, but it has
low accuracy for intermediate values, ({\it e.g.\@} $z=5$).
(TODO: Seems quite bad for complex values)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_airy\_bi z} = \Bi(z)$}
\begin{code}
sf_airy_bi :: (Value v) => v -> v
sf_airy_bi !z
  | (rabs z)>=9 && (re z)>=0 = airy_bi__asympt_pos z
  | (rabs z)>=9 && (re z)< 0 = airy_bi__asympt_neg z
  | otherwise                = airy_bi__series z
\end{code}
\end{titled-frame}

Initial conditions
\begin{eqnarray*}
\Bi(0)  &=& \frac{1}{3^{1/6} \Gamma(2/3)} \\
\Bi'(0) &=& \frac{3^{1/6}}{\Gamma(1/3)}
\end{eqnarray*}
\begin{code}
bi0 :: (Value v) => v
bi0 = 3**(-1/6)/sf_gamma(2/3)
bi'0 :: (Value v) => v
bi'0 = 3**(1/6)/sf_gamma(1/3)
\end{code}

Series expansion, where $n!!!=\max(n,1)$ for $n\leq2$ and otherwise $n!!!=n\cdot(n-3)!!!$:
\[ \Bi(z) = \Bi(0)\left(\sum_{n=0}^\infty\frac{(3n-2)!!!}{(3n)!}z^{3n}\right)
    + \Bi'(0)\left(\frac{(3n-1)!!!}{(3n+1)!}z^{3n+1}\right) \]
\begin{code}
airy_bi__series z =
    let !z3 = z^3
        !biterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        !bi'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in bi0 * (ksum biterms) + bi'0 * (ksum bi'terms)
\end{code}

The asymptotic expansion for $\Bi(z)$ when $z\to\infty$ with $|\ph z|\leq\pi-\delta$ is given by
\[ \Bi(z) \sim \frac{e^{\zeta}}{\sqrt{\pi} z^{1/4}}\sum_{k=0}^\infty \frac{u_k}{\zeta^k} \]
where $\zeta=(2/3)z^{3/2}$ and where (with $u_0=1$)
\[ u_{k} = \frac{(2k+1)(2k+3)\cdots(6k-1)}{216^k k!} = \frac{(6k-5)(6k-3)(6k-1)}{(2k-1)216 k} u_{k-1} \]
\begin{code}
airy_bi__asympt_pos :: (Value v) => v -> v
airy_bi__asympt_pos z =
  let !zeta = z**(3/2)*2/3
      !uk = ixiter 1 1 $ \ k u -> let k'=(#)k in u*(6*k'-5)*(6*k'-3)*(6*k'-1)/(2*k'-1)/216/k'
      !zn = iterate (/zeta) 1
      !tterms = zipWith (*) uk zn
      !terms = tk tterms
  in (sf_exp(zeta))/((sf_sqrt pi)*(z**(1/4))) * (ksum terms)
  where tk !(a:b:c:ts) = if (rabs b)<(rabs c) then [a] else a:(tk$b:c:ts)
\end{code}

We also have the asymptotic expansion
\[ \Bi(-z) \sim \frac{1}{\sqrt\pi z^{1/4}} \left(
    \cos(\zeta-\frac\pi4) \sum_{k=0}^\infty (-)^k\frac{u_{2k}}{\zeta^{2k}}
   + \sin(\zeta-\frac\pi4) \sum_{k=0}^\infty (-)^k\frac{u_{2k=1}}{\zeta^{2k+1}}
  \right) \]
\begin{code}
airy_bi__asympt_neg :: (Value v) => v -> v
airy_bi__asympt_neg z' =
  let !z = -z'
      !zeta = z**(3/2)*2/3
      !zp4 = zeta - pi/4
      !uk = ixiter 1 1 $ \ k u -> let k'=(#)k in u*(6*k'-5)*(6*k'-3)*(6*k'-1)/(2*k'-1)/216/k'
      !uke = evel uk
      !uko = evel (tail uk)
      !eterms = tk $ zipWith (*) uke (iterate (/(-zeta^2)) 1)
      !oterms = tk $ zipWith (*) uko (iterate (/(-zeta^2)) (1/zeta))
  in (-(sf_sin zp4)*(ksum eterms) + (sf_cos zp4)*(ksum oterms))/((sf_sqrt pi)*z**(1/4))
  where tk (a:b:c:ts) = if (rabs b)<(rabs c) then [a] else a:(tk$b:c:ts)
        evel (a:b:cs) = a:(evel cs)
\end{code}


