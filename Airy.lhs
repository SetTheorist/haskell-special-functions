\section{Airy}

The Airy functions $Ai$ and $Bi$,
standard solutions of the ode $y''-zy=0$.

\subsection{Preamble}
A basic preamble.
\begin{code}
module Airy (sf_airy_ai, sf_airy_bi) where
import Gamma
import Util
\end{code}

\subsection{Ai}

For now, just use a simple series expansion.
\begin{code}
sf_airy_ai :: (Value v) => v -> v
sf_airy_ai z = airy_ai_series z
\end{code}

Initial conditions
$Ai(0) = 3^{-2/3}\frac{1}{\Gamma(2/3)}$
and 
$Ai'(0) = -3^{-1/3}\frac{1}{\Gamma(1/3)}$
\begin{code}
ai0 :: (Value v) => v
ai0 = 3**(-2/3)/sf_gamma(2/3)

ai'0 :: (Value v) => v
ai'0 = -3**(-1/3)/sf_gamma(1/3)
\end{code}

Series expansion, where $n!!!=\max(n,1)$ for $n\leq2$ and otherwise $n!!!=n\cdot(n-3)!!!$:
\[ Ai(z) = Ai(0)\left(\sum_{n=0}^\infty\frac{(3n-2)!!!}{(3n)!}z^{3n}\right) + Ai'(0)\left(\frac{(3n-1)!!!}{(3n+1)!}z^{3n+1}\right) \]
\begin{code}
airy_ai_series z =
    let z3 = z^3
        aiterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        ai'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in ai0 * (ksum aiterms) + ai'0 * (ksum ai'terms)
\end{code}

\subsection{Bi}

For now, just use a simple series expansion.
\begin{code}
sf_airy_bi :: (Value v) => v -> v
sf_airy_bi z = airy_bi_series z
\end{code}

Initial conditions
$Bi(0) = 3^{-1/6}\frac{1}{\Gamma(2/3)}$
and 
$Bi'(0) = 3^{1/6}\frac{1}{\Gamma(1/3)}$
\begin{code}
bi0 :: (Value v) => v
bi0 = 3**(-1/6)/sf_gamma(2/3)

bi'0 :: (Value v) => v
bi'0 = 3**(1/6)/sf_gamma(1/3)
\end{code}

Series expansion, where $n!!!=\max(n,1)$ for $n\leq2$ and otherwise $n!!!=n\cdot(n-3)!!!$:
\[ Bi(z) = Bi(0)\left(\sum_{n=0}^\infty\frac{(3n-2)!!!}{(3n)!}z^{3n}\right) + Bi'(0)\left(\frac{(3n-1)!!!}{(3n+1)!}z^{3n+1}\right) \]
\begin{code}
airy_bi_series z =
    let z3 = z^3
        biterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        bi'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in bi0 * (ksum biterms) + bi'0 * (ksum bi'terms)
\end{code}

