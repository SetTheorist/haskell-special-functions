\section{Riemann zeta function}

\subsection{Preamble}
\begin{code}
{-# Language BangPatterns #-}
module Zeta (
  sf_zeta,
  sf_zeta_m1,
) where
import Gamma
import Trig
import Util
\end{code}

\subsection{Zeta}

\subsubsection{\tt sf\_zeta z}
Compute the Riemann zeta function $\verb|sf_zeta z| = \zeta(z)$ where
\[ \zeta(z) = \sum_{n=1}^\infty n^{-z} \]
(for $\Re z>1$ and defined by analytic continuation elsewhere).
\begin{code}
sf_zeta :: (Value v) => v -> v
sf_zeta z
  | z==1      = (1/0)
  | (re z)<0  = 2 * (2*pi)**(z-1) * (sf_sin$pi*z/2) * (sf_gamma$1-z) * (sf_zeta$1-z)
  | otherwise = zeta_series 1.0 z

\end{code}

\subsubsection{\tt sf\_zeta\_m1 z}
For numerical purposes, it is useful to have $\verb|sf_zeta_m1 z| = \zeta(z)-1$.
\begin{code}
sf_zeta_m1 :: (Value v) => v -> v
sf_zeta_m1 z
  | z==1      = (1/0)
  | (re z)<0  = 2 * (2*pi)**(z-1) * (sf_sin$pi*z/2) * (sf_gamma$1-z) * (sf_zeta$1-z) - 1  -- TODO:
  | otherwise = zeta_series 0.0 z
\end{code}

\subsubsection*{*\tt zeta\_series i z}
We use the simple series expansion for $\zeta(z)$ with an
Euler-Maclaurin correction:
\[ \zeta(z) = \sum_{n=1}^{N}\frac{1}{n^z} + \sum_{k=1}^{p}\cdots \]
\begin{code}
zeta_series :: (Value v) => v -> v -> v
zeta_series !init !z = 
  let terms = map (\n->((#)n)**(-z)) [2..]
      corrs = map correction [2..]
  in summer terms corrs init 0.0 0.0
  where
    --TODO: make general "corrected" kahan_sum!
    summer !(t:ts) !(c:cs) !s !e !r = 
      let !y = t + e
          !s' = s + y
          !e' = (s - s') + y
          !r' = s' + c + e'
      in if r==r' then r'
         else summer ts cs s' e' r'
    !zz1 = z/12
    !zz2 = z*(z+1)*(z+2)/720
    !zz3 = z*(z+1)*(z+2)*(z+3)*(z+4)/30240
    !zz4 = z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6)/1209600
    !zz5 = z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6)*(z+7)*(z+8)/239500800
    correction !n' =
      let n=(#)n'
      in n**(1-z)/(z-1) - n**(-z)/2
         + n**(-z-1)*zz1 - n**(-z-3)*zz2 + n**(-z-5)*zz3
         - n**(-z-7)*zz4 + n**(-z-9)*zz5
\end{code}


