\section{Exponential Integral}

\subsection{Preamble}
\begin{titled-frame}{\color{blue}\tt module ExpInt}
\begin{code}
module ExpInt(sf_expint_ei, sf_expint_en)
where
import Exp
import Gamma
import Util
\end{code}
\end{titled-frame}

\subsection{Exponential integral $\Ei$}
The exponential integral $\Ei z$ is defined for $x<0$ by
\[ \Ei(z) = -\int_{-x}^\infty \frac{e^{-t}}{t}\,dt \]
It can be defined 

\subsubsection{\tt sf\_expint\_ei z}
We give only an implementation for $\Re z\geq0$.
We use a series expansion for $|z|<40$ and an asymptotic expansion otherwise.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_expint\_ei z} = \Ei(z)$\marginnote{\tt sf\_expint\_ei}}
\begin{code}
sf_expint_ei :: (Value v) => v -> v
sf_expint_ei z
  | (re z) < 0.0  = (0/0)  -- (NaN)
  | z == 0.0      = (-1/0) -- (-Inf)
  | (rabs z) < 40 = expint_ei__series z
  | otherwise     = expint_ei__asymp z
\end{code}
\end{titled-frame}

The series expansion is given (for $x>0$)
\[ \Ei(x) = \gamma + \ln x + \sum_{n=1}^\infty \frac{x^n}{n! n} \]
We evaluate the addition of the two terms with the sum slightly differently
when $\Re z<1/2$ to reduce floating-point cancellation error slightly.
\begin{titled-frame}{{\color{blue}\tt expint\_ei\_\_series z}\marginnote{\tt expint\_ei\_\_series}}
\begin{code}
expint_ei__series :: (Value v) => v -> v
expint_ei__series z =
  let tterms = ixiter 2 z $ \n t -> t*z/(#)n
      terms = zipWith (\ t n ->t/(#)n) tterms [1..]
      res = ksum terms
  in if (re z)<0.5
     then sf_log(z * sf_exp(euler_gamma + res))
     else res + sf_log(z) + euler_gamma
\end{code}
\end{titled-frame}

The asymptotic expansion as $x\to+\infty$ is
\[ \Ei(x) \sim \frac{e^x}{x}\sum_{n=0}^\infty \frac{n!}{x^n} \]
\begin{titled-frame}{{\color{blue}\tt expint\_ei\_\_asymp z}\marginnote{\tt expint\_ei\_\_asymp}}
\begin{code}
expint_ei__asymp :: (Value v) => v -> v
expint_ei__asymp z =
  let terms = tk $ ixiter 1 1.0 $ \n t -> t/z*(#)n
      res = ksum terms
  in res * (sf_exp z) / z
  where tk (a:b:cs) = if (rabs a)<(rabs b) then [a] else a:(tk$b:cs)
\end{code}
\end{titled-frame}

\subsection{Exponential integral $E_n$}
The exponential integrals $E_n(z)$ are defined as
\[ E_n(z) = z^{n-1}\int_{z}^\infty \frac{e^{-t}}{t^n}\,dt \]
They satisfy the following relations:
\begin{eqnarray*}
E_0(z) &=& \frac{e^{-z}}{z} \\
E_{n+1}(z) &=& \int_z^\infty E_{n}(t)\,dt \\
\end{eqnarray*}
And they can be expressed in terms of incomplete gamma functions:
\[ E_n(z) = z^{n-1}\Gamma(1-n,z) \]
(which also gives a generalization for non-integer $n$).

\subsubsection{\tt sf\_expint\_en n z}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_expint\_en n z} = E_n(z)$\marginnote{\tt sf\_expint\_en}}
\begin{code}
sf_expint_en :: (Value v) => Int -> v -> v
sf_expint_en n z | (re z)<0 = (0/0) -- (NaN) TODO: confirm this
                 | z == 0   = (1/(#)(n-1)) -- TODO: confirm this
sf_expint_en 0 z = sf_exp(-z) / z
sf_expint_en 1 z = expint_en__1 z
sf_expint_en n z | (rabs z) <= 1.0 = expint_en__series n z
                 | otherwise = expint_en__contfrac n z
\end{code}
\end{titled-frame}

We use this series expansion for $E_1(z)$:
\[ E_1(z) = -\gamma - \ln z + \sum_{k=1}^\infty(-)^k\frac{z^k}{k! k}\]
(Note that this will not be good for large values of $z$.)
\begin{code}
expint_en__1 :: (Value v) => v -> v
expint_en__1 z =
  let r0 = -euler_gamma - (sf_log z)
      tterms = ixiter 2 (z) $ \k t -> -t*z/(#)k
      terms = zipWith (\ t k -> t/(#)k) tterms [1..]
  in ksum (r0:terms)
\end{code}

\begin{code}
-- assume n>=2, z<=1
expint_en__series :: (Value v) => Int -> v -> v
expint_en__series n z =
  let n' = (#)n
      res = (-(sf_log z) + (sf_digamma n')) * (-z)^(n-1)/(#)(factorial$n-1) + 1/(n'-1)
      terms' = ixiter 2 (-z) (\m t -> -t*z/(#)m)
      terms = map (\(m,t)->(-t)/(#)(m-(n-1))) $ filter ((/=(n-1)) . fst) $ zip [1..] terms'
  in ksum (res:terms)

-- assume n>=2, z>1
-- modified Lentz algorithm
expint_en__contfrac :: (Value v) => Int -> v -> v
expint_en__contfrac n z =
  let fj = zeta
      cj = fj
      dj = 0
      j = 1
      n' = (#)n
  in lentz j cj dj fj
  where
    zeta = 1e-100
    eps = 5e-16
    nz x = if x==0 then zeta else x
    lentz j cj dj fj =
      let aj = (#) $ if j==1 then 1 else -(j-1)*(n+j-2)
          bj = z + (#)(n + 2*(j-1))
          dj' = nz $ bj + aj*dj
          cj' = nz $ bj + aj/cj
          dji = 1/dj'
          delta = cj'*dji
          fj' = fj*delta
      in if (rabs$delta-1)<eps
         then fj' * sf_exp(-z)
         else lentz (j+1) cj' dji fj'
\end{code}
