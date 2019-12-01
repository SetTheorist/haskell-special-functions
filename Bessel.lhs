\section{Bessel Functions}

Bessel's differential equation is:
\begin{equation} z^2 w'' + z w' + (z^2-\nu^2)w = 0 \label{bessel:ode}\end{equation}

When $\nu$ is not an integer, then this has two
linearly independent solutions~$J_{\pm}\nu(z)$.
If $\nu=n$ is an integer, then $J_n(z)$ is still a solution,
but $J_{-n}(z)=(-)^nJ_n(z)$ so it is not a second linearly
independent solution of Eqn.~\ref{bessel:ode}.

\subsection{Preamble}
\begin{titled-frame}{\color{blue}\tt module Bessel}
\begin{code}
{-# Language BangPatterns #-}
{-# Language ScopedTypeVariables #-}
module Bessel where
import Gamma
import Trig
import Util
\end{code}
\end{titled-frame}

\subsection{Bessel function $J$ of the first kind}

The Bessel functions $J_\nu(z)$ are defined as 

\subsubsection{\tt sf\_bessel\_j nu z}
Compute Bessel $J\_\nu(z)$ function
\begin{titled-frame}{$\text{\color{blue}\tt sf\_bessel\_j nu z} = J_\nu(z)$}
\begin{code}
sf_bessel_j :: (Value v) => v -> v -> v
sf_bessel_j nu z 
  | (rabs z) > (15+rabs(nu)) = bessel_j__asympt_z nu z
  | otherwise                = bessel_j__series nu z
  --rec = recur_back(z, nu);
  --ref = recur_fore(z, nu);
  --re2 = recur_backwards(nu, z, round(abs(max(z, nu)))+21);
\end{code}
\end{titled-frame}

\subsubsection*{*\tt bessel\_j\_\_series nu z}
The power-series expansion given by
\[ J_\nu(z) = \left(\frac{z}{2}\right)^\nu \frac{1}{1+\nu} \sum_{k=0}^\infty(-)^k\frac{z^{2k}}{2^{2k}k!\Gamma(\nu+k+1)} \]
\begin{titled-frame}{$\text{\color{blue}\tt bessel\_j\_\_series nu z}$}
\begin{code}
bessel_j__series :: (Value v) => v -> v -> v
bessel_j__series !nu !z = 
  let !z2 = -(z/2)^2
      !terms = ixiter 1 1 $ \ !n !t -> t*z2/((#)n)/(nu+(#)n)
      !res = ksum terms
  in res * (z/2)**nu / sf_gamma(1+nu)
\end{code}
\end{titled-frame}

\subsubsection*{*\tt bessel\_j\_\_asympt nu z}
Asymptotic expansion for $|z|>>\nu$ with $|arg z|<\pi$. is given by
\[ J_\nu(z) \sim \left(\frac{2}{\pi z}\right)^{1/2}\
    \left( \cos\omega \sum_{k=0}^\infty (-)^k \frac{a_{2k}(\nu)}{z^{2k}}
        -  \sin\omega \sum_{k=0}^\infty (-)^k \frac{a_{2k+1}(\nu)}{z^{2k+1}} \right) \]
where $\omega = z - \frac{\pi\nu}{2} - \frac{\pi}{4}$ and
\[ a_k(\nu) = \frac{(4\nu^2-1^2)(4\nu^2-3^2)\cdots(4\nu^2-(2k-1)^2)}{k! 8^k} \]
(with $a_0(\nu) = 1$).
\begin{code}
bessel_j__asympt_z :: forall v.(Value v) => v -> v -> v
bessel_j__asympt_z !nu !z =
  let !om = z - (nu/2 + 1/4)*pi
      !nu2 = nu^2
      !aks = ixiter 1 1 $ \k t -> t*(4*nu2 - ((#)$((2*k-1)^2)))/((#)$8*k)/z
      !akse = tk $ zipWith (\k t -> (-1)^k*t) [0..] (evel aks)
      !akso = tk $ zipWith (\k t -> (-1)^k*t) [0..] (evel (tail aks))
  in (sf_sqrt(2/pi/z))*((sf_cos om)*(ksum akse) - (sf_sin om)*(ksum akso))
  where
    tk :: [v] -> [v]
    tk (a:b:c:xs) = if (rabs b) < (rabs c) then [a] else a:(tk (b:c:xs))
    evel (a:b:cs) = a:(evel cs)
\end{code}

This approach uses the recursion in order (for large order) in a backward direction
\[ J_{\nu-1}(z) = \frac{2\nu}{z} J_{\nu}(z) - J_{\nu+1}(z) \]
(largest to smallest).
We start by iterating downward from 20 terms above the largest order we'd like with
initial values~0 and~1.  We then compute the initial (smallest order) term and
scale the whole series with the iterated value and the computed value.
\begin{code}
--bessel_j_recur_back :: (Value v) => Double -> v -> v
bessel_j_recur_back :: forall v.(Value v) => (RealKind v) -> v -> [v]
bessel_j_recur_back !nu !z =
  let !jjs = runback (nnx-2) [1.0,0.0]
      !scale = if (rabs z)<10 then (bessel_j__series nuf z) else (bessel_j__asympt_z nuf z)
      --scale2 = ((head jjs)^2) + 2*(ksum (map (^2) $ tail jjs)) -- only integral nu
  --in jjs!!(nnn) * scale / (jjs!!0)
  in map (\j -> j * scale / (jjs!!0)) (take (nnn+1) jjs)
  --in map (\j -> j/scale2) (take (nnn+1) jjs)
  where
    !nnn = truncate nu
    !nuf = fromReal $ nu - (#)nnn
    !nnx = nnn + 20
    runback :: Int -> [v] -> [v]
    runback !0 !j = j
    runback !nx !j@(jj1:jj2:jjs) =
      let !jj = jj1*2*(nuf+(#)nx)/z - jj2
      in runback (nx-1) (jj:j)
\end{code}

