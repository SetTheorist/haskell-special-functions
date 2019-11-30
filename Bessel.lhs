\section{Bessel Functions}

Bessel's differential equation is:
\[ z^2 w'' + z w' + (z^2-\nu^2)w = 0 \]

\subsection{Preamble}
\begin{code}
{-# Language BangPatterns #-}
{-# Language ScopedTypeVariables #-}
module Bessel where
import Gamma
import Trig
import Util
\end{code}

\subsection{Bessel function $J$ of the first kind}

The Bessel functions $J_\nu(z)$ are defined as

\subsubsection{\tt sf\_bessel\_j nu z}
Compute Bessel $J\_\nu(z)$ function
\begin{titled-frame}{$\text{\color{blue}\tt sf\_bessel\_j nu z} = J_\nu(z)$}
\begin{code}
sf_bessel_j :: (Value v) => v -> v -> v
sf_bessel_j nu z 
  | (rabs z) < 5 = bessel_j__series nu z
  | otherwise    = bessel_j__asympt_z nu z
  --rec = recur_back(z, nu);
  --ref = recur_fore(z, nu);
  --re2 = recur_backwards(nu, z, round(abs(max(z, nu)))+21);
  --res = sys;
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
      !terms = ixiter 1 1 $ \n t -> t*z2/((#)n)/(nu+(#)n)
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
(largest to smallest) wtih ...
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


\begin{code}
{--
-- recursion in order (forewards)
bessel_j_recur_fore :: (Value v) => Double -> v -> v
bessel_j_recur_fore !nu !z =
  let !jj1 = bessel_j_series nuf z
      !jj2 = bessel_j_series (nuf+1) z
  in loop 3 jj1 jj2
  where
    !nnn = truncate nu
    !nuf = nu - (#)nnn
    !nnx = nnn + 10
    loop :: Int -> v -> v -> v
    loop j jjm2 jjm1
      | j==(nnx+1) = jjm1
      | otherwise =
          let jjj = jjm1*2*(fromDouble(nuf+(#)j))/z - jjm2
          in loop (j+1) jjm1 jjj

{--
function res = recur_backwards(n, z, topper)
  jjp2 = zeros(size(z));
  jjp1 = ones(size(z));
  jjp2_e_ = 1e-40 * ones(size(z));
  jjp1_e_ = 1e-20 * ones(size(z));
  scale = 2*ones(size(z));
  res = zeros(size(z));
  for m = (topper-2):(-1):1
    #jj(m) = (2*nu/z)*jj(m+1) - jj(m+2);
    s_ = -jjp2;
    e_ = -jjp2_e_;
    # add high
      t_ = s_;
      y_ = ((2*m./z).*jjp1) + e_;
      s_ = t_ + y_;
      e_ = (t_ - s_) + y_;
    # add low
      t_ = s_;
      y_ = ((2*m./z).*jjp1_e_) + e_;
      s_ = t_ + y_;
      e_ = (t_ - s_) + y_;
    jjp2 = jjp1;
    jjp2_e_ = jjp1_e_;
    jjp1 = s_;
    jjp1_e_ = e_;

    if (m==n+1)
      # store the desired result,
      # but keep recursing to get scale factor
      res = jjp1;
    endif

    if (m!=1)
      scale += 2 * (s_.^2 + e_.^2 + 2*s_.*e_);
    else
      scale += 1 * (s_.^2 + e_.^2 + 2*s_.*e_);
    endif

    if (scale>1e20)
      jjp2 /= 1024;
      jjp2_e_ /= 1024;
      jjp1 /= 1024;
      jjp1_e_ /= 1024;
      res /= 1024;
      scale /= 1024^2;
    endif
  endfor
  res ./= sqrt(scale);
endfunction
--}
--}
\end{code}
