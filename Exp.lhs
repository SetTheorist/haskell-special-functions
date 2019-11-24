\section{Exponential \& Logarithm}

In this section, we implement the exponential function
and logarithm function, as well as useful variations.

\subsection{Preamble}
We begin with a typical preamble.

\begin{code}
{-# Language BangPatterns #-}
{-# Language FlexibleInstances #-}
module Exp (
    sf_exp, sf_expn, sf_exp_m1, sf_exp_m1vx, sf_exp_men, sf_exp_menx,
    sf_log, sf_log_p1,
) where
import Numbers
import Util
\end{code}

\subsection{Exponential}

We start with implementation of the most basic special function, $exp(x)$ or $e^x$
and variations thereof.

\subsubsection{\tt sf\_exp x}
For the exponential $\verb|sf_exp x|=\exp(x)$ we use a simple series expansions
\[ e^x = \sum_{n=0}^\infty\frac{x^n}{n!} \]
after first using the identity $e^{-x}=1/e^x$ to ensure
that the real part of the argument is positive.
This avoids disastrous cancellation for negative arguments,
(though note that for complex arguments this is not sufficient.)
TODO: should do range-reduction first...
TODO: maybe for complex, use explicit cis?
%\marginpar{\tt sf\_exp} % ugly placement
\begin{code}
sf_exp :: (Value v) => v -> v 
sf_exp !x
  | is_inf x  = if (re x)<0 then 0 else (1/0)
  | is_nan x  = x
  | (re x)<0  = 1/(sf_exp (-x))
  | otherwise = kahan_sum $ ixiter 1 1.0 $ \n t -> t*x/((#)n)
\end{code}

\subsubsection{\tt sf\_exp\_m1 x}
For numerical calculations, it is useful to have $\verb|sf_exp_m1 x| = e^x-1$
as explicitly calculating this expression will give poor results for $x$ near~1.
We use a series expansion for the calculation.
Again for negative real part we reflect using $e^{-x}-1 = -e^{-x}(e^x-1)$.
TODO: should do range-reduction first...
TODO: maybe for complex, use explicit cis?
\begin{code}
sf_exp_m1 :: (Value v) => v -> v
sf_exp_m1 !x
  | is_inf x  = if (re x)<0 then -1 else (1/0)
  | is_nan x  = x
  | (re x)<0  = -sf_exp x * sf_exp_m1 (-x)
  | otherwise = kahan_sum $ ixiter 2 x $ \n t -> t*x/((#)n)
\end{code}

\subsubsection{\tt sf\_exp\_m1vx x}
Similarly, it is useful to have the scaled variant $\verb|sf_exp_m1vx x| = \frac{e^x-1}{x}$.
In this case, we use a continued-fraction expansion
\[ \frac{e^x-1}{x} = \frac{2}{2-x+} \frac{x^2/6}{1+}
    \frac{x^2/4\cdot3\cdot5}{1+}
    \frac{x^2/4\cdot5\cdot7}{1+}
    \frac{x^2/4\cdot7\cdot9}{1+}
    \cdots \]
For complex values, simple calculation is inaccurate (when $\Re z\sim 1$).
\begin{code}
sf_exp_m1vx :: (Value v) => v -> v
sf_exp_m1vx !x
  | is_inf x = if (re x)<0 then 0 else (1/0)
  | is_nan x = x
  | rabs(x)>(1/2) = (sf_exp x - 1)/x -- inaccurate for some complex points
  | otherwise =
      let x2 = x^2
      in 2/(2 - x + x2/6/(1
          + x2/(4*(2*3-3)*(2*3-1))/(1
          + x2/(4*(2*4-3)*(2*4-1))/(1
          + x2/(4*(2*5-3)*(2*5-1))/(1
          + x2/(4*(2*6-3)*(2*6-1))/(1
          + x2/(4*(2*7-3)*(2*7-1))/(1
          + x2/(4*(2*8-3)*(2*8-1))/(1
          ))))))));
\end{code}

\subsubsection{\tt sf\_exp\_menx n x}
\begin{code}
-- Compute scaled tail of series expansion of exponential $(e^z - \sum_(k=0)^(n-1) x^k/k!) / (x^n/n!)$
-- ($n=0, 1, 2, ...$)
sf_exp_menx :: (Value v) => Int -> v -> v
sf_exp_menx 0 z = sf_exp z
sf_exp_menx 1 z = sf_exp_m1vx z
sf_exp_menx n z
  | is_inf z  = if (re z)>0 then (1/0) else (0) -- TODO: verify
  | is_nan z  = z
  | otherwise = exp_menx__contfrac n z
  where
    !zeta = 1e-150
    !eps = 1e-16
    nz !z = if z==0 then zeta else z
    -- modified Lentz for continued fraction
    exp_menx__contfrac n z =
      let !fj = (#)$ n+1
          !cj = fj
          !dj = 0
          !j  = 1
      in lentz j dj cj fj
    lentz !j !dj !cj !fj =
      let !aj = if (odd j)
                then z*((#)$(j+1)`div`2)
                else -z*((#)$(n+(j`div`2)))
          bj = (#)$n+1+j
          !dj' = nz$ bj + aj*dj
          !cj' = nz$ bj + aj/cj
          !dji = 1/dj'
          !deltaj = cj'*dji
          !fj' = fj*deltaj
      in if (rabs(deltaj-1)<eps)
         then 1/(1-z/fj')
         else lentz (j+1) dji cj' fj'
\end{code}

\subsubsection{\tt sf\_exp\_men n x}
This is the generalization of \verb|sf_exp_m1 x|, giving the tail of the
series expansion of the exponential function, for $n=0, 1, \dots$.
\[ \verb|sf_exp_men n z| = e^z - \sum_{k=0}^{n-1}\frac{z^k}{k!} = \sum_{k=n}^\infty\frac{z^k}{k!} \]
The special cases are: $n=0$ gives $e^x=\verb|sf_exp x|$ and $n=1$ gives $e^x-1=\verb|sf_exp_m1 x|$.
We compute this by calling the scaled version \verb|sf_exp_menx| and rescaling back.
\begin{code}
-- ($n=0, 1, 2, ...$)
sf_exp_men :: (Value v) => Int -> v -> v
sf_exp_men !n !x = (sf_exp_menx n x) * x^n / ((#)$factorial n)
\end{code}

\begin{code}
-- Compute initial part of series for exponential, $\sum_(k=0)^n z^k/k!$ 
-- ($n=0,1,2,...$)
sf_expn :: (Value v) => Int -> v -> v
sf_expn n z 
  | is_inf z  = if (re z)>0 then (1/0) else (if (odd n) then (-1/0) else (1/0))
  | is_nan z  = z
  | otherwise = expn__series n z
  where
    -- TODO: just call sf_exp when possible
    -- TODO: better handle large -ve values!
    expn__series :: (Value v) => Int -> v -> v
    expn__series n z = kahan_sum $ take (n+1) $ ixiter 1 1.0 $ \k t -> t*z/(#)k

----------------------------------------

sf_log :: (Value v) => v -> v
sf_log = log

-- log(1+x)
sf_log_p1 :: (Value v) => v -> v
sf_log_p1 !z
  | is_nan z = z
  | (rabs z)>0.25 = sf_log (1+z)
  | otherwise = log_p1__series z
  where
    -- ln(1+x) = 2\sum_{n=0}^\infty (x / x+2)^(2n+1) / 2n+1
    log_p1__series z =
      let !r = z/(z+2)
          !zr2 = r^2
          !tterms = iterate (*zr2) (r*zr2)
          !terms = map (\(n,t)->t/((#)$2*n+1)) (zip [1..] tterms)
      in 2*(kahan_sum (r:terms))

----------------------------------------

-- $$\ln(1+z) = z/(1+ z/(2+ z/(3+ 4z/(4+ 4z/(5+ 9z/(6+ 9z/(7+ ...)))))))$$
-- NB not bad convergence
--ln_1_z_cf z = steeds (z:(ts 1)) [0..]
--  where ts n = (n^2*z):(n^2*z):(ts (n+1))
\end{code}

