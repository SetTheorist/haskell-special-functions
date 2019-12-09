\chapter{Exponential \& Logarithm}

In this section, we implement the exponential function
and logarithm function, as well as useful variations.

\section{Preamble}
We begin with a typical preamble.

\begin{titled-frame}{\color{blue}\tt module Exp}
\begin{code}
{-# Language BangPatterns #-}
{-# Language ScopedTypeVariables #-}
module Exp (
    sf_exp, sf_expn, sf_exp_m1, sf_exp_m1vx, sf_exp_men, sf_exp_menx,
    sf_log, sf_log_p1,
) where
import Numbers
import Util
import System.IO.Unsafe
\end{code}
\end{titled-frame}

\section{Exponential}

We start with implementation of the most basic special function, $exp(x)$ or $e^x$
and variations thereof.

\subsection{\tt sf\_exp x}
For the exponential $\verb|sf_exp x|=\exp(x)$ we use a simple series expansions
\[ e^x = \sum_{n=0}^\infty\frac{x^n}{n!} \]
after first using the identity $e^{-x}=1/e^x$ to ensure
that the real part of the argument is positive.
This avoids disastrous cancellation for negative arguments,
(though note that for complex arguments this is not sufficient.)

We also do a range-reduction so that we require fewer terms in the series.
We write $x = n\ln 2 + r$ where $|r|<\ln 2$ and then
\[ e^x = e^{n\ln 2 + r} = 2^n e^r \]
TODO: This needs to be done with enhanced precision; currently loses accuracy.
(TODO: maybe for complex, use explicit cis?)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_exp x} = e^x$}
\begin{code}
{-# SPECIALISE sf_exp :: Double -> Double #-}
{-# SPECIALISE sf_exp :: CDouble -> CDouble #-}
sf_exp :: (Value v) => v -> v 
sf_exp !x
  | is_inf x  = if (re x)<0 then 0 else pos_infty
  | is_nan x  = x
  | (re x)<0  = 1/(sf_exp (-x))
  -- | otherwise = ksum $ ixiter 1 1.0 $ \n t -> t*x/(#)n
  | otherwise = -- TODO: split real/imaginary!!!
      let !ln2 = 0.69314718055994530941723212145817656807
          !n   = floor $ (rabs x)/ln2
          !r   = x - (fromReal$ln2*((#)n))
          !sm  = ksum $ ixiter 1 1.0 $ \n t -> t*r/(#)n
      in sm * 2^n
\end{code}
\end{titled-frame}

\subsection{\tt sf\_exp\_m1 x}
For numerical calculations, it is useful to have $\verb|sf_exp_m1 x| = e^x-1$
as explicitly calculating this expression will give poor results for $x$ near~1.
We use a series expansion for the calculation.
Again for negative real part we reflect using $e^{-x}-1 = -e^{-x}(e^x-1)$.
TODO: should do range-reduction first...
TODO: maybe for complex, use explicit cis?
\begin{titled-frame}{$\text{\color{blue}\tt sf\_exp\_m1 x} = e^x-1$}
\begin{code}
{-# SPECIALISE sf_exp_m1 :: Double -> Double #-}
sf_exp_m1 :: (Value v) => v -> v
sf_exp_m1 !x
  | is_inf x  = if (re x)<0 then -1 else pos_infty
  | is_nan x  = x
  | (re x)<0  = -sf_exp x * sf_exp_m1 (-x)
  | otherwise = ksum $ ixiter 2 x $ \n t -> t*x/((#)n)
\end{code}
\end{titled-frame}

\subsection{\tt sf\_exp\_m1vx x}
Similarly, it is useful to have the scaled variant $\verb|sf_exp_m1vx x| = \frac{e^x-1}{x}$.
In this case, we use a continued-fraction expansion
\[ \frac{e^x-1}{x} = \frac{2}{2-x+{}} \frac{x^2/6}{1+{}}
    \frac{x^2/4\cdot3\cdot5}{1+{}}
    \frac{x^2/4\cdot5\cdot7}{1+{}}
    \frac{x^2/4\cdot7\cdot9}{1+{}}
    \cdots \]
For complex values, simple calculation is inaccurate (when $\Re z\sim 1$).
\begin{titled-frame}{$\text{\color{blue}\tt sf\_exp\_m1vx x} = \frac{e^x-1}{x}$}
\begin{code}
{-# SPECIALISE sf_exp_m1vx :: Double -> Double #-}
sf_exp_m1vx :: (Value v) => v -> v
sf_exp_m1vx !x
  | is_inf x = if (re x)<0 then 0 else pos_infty
  | is_nan x = x
  | rabs(x)>(1/2) = (sf_exp x - 1)/x -- inaccurate for some complex points
  | otherwise =
      let x2 = x^2
      in 2/(2 - x + x2/6/(1
          + x2/(4*(2*3-3)*(2*3-1))/(1 + x2/(4*(2*4-3)*(2*4-1))/(1
          + x2/(4*(2*5-3)*(2*5-1))/(1 + x2/(4*(2*6-3)*(2*6-1))/(1
          + x2/(4*(2*7-3)*(2*7-1))/(1 + x2/(4*(2*8-3)*(2*8-1))/(1
          ))))))));
\end{code}
\end{titled-frame}

\subsection{\tt sf\_exp\_menx n x}
Compute the scaled tail of series expansion of the exponential function, $exd_n(x)$:
\begin{eqnarray*}
\verb|sf_exp_menx n x|
    &=& \frac{n!}{x^n} \left(e^z - \sum_{k=0}^{n-1}\frac{x^k}{k!}\right) \\
    &=& \frac{n!}{x^n} \sum_{k=n}^{\infty}\frac{x^k}{k!} \\
    &=& n!\sum_{k=0}^{\infty}\frac{x^{k}}{(k+n)!}
\end{eqnarray*}
We use a continued fraction expansion
\[ exd_n(z) = \frac{1}{1 - {}}\ \frac{z}{(n+1)+{}}\ \frac{z}{(n+2)-{}}%
    \ \frac{(n+1)z}{(n+3)+{}}\ \frac{2 z}{(n+4)-{}}\ \frac{(n+2)z}{(n+5)+{}}%
    \ \frac{3 z}{(n+6)-\cdots} 
\]
which is evaluated with the modified Lentz algorithm.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_exp\_menx n z} = exd_n(x)$}
\begin{code}
{-# SPECIALISE sf_exp_menx :: Int -> Double -> Double #-}
sf_exp_menx :: (Value v) => Int -> v -> v
sf_exp_menx 0 z = sf_exp z
sf_exp_menx 1 z = sf_exp_m1vx z
sf_exp_menx n z
  | is_inf z  = if (re z)>0 then pos_infty else 0
  | is_nan z  = z
  | otherwise =
      let !an = 1:(-z):(map aterm [2..])
          !bn = 0:1:(map (#) [(n+1)..])
      in sf_cf_lentz an bn
      where
        aterm k | even k    = z*((#)(k`div`2))
                | otherwise = -z*((#)(n+(k`div`2)))
\end{code}
\end{titled-frame}

\subsection{\tt sf\_exp\_men n x}
This is the generalization of \verb|sf_exp_m1 x|, giving the tail of the
series expansion of the exponential function, for $n=0, 1, \dots$.
\[ \verb|sf_exp_men n z| = e^z - \sum_{k=0}^{n-1}\frac{z^k}{k!} = \sum_{k=n}^\infty\frac{z^k}{k!} \]
The special cases are: $n=0$ gives $e^x=\verb|sf_exp x|$ and $n=1$ gives $e^x-1=\verb|sf_exp_m1 x|$.
We compute this by calling the scaled version \verb|sf_exp_menx| and rescaling back.
Though note that it this, of course, has the continued fraction expansion:
\[ ex_n(z) = \frac{z^n}{n! - {}}\ \frac{n! z}{(n+1)+{}}\ \frac{z}{(n+2)-{}}%
    \ \frac{(n+1)z}{(n+3)+{}}\ \frac{2 z}{(n+4)-{}}\ \frac{(n+2)z}{(n+5)+{}}%
    \ \frac{3 z}{(n+6)-\cdots} 
\]
\begin{titled-frame}{$\text{\color{blue}\tt sf\_exp\_men n z} = ex_n(x)$}
\begin{code}
sf_exp_men :: (Value v) => Int -> v -> v
sf_exp_men !n !x = (sf_exp_menx n x) * x^n / ((#)$factorial n)
\end{code}
\end{titled-frame}

\subsection{\tt sf\_expn n x}
We compute the initial part of the series for the exponential function
\[ e_n(z) = \sum_{k=0}^n \frac{z^k}{k!} \]
This implementation simply computes the series directly.
Note that this will suffer from catastrophic cancellation for non-small -ve values.
(TODO: just call {\tt sf\_exp} when possible \& handle large -ve values better!)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_exp\_men n z} = ex_n(x)$}
\begin{code}
sf_expn :: forall v.(Value v) => Int -> v -> v
sf_expn n z 
  | is_inf z  = z^n
  | is_nan z  = z
  | otherwise = ksum $ take (n+1) $ ixiter 1 1.0 $ \k t -> t*z/(#)k
\end{code}
\end{titled-frame}


\section{Logarithm}

\subsection{\tt sf\_log x}
We simply use the built-in implementation (from the \verb|Floating| typeclass).
\begin{code}
sf_log :: (Value v) => v -> v
sf_log = log
\end{code}

\subsection{\tt sf\_log\_p1 x}
The accuracy preserving $\verb|sf_log_p1 x|=\ln 1+x$.
For values close to zero, we use a power series expansion
\[ \ln(1+x) = 2\sum_{n=0}^\infty \frac{(\frac{x}{x+2})^{2n+1}}{2n+1} \]
and otherwise just compute it directly.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_log\_p1 z} = \ln z+1$}
\begin{code}
sf_log_p1 :: (Value v) => v -> v
sf_log_p1 !z
  | is_nan z = z
  | (rabs z)>0.25 = sf_log (1+z)
  | otherwise = ser z
  where
    ser z =
      let !r = z/(z+2)
          !r2 = r^2
          !zterms = iterate (*r2) (r*r2)
          !terms = zipWith (\n t -> t/((#)$2*n+1)) [1..] zterms
      in 2*(ksum (r:terms))
\end{code}
\end{titled-frame}

A simple continued fraction implementation for $\ln 1+z$
\[\ln(1+z) = z/(1+ z/(2+ z/(3+ 4z/(4+ 4z/(5+ 9z/(6+ 9z/(7+ \cdots)))))))\]
Though unused for now, it seems to have decent convergence properties.
Steeds may give better results that modified Lentz here.
\begin{code}
ln_1_z_cf z = sf_cf_steeds (z:(ts 1)) (map (#) [0..])
  where ts n = (n^2*z):(n^2*z):(ts (n+1))
ln_1_z_cf' z = sf_cf_lentz (z:(ts 1)) (map (#) [0..])
  where ts n = (n^2*z):(n^2*z):(ts (n+1))
\end{code}

