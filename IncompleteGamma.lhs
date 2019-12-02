\chapter{Incomplete Gamma}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Preamble}
A basic preamble.
\begin{titled-frame}{\color{blue}\tt module IncompleteGamma}
\begin{code}
{-# Language BangPatterns #-}
module IncompleteGamma (sf_incomplete_gamma, sf_incomplete_gamma_co) where
import Exp
import Gamma
import Numbers(factorial)
import Trig
import Util
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Incomplete Gamma functions}

We define the two basic incomplete Gamma functions
(the incomplete Gamma function and the complementary incomplete Gamma function, resp.\@)
via
\begin{eqnarray*}
\Gamma(a,z) &=& \int_z^\infty e^{-t}t^{a} \,\frac{dt}{t} \marginnote{$\Gamma(a,z)$} \marginnote{$\Gamma(a,z)$}\\
\gamma(a,z) &=& \int_0^z e^{-t}t^{a} \,\frac{dt}{t}      \marginnote{$\gamma(a,z)$} \marginnote{$\gamma(a,z)$}
\end{eqnarray*}
where we clearly have $\Gamma(a,z) + \gamma(a,z) = \Gamma(a)$.

\subsection{\tt sf\_incomplete\_gamma a z}
The incomplete gamma function implemented via ...
Seems to work okay for $z>0$, not great for complex values.
Untested for $z<0$.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_incomplete\_gamma a z} = \Gamma(a,z)$}
\begin{code}
sf_incomplete_gamma :: (Value v) => v -> v -> v
sf_incomplete_gamma a z 
  | (rabs z)>(rabs a) && (re z)<5 = incgam__contfrac a z
  | (rabs z)>(rabs a)             = incgam__asympt_z a z
  | otherwise = (sf_gamma a) - (incgamco__series a z)
\end{code}
\end{titled-frame}

\subsubsection{\tt incgam\_\_contfrac}
This continued fraction expansion converges for $\Re z>0$ where $v=1/z$:
(Perhaps even for $|\ph z|<\pi$.)
\[ e^{z} z^{1-a} \Gamma(a,z) = \frac{1}{1+{}}%
    \ \frac{(1-a)v}{1+{}}\ \frac{v}{1+{}}\ \frac{(2-a)v}{1+{}}\ \frac{2v}{1+{}}%
    \ \frac{(3-a)v}{1+{}}\ \frac{3v}{1+\cdots} \]
Seems to work best for $z>a$
\begin{code}
incgam__contfrac !a !z =
  let !v = 1/z
      !ane = map (\k->v*((#)k-a)) [1..]
      !ano = map (\k->v*((#)k)) [1..]
      !an = 1:(merge ane ano)
      !bn = 0:(repeat 1)
  in (sf_exp(-z))*(z**(a-1))*(sf_cf_lentz an bn)
  where merge (a:as) bs = a:(merge bs as)
\end{code}
Can be written in an equivalent form
\[ \Gamma(a,z) = \frac{e^{-z}z^a}{z+{}}\ \frac{1-a}{1+{}}\ \frac{1}{z+{}}%
    \ \frac{2-a}{1+{}}\ \frac{2}{z+\cdots} \]
\begin{code}
incgam__contfrac' !a !z =
  let !ane = map (\k->((#)k-a)) [1..]
      !ano = map (\k->((#)k)) [1..]
      !an = ((sf_exp(-z))*(z**a)):(merge ane ano)
      !bn = 0:(merge (repeat z) (repeat 1))
  in (sf_cf_lentz an bn)
  where merge (a:as) bs = a:(merge bs as)
\end{code}

\subsubsection{\tt incgam\_\_asympt\_z}
We have the asymptotic expansion as $z\to\infty$
\[ \Gamma(a,z) \sim z^{a-1}e^{-z} \sum_{n=0}^\infty \frac{\Gamma(a)}{\Gamma(a-n)}z^{-n} \]
This seems to give good results for $z>a$.
\begin{code}
incgam__asympt_z !a !z =
  let tterms = ixiter 1 1 $ \n t -> t*(a-(#)n)/z
      terms = tk tterms
  in z**(a-1) * (sf_exp(-z)) * (ksum terms)
  where
    tk (a:b:c:ts) = if (rabs b)<(rabs c) then [a] else a:(tk$b:c:ts)
\end{code}

\subsection{\tt sf\_incomplete\_gamma\_co a z}
The complementary incomplete Gamma function implemented via
TODO: this is just a quick hack implementation!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_incomplete\_gamma\_co a z} = \gamma(a,z)$}
\begin{code}
sf_incomplete_gamma_co :: (Value v) => v -> v -> v
sf_incomplete_gamma_co a z
  | (rabs z)>(rabs a) && (re z)<5 = (sf_gamma a) - (incgam__contfrac a z)
  | (rabs z)>(rabs a)             = (sf_gamma a) - (incgam__asympt_z a z)
  | otherwise = (incgamco__series a z)
\end{code}
\end{titled-frame}

\subsubsection{\tt incgamco\_\_series}
A series expansion for the complementary incomplete Gamma function
where $a\neq0,-1,-2,\dots$
\[ \gamma(a,z) = e^{-x}\sum_{k=0}^\infty \frac{x^{a+k}}{(a)_{k+1}} \]
This should converge well for $a\geq z$.
\begin{code}
incgamco__series a z =
  let terms = ixiter 1 (z**a / a) $ \k t -> t*z/(a+(#)k)
  in (sf_exp(-z)) * (ksum terms)
\end{code}



