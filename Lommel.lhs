\chapter{Lommel functions}

\section{Preamble}
\begin{titled-frame}{\color{blue}\tt module Lommel}
\begin{code}
{-# Language BangPatterns #-}
module Lommel (sf_lommel_s, sf_lommel_s2) where
import Bessel
import Trig
import Util
\end{code}
\end{titled-frame}

--TODO: These are completely untested!

\section{First Lommel function}
For $\mu\pm\nu \neq \pm1, \pm3, \pm5, \cdots$ we define
the first Lommel function $\verb|sf_lommel_s mu nu z| = S_{\mu,\nu}(z)$ 
via series-expansion:
\[ S_{\mu,\nu}(z) = \frac{z^{mu+1}}{(\mu+1)^2-\nu^2} \sum_{k=0}^\infty t_k \]
where
\[ t_0=1 \qquad t_{k}=t_{k-1}\frac{-z^2}{(\mu+2k+1)^2-\nu^2} \]

\subsection{\tt sf\_lommel\_s mu nu z}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_lommel\_s mu nu z} = S_{\mu,\nu}(z)$}
\begin{code}
sf_lommel_s mu nu z = 
  let terms = ixiter 1 1.0 $ \ k t -> -t*z^2 / ((mu+((#)$2*k+1))^2 - nu^2)
      res = ksum terms
  in res * z**(mu+1) / ((mu+1)^2 - nu^2)
\end{code}
\end{titled-frame}

\section{Second Lommel function}
For $\mu\pm\nu \neq \pm1, \pm3, \pm5, \cdots$ 
the second Lommel function $\verb|sf_lommel_s2 mu nu z| = s_{\mu,\nu}(z)$ 
is given via an asymptotic expansion:
\[ s_{\mu,\nu}(z) \sim \sum_{k=0}^\infty u_k \]
where
\[ u_0=1 \qquad u_{k}=u_{k-1}\frac{-(\mu-2k+1)^2-\nu^2}{z^2} \]

\subsection{\tt sf\_lommel\_s2 mu nu z}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_lommel\_s2 mu nu z} = s_{\mu,\nu}(z)$}
\begin{code}
sf_lommel_s2 mu nu z =
  let tterms = ixiter 1 1.0 $ \ k t -> -t*((mu-((#)$2*k+1))^2 - nu^2) / z^2
      terms = tk tterms
      res = ksum terms
  in res
  where tk (a:b:cs) = if (rabs a)<(rabs b) then [a] else a:(tk$b:cs)
\end{code}
\end{titled-frame}


\section{Weber E function}
Compute Weber's function $E_\nu(z)$

TODO: untested, not sure where it even converges...

\begin{code}
sf_weber_e :: (Value v) => v -> v -> v
sf_weber_e !nu !z
  | (rabs z)<10 || (rabs (z/nu))<1 = weber__series nu z
  | otherwise = weber__asympt nu z
\end{code}

Series expansions
\begin{code}
weber__series :: (Value v) => v -> v -> v
weber__series !nu !z =
  let !sinnupi = (sf_sin $ nu*pi)/pi
      !cosnupi = sf_cos $ nu*pi
      !lom_0   = sf_lommel_s   0  nu z
      !lom_1   = sf_lommel_s (-1) nu z
  in -(1 + cosnupi) * lom_0/pi -nu*(1 - cosnupi) * lom_1/pi
\end{code}

Asymptotic expansions
\begin{code}
weber__asympt :: (Value v) => v -> v -> v
weber__asympt !nu !z =
  let !cosnupi = sf_cos $ nu*pi
      !lom2_0   = sf_lommel_s2   0  nu z
      !lom2_1   = sf_lommel_s2 (-1) nu z
  in -(sf_bessel_y z  nu) - (1 + cosnupi) * lom2_0 / (pi*z) -nu * (1 - cosnupi) * lom2_1 / (pi*z^2)
\end{code}

\begin{code}
{--
# quadrature --- very slow, not a good approach
function res = anger_quad(nu, z)
  if (imag(z)!=0 || imag(nu)!=0)
    qre = quad(@(th)(real( sf_cos(nu*th - z*sf_sin(th))/pi )), 0, pi, 1e-8);
    qim = quad(@(th)(imag( sf_cos(nu*th - z*sf_sin(th))/pi )), 0, pi, 1e-8);
    res = qre + I*qim
  else
    res = quad(@(th)( sf_cos(nu*th - z*sf_sin(th))/pi ), 0, pi, 1e-8)
  endif
endfunction
--}
\end{code}

