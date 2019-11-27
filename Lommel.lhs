\section{Lommel functions}

\subsection{Preamble}
\begin{code}
module Lommel (
  sf_lommel_s,
  sf_lommel_s2,
) where
import Util
\end{code}

--TODO: These are completely untested!

\subsection{First Lommel function}
For $\mu\pm\nu \neq \pm1, \pm3, \pm5, \cdots$ we define
the first Lommel function $\verb|sf_lommel_s mu nu z| = S_{\mu,\nu}(z)$ 
via series-expansion:
\[ S_{\mu,\nu}(z) = \frac{z^{mu+1}}{(\mu+1)^2-\nu^2} \sum_{k=0}^\infty t_k \]
where
\[ t_0=1 \qquad t_{k}=t_{k-1}\frac{-z^2}{(\mu+2k+1)^2-\nu^2} \]

\subsubsection{\tt sf\_lommel\_s mu nu z}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_lommel\_s mu nu z} = S_{\mu,\nu}(z)$}
\begin{code}
sf_lommel_s mu nu z = 
  let terms = ixiter 1 1.0 $ \ k t -> -t*z^2 / ((mu+((#)$2*k+1))^2 - nu^2)
      res = ksum terms
  in res * z**(mu+1) / ((mu+1)^2 - nu^2)
\end{code}
\end{titled-frame}

\subsection{Second Lommel function}
For $\mu\pm\nu \neq \pm1, \pm3, \pm5, \cdots$ 
the second Lommel function $\verb|sf_lommel_s2 mu nu z| = s_{\mu,\nu}(z)$ 
is given via an asymptotic expansion:
\[ s_{\mu,\nu}(z) \sim \sum_{k=0}^\infty u_k \]
where
\[ u_0=1 \qquad u_{k}=u_{k-1}\frac{-(\mu-2k+1)^2-\nu^2}{z^2} \]
\subsubsection{\tt sf\_lommel\_s2 mu nu z}
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


