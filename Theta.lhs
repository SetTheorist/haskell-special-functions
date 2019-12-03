\chapter{Jacobian Theta functions}

General notation: we assume $\Im\tau > 0$ and $0<|q|<1$ where $q = e^{\ii \pi \tau}$.

\section{Preamble}
\begin{code}
{-# Language BangPatterns #-}
module Theta where
import Exp
import Trig
import Util
\end{code}


\subsection{Theta1}
\[ \theta_1(z\mid\tau) = \theta_1(z,q) = 2\sum_{n=0}^\infty (-)^n q^{(n+\frac12)^2} \sin((2n+1)z) \]

\subsubsection{\tt sf\_theta\_1 z q}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_theta\_1 z q} = \theta_1(z,q)$}
\begin{code}
sf_theta_1 :: (Value v) => v -> v -> v
sf_theta_1 !z !q =
  let !qpows = map (\n -> q**(((#)n+1/2)^2) * (-1)^n) [0..]
      !sins  = map (\n -> sf_sin $ z*(#)(2*n+1)) [0..]
      !terms = zipWith (*) qpows sins
  in 2 * (ksum terms)
\end{code}
\end{titled-frame}

\subsection{Theta2}
\[ \theta_2(z\mid\tau) = \theta_2(z,q) = 2\sum_{n=0}^\infty (-)^n q^{(n+\frac12)^2} \cos((2n+1)z) \]

\subsubsection{\tt sf\_theta\_2 z q}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_theta\_2 z q} = \theta_2(z,q)$}
\begin{code}
sf_theta_2 :: (Value v) => v -> v -> v
sf_theta_2 !z !q =
  let !qpows = map (\n -> q**(((#)n+1/2)^2)) [0..]
      !coss  = map (\n -> sf_cos $ z*(#)(2*n+1)) [0..]
      !terms = zipWith (*) qpows coss
  in 2 * (ksum terms)
\end{code}
\end{titled-frame}

\subsection{Theta4}
\[ \theta_4(z\mid\tau) = \theta_4(z,q) = 1 + 2\sum_{n=1}^\infty (-)^n q^{n^2} \cos(2nz) \]

\subsubsection{\tt sf\_theta\_4 z q}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_theta\_4 z q} = \theta_4(z,q)$}
\begin{code}
sf_theta_4 :: (Value v) => v -> v -> v
sf_theta_4 !z !q =
  let !qpows = map (\n -> q^(n^2) * (-1)^n) [1..]
      !coss  = map (\n -> sf_cos $ z*(#)(2*n)) [1..]
      !terms = zipWith (*) qpows coss
  in 1 + 2 * (ksum terms)
\end{code}
\end{titled-frame}

\subsection{Theta3}
\[ \theta_3(z\mid\tau) = \theta_3(z,q) = 1 + 2\sum_{n=1}^\infty q^{n^2} \cos(2nz) \]

\subsubsection{\tt sf\_theta\_3 z q}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_theta\_3 z q} = \theta_3(z,q)$}
\begin{code}
sf_theta_3 :: (Value v) => v -> v -> v
sf_theta_3 !z !q =
  let !terms = map (\n-> q^(n^2) * (sf_cos $ z*(#)(2*n))) [1..]
  in 1 + 2 * (ksum terms)
\end{code}
\end{titled-frame}

TODO: this looks incorrect, need to fix
\begin{code}
sf_theta_3' :: (Value v) => v -> v -> v
sf_theta_3' !z !q =
  let !phi   = -(sf_log q)/pi
      !q'    = sf_exp(-pi/phi)
      !qpows = map (\n->q'^(n^2)) [1..]
      !ees   = map (\n -> (sf_exp$ -((#)n)^2*pi/phi + 2*(#)n*z/phi)
                * (1 + (sf_exp$ -4*(#)n*z/phi)) * 0.5) [1..]
      !terms = zipWith (*) qpows ees
      -- terms = ees
      !res   = ksum terms
  in ((sf_exp$ -z^2/(pi*phi)) + (sf_log_p1$2*res)) / (sf_sqrt phi)
\end{code}
