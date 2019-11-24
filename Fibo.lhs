\section{Fibonacci Numbers}

A silly approach to efficient computation of Fibonacci numbers
\[ f_n = f_{n-1} + f_{n-2} \qquad f_0=0 \qquad f_1=1 \]

The idea is to use the closed-form solution:
\[ f_n = \frac{1}{\sqrt5}\left(\frac{1+\sqrt5}{2}\right)^n + \frac{-1}{\sqrt5}\left(\frac{1-\sqrt5}{2}\right)^n \]
and note that we can work in $\mathbb{Q}[\sqrt5]$ with terms of the form $a+b\sqrt5$ with $a,b\in\mathbb{Q}$
(notice that $\frac{1}{\sqrt5}=\frac{\sqrt5}{5}$.)
\begin{eqnarray*}\
(a+b\sqrt5) + (c+d\sqrt5) &=& (a+c)+(b+d)\sqrt5 \\
(a+b\sqrt5) * (c+d\sqrt5) &=& (ac+5bd)+(ad+bc)\sqrt5
\end{eqnarray*}

We use the \verb|Rational| type to represent elements of $\mathbb{Q}$, which is a bit more than we actually need,
as in the computations above the denominator of $\left(\frac{1\pm\sqrt5}{2}\right)^n$ is always, in fact,~1~or~2.
\begin{code}
module Fibo (fibonacci) where
import Data.Ratio

data Q5 = Q5 Rational Rational
  deriving (Eq)
\end{code}
The number-theoretic norm, unused for our application.
\begin{code}
norm (Q5 ra qa) = ra^2-5*qa^2
\end{code}

Human-friendly \verb|Show| instantiation.
\begin{code}
instance Show Q5 where
  show (Q5 ra qa) = (show ra)++"+"++(show qa)++"*sqrt(5)"
\end{code}

Implementation of the operations.
The \verb|abs| and \verb|signum| functions are unused, so we just give placeholder values.
\begin{code}
instance Num Q5 where
  (Q5 ra qa) + (Q5 rb qb) = Q5 (ra+rb) (qa+qb)
  (Q5 ra qa) - (Q5 rb qb) = Q5 (ra-rb) (qa-qb)
  (Q5 ra qa) * (Q5 rb qb) = Q5 (ra*rb+5*qa*qb) (ra*qb+rb*qa)
  negate (Q5 ra qa) = Q5 (-ra) (-qa)
  abs a = Q5 (norm a) 0
  signum a@(Q5 ra qa) = if a==0 then 0 else Q5 (ra/(norm a)) (qa/(norm a))
  fromInteger n = Q5 (fromInteger n) 0

instance Fractional Q5 where
  recip a@(Q5 ra qa) = Q5 (ra/(norm a)) (-qa/(norm a))
  fromRational r = (Q5 r 0)
\end{code}

Finally, we define $\phi_{\pm}=\frac12(1\pm\sqrt5)$ and $c_{\pm}=\pm\frac15\sqrt5$ so that
$f_n = c_+\phi_+^n+c_-\phi_-^n$.
\begin{code}
phip = Q5 (1%2) (1%2)
cp   = Q5 0     (1%5)
phim = Q5 (1%2) (-1%2)
cm   = Q5 0     (-1%5)
fibonacci' n = let (Q5 r q) = cp*phip^n + cm*phim^n in numerator r
fibonacci n = let (Q5 _ q) = phip^^n in numerator (2*q)
\end{code}

