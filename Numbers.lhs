\chapter{Numbers}

\section{Preamble}
\begin{titled-frame}{\color{blue}\tt module Numbers}
\begin{code}
{-# Language BangPatterns #-}
module Numbers where
import Data.Ratio
import qualified Fibo
import Util
\end{code}
\end{titled-frame}


\section{misc}
\begin{code}
fibonacci_number :: Int -> Integer
fibonacci_number n = Fibo.fibonacci n
\end{code}

\begin{code}
lucas_number :: Int -> Integer
lucas_number = undefined

catalan_number :: Integer -> Integer
catalan_number 0 = 1
catalan_number n = 2*(2*n-1)*(catalan_number (n-1))`div`(n+1)
\end{code}

\section{Bernoulli numbers}

The Bernoulli numbers, $B_n$, are defined via their exponential generating function
\[ \frac{t}{e^t-1} = \sum_{n=1}^\infty B_n \frac{t^n}{n!} \]

\subsection{\tt sf\_bernoulli\_b}
To compute the Bernoulli numbers $B_n$, we use the relation
\[ \sum_{k=0}^{n}\binom{n+1}{k}B_k = 0 \]
we compute with rational numbers, so the result will be exact.
(Note that this is not the most efficient approach to computing
the Bernoulli numbers, but it suffices for now.)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_bernoulli\_b !! n} = B_n$}
\begin{code}
sf_bernoulli_b :: [Rational]
sf_bernoulli_b = map _bernoulli_number_computation [0..]
_bernoulli_number_computation :: Int -> Rational
_bernoulli_number_computation n
  | n == 0    = 1
  | n == 1    = -1%2
  | (odd n)   = 0
  | otherwise =
      let !terms = map (\k -> ((#)$binomial (n+1) k)*(sf_bernoulli_b!!k)) [k|k<-(0:1:[2..(n-1)]),k<=n-1]
      in -(sum terms)/((#)n+1)
\end{code}
\end{titled-frame}

\subsection{\tt sf\_bernoulli\_b\_scaled}
To compute the scaled Bernoulli numbers $\widetilde{B}_n = \frac{B_n}{n!}$, we simply divide
the (unscaled) Bernoulli number by~$n!$.
Again, this is not the most efficient approach, but it suffices for now.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_bernoulli\_b\_scaled !! n} = \widetilde{B}_n = B_n/n!$}
\begin{code}
sf_bernoulli_b_scaled :: [Rational]
sf_bernoulli_b_scaled = zipWith (/) sf_bernoulli_b (map (fromIntegral.factorial) [0..])
\end{code}
\end{titled-frame}

\section{Euler numbers}

The Euler numbers, $E_n$, are defined via their exponential generating function
\[ \frac{2t}{e^{2t}-1} = \sum_{n=1}^\infty E_n \frac{t^n}{n!} \]

\subsection{\tt sf\_euler\_e}
To compute the Euler numbers $E_n$, we use the relation
\[ \sum_{k=0}^{n}\binom{2n}{2k}E_{2k} = 0 \]
as Euler numbers are all integers, we compute with \verb|Integer| type
to get exact results.
(Note that this is not the most efficient approach to computing
the Euler numbers, but it suffices for now.)
\begin{titled-frame}{$\text{\color{blue}\tt sf\_euler\_e !! n} = E_n$}
\begin{code}
sf_euler_e :: [Integer]
sf_euler_e = map _euler_number_computation [0..]
_euler_number_computation :: Int -> Integer
_euler_number_computation n
  | n == 0    = 1
  | (odd n)   = 0
  | otherwise =
      let !n' = n`div`2
          !terms = map (\k -> ((#)$binomial (2*n') (2*k))*(sf_euler_e!!(2*k))) [0..(n'-1)]
      in -(sum terms)
\end{code}
\end{titled-frame}

\subsection{\tt sf\_euler\_e\_scaled}
To compute the scaled Euler numbers $\widetilde{E}_n = \frac{E_n}{n!}$, we simply divide
the (unscaled) Euler number by~$n!$.
Again, this is not the most efficient approach, but it suffices for now.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_euler\_e\_scaled !! n} = \widetilde{E}_n = E_n/n!$}
\begin{code}
sf_euler_e_scaled :: [Rational]
sf_euler_e_scaled = zipWith (\a b->(#)a/(#)b) sf_euler_e (map factorial [0..])
\end{code}
\end{titled-frame}

\section{misc}

\begin{code}
tangent_number :: Int -> Integer
tangent_number = undefined

triangular_number :: Integer -> Integer
triangular_number n = n*(n+1)`div`2

factorial :: (Integral a) => a -> a
factorial 0 = 1
factorial 1 = 1
factorial n = product [1..n]

binomial :: (Integral a) => a -> a -> a
binomial n k
    | k<0 = 0
    | n<0 = 0
    | k>n = 0
    | k==0 = 1
    | k==n = 1
    | k>n`div`2 = binomial n (n-k)
    | otherwise = (product [n-(k-1)..n]) `div` (product [1..k])
\end{code}

\section{Stirling numbers}
\begin{code}
-- TODO: this is extremely inefficient approach
stirling_number_first_kind n k = s n k
  where s n k | k<=0 || n<=0 = 0
        s n 1 = (-1)^(n-1)*(factorial (n-1))
        s n k = (s (n-1) (k-1)) - (n-1)*(s (n-1) k)

-- TODO: this is extremely inefficient approach
stirling_number_second_kind n k = s n k
  where s n k | k<=0 || n<=0 = 0
        s n 1 = 1
        s n k = k*(s (n-1) k) + (s (n-1) (k-1))
\end{code}
