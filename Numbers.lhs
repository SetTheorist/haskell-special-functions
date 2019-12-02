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

\begin{code}
fibonacci_number :: Int -> Integer
fibonacci_number n = Fibo.fibonacci n
\end{code}

\begin{code}
lucas_number :: Int -> Integer
lucas_number = undefined

euler_number :: Int -> Integer
euler_number = undefined

catalan_number :: Integer -> Integer
catalan_number 0 = 1
catalan_number n = 2*(2*n-1)*(catalan_number (n-1))`div`(n+1)
\end{code}

To compute the Bernoulli numbers $B_n$, we use the relation
\[ \sum_{k=0}^{n}\binom{n+1}{k}B_k = 0 \]
we compute with rational numbers, so the result will be exact.
(Note that this is not the most efficient approach to computing
the Bernoulli numbers, but it suffices for now.)
\begin{code}
sf_bernoulli_b :: [Rational]
sf_bernoulli_b = map _bernoulli_number_computation [0..]
_bernoulli_number_computation :: Int -> Rational
_bernoulli_number_computation n
  | n == 0    = 1
  | n == 1    = -1%2
  | (odd n)   = 0
  | otherwise =
      let !n' = (#)n
          !terms = map (\k -> let k'=(#)k in ((#)$binomial (n+1) k)*(sf_bernoulli_b!!k)) [0..(n-1)]
      in -(sum terms)/(n'+1)
\end{code}

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
