\section{Utility}

\subsection{Preamble}
We start with the basic preamble.
\begin{code}
{-# Language BangPatterns #-}
{-# Language FlexibleContexts #-}
{-# Language FlexibleInstances #-}
{-# Language TypeFamilies #-}
-- {-# Language UndecidableInstances #-}
module Util where
import Data.Complex
\end{code}

\subsection{Data Types}

We start by defining a convenient type synonym for complex numbers over \verb|Double|.
\begin{code}
type CDouble = Complex Double
\end{code}

Next, we define the \verb|Value| typeclass which is useful for defining
our special functions to work over both real (\verb|Double|) values and over
complex (\verb|CDouble|) values with uniform implementations.
This will also make it convenient for handling \verb|Quad| values (later).
\begin{code}
class (Eq t, Floating t, Fractional t, Num t,
       Enum (RealKind t), Eq (RealKind t), Floating (RealKind t),
         Fractional (RealKind t), Num (RealKind t), Ord (RealKind t),
       Eq (ComplexKind t), Floating (ComplexKind t), Fractional (ComplexKind t),
         Num (ComplexKind t)
      ) => Value t where
  type RealKind t :: *
  type ComplexKind t :: *
  re :: t -> (RealKind t)
  im :: t -> (RealKind t)
  rabs :: t -> (RealKind t)
  is_inf :: t -> Bool
  is_nan :: t -> Bool
  fromDouble :: Double -> t
  fromReal :: (RealKind t) -> t
  toComplex :: t -> (ComplexKind t)
\end{code}
Both \verb|Double| and \verb|CDouble| are instances of the \verb|Value| typeclass
in the obvious ways.
\begin{code}
instance Value Double where
  type RealKind Double = Double
  type ComplexKind Double = CDouble
  re = id
  im = const 0
  rabs = abs
  is_inf = isInfinite
  is_nan = isNaN
  fromDouble = id
  fromReal = id
  toComplex x = x :+ 0

instance Value CDouble where
  type RealKind CDouble = Double
  type ComplexKind CDouble = CDouble
  re = realPart
  im = imagPart
  rabs = realPart.abs
  is_inf z = (is_inf.re$z) || (is_inf.im$z)
  is_nan z = (is_nan.re$z) || (is_nan.im$z)
  fromDouble x = x :+ 0
  fromReal x = x :+ 0
  toComplex = id
\end{code}

TODO: add quad versions also


\subsection{Helper functions}

A convenient shortcut, as we often find ourselves converting
indices (or other integral values) to our computation type.
\begin{code}
(#) :: (Integral a, Num b) => a -> b
(#) = fromIntegral
\end{code}

A version of \verb|iterate| which passes along an index also
(very useful for computing terms of a power-series, for example.)
\begin{code}
ixiter :: (Enum ix) => ix -> a -> (ix->a->a) -> [a]
ixiter i x f = x:(ixiter (succ i) (f i x) f)
\end{code}

Computes the relative error in terms of decimal digits, handy for testing.
Note that this fails when the exact value is zero.
\[ \verb|relerr e a| = \log_{10}\left|\frac{a-e}{e}\right|\]
\begin{code}
relerr :: (Value v) => v -> v -> (RealKind v)
relerr !exact !approx = re $! logBase 10 (abs ((approx-exact)/exact))
\end{code}

\subsection{Kahan summation}

A useful tool is so-called Kahan summation, based on the observation that
in floating-point arithmetic, one can \dots

Here \verb|kadd t s e k| is a single step of addition, adding
a term to a sum+error and passing the updated sum+error to the continuation.
\begin{code}
-- kadd value oldsum olderr ---> newsum newerr
kadd :: (Value v) => v -> v -> v -> (v -> v -> a) -> a
kadd t s e k =
  let y = t - e
      s' = s + y
      e' = (s' - s) - y
  in k s' e'
\end{code}

Here \verb|ksum terms| sums a list with Kahan summation.
The list is assumed to be (eventually) decreasing and the
summation is terminated as soon as adding a term doesn't change the value.
(Thus any zeros in the list will immediately terminate the sum.)
This is typically used for power-series or asymptotic expansions.
\begin{code}
ksum :: (Value v) => [v] -> v
ksum terms = k 0 0 terms
  where
    k !s !e [] = s
    k !s !e (t:terms) =
      let !y = t - e
          !s' = s + y
          !e' = (s' - s) - y
      in if s' == s
         then s
         else k s' e' terms
\end{code}


\subsection{Continued fraction evaluation}
This is Steed's algorithm for evaluation of a continued fraction
\[ C = b_0 + a_1/(b_1 + a_2/(b_2 + a_3/(b_3 + \cdots))) \]
where $C_n=A_n/B_n$ is the partial evaluation up to $\dots a_n/b_n$.
Here \verb|steeds as bs| evaluates until $C_n=C_{n+1}$.
TODO: describe the algorithm.
\begin{code}
steeds :: (Value v) => [v] -> [v] -> v
steeds (a1:as) (b0:b1:bs) =
    let !c0 = b0
        !d1 = 1/b1
        !delc1 = a1*d1
        !c1 = c0 + delc1
    in recur c1 delc1 d1 as bs
    where recur !cn_1 !delcn_1 !dn_1 !(an:as) !(bn:bs) = 
            let !dn = 1/(dn_1*an+bn)
                !delcn = (bn*dn - 1)*delcn_1
                !cn = cn_1 + delcn
            in if (cn == cn_1) || is_nan cn then cn else (recur cn delcn dn as bs)
\end{code}

\subsection{TO BE MOVED}
\begin{code}
sf_sqrt :: (Value v) => v -> v
sf_sqrt = sqrt
\end{code}

