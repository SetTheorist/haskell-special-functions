\section{Utility}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Preamble}
We start with the basic preamble.
\begin{titled-frame}{\color{blue}\tt module Util}
\begin{code}
{-# Language BangPatterns #-}
{-# Language FlexibleContexts #-}
{-# Language FlexibleInstances #-}
{-# Language ScopedTypeVariables #-}
{-# Language TypeFamilies #-}
-- {-# Language UndecidableSuperClasses #-}
-- {-# Language UndecidableInstances #-}
module Util where
import Data.Complex
import Data.List(zipWith5)
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Data Types}

We start by defining a convenient type synonym for complex numbers over \verb|Double|.
\begin{code}
type CDouble = Complex Double
\end{code}

Next, we define the \verb|Value| typeclass which is useful for defining
our special functions to work over both real (\verb|Double|) values and over
complex (\verb|CDouble|) values with uniform implementations.
This will also make it convenient for handling \verb|Quad| values (later).
\begin{titled-frame}{{\color{blue}\tt class Value v}\marginnote{\tt Value}}
\begin{code}
class (Eq v, Floating v, Fractional v, Num v,
       Enum (RealKind v), Eq (RealKind v), Floating (RealKind v),
         Fractional (RealKind v), Num (RealKind v), Ord (RealKind v),
         RealFrac (RealKind v),
       Eq (ComplexKind v), Floating (ComplexKind v), Fractional (ComplexKind v),
         Num (ComplexKind v)
      ) => Value v where
  type RealKind v :: *
  type ComplexKind v :: *
  pos_infty :: v
  neg_infty :: v
  nan :: v
  re :: v -> (RealKind v)
  im :: v -> (RealKind v)
  rabs :: v -> (RealKind v)
  is_inf :: v -> Bool
  is_nan :: v -> Bool
  is_real :: v -> Bool
  fromDouble :: Double -> v
  fromReal :: (RealKind v) -> v
  toComplex :: v -> (ComplexKind v)
\end{code}
\end{titled-frame}
Both \verb|Double| and \verb|CDouble| are instances of the \verb|Value| typeclass
in the obvious ways.
\begin{titled-frame}{{\color{blue}\tt instance Value Double}\marginnote{\tt Value Double}}
\begin{code}
instance Value Double where
  type RealKind Double = Double
  type ComplexKind Double = CDouble
  pos_infty = 1.0/0.0
  neg_infty = -1.0/0.0
  nan = 0.0/0.0
  re = id
  im = const 0
  rabs = abs
  is_inf = isInfinite
  is_nan = isNaN
  is_real _ = True
  fromDouble = id
  fromReal = id
  toComplex x = x :+ 0
\end{code}
\end{titled-frame}
\begin{titled-frame}{{\color{blue}\tt instance Value CDouble}\marginnote{\tt Value CDouble}}
\begin{code}
instance Value CDouble where
  type RealKind CDouble = Double
  type ComplexKind CDouble = CDouble
  pos_infty = (1.0/0.0) :+ 0
  neg_infty = (-1.0/0.0) :+ 0
  nan = (0.0/0.0) :+ 0
  re = realPart
  im = imagPart
  rabs = realPart.abs
  is_inf z = (is_inf.re$z) || (is_inf.im$z)
  is_nan z = (is_nan.re$z) || (is_nan.im$z)
  is_real _ = False
  fromDouble x = x :+ 0
  fromReal x = x :+ 0
  toComplex = id
\end{code}
\end{titled-frame}

TODO: add quad versions also


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Helper functions}

A convenient shortcut, as we often find ourselves converting
indices (or other integral values) to our computation type.
\begin{code}
{-# INLINE (#) #-}
(#) :: (Integral a, Num b) => a -> b
(#) = fromIntegral
\end{code}

A version of \verb|iterate| which passes along an index also
(very useful for computing terms of a power-series, for example.)
\begin{titled-frame}{{\color{blue}\tt ixiter i x f}\marginnote{\tt ixiter}}
\begin{code}
{-# INLINE ixiter #-}
ixiter :: (Enum ix) => ix -> a -> (ix->a->a) -> [a]
ixiter i x f = x:(ixiter (succ i) (f i x) f)
\end{code}
\end{titled-frame}

Computes the relative error in terms of decimal digits, handy for testing.
Note that this fails when the exact value is zero.
\[ \verb|relerr e a| = \log_{10}\left|\frac{a-e}{e}\right|\]
\begin{code}
relerr :: forall v.(Value v) => v -> v -> (RealKind v)
relerr !exact !approx = re $! logBase 10 (abs ((approx-exact)/exact))
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Kahan summation}

A useful tool is so-called Kahan summation, based on the observation that
in floating-point arithmetic, one can \dots

Here \verb|kadd t s e k| is a single step of addition, adding
a term to a sum+error and passing the updated sum+error to the continuation.
\begin{code}
-- kadd value oldsum olderr ---> newsum newerr
{-# INLINE kadd #-}
{-# SPECIALISE kadd :: Double -> Double -> Double -> (Double -> Double -> a) -> a #-}
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
(TODO: make generic over stopping condition)
\begin{titled-frame}{{\color{blue}\tt ksum terms}\marginnote{\tt ksum}}
\begin{code}
{-# SPECIALISE ksum :: [Double] -> Double #-}
{-# SPECIALISE ksum' :: [Double] -> (Double -> Double -> a) -> a #-}
ksum :: (Value v) => [v] -> v
ksum terms = ksum' terms const

ksum' :: (Value v) => [v] -> (v -> v -> a) -> a
ksum' terms k = f 0 0 terms
  where
    f !s !e [] = k s e
    f !s !e (t:terms) =
      let !y  = t - e
          !s' = s + y
          !e' = (s' - s) - y
      in if s' == s
         then k s' e'
         else f s' e' terms
\end{code}
\end{titled-frame}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Continued fraction evaluation}

Given two sequences $\{a_n\}_{n=1}^\infty$ and $\{b_n\}_{n=0}^\infty$ we have the continued fraction
\[ b_0 + \frac{a_1}{b_1 + \frac{a_2}{b_2 + \frac{a_3}{b_3 + \frac{a_4}{b_4 + \cdots}}}} \]
or 
\[ b_0 + a_1/(b_1 + a_2/(b_2 + a_3/(b_3 + a_4/(b_4 + \cdots)))) \]
though for typesetting purposes this is often written
\[ b_0 + \frac{a_1}{b_1 + {}}\ \frac{a_2}{b_2 + {}}\ \frac{a_3}{b_3 + {}}\ \frac{a_4}{b_4 + \cdots} \]

We conventionally notate the $n$'th approximant or convergent as
\[ C_n = b_0 + \frac{a_1}{b_1 + {}}\ \frac{a_2}{b_2 + {}}\ \frac{a_3}{b_3 + \dots}\ \frac{a_n}{b_n} \]

\subsubsection{Backwards recurrence algorithm}

We can compute the $n$'th convergent $C_n$ for a predetermined $n$ by evaluating
\[ u_k = b_k + \frac{a_{k+1}}{u_{k+1}} \]
for $k=n-1, n-2, \dots, 0$, with $u_n = b_n$.  Then $u_0 = C_n$.

\begin{titled-frame}{\color{blue}\tt sf\_cf\_back}
\begin{code}
sf_cf_back :: forall v.(Value v) => Int -> [v] -> [v] -> v
sf_cf_back !n !as !bs =
  let !an = reverse $ take n as
      !(un:bn) = reverse $ take (n+1) bs
  in go un an bn
  where 
    go :: v -> [v] -> [v] -> v
    go !ukp1 ![] ![] = ukp1
    go !ukp1 !(a:an) !(b:bn) =
        let uk = b + a/ukp1
        in go uk an bn
\end{code}
\end{titled-frame}

\subsubsection{Steed's algorithm}
This is Steed's algorithm for evaluation of a continued fraction
It evaluates the partial convergents $C_n$ in a forward direction.
This implementation will evaluate until $C_n=C_{n+1}$.
TODO: describe algorithm.
\begin{titled-frame}{\color{blue}\tt sf\_cf\_steeds}
\begin{code}
sf_cf_steeds :: (Value v) => [v] -> [v] -> v
sf_cf_steeds (a1:as) (b0:b1:bs) =
    let !c0 = b0
        !d1 = 1/b1
        !delc1 = a1*d1
        !c1 = c0 + delc1
    in recur c1 delc1 d1 as bs
    where
      !eps = 5e-16
      recur !cn' !delcn' !dn' !(an:as) !(bn:bs) = 
        let !dn = 1/(dn'*an+bn)
            !delcn = (bn*dn - 1)*delcn'
            !cn = cn' + delcn
        in if cn == cn' || (rabs delcn)<eps || is_nan cn
           then cn
           else (recur cn delcn dn as bs)
\end{code}
\end{titled-frame}

\subsubsection{Modified Lentz algorithm}
An alternative algorithm for evaluating a continued fraction in a forward directions.
This algorithm can be less susceptible to contamination from rounding errors.
TODO: describe algorithm
\begin{titled-frame}{$\text{\color{blue}\tt sf\_cf\_lentz}$}
\begin{code}
sf_cf_lentz :: (Value v) => [v] -> [v] -> v
sf_cf_lentz as (b0:bs) =
  let !c0 = nz b0
      !e0 = c0
      !d0 = 0
  in iter c0 d0 e0 as bs
  where
    !eps = 5e-16
    !zeta = 1e-100
    nz !x = if x==0 then zeta else x
    iter cn dn en (an:as) (bn:bs) = 
      let !idn = nz $ bn + an*dn
          !en' = nz $ bn + an/en
          !dn' = 1 / idn
          !hn  = en' * dn'
          !cn' = cn * hn
          !delta = rabs(hn - 1)
      in if cn==cn' || delta<eps || is_nan cn'
         then cn
         else iter cn' dn' en' as bs
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Solving ODEs}

\subsubsection{Runge-Kutta IV}
Solve a system of first-order ODEs using the Runge-Kutta IV method.
To solve ${\bf y}' = {\bf f}(t,{\bf y})$ from $t=t_0$ to $t=t_n$ with initial condition ${\bf y}(t_0)={\bf y}_0$,
first choose a step-size $h>0$.  Then iteratively proceed by letting
\begin{eqnarray*}
  {\bf k}_1 &=& h {\bf f}(t_i, {\bf y}_i) \\
  {\bf k}_2 &=& h {\bf f}(t_i+\frac{h}{2}, {\bf y}_i+\frac{1}{2}{\bf k}_1) \\
  {\bf k}_3 &=& h {\bf f}(t_i+\frac{h}{2}, {\bf y}_i+\frac{1}{2}{\bf k}_2) \\
  {\bf k}_4 &=& h {\bf f}(t_i+h, {\bf y}_i+{\bf k}_3)
\end{eqnarray*}
and then
\begin{eqnarray*}
  t_{i+1} &=& t_i + h \\
  {\bf y}_{i+1} &=& {\bf y}_{i} + \frac16({\bf k}_1 + 2{\bf k}_2 + 2{\bf k}_3 + {\bf k}_4)
\end{eqnarray*}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_runge\_kutta\_4}$}
\begin{code}
sf_runge_kutta_4 :: forall v.(Value v) =>
    (RealKind v) -> (RealKind v) -> (RealKind v) -> [v] -> ((RealKind v)->[v]->[v])
      -> [(RealKind v,[v])]
sf_runge_kutta_4 !h !t0 !tn !x0 !f = iter t0 x0 [(t0,x0)]
  where
    iter :: (RealKind v) -> [v] -> [(RealKind v,[v])] -> [(RealKind v,[v])]
    iter !ti !xi !path
      | ti>=tn    = path
      | otherwise =
          let !h'  = (min h (tn-ti))
              !h'2 = h'/2
              !h'' = fromReal h'
              !k1  = fmap (h''*) (f ti xi)
              !k2  = fmap (h''*) (f (ti+h'2) (zipWith (\x k->x+k/2) xi k1))
              !k3  = fmap (h''*) (f (ti+h'2) (zipWith (\x k->x+k/2) xi k2))
              !k4  = fmap (h''*) (f (ti+h' ) (zipWith (\x k->x+k  ) xi k3))
              !ti1 = ti + h'
              !xi1 = zipWith5 (\x k1 k2 k3 k4 -> x + (k1+2*k2+2*k3+k4)/6) xi k1 k2 k3 k4
          in iter ti1 xi1 ((ti1,xi1):path)
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{TO BE MOVED}
\begin{code}
sf_sqrt :: (Value v) => v -> v
sf_sqrt = sqrt
\end{code}


