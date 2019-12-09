\chapter{Orthogonal Polynomials}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Preamble}

\begin{code}
{-# Language BangPatterns #-}
{-# Language ScopedTypeVariables #-}
module OrthoPoly where
import qualified Polynomial as P
import Trig
import Util
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Typeclass}
\begin{code}
-- TODO: make the general (Value, etc.)
class OrthogonalPolynomial a where
  support :: a -> (Double,Double) -- TODO: inclusive/exclusive?!
  -- in_support :: a -> Double -> Bool
  weight :: a -> Double -> Double

  coeffs :: a -> Int -> [Double]
  scale :: a -> Int -> Double
  value :: a -> Int -> Double -> Double
  weights :: a -> Int -> [Double]
  zeros :: a -> Int -> [Double]

  value_arr :: a -> Int -> Double -> [Double]
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Chebyshev Polynomials of the First Kind}

\begin{code}
data ChebyshevTPolynomial = ChebyshevTPolynomial

instance OrthogonalPolynomial ChebyshevTPolynomial where
  support _ = (-1,1) 
  -- in_support _ x = -1<x && x<1
  weight _ x = if -1<x && x<1 then 1/(sf_sqrt(1-x^2))

  coeffs _ = sf_orthopoly_chebyshev_t_coeffs
  scale _ = sf_orthopoly_chebyshev_t_scale
  value _ = sf_orthopoly_chebyshev_t_value
  weights _ = sf_orthopoly_chebyshev_t_weights
  zeros _ = sf_orthopoly_chebyshev_t_zeros

  value_arr _ = sf_orthopoly_chebyshev_t_values
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}
Compute the coefficients of the $n$'th Chebyshev polynomial of the first kind: $T_n(z) = a_1 + a_2*x + \cdots + a_{n+1}*x^n$,
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_coeffs !n
  | n==0 = [1]
  | n==1 = [0,1]
  | otherwise = P.coeffs $ iter 2 1 P.px
  where
    iter !k !rm2 !rm1 =
      let !rm0 = 2 P.<* (P.px*rm1) - rm2
      in if k==n then rm0
         else iter (k+1) rm1 rm0
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th Chebyshev polynomial of the first kind, $T_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_scale !n
  | n==0      = sf_sqrt (1/pi)
  | otherwise = sf_sqrt (2/pi)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Value}
Compute the value of the $n$'th Chebyshev polynomial at $z$ of the first kind: $T_n(z)$,
$n=0, 1, 2, \dots$, typically $z\in[-1,1]$
\begin{code}
sf_orthopoly_chebyshev_t_value_accum !n !z !k1 !k0
  | n==0 = k1 1 k0
  | n==1 = k1 z $ k1 1 k0
  | otherwise = iter 2 1 z (k1 z $ k1 1 k0)
  where
    iter !k !rm2 !rm1 !acc = 
      let !rm0 = 2*z*rm1 - rm2
      in if k==n then (k1 rm0 acc)
         else iter (k+1) rm1 rm0 (k1 rm0 acc)
\end{code}

Returns values in reverse order: $[T_n(z),T_{n-1}(z),\dots,T_0(z)]$
\begin{code}
sf_orthopoly_chebyshev_t_value  !n !z = sf_orthopoly_chebyshev_t_value_accum n z const 0
sf_orthopoly_chebyshev_t_values !n !z = sf_orthopoly_chebyshev_t_value_accum n z (:)   []
\end{code}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gauss quadrature weights for the $n$'th Chebyshev polynomial of the first kind: $T_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_weights !n
  | n==0      = []
  | otherwise = replicate n (pi/(#)n)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Zeros}
Compute the zeros of the $n$'th Chebyshev polynomial of the first kind: $T_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_zeros !n
  | n==0      = []
  | otherwise =
      let !ths = [ ((#)$2*k-1)/(#)n | k<-[n,(n-1)..((n+3)`div`2)] ]
          !cos = map (\t->sf_cos $ (pi/2)*t) ths
      in cos ++ (if n`mod`2==0 then [] else [0]) ++ (map negate $ reverse cos)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Legendre Polynomials}

Legendre polynomials $P_n(x)$ are orthogonal on the interval $(-1,1)$ 
for the weight function $w(x)=1$.

They satisfy the following recurrence, with $P_0(x)=1$, $P_1(x)=x$,
\[ P_n(x) = \frac{(2n-1) x P_{n-1}(x) - (n-1)P_{n-2}}{n} \]

\begin{code}
data LegendrePolynomial = LegendrePolynomial

instance OrthogonalPolynomial LegendrePolynomial where
  support _ = (-1,1)
  -- in_support _ x = -1<x && x<1
  weight l x = if -1<x && x<1 then 1 else 0

  coeffs _ = sf_orthopoly_legendre_coeffs
  scale _ = sf_orthopoly_legendre_scale
  value _ = sf_orthopoly_legendre_value
  weights = undefined
  zeros = undefined

  value_arr _ = sf_orthopoly_legendre_values
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}

Compute the coefficients of the $n$'th Legendre polynomial:
$P_n(z) = a_1 + a_2*x + \cdots + a_{n+1}*x^n$,
$n=0, 1, 2, \dots$
(TODO: cache?)
\begin{code}
sf_orthopoly_legendre_coeffs :: forall v.(Eq v,Num v,Fractional v) => Int -> [v]
sf_orthopoly_legendre_coeffs !n
  | n == 0    = [1]
  | n == 1    = [0,1]
  | otherwise = P.coeffs $ iter 2 1 P.px
  where
    iter :: Int -> P.Poly v -> P.Poly v -> P.Poly v
    iter !k !rm2 !rm1 =
      let rm0 = ((((#)$2*k-1)P.<*(P.px*rm1)) - (((#)$k-1)P.<*rm2)) P./> (#)k
      in if k==n then rm0
         else iter (k+1) rm1 rm0
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th Legendre polynomial:
$P_n(z) = a_1 + a_2*x + \cdots + a_{n+1}*x^n$,
$n=0, 1, 2, \dots$

\begin{code}
--TODO: validate non-negative integer input
sf_orthopoly_legendre_scale :: (Integral a, Value v) => a -> v
sf_orthopoly_legendre_scale !n = sf_sqrt $ (2*(#)n+1)/2
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Value}
Compute the value of the $n$'th Legendre polynomial: $P_n(z)$,
$n=0, 1, 2, \dots$, typically $z\in[-1,1]$
Uses recursion (which is more robust than direct evaluation of polynomial.)

\begin{code}
sf_orthopoly_legendre_value_accum !n !z !k1 !k0
  | n==0      = k1 1 k0
  | n==1      = k1 z $ k1 1 k0
  | otherwise = iter 1 z 2 (k1 z $ k1 1 k0)
  where
    iter !rm2 !rm1 !k !acc =
      let !rm0 = (((#)$2*k-1)*z*rm1 - ((#)k-1)*rm2)/(#)k
      in if k==n then (k1 rm0 acc)
         else iter rm1 rm0 (k+1) (k1 rm0 acc)
\end{code}

Returns values in reverse order: $[P_n(z),P_{n-1}(z),\dots,P_0(z)]$
\begin{code}
sf_orthopoly_legendre_value  !n !z = sf_orthopoly_legendre_value_accum n z const 0
sf_orthopoly_legendre_values !n !z = sf_orthopoly_legendre_value_accum n z (:)   []
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gaussian quadrature weights of the $n$'th Legendre polynomial: $P_n(z)$
$n=0, 1, 2, \dots$

\begin{code}
-- cache?
--if (!sf_is_nonnegint(n)) print_usage; endif
sf_orthopoly_legendre_weights n 
  | n==0 = []
  | n==1 = [1]
  | otherwise =
      let !zs = sf_orthopoly_legendre_zeros n
          --  need ortho_normal_ polynomials here
          jth j =
            let !jvs = map (sf_orthopoly_legendre_value j) zs
                !nj = sf_orthopoly_legendre_scale j
            in map (\x->(x*nj)^2) jvs
          !ires = foldl (P..+.) (replicate n 0) (map jth [0..(n-1)])
          !res = map (1/) ires
      in res
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Zeros}
Compute the zeros of the $n$'th Legendre polynomial: $P_n(z)$
$n=0, 1, 2, ...$
not good for large
TODO: domain check, cache
\begin{code}
sf_orthopoly_legendre_zeros n
  | n==0 = []
  | n==1 = [0]
  | otherwise = undefined
  {--
    m = zeros(n);
    for k=1:n-1
      m(k,k+1) = m(k+1,k) = k/sqrt(4*k^2-1);
    endfor
    res = sort(eig(m));

    # "polish" the results
    fx = sf_orthopoly_legendre_value(n, res);
    dfx = (-n*res.*fx + n*sf_orthopoly_legendre_value(n-1, res)) ./ (1-res.^2);
    nwt = res - fx./dfx;
    #nrs = sf_orthopoly_legendre_value(n,nwt)
    res = nwt;

    # "polish" the results
    fx = sf_orthopoly_legendre_value(n, res);
    dfx = (-n*res.*fx + n*sf_orthopoly_legendre_value(n-1, res)) ./ (1-res.^2);
    nwt = res - fx./dfx;
    #nrs = sf_orthopoly_legendre_value(n,nwt)
    res = nwt;
    --}
\end{code}

