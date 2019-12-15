\chapter{Orthogonal Polynomials}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Preamble}

\begin{code}
{-# Language BangPatterns #-}
{-# Language ScopedTypeVariables #-}
module OrthoPoly where
import Exp
import Gamma
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Utility}

Compute the polynomial $p_n(x)$ satisfying the second-order recurrence
\[ p_{n+1}(x) = (A_n x + B_n)p_n(x) - C_n p_{n-1}(x) \]
with $p_0(x)$ and $p_1(x)$ given.
\begin{code}
{-# INLINABLE recur_coeffs #-}
recur_coeffs :: P.Poly Double -> P.Poly Double -> (Int->Double) -> (Int->Double) -> (Int->Double) -> Int -> P.Poly Double
recur_coeffs p0 p1 an bn cn n
  | n == 0    = p0
  | n == 1    = p1
  | otherwise = iter 2 p0 p1
  where
    iter :: Int -> P.Poly Double -> P.Poly Double -> P.Poly Double
    iter !k !rm2 !rm1 =
      let !ak = an (k-1)
          !bk = bn (k-1)
          !ck = cn (k-1)
          !rm0 = (P.Poly [bk,ak])*rm1 - (ck P.<* rm2)
      in if k==n then rm0
         else iter (k+1) rm1 rm0
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Chebyshev Polynomials of the First Kind}

\begin{code}
data ChebyshevTPolynomial = ChebyshevTPolynomial

instance OrthogonalPolynomial ChebyshevTPolynomial where
  support _ = (-1,1) 
  -- in_support _ x = -1<x && x<1
  weight _ x  = if -1<x && x<1 then 1/(sf_sqrt(1-x^2)) else 0
  coeffs _    = sf_orthopoly_chebyshev_t_coeffs
  scale _     = sf_orthopoly_chebyshev_t_scale
  value _     = sf_orthopoly_chebyshev_t_value
  value_arr _ = sf_orthopoly_chebyshev_t_values
  weights _   = sf_orthopoly_chebyshev_t_weights
  zeros _     = sf_orthopoly_chebyshev_t_zeros
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}
Compute the coefficients of the $n$'th Chebyshev polynomial of the first kind: $T_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_coeffs :: Int -> [Double]
sf_orthopoly_chebyshev_t_coeffs !n = P.coeffs $ recur_coeffs (P.Poly[1]) (P.Poly[0,1]) (const 2) (const 0) (const 1) n
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th Chebyshev polynomial of the first kind, $T_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_scale :: Int -> Double
sf_orthopoly_chebyshev_t_scale !n
  | n==0      = sf_sqrt (1/pi)
  | otherwise = sf_sqrt (2/pi)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Value}
Compute the value of the $n$'th Chebyshev polynomial at $z$ of the first kind: $T_n(z)$,
$n=0, 1, 2, \dots$, typically $z\in[-1,1]$
\begin{code}
sf_orthopoly_chebyshev_t_value_accum !k1 !k0 !n !z
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
sf_orthopoly_chebyshev_t_value  :: Int -> Double -> Double
sf_orthopoly_chebyshev_t_value  = sf_orthopoly_chebyshev_t_value_accum const 0
sf_orthopoly_chebyshev_t_values :: Int -> Double -> [Double]
sf_orthopoly_chebyshev_t_values = sf_orthopoly_chebyshev_t_value_accum (:)   []
\end{code}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gauss quadrature weights for the $n$'th Chebyshev polynomial of the first kind: $T_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_t_weights :: Int -> [Double]
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
\section{Chebyshev Polynomials of the Second kind}

\begin{code}
data ChebyshevUPolynomial = ChebyshevUPolynomial

instance OrthogonalPolynomial ChebyshevUPolynomial where
  support _ = (-1,1) 
  -- in_support _ x = -1<x && x<1
  weight _ x  = if -1<x && x<1 then (sf_sqrt(1-x^2)) else 0
  coeffs _    = sf_orthopoly_chebyshev_u_coeffs
  scale _     = sf_orthopoly_chebyshev_u_scale
  value _     = sf_orthopoly_chebyshev_u_value
  value_arr _ = sf_orthopoly_chebyshev_u_values
  weights _   = sf_orthopoly_chebyshev_u_weights
  zeros _     = sf_orthopoly_chebyshev_u_zeros
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}
Compute the coefficients of the $n$'th Chebyshev polynomial of the second kind:
$U_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_u_coeffs :: Int -> [Double]
sf_orthopoly_chebyshev_u_coeffs !n = P.coeffs $ recur_coeffs (P.Poly[1]) (P.Poly[0,2]) (const 2) (const 0) (const 1) n
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th Chebyshev polynomial of the second kind, $U_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_u_scale :: Int -> Double
sf_orthopoly_chebyshev_u_scale _ = sf_sqrt $ 2/pi
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Value}
Compute the value of the $n$'th Chebyshev polynomial of the second kind at $z$: $U_n(z)$,
$n=0, 1, 2, \dots$, typically $z\in[-1,1]$
\begin{code}
sf_orthopoly_chebyshev_u_value_accum !k1 !k0 !n !z
  | n==0 = k1 1 k0
  | n==1 = k1 (2*z) $ k1 1 k0
  | otherwise = iter 2 1 (2*z) (k1 (2*z) $ k1 1 k0)
  where
    iter !k !rm2 !rm1 !acc = 
      let !rm0 = 2*z*rm1 - rm2
      in if k==n then (k1 rm0 acc)
         else iter (k+1) rm1 rm0 (k1 rm0 acc)
\end{code}

Returns values in reverse order: $[U_n(z),U_{n-1}(z),\dots,U_0(z)]$
\begin{code}
sf_orthopoly_chebyshev_u_value  :: Int -> Double -> Double
sf_orthopoly_chebyshev_u_value  = sf_orthopoly_chebyshev_u_value_accum const 0
sf_orthopoly_chebyshev_u_values :: Int -> Double -> [Double]
sf_orthopoly_chebyshev_u_values = sf_orthopoly_chebyshev_u_value_accum (:)   []
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gauss quadrature weights for the $n$'th Chebyshev polynomial of the second kind: $U_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_u_weights :: Int -> [Double]
sf_orthopoly_chebyshev_u_weights !n
  | n==0 = []
  | otherwise = map (\j -> (pi/(#)(n+1)) * (sf_sin $ (pi*(#)j/(#)(n+1)))^2) [1..n]
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Zeros}

Compute the zeros of the $n$'th Chebyshev polynomial of the second kind: $U_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_chebyshev_u_zeros :: Int -> [Double]
sf_orthopoly_chebyshev_u_zeros !n
  | n==0      = []
  | otherwise =
      let !ths = [ (#)k/(#)(n+1) | k<-[n,(n-1)..((n+3)`div`2)] ]
          !cos = map (\t->sf_cos $ pi*t) ths
      in cos ++ (if n`mod`2==0 then [] else [0]) ++ (map negate $ reverse cos)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Jacobi Polynomials}

\begin{code}
data JacobiPolynomial = JacobiPolynomial Double Double

instance OrthogonalPolynomial JacobiPolynomial where
  support _ = (-1,1)
  weight    (JacobiPolynomial a b) x = (1-x)**a * (1+x)**b
  coeffs    (JacobiPolynomial a b) = sf_orthopoly_jacobi_coeffs a b
  scale     (JacobiPolynomial a b) = sf_orthopoly_jacobi_scale a b
  value     (JacobiPolynomial a b) = sf_orthopoly_jacobi_value a b
  value_arr (JacobiPolynomial a b) = sf_orthopoly_jacobi_values a b
  weights   (JacobiPolynomial a b) = sf_orthopoly_jacobi_weights a b
  zeros     (JacobiPolynomial a b) = sf_orthopoly_jacobi_zeros a b
\end{code}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}

Compute the coefficients of the $n$'th Jacobi polynomial:
$J^{(\alpha,\beta)}_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
$\alpha>-1$, $\beta>-1$, $n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_jacobi_coeffs :: Double -> Double -> Int -> [Double]
sf_orthopoly_jacobi_coeffs a b n
--  if (!sf_is_nonnegint(n) || a<-1 || b<-1) print_usage; endif
  | n==0 = P.coeffs $ p0
  | n==1 = P.coeffs $ p1
  | otherwise = P.coeffs $ go 2 p0 p1
  where
    p0 = 1
    p1 = P.Poly [(a-b)/2, (2+a+b)/2]
    go !k !rm2 !rm1 = 
      let !ax = (2*#k+a+b-1)*(a^2-b^2)
          !bx = (2*#k+a+b-2)*(2*k#+a+b-1)*(2*k#+a+b)
          !cx = 2*(k#+a-1)*(k#+b-1)*(2*k#+a+b)
          !dx = 2*k#*(k#+a+b)*(2*k#+a+b-2)
          !rm0 = (ax P.<* rm1 + (bx P.<*P.px)*rm1 - cx P.<* rm2) P./> dx
      in if k==n then rm0 else go (k+1) rm1 rm0
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th Jacobi polynomial:
$J^{(\alpha,\beta)}_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
$\alpha>-1$, $\beta>-1$, $n=0, 1, 2, \dots$
\begin{code}
-- if (any(!sf_is_nonnegint(n)) || a<-1 || b<-1) print_usage; endif
sf_orthopoly_jacobi_scale a b n = sf_sqrt $
    (2*n#+a+b+1)/(2**(a+b+1))
    * (sf_gamma$(#)n+1)/(sf_gamma$n#+a+1)
    * (sf_gamma$n#+a+b+1)/(sf_gamma$n#+b+1)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Value}
Compute the value of the $n$'th Jacobi polynomial: $J^{(\alpha,\beta)}_n(z)$,
$\alpha>-1$, $\beta>-1$, $n=0, 1, 2, \dots$, typically $z\in(-1,1)$
\begin{code}
-- if (!sf_is_nonnegint(n) || a<=-1 || b<=-1) print_usage; endif
sf_orthopoly_jacobi_value_accum !k1 !k0 !a !b !n !z
  | n==0 = k1 v0 k0
  | n==1 = k1 v1 $ k1 v0 k0
  | otherwise = iter 2 v0 v1 (k1 v1 $ k1 v0 k0)
  where
    v0 = 1
    v1 = (a-b)/2 + z*(2+a+b)/2
    iter !k !rm2 !rm1 !acc =
      let !ax = (2*k#+a+b-1)*(a^2-b^2)
          !bx = (2*k#+a+b-2)*(2*k#+a+b-1)*(2*k#+a+b)
          !cx = 2*(k#+a-1)*(k#+b-1)*(2*k#+a+b)
          !dx = 2*((#)k)*(k#+a+b)*(2*k#+a+b-2)
          !rm0 = ((ax+bx*z)*rm1 - cx*rm2)/dx
      in if k==n then (k1 rm0 acc)
         else iter (k+1) rm1 rm0 (k1 rm0 acc)

sf_orthopoly_jacobi_value = sf_orthopoly_jacobi_value_accum  const 0
sf_orthopoly_jacobi_values = sf_orthopoly_jacobi_value_accum (:)   []
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gaussian quadrature weights of the $n$'th Jacobi polynomial: $P^{(\alpha,\beta)}_n(z)$
$n=0, 1, 2, \dots$, $\alpha>-1$, $\beta>-1$
\begin{code}
sf_orthopoly_jacobi_weights :: Double -> Double -> Int -> [Double]
sf_orthopoly_jacobi_weights a b n
  -- if (!sf_is_nonnegint(n) || a<=-1 || b<=-1) print_usage; endif
  | n==0 = []
  | n==1 = [1]
  | otherwise =
    let !zs = sf_orthopoly_jacobi_zeros a b n
        -- need ortho_normal_ polynomials here
        jth j =
          let !jvs = map (sf_orthopoly_jacobi_value a b j) zs
              !nj = sf_orthopoly_jacobi_scale a b j
          in map (\x->(x*nj)^2) jvs
        !ires = foldl (P..+.) (replicate n 0) (map jth [0..(n-1)])
        !res = map (1/) ires
    in res
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Zeros}
Compute the zeros of the $n$'th Jacobi polynomial: $P^{(\alpha,\beta)}_n(z)$
$n=0, 1, 2, \dots$, $\alpha>-1$, $\beta>-1$
\begin{code}
sf_orthopoly_jacobi_zeros :: Double -> Double -> Int -> [Double]
sf_orthopoly_jacobi_zeros a b n = undefined
{--
## -*- texinfo -*-
## @deftypefn {Function File} {@var{res} =} sf_orthopoly_jacobi_zeros (@var{n}, @var{a}, @var{b})
## @end deftypefn

function res = sf_orthopoly_jacobi_zeros(n, a, b)
  if (nargin < 3) print_usage; endif
  if (!sf_is_nonnegint(n) || a<=-1 || b<=-1) print_usage; endif
  if (n==0)
    res = [];
  elseif (n==1)
    res = [(b-a)/(2+a+b)];
  else
    m = zeros(n);
    for k=1:n-1
      m(k+1,k+1) = (b^2-a^2) / ((2*k+a+b)*(2*k+a+b+2));
      m(k,k+1) = m(k+1,k) = 2/(2*k+a+b) * sf_sqrt((k*(k+a)*(k+b)*(k+a+b)) / ((2*k+a+b+1)*(2*k+a+b-1)));
    endfor
    m(1,1) = (b-a)/(2+a+b);
    res = sort(eig(m));

    # "polish" the results
    fx = sf_orthopoly_jacobi_value(n,a,b,res);
    dfx = (n+a+b+1)/2 * sf_orthopoly_jacobi_value(n-1, a+1, b+1, res);
    nwt = res - fx./dfx;
    #nrs = sf_orthopoly_jacobi_value(n,a,b,nwt)
    res = nwt;
  endif
endfunction
--}
\end{code}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Laguerre Polynomials}

TODO: how to make convenient default a=0?

\begin{code}
data LaguerrePolynomial = LaguerrePolynomial Double

instance OrthogonalPolynomial LaguerrePolynomial where
  support _ = (0,1.0/0.0)
  weight    (LaguerrePolynomial a) x = if x>0 then (sf_exp(-x))*(x**a) else 0
  coeffs    (LaguerrePolynomial a) = sf_orthopoly_laguerre_coeffs a
  scale     (LaguerrePolynomial a) = sf_orthopoly_laguerre_scale a
  value     (LaguerrePolynomial a) = sf_orthopoly_laguerre_value a
  value_arr (LaguerrePolynomial a) = sf_orthopoly_laguerre_values a
  weights   (LaguerrePolynomial a) = sf_orthopoly_laguerre_weights a
  zeros     (LaguerrePolynomial a) = sf_orthopoly_laguerre_zeros a
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}
Compute the coefficients of the $n$'th (generalized) Laguerre polynomial:
$L^{\alpha}_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
(or its $k$'th derivative)
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_laguerre_coeffs_dx :: Double -> Int -> Int -> [Double]
sf_orthopoly_laguerre_coeffs_dx !a !n !k
  | k>n = [0]
  | otherwise = ((-1)^(k`mod`2)) P.*. (sf_orthopoly_laguerre_coeffs (a+#k) (n-k))
  
sf_orthopoly_laguerre_coeffs :: Double -> Int -> [Double]
sf_orthopoly_laguerre_coeffs !a !n = P.coeffs $ 
  recur_coeffs (P.Poly[1]) (P.Poly[1+a,-1])
    (\j -> -1        /#(j+1))
    (\j -> (2*j#+a+1)/#(j+1))
    (\j -> (j#+a)    /#(j+1))
    n
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th (generalized) Laguerre polynomial:
$L^{\alpha}_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_laguerre_scale :: Double -> Int -> Double
sf_orthopoly_laguerre_scale !a !n
  | a==0      = 1
  | otherwise = (sf_sqrt.(#) $ factorial n) / (sf_gamma $ n#+a+1)
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Value}
Compute the value of the $n$'th (generalized) Laguerre polynomial: $L^\alpha_n(z)$,
(or its $k$'th derivative),
$n=0, 1, 2, \dots$, typically $z\in[0,\infty]$
\begin{code}
sf_orthopoly_laguerre_value_dx !a !n !k !z
  | k>n = 0
  | otherwise = (-1)^(k`mod`2) * (sf_orthopoly_laguerre_value (a+#k) (n-k) z)

sf_orthopoly_laguerre_value_accum !k1 !k0 !a !n !z
  | n==0 = k1 1 k0
  | n==1 = k1 (1 + a - z) $ k1 1 k0
  | otherwise = iter 2 1 (1+a-z) (k1 (1+a-z) $ k1 1 k0)
  where
    iter !k !rm2 !rm1 !acc =
      let !rm0 = (((2*k)#+a-1-z)*rm1 - (k#+a-1)*rm2) / (#)k
      in if k==n then (k1 rm0 acc)
         else iter (k+1) rm1 rm0 (k1 rm0 acc)

sf_orthopoly_laguerre_value :: Double -> Int -> Double -> Double
sf_orthopoly_laguerre_value = sf_orthopoly_laguerre_value_accum  const 0

sf_orthopoly_laguerre_values :: Double -> Int -> Double -> [Double]
sf_orthopoly_laguerre_values = sf_orthopoly_laguerre_value_accum (:)   []

-- TODO: tests
{--
%!test assert( all(sf_orthopoly_laguerre_value(0,0:20)==sf_polynomial_value(sf_orthopoly_laguerre_coeffs(0),0:20)) )
%!test assert( all(sf_orthopoly_laguerre_value(1,0:20)==sf_polynomial_value(sf_orthopoly_laguerre_coeffs(1),0:20)) )
%!test assert( all(abs(1 - sf_orthopoly_laguerre_value(4,0:20)./sf_polynomial_value(sf_orthopoly_laguerre_coeffs(4),0:20))<1e-12) )
%!test assert( all(abs(1 - sf_orthopoly_laguerre_value(11,0:20)./sf_polynomial_value(sf_orthopoly_laguerre_coeffs(11),0:20))<1e-10) )
%!test assert( all(abs(1 - sf_orthopoly_laguerre_value(12,0:20)./sf_polynomial_value(sf_orthopoly_laguerre_coeffs(12),0:20))<1e-8) )

%!test assert( all(sf_orthopoly_laguerre_value(0,2,0:20)==sf_polynomial_value(sf_orthopoly_laguerre_coeffs(0,2),0:20)) )
%!test assert( all(sf_orthopoly_laguerre_value(1,2,0:20)==sf_polynomial_value(sf_orthopoly_laguerre_coeffs(1,2),0:20)) )
%!test assert( all(abs(1 - sf_orthopoly_laguerre_value(4,2,0:20)./sf_polynomial_value(sf_orthopoly_laguerre_coeffs(4,2),0:20))<1e-12) )
%!test assert( all(abs(1 - sf_orthopoly_laguerre_value(11,2,0:20)./sf_polynomial_value(sf_orthopoly_laguerre_coeffs(11,2),0:20))<1e-8) )
%!test assert( all(abs(1 - sf_orthopoly_laguerre_value(12,2,0:20)./sf_polynomial_value(sf_orthopoly_laguerre_coeffs(12,2),0:20))<1e-8) )
--}
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gaussian quadrature weights of the $n$'th (generalized) Laguerre polynomial: $L^{a}_n(z)$
$n=0, 1, 2, \dots$
\begin{code}
sf_orthopoly_laguerre_weights :: Double -> Int -> [Double]
sf_orthopoly_laguerre_weights !a !n
  | n==0 = []
  | n==1 = [1]
  | otherwise =
      let !zs = sf_orthopoly_laguerre_zeros a n
          --  need ortho_normal_ polynomials here
          jth j =
            let !jvs = map (sf_orthopoly_laguerre_value a j) zs
                !nj = sf_orthopoly_laguerre_scale a j
            in map (\x->(x*nj)^2) jvs
          !ires = foldl (P..+.) (replicate n 0) (map jth [0..(n-1)])
          !res = map (1/) ires
      in res
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Zeros}
Compute the zeros of the $n$'th (generalized) Laguerre polynomial: $L^\alpha_n(z)$
(or its $k$'th derivative)
$n=0, 1, 2, \dots$
not good for large
\begin{code}
sf_orthopoly_laguerre_zeros_dx !a !n !k 
  | k>n       = []  -- degenerate case... technically every point is a zero...
  | otherwise = sf_orthopoly_laguerre_zeros (a+#k) (n-k)

sf_orthopoly_laguerre_zeros !a !n
  | n==0 = []
  | n==1 = [1+a]
  | otherwise = undefined
{--
function res = sf_orthopoly_laguerre_zeros(n, a, dum, k)
  if (n==0)
    res = [];
  elseif (n==1)
    res = [1+a];
  else
    m = zeros(n);
    for k=1:n-1
      m(k,k+1) = -k;
      m(k,k) = 2*k-1 + a;
      m(k+1,k) = -k - a;
    endfor
    m(n,n) = 2*n-1 + a;
    res = sort(eig(m));

    # "polish" the results
    #fx = sf_orthopoly_laguerre_value(n, a, res)
    #dfx = sf_orthopoly_laguerre_value(n, a, res, [], 1);
    #nwt = res - fx./dfx;
    #nrs = sf_orthopoly_laguerre_value(n,nwt)
    #res = nwt;
  endif
endfunction
--}
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
  weight l x  = if -1<x && x<1 then 1 else 0
  coeffs _    = sf_orthopoly_legendre_coeffs
  scale _     = sf_orthopoly_legendre_scale
  value _     = sf_orthopoly_legendre_value
  value_arr _ = sf_orthopoly_legendre_values
  weights _   = sf_orthopoly_legendre_weights
  zeros _     = sf_orthopoly_legendre_zeros
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Coefficients}

Compute the coefficients of the $n$'th Legendre polynomial:
$P_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
$n=0, 1, 2, \dots$
(TODO: cache?)
\begin{code}
sf_orthopoly_legendre_coeffs :: Int -> [Double]
sf_orthopoly_legendre_coeffs !n = P.coeffs $
  recur_coeffs (P.Poly[1]) (P.Poly[0,1])
    (\k -> ((#)$2*k+1)/((#)$k+1))
    (const 0)
    (\k -> (#)k       /((#)$k+1))
    n
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Scale-factor}
Compute the scale-factor to normalize the $n$'th Legendre polynomial:
$P_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
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
sf_orthopoly_legendre_value_accum !k1 !k0 !n !z
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
sf_orthopoly_legendre_value  = sf_orthopoly_legendre_value_accum const 0
sf_orthopoly_legendre_values = sf_orthopoly_legendre_value_accum (:)   []
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Weights}
Compute the Gaussian quadrature weights of the $n$'th Legendre polynomial: $P_n(z)$
$n=0, 1, 2, \dots$

\begin{code}
-- cache?
--if (!sf_is_nonnegint(n)) print_usage; endif
sf_orthopoly_legendre_weights :: Int -> [Double]
sf_orthopoly_legendre_weights !n
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
sf_orthopoly_legendre_zeros !n
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Gegenbauer (ultraspherical) Polynomials}


Compute the coefficients of the $n$'th Gegenbauer (ultraspherical) polynomial
(or its $k$'th derivative):
$C^{(\lambda)}_n(z) = a_1 + a_2x + \cdots + a_{n+1}x^n$,
$n=0, 1, 2, \dots$, $\lambda>-1/2, \lambda\neq0$
\begin{code}
sf_orthopoly_gegenbauer_coeffs_dx :: Double -> Int -> Int -> [Double]
sf_orthopoly_gegenbauer_coeffs_dx !lambda !n !k
  -- if (!sf_is_nonnegint(k)) print_usage; endif
  | k==0      = sf_orthopoly_gegenbauer_coeffs lambda n
  | k>n       = [0]
  | otherwise = (2*lambda) P.*. (sf_orthopoly_gegenbauer_coeffs_dx (lambda+1) (n-1) (k-1))

sf_orthopoly_gegenbauer_coeffs :: Double -> Int -> [Double]
sf_orthopoly_gegenbauer_coeffs !lambda !n = P.coeffs $ recur_coeffs
  -- if (!sf_is_nonnegint(n) || !(a>-1/2) || a==0) print_usage; endif
  (P.Poly [1]) (P.Poly [0, 2*lambda])
  (\j -> 2*(j#+lambda)/((#)$j+1))
  (const 0)
  (\j -> (j#+2*lambda-1)/((#)$j+1))
  n
\end{code}

Compute the scale-factor to normalize the $n$'th Gegenbauer (ultraspherical) polynomial: $C^{(\lambda)}_n(z)$,
$n=0, 1, 2, \dots$, $\lambda>-1/2, \lambda\neq0$
\begin{code}
sf_orthopoly_gegenbauer_scale a n = undefined
{--
  -- if (!sf_is_nonnegint(n) || !(a>-1/2) || a==0) print_usage; endif
  (sf_pochhammer (a+1/2) n)*(sf_sqrt(2*(n+a)*(sf_factorial n)*(sf_gamma$n+2*a)))/(sf_pochhammer (2*a) n)/2^a/(sf_gamma $ n+a+1/2)
--}
\end{code}

Compute the value of the $n$'th Gegenbauer (ultraspherical) polynomial: $C^{(\lambda)}_n(z)$,
(or its $k$'th derivative)
$n=0, 1, 2, \dots$, $\lambda>-1/2, \lambda\neq0$, typically $z\in[-1,1]$
\begin{code}
sf_orthopoly_gegenbauer_value_dx :: Double -> Int -> Int -> Double -> Double
sf_orthopoly_gegenbauer_value_dx !lambda !n !k !z
  -- if (!sf_is_nonnegint(k)) print_usage; endif
  | k==0 = sf_orthopoly_gegenbauer_value lambda n z
  | k>n = 0
  | otherwise = (2*lambda) * (sf_orthopoly_gegenbauer_value_dx (lambda+1) (n-1) (k-1) z)

-- if (!sf_is_nonnegint(n) || !(a>-1/2) || a==0) print_usage; endif
sf_orthopoly_gegenbauer_value_accum :: (Double -> t -> t) -> t -> Double -> Int -> Double -> t
sf_orthopoly_gegenbauer_value_accum !k1 !k0 !lambda !n !z
  | n==0 = 1 `k1` k0
  | n==1 = (2*lambda*z) `k1` (1 `k1` k0)
  | otherwise = iter 2 1 (2*lambda*z) ((2*lambda*z) `k1` (1 `k1` k0))
  where
    iter !k !rm2 !rm1 !acc =
      let !rm0 = (2*z*((#)k+lambda-1)*rm1 - ((#)k+2*lambda-2)*rm2)/#k
      in if k==n then (rm0 `k1` acc)
         else iter (k+1) rm1 rm0 (rm0 `k1` acc)

sf_orthopoly_gegenbauer_value  :: Double -> Int -> Double -> Double
sf_orthopoly_gegenbauer_value  = sf_orthopoly_gegenbauer_value_accum const 0
sf_orthopoly_gegenbauer_values :: Double -> Int -> Double -> [Double]
sf_orthopoly_gegenbauer_values = sf_orthopoly_gegenbauer_value_accum (:)   []
\end{code}

Compute the Gaussian quadrature weights of the $n$'th Gegenbauer (ultraspherical) polynomial: $C^{(\lambda)}_n(z)$
$n=0, 1, 2, \dots$, $\lambda>-1/2, \lambda\neq0$
\begin{code}
-- if (!sf_is_nonnegint(n) || !(a>-1/2) || a==0) print_usage; endif
sf_orthopoly_gegenbauer_weights lambda n
  | n==0 = []
  | n==1 = [1]
  | otherwise =
    let !zs = sf_orthopoly_gegenbauer_zeros lambda n
        -- need ortho_normal_ polynomials here
        jth j =
          let !jvs = map (sf_orthopoly_gegenbauer_value lambda j) zs
              !nj = sf_orthopoly_gegenbauer_scale lambda j
          in map (\x->(x*nj)^2) jvs
        !ires = foldl (P..+.) (replicate n 0) (map jth [0..(n-1)])
        !res = map (1/) ires
    in res
\end{code}

\begin{code}
sf_orthopoly_gegenbauer_zeros :: Double -> Int -> [Double]
sf_orthopoly_gegenbauer_zeros lambda n = undefined
{--
## -*- texinfo -*-
## @deftypefn {Function File} {@var{res} =} sf_orthpoly_gegenbauer_zeros (@var{n}, @var{alpha})
## Compute the zeros of the $n$'th Gegenbauer (ultraspherical) polynomial: $C^(\alpha)_n(z)$
## $n=0, 1, 2, ...$, $\alpha>-1/2, \alpha!=0$
## @end deftypefn

function res = sf_orthpoly_gegenbauer_zeros(n, a)
  if (nargin < 2) print_usage; endif
  if (!sf_is_nonnegint(n) || !(a>-1/2) || a==0) print_usage; endif
  if (n==0)
    res = [];
  elseif (n==1)
    res = [0];
  else
    m = zeros(n);
    m(1,1+1) = m(1+1,1) = sf_sqrt(1/(2+2*a));
    for k=2:n-1
      m(k,k+1) = m(k+1,k) ...
               = 2/(2*k+2*a-1) * sf_sqrt((k*(k+a-1/2)^2*(k+2*a-1))/((2*k+2*a)*(2*k+2*a-2)));
    endfor
    res = sort(eig(m));

    # "polish" the results
    fx = sf_orthpoly_gegenbauer_value(n, a, res);
    dfx = sf_orthpoly_gegenbauer_value(n, a, res, [], 1);
    nwt = res - fx./dfx;
    #nrs = sf_orthpoly_gegenbauer_value(n, a, nwt)
    res = nwt;
  endif
endfunction
--}
\end{code}
