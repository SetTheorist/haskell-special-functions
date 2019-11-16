module ExpInt(
 sf_expint_ei,
 sf_expint_en,
 )
where

import Exp
import Gamma
import Util

----------------------------------------

sf_expint_ei :: Value -> Value
sf_expint_ei z
  | z<0.0  = (0/0)
  | z==0.0 = (-1/0)
  | z<40   = expint_ei__series z
  | otherwise = expint_ei__asymp z

expint_ei__series :: Value -> Value
expint_ei__series z =
  let tterms = ixiter 2 z $ \n t -> t*z/(#)n
      terms = map (\(t,n)->t/(#)n) $ zip tterms [1..]
      res = kahan_sum terms
  in if z<0.5
     then sf_log(z * sf_exp(euler_gamma + res))
     else res + sf_log(z) + euler_gamma

expint_ei__asymp :: Value -> Value
expint_ei__asymp z =
  let terms = tk $ ixiter 1 1.0 $ \n t -> t/z*(#)n
      res = kahan_sum terms
  in res * (sf_exp z) / z
  where tk (a:b:cs) = if (abs a)<(abs b) then [a] else a:(tk$b:cs)

----------------------------------------


sf_expint_en :: Int -> Value -> Value
sf_expint_en 0 z = sf_exp(-z) / z
sf_expint_en 1 z = expint_en__1 z
sf_expint_en n z | z<=1.0    = expint_en__series n z
                 | otherwise = expint_en__contfrac n z

expint_en__1 :: Value -> Value
expint_en__1 z =
    let r0 = -euler_gamma - (sf_log z)
        tterms = ixiter 2 (z) $ \k t -> -t*z/(#)k
        terms = map (\(t,k)->t/(#)k) $ zip tterms [1..]
    in kahan_sum (r0:terms)

expint_en__series :: Int -> Value -> Value
expint_en__series = undefined

expint_en__contfrac :: Int -> Value -> Value
expint_en__contfrac = undefined
{--
# assume n>=2, x<=1
function res = sf_expint_en_series(z,n)
  res = (-log(z) + sf_digamma(n)) * (-z)^(n-1)/factorial(n-1) + 1.0/(n-1);
  m = 1;
  term = 1.0;
  old_res = 0.0;
  do
    term *= -z/m;
    if (m==n-1)
      ++m;
      continue;
    endif
    old_res = res;
    res -= term / (m - (n-1));
    ++m;
    if (m>999) break; endif
  until (res == old_res);
endfunction

# assume n>=2, z>1
function res = sf_expint_en_contfrac(z,n)
  zeta = 1e-100;
  eps = 5e-16;
  # modified Lentz algorithm
  fj = zeta;
  Cj = fj;
  Dj = 0;
  j = 1;
  do
    if (j==1) aj=1; else aj=-(j-1)*(n+j-2); endif
    bj = z + n + 2*(j-1);
    Dj = bj + aj*Dj; if (Dj==0.0) Dj=zeta; endif
    Cj = bj + aj/Cj; if (Cj==0.0) Cj=zeta; endif
    Dj = 1/Dj;
    Delta = Cj*Dj;
    fj = fj*Delta;
    ++j;
    if (j>999) break; endif
  until (abs(Delta-1)<eps)
  res = fj * sf_exp(-z);
endfunction
--}


