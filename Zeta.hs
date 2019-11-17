module Zeta (
  sf_zeta,
) where

import Gamma
import Trig
import Util

-- Compute the Riemann zeta function
sf_zeta :: Value -> Value
sf_zeta z
  | z==1      = (1/0)
  | z<0       = 2 * (2*pi)**(z-1) * (sf_sin$pi*z/2) * (sf_gamma$1-z) * (sf_zeta$1-z)
  | otherwise = zeta_series z


zeta_series :: Value -> Value
zeta_series z = 
  let terms = map (\n->((#)n)**(-z)) [2..]
      corrs = map correction [2..]
  in summer terms corrs 1.0 0.0 0.0
  where
    --TODO: make general "corrected" kahan_sum!
    summer (t:ts) (c:cs) s e r = 
      let y = t + e
          s' = s + y
          e' = (s - s') + y
          r' = s' + c + e'
      in if r==r' then r'
         else summer ts cs s' e' r'
    zz1 = z/12
    zz2 = z*(z+1)*(z+2)/720
    zz3 = z*(z+1)*(z+2)*(z+3)*(z+4)/30240
    zz4 = z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6)/1209600
    zz5 = z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6)*(z+7)*(z+8)/239500800
    correction n' =
      let n=(#)n'
      in n**(1-z)/(z-1) - n**(-z)/2 + n**(-z-1)*zz1 - n**(-z-3)*zz2 + n**(-z-5)*zz3 - n**(-z-7)*zz4 + n**(-z-9)*zz5

{--
function res = sf_zeta_1(z)
  oldold_res = old_res = res = 0.0;
  smm = 1.0;
  e_ = 0;
  n = 2;
  do
    #smm += n^(-z);
    t_ = smm;
    y_ = n^(-z) + e_;
    smm = t_ + y_;
    e_ = (t_ - smm) + y_;
    oldold_res = old_res;
    old_res = res;
    res = smm ...
        + n^(1-z)/(z-1) - n^(-z)/2   ...
        + n^(-z-1)*zz1 - n^(-z-3)*zz2 + n^(-z-5)*zz3 - n^(-z-7)*zz4 + n^(-z-9)*zz5 ...
        + e_;
    ++n;
    if (n>999) break; endif
  until ((res == old_res) && (oldold_res == res))
endfunction
--}

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

{--
## -*- texinfo -*-
## @deftypefn {Function File} {@var{res} =} sf_zeta_m1 (@var{z})
## Compute Riemann zeta function minus 1.
## @end deftypefn
function res = sf_zeta_m1(z)
  if (nargin < 1) print_usage; endif
  res = ones(size(z));
  for kk = 1:prod(size(z))
    res(kk) = sf_zeta_m1_1(z(kk));
  endfor
endfunction

function res = sf_zeta_m1_1(z)
  if (z==1.0) res = Inf; return endif
# TODO: fixup reflection here
#  if (z < 1.0)
#    res = nan;
#    return;
#  endif
  oldold_res = old_res = res = 0.0;
  sum = 0.0;
  n = 2;
  em1 = z/12.0;
  em2 = z*(z+1)*(z+2)/720.0;
  em3 = z*(z+1)*(z+2)*(z+3)*(z+4)/30240.0;
  em4 = z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6)/1209600.0;
  em5 = z*(z+1)*(z+2)*(z+3)*(z+4)*(z+5)*(z+6)*(z+7)*(z+8)/239500800.0;
  do
    sum += n^(-z);
    oldold_res = old_res;
    old_res = res;
    res = sum + n^(1-z)/(z-1) - n^(-z)/2   ...
        + n^(-z-1)*em1 - n^(-z-3)*em2 ...
        + n^(-z-5)*em3 - n^(-z-7)*em4 ...
        + n^(-z-9)*em5 ;
    ++n;
    if (n>999) break; endif
  until ((res == old_res) && (oldold_res == res))
endfunction
--}
