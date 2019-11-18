module Spence (
    sf_spence,
) where

import Exp
import Util

-- Compute Spence's integral.

pi2_6 = pi^2/6
sf_spence :: Value -> Value
sf_spence z
  | z/=z  = z
  | z<0   = 0/0
  | z==0  = pi2_6
  | z<0.5 =
      (spence__series z) + (pi2_6 - (sf_log z)*(sf_log (1-z)))
  | z<1.0 =
      -(spence__series (1-z))
  | z<2.5 =
      (spence__series ((z-1)/z)) - (sf_log z)^2/2
  | otherwise = 
      (spence__series (1/(1-z))) - pi2_6 - (sf_log (z-1))^2/2

spence__series z = 
  let zk = iterate (*z) z
      terms = map (\(t,k)-> -t/(#)k^2) (zip zk [1..])
  in kahan_sum terms

{--
function res = spence_single(z)
  persistent pi2_6 = pi^2/6;

  if (isnan(z)) res = z; return; endif
  if (z < 0) res = nan; return; endif
  if (z == 0) res = pi2_6; return; endif

  if (z < 0.5) method = 0;
  elseif (z < 1.0) method = 1;
  elseif (z < 2.5) method = 3;
  else method = 2;
  endif

  # different reflections
  switch (method)
  case 0
    z_mult = z;
    [res,e_] = series(z);
    term = pi2_6 - sf_log(z)*sf_log(1-z);
  case 1
    [res,e_] = series(1 - z);
     res = -res; e_ = -e_;
     term = 0.0;
  case 2
    [res,e_] = series(1/(1-z));
    term = -pi2_6 - sf_log(z-1)^2/2;
  case 3
    #if (abs((z-1)/z)<1)
      [res,e_] = series((z - 1)/z);
      #[res,e_] = series(1 - 1/z);
    #else
    #  [res,e_] = series(z/(z-1));
    # TODO: - finish this!  need this inversion
    term = -sf_log(z)^2/2;
  endswitch
  res = (res + term) + e_;
endfunction
--}
