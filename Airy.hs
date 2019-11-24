module Airy (
    sf_airy_ai,
    sf_airy_bi,
) where

import Gamma
import Util

----------------------------------------

sf_airy_ai :: (Value v) => v -> v
sf_airy_ai z = airy_ai_series z

sf_airy_bi :: (Value v) => v -> v
sf_airy_bi z = airy_bi_series z

ai0 :: (Value v) => v
ai0 = 3**(-2/3)/sf_gamma(2/3)

ai'0 :: (Value v) => v
ai'0 = -3**(-1/3)/sf_gamma(1/3)

airy_ai_series z =
    let z3 = z^3
        aiterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        ai'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in ai0 * (kahan_sum aiterms) + ai'0 * (kahan_sum ai'terms)

bi0 :: (Value v) => v
bi0 = 3**(-1/6)/sf_gamma(2/3)

bi'0 :: (Value v) => v
bi'0 = 3**(1/6)/sf_gamma(1/3)

airy_bi_series z =
    let z3 = z^3
        biterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        bi'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in bi0 * (kahan_sum biterms) + bi'0 * (kahan_sum bi'terms)

