module Erf where

import Exp
import Util

----------------------------------------
-- not perfect, but workable for now

sf_erf :: Value -> Value
sf_erf z 
    | z<(-1)    = -sf_erf(-z)
    | z<1       = erf_series z
    | otherwise = 1 - sf_erfc z

sf_erfc :: Value -> Value
sf_erfc z 
    | z<(-1)    = 2-(sf_erfc (-z))
    | z<1       = 1-(sf_erf z)
    | z<10      = erfc_cf_pos1 z
    | otherwise = erfc_asymp_pos z

erf_series z =
    let z2 = z^2
        rts = ixiter 1 z $ \n t -> (-t)*z2/(#)n
        terms = map (\(n,t)->t/(#)(2*n+1)) $ zip [0..] rts
    in (2/sf_sqrt pi)  * kahan_sum terms

erfc_asymp_pos z =
    let z2 = z^2
        iz2 = 1/2/z2
        terms = ixiter 1 (1/z) $ \n t -> (-t*iz2)*(#)(2*n-1)
        tterms = tk terms
    in (sf_exp (-z2))/(sqrt pi) * kahan_sum tterms
    where tk (a:b:cs) = if (abs a)<(abs b) then [a] else a:(tk$b:cs)

erfc_cf_pos1 z = 
    let z2 = z^2
        as = z:[1/2,1..]
        bs = 0:cycle [z2,1]
        cf = steeds as bs
    in sf_exp(-z2) / (sqrt pi) * cf

erfc_cf_pos2 z = 
    let z2 = z^2
        as = (2*z):(map (\n->(#)$ -(2*n+1)*(2*n+2)) [0..])
        bs = 0:(map (\n->2*z2+(#)4*n+1) [0..])
        cf = steeds as bs
    in sf_exp(-z2) / (sqrt pi) * cf

