module Trig where

import Util

----------------------------------------
-- cos(x)
-- TODO: range reduce
sf_cos :: Value -> Value
sf_cos x = kahan_sum $ ixiter 0 1 $ \n t -> -t*x^2/((#)$2*n+1)/((#)$2*n+2)

sf_cosh = cosh

-- vercosine function
sf_vcos :: Value -> Value
sf_vcos x = 2 * (sf_cos $ x/2)^2

-- havercosine function
sf_hvcos :: Value -> Value
sf_hvcos x = (sf_cos $ x/2)^2

----------------------------------------
-- sin(x)
-- TODO: range reduce
sf_sin :: Value -> Value
sf_sin x = kahan_sum $ ixiter 1 x $ \n t -> -t*x^2/((#)$2*n)/((#)$2*n+1)

sf_sinh = sinh

sf_asin :: Value -> Value
sf_asin = asin

-- versine function
sf_vsin :: Value -> Value
sf_vsin x = 2 * (sf_sin $ x/2)^2

-- haversine function
sf_hvsin :: Value -> Value
sf_hvsin x = (sf_sin $ x/2)^2

----------------------------------------
-- secant function
sf_sec x = 1 / (cos x)

-- exsecant function
sf_exsec x = 2 * (sf_sin $ x/2)^2 * (sf_sec x)

----------------------------------------
-- tangent function
sf_tan = tan

-- hyperbolic tangent function
sf_tanh = tanh

-- cotangent
sf_cot x = sf_cos(x) / sf_sin(x)
-- = 1/sf_tan x


----------------------------------------
-- Gudermannian function
sf_gud x = sf_asin(sf_tanh x)
-- = 2*sf_atan(sf_exp(x)) - pi/2
-- = 2*sf_atan(sf_tanh(x/2))
-- = sf_atan(sf_sinh(x))

-- $$\tan z = z/(1- z^2/(3- z^2/(5- z^2/(7- ...))))$$
-- NB terrible convergence
tan_cf z = steeds (z:(cycle [-z^2,z^2])) (0:[1,3..])


