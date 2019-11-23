module Trig where

import Util

----------------------------------------
-- cos(x)
-- TODO: range reduce
sf_cos :: (Value v) => v -> v
sf_cos x = kahan_sum $ ixiter 0 1 $ \n t -> -t*x^2/((#)$2*n+1)/((#)$2*n+2)

sf_acos :: (Value v) => v -> v
sf_acos = acos

-- vercosine function
sf_vcos :: (Value v) => v -> v
sf_vcos x = 2 * (sf_cos $ x/2)^2

-- havercosine function
sf_hvcos :: (Value v) => v -> v
sf_hvcos x = (sf_cos $ x/2)^2

----------------------------------------
-- sin(x)
-- TODO: range reduce
sf_sin :: (Value v) => v -> v
sf_sin x = kahan_sum $ ixiter 1 x $ \n t -> -t*x^2/((#)$2*n)/((#)$2*n+1)

sf_asin :: (Value v) => v -> v
sf_asin = asin

-- versine function
sf_vsin :: (Value v) => v -> v
sf_vsin x = 2 * (sf_sin $ x/2)^2

-- haversine function
sf_hvsin :: (Value v) => v -> v
sf_hvsin x = (sf_sin $ x/2)^2

----------------------------------------
-- secant function
sf_sec :: (Value v) => v -> v
sf_sec x = 1 / (cos x)

-- exsecant function
sf_exsec :: (Value v) => v -> v
sf_exsec x = 2 * (sf_sin $ x/2)^2 * (sf_sec x)

-- cosecant function
sf_csc :: (Value v) => v -> v
sf_csc x = 1 / (sin x)

----------------------------------------
-- tangent function
sf_tan :: (Value v) => v -> v
sf_tan = tan

-- arctangent function
sf_atan :: (Value v) => v -> v
sf_atan = atan

-- $$\tan z = z/(1- z^2/(3- z^2/(5- z^2/(7- ...))))$$
-- NB terrible convergence
tan_cf :: (Value v) => v -> v
tan_cf z = steeds (z:(cycle [-z^2,z^2])) (map (#) (0:[1,3..]))

-- cotangent
sf_cot :: (Value v) => v -> v
sf_cot x = sf_cos(x) / sf_sin(x)
-- = 1/sf_tan x

----------------------------------------
-- hyperbolic cosine function
sf_cosh :: (Value v) => v -> v
sf_cosh = cosh

sf_acosh :: (Value v) => v -> v
sf_acosh = acosh

-- hyperbolic sine function
sf_sinh :: (Value v) => v -> v
sf_sinh = sinh

sf_asinh :: (Value v) => v -> v
sf_asinh = asinh

-- hyperbolic secant function
sf_sech :: (Value v) => v -> v
sf_sech z = 1/sf_cosh z

-- hyperbolic tangent function
sf_tanh :: (Value v) => v -> v
sf_tanh = tanh

sf_atanh :: (Value v) => v -> v
sf_atanh = atanh

----------------------------------------
-- Gudermannian function
sf_gud :: (Value v) => v -> v
sf_gud x = sf_asin(sf_tanh x)
-- = 2*sf_atan(sf_exp(x)) - pi/2
-- = 2*sf_atan(sf_tanh(x/2))
-- = sf_atan(sf_sinh(x))

-- Compute the inverse Gudermannian function
sf_agud :: (Value v) => v -> v
sf_agud z = sf_asinh (sf_tan z)
-- = sf_log(abs((1+sf_sin(z))/sf_cos(z)))
-- = sf_log(abs((1+sf_sin(z))/(1-sf_sin(z))))/2
-- = sf_log(abs(sf_tan(z) + sf_sec(z))
-- = sf_log(abs(sf_tan(pi/4 + z/2)))
-- = sf_atanh(sf_sin(z))

