module Probability where

import Exp
import Util

class ProbabilityDistribution d where
    in_support :: d -> Double -> Bool
    pdf :: d -> Double -> Double
    cdf :: d -> Double -> Double
    mean :: d -> Double
    variance :: d -> Double
    stddev :: d -> Double
    skew :: d -> Double
    kurtosis :: d -> Double
    moment :: d -> Int -> Double

----------------------------------------

data Uniform = Uniform Double Double

instance ProbabilityDistribution Uniform where
    in_support (Uniform a b)   x = a<=x && x<=b
    pdf        u@(Uniform a b) x = if in_support u x then 1/(b-a) else 0
    cdf        (Uniform a b)   x = if x<=a then 0 else if x>=b then 1 else (x-a)/(b-a)
    mean       (Uniform a b)     = a+(b-a)/2
    variance   (Uniform a b)     = (b-a)^2/12
    stddev     (Uniform a b)     = (b-a)/(sf_sqrt 12)
    skew       (Uniform a b)     = 0
    kurtosis   (Uniform a b)     = -6/5
    moment     (Uniform a b)   n = undefined

----------------------------------------

data Normal = Normal Double Double

instance ProbabilityDistribution Normal where
    in_support (Normal m v)   _ = True
    pdf        (Normal m v)   x = sf_exp (-(x-m)^2/2/v) / sqrt(2*pi*v)
    cdf        (Normal m v)   x = undefined
    mean       (Normal m v)     = m
    variance   (Normal m v)     = v
    stddev     (Normal m v)     = sf_sqrt v
    skew       (Normal m v)     = 0
    kurtosis   (Normal m v)     = 0
    moment     (Normal m v)   n = undefined

----------------------------------------




