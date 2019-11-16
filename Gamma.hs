module Gamma where

import Exp
import Util

euler_gamma :: (Floating a) => a
euler_gamma = 0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749

factorial :: (Integral a) => a -> a
factorial 0 = 1
factorial 1 = 1
factorial n = n * (factorial $ n-1)

----------------------------------------
-- gamma function
-- $$\Gamma(z) = \int_0^\infty e^{-t}t^{z}\frac{dz}{z}$$
sf_gamma :: Value -> Value
sf_gamma x =
    let (x',t) = redup x 1
        lngx = lngamma_asymp x'
    in t * (sf_exp lngx)
    where redup x t
            | x>15 = (x,t)
            | otherwise = redup (x+1) (t/x)

-- $$\frac{1}{\Gamma(z)}$$
sf_invgamma :: Value -> Value
sf_invgamma x =
    let (x',t) = redup x 1
        lngx = lngamma_asymp x'
    in t * (sf_exp$ -lngx)
    where redup x t
            | x>15 = (x,t)
            | otherwise = redup (x+1) (t*x)

-- log-gamma function
-- $$\ln\Gamma(z)$$
sf_lngamma :: Value -> Value
sf_lngamma x =
    let (x',t) = redup x 0
        lngx = lngamma_asymp x'
    in t + lngx
    where redup x t
            | x>15 = (x,t)
            | otherwise = redup (x+1) (t-sf_log x)

lngamma_asymp z = (z - 1/2)*(sf_log z) - z + (1/2)*sf_log(2*pi) + (kahan_sum terms)
    where terms = [b2k/(2*k*(2*k-1)*z^(2*k'-1)) | k'<-[1..10], let k=(#)k', let b2k=bernoulli_b$2*k']

bernoulli_b 1 = -1/2
bernoulli_b k | k`mod`2==1 = 0
bernoulli_b 0 = 1
bernoulli_b 2 = 1/6
bernoulli_b 4 = -1/30
bernoulli_b 6 = 1/42
bernoulli_b 8 = -1/30
bernoulli_b 10 = 5/66
bernoulli_b 12 = -691/2730
bernoulli_b 14 = 7/6
bernoulli_b 16 = -3617/510
bernoulli_b 18 = 43867/798
bernoulli_b 20 = -174611/330

{--
sf_gamma :: Value -> Value
sf_gamma x = spouge_approx 17 x
-- Spouge's approximation
spouge_approx a z' =
  let z = z' - 1
      a' = (#) a
      res = (z+a')**(z+(1/2)) * sf_exp (-(z+a'))
      sm = sqrt(2*pi)
      terms = [(spouge_c k a') / (z+k') | k<-[1..(a-1)], let k' = (#) k]
      sum = kahan_sum terms
  in res*sum
  where
    spouge_c k a = ((if k`mod`2==0 then -1 else 1) / ((#) $ factorial (k-1)))
                    * (a-((#) k))**(((#) k)-1/2) * sf_exp(a-((#) k))

--}
