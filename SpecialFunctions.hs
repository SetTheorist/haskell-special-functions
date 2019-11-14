--module SpecialFunctions where

import Control.Monad

type Value = Double

(#) :: (Integral a, Num b) => a -> b
(#) = fromIntegral

relerr :: Value -> Value -> Value
relerr ex ap = logBase 10 (abs ((ex-ap)/ex))

euler_gamma :: (Floating a) => a
euler_gamma = 0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749

factorial :: (Integral a) => a -> a
factorial 0 = 1
factorial 1 = 1
factorial n = n * (factorial $ n-1)

kahan_sum :: [Value] -> Value
kahan_sum terms = k 0.0 0.0 terms
  where
    k sum err [] = sum
    k sum err (t:terms) =
      let y = t - err
          x = sum + y
          err' = (x - sum) - y
      in if x==sum
         then sum
         else (k x err' terms)

ixiter :: (Enum ix) => ix -> a -> (ix->a->a) -> [a]
ixiter i x f = x:(ixiter (succ i) (f i x) f)

----------------------------------------
-- exp(x)
sf_exp :: Value -> Value 
-- exp(-x) = 1/exp(x)
sf_exp x | isInfinite x = if x<0 then 0 else (1.0/0.0)
sf_exp x | x<0 = 1/(sf_exp (-x))
sf_exp x = kahan_sum $ ixiter 1 1.0 $ \n t -> t*x/((#)n)

test__sf_exp :: IO ()
test__sf_exp = do
  cs <- readFile "exp_1.dat"
  let ls = lines cs
  let ws = map words ls
  let ns = map (map read) ws :: [[Value]]
  let re = map (\(x:fx:_)->relerr fx (sf_exp x)) ns
  mapM_ (putStrLn.show) re

----------------------------------------
-- exp(x)-1
sf_expm1 :: Value -> Value 
-- exp(-x)-1 = -exp(-x)*(exp(x)-1)
sf_expm1 x | isInfinite x = if x<0 then -1 else (1.0/0.0)
sf_expm1 x | x<0 = -sf_exp x * sf_expm1 (-x)
sf_expm1 x = kahan_sum $ ixiter 2 x $ \n t -> t*x/((#)n)

----------------------------------------
sf_log = log

sf_sqrt = sqrt

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

----------------------------------------
-- Bessel J(nu,x)
bessel_j_series :: Value -> Value -> Value
bessel_j_series nu z = 
  let z2 = -(z/2)^2
      terms = ixiter 1 1 $ \n t -> t*z2/((#) n)/(nu+(#) n)
      res = kahan_sum terms
  in res * (z/2)**nu / sf_gamma (1+nu)

----------------------------------------

-- $$\ln(1+z) = z/(1+ z/(2+ z/(3+ 4z/(4+ 4z/(5+ 9z/(6+ 9z/(7+ ...)))))))$$
-- NB not bad convergence
ln_1_z_cf z = steeds (z:(ts 1)) [0..]
    where ts n = (n^2*z):(n^2*z):(ts (n+1))

-- $$\tan z = z/(1- z^2/(3- z^2/(5- z^2/(7- ...))))$$
-- NB terrible convergence
tan_cf z = steeds (z:(cycle [-z^2,z^2])) (0:[1,3..])

-- Steed's algorithm for evaluation of continued fraction
-- C = b0 + a1/(b1+ a2/(b2+ ...))
-- Cn = ...an/bn = An/Bn
--steeds :: [Value] -> [Value] -> [Value]
steeds :: [Value] -> [Value] -> Value
steeds (a1:as) (b0:b1:bs) =
    let c0 = b0
        d1 = 1/b1
        delc1 = a1*d1
        c1 = c0 + delc1
    in recur c1 delc1 d1 as bs
    --in c0:c1:(recur c1 delc1 d1 as bs)
    where recur cn_1 delcn_1 dn_1 (an:as) (bn:bs) = 
            let dn = 1/(dn_1*an+bn)
                delcn = (bn*dn - 1)*delcn_1
                cn = cn_1 + delcn
            in if (cn==cn_1) then cn else (recur cn delcn dn as bs)
            --in cn:(recur cn delcn dn as bs)

----------------------------------------

sf_airy_ai z = airy_ai_series z
sf_airy_bi z = airy_bi_series z

ai0 = 3**(-2/3)/sf_gamma(2/3)
ai'0 = -3**(-1/3)/sf_gamma(1/3)
airy_ai_series z =
    let z3 = z^3
        aiterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        ai'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in ai0 * (kahan_sum aiterms) + ai'0 * (kahan_sum ai'terms)

bi0 = 3**(-1/6)/sf_gamma(2/3)
bi'0 = 3**(1/6)/sf_gamma(1/3)
airy_bi_series z =
    let z3 = z^3
        biterms  = ixiter 0 1 $ \n t -> t*z3*((#)$3*n+1)/((#)$(3*n+1)*(3*n+2)*(3*n+3))
        bi'terms = ixiter 0 z $ \n t -> t*z3*((#)$3*n+2)/((#)$(3*n+2)*(3*n+3)*(3*n+4))
    in bi0 * (kahan_sum biterms) + bi'0 * (kahan_sum bi'terms)

main :: IO ()
main = do
    forM_ [-10,-9.99..5] $ \x -> do
        let ax = sf_airy_ai x
        let bx = sf_airy_bi x
        putStr $ (show x)++" "++(show ax)++" "++(show bx)++"\n"

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

{--
main :: IO ()
main = do
    forM_ [-20,-19.99..20] $ \x -> do
        let ex = sf_erf x
        let cx = sf_erfc x
        putStr $ (show x)++" "++(show ex)++" "++(show cx)++"\n"
--}
