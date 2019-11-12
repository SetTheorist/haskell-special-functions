module SpecialFunctions where

type Value = Double

relerr :: Value -> Value -> Value
relerr ex ap = logBase 10 (abs ((ex-ap)/ex))

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

----------------------------------------
-- exp(x)
sf_exp :: Value -> Value 
-- exp(-x) = 1/exp(x)
sf_exp x | isInfinite x = if x<0 then 0 else (1.0/0.0)
sf_exp x | x<0 = 1/(sf_exp (-x))
sf_exp x = kahan_sum $ map snd $ iterate (\(n,t)->(n+1,t*x/(fromIntegral n))) (1,1.0)

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
sf_expm1 x = kahan_sum $ map snd $ iterate (\(n,t)->(n+1,t*x/(fromIntegral n))) (2,x)

----------------------------------------
-- cos(x)
-- TODO: range reduce
sf_cos :: Value -> Value
sf_cos x = kahan_sum $ map snd $ iterate (\(n,t)->(n+2,(-t*x*x/(fromIntegral n)/(fromIntegral (n+1))))) (1,1)

-- vercosine function
sf_vcos :: Value -> Value
sf_vcos x = 2 * (sf_cos $ x/2)^2

----------------------------------------
-- sin(x)
-- TODO: range reduce
sf_sin :: Value -> Value
sf_sin x = kahan_sum $ map snd $ iterate (\(n,t)->(n+2,(-t*x*x/(fromIntegral n)/(fromIntegral (n+1))))) (2,x)

-- versine function
sf_vsin :: Value -> Value
sf_vsin x = 2 * (sf_sin $ x/2)^2

----------------------------------------
-- Compute the exsecant function
sf_exsec x = 2 * (sf_sin $ x/2)^2 * (sf_sec x)


----------------------------------------
sf_gamma :: Value -> Value
sf_gamma x = spouge_approx 17 x

-- Spouge's approximation
spouge_approx a z' =
  let z = z' - 1
      a' = fromIntegral a
      res = (z+a')**(z+(1/2)) * sf_exp (-(z+a'))
      sm = sqrt(2*pi)
      terms = [(spouge_c k a') / (z+k') | k<-[1..(a-1)], let k' = fromIntegral k]
      sum = kahan_sum terms
  in res*sum
  where
    spouge_c k a = ((if k`mod`2==0 then -1 else 1) / (fromIntegral $ factorial (k-1)))
                    * (a-(fromIntegral k))**((fromIntegral k)-1/2) * sf_exp(a-(fromIntegral k))


----------------------------------------
-- Bessel J(nu,x)
bessel_j_series :: Value -> Value -> Value
bessel_j_series nu z = 
  let z2 = -(z/2)^2
      terms = map snd $ iterate (\(n,t)->(n+1,t*z2/(fromIntegral n)/(nu+fromIntegral n))) (1,1)
      res = kahan_sum terms
  in res * (z/2)**nu / sf_gamma (1+nu)



----------------------------------------

