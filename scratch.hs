
fac :: Integer -> Integer
fac 0 = 1
fac 1 = 1
fac n = n*(fac (n-1))

sums :: [Double] -> [Double]
sums = scanl (+) 0

converge (s0:s1:ss) = if s0==s1 then [s0] else s0:(converge (s1:ss))

aitken (s0:s1:s2:ss) = (s2 - (s2-s1)^2/(s2-2*s1+s0)):(aitken (s1:s2:ss))
aitken' (s0:s1:s2:ss) = (s2 - (s2-s1)^2/(s2-2*s1+s0)):(aitken ss)

expts :: Double -> [Double]
expts x = map (\n->x^n/(fromIntegral$fac n)) [0..]

ln1xs :: Double -> [Double]
ln1xs x = map (\n->x^n/(fromIntegral$n)) [1..]

