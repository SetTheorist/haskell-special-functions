module Fibo where

-- A silly way to compute Fibonacci numbers

import Data.Ratio

data Q5 = Q5 Rational Rational
  deriving (Eq)

norm (Q5 ra qa) = ra^2-5*qa^2

instance Show Q5 where
  --show (Q5 0  0) = "0"
  --show (Q5 ra 0) = (show ra)
  --show (Q5 0  qa) = (show qa)++"*sqrt(5)"
  show (Q5 ra qa) = (show ra)++"+"++(show qa)++"*sqrt(5)"

instance Num Q5 where
  (Q5 ra qa) + (Q5 rb qb) = Q5 (ra+rb) (qa+qb)
  (Q5 ra qa) - (Q5 rb qb) = Q5 (ra-rb) (qa-qb)
  (Q5 ra qa) * (Q5 rb qb) = Q5 (ra*rb+5*qa*qb) (ra*qb+rb*qa)
  negate (Q5 ra qa) = Q5 (-ra) (-qa)
  abs a = Q5 (norm a) 0
  signum a@(Q5 ra qa) = if a==0 then 0 else Q5 (ra/(norm a)) (qa/(norm a))
  fromInteger n = Q5 (fromInteger n) 0

instance Fractional Q5 where
  --(Q5 ra qa) / (Q5 rb qb) =
  recip a@(Q5 ra qa) = Q5 (ra/(norm a)) (-qa/(norm a))
  fromRational r = (Q5 r 0)

phip = Q5 (1%2) (1%2)
cp   = Q5 0     (1%5)
phim = Q5 (1%2) (-1%2)
cm   = Q5 0     (-1%5)
--fibonacci n = let (Q5 r q) = cp*phip^n + cm*phim^n in numerator r
--fibonacci n = let (Q5 r _) = cp*phip^n in numerator (2*r)
fibonacci n = let (Q5 _ q) = phip^^n in numerator (2*q)

