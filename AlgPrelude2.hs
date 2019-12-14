{-# Language AllowAmbiguousTypes #-}
{-# Language BangPatterns #-}
{-# Language FlexibleContexts #-}
{-# Language FlexibleInstances #-}
{-# Language MagicHash #-}
{-# Language MultiParamTypeClasses #-}
{-# Language NoImplicitPrelude #-}
{-# Language RebindableSyntax #-}

module AlgPrelude (
  ) where

import Prelude(Int,Integer,Float,Double,Ord,Eq,Bool,Rational,(<),(==),(/=),error)
import qualified Prelude as P

----------------------------------------

infixl 6 .+
infixl 6 .-
infixl 7 .*
infixl 7 ./
infixr 8 .^
infixr 8 .^^
infixr 8 .**

(.+) :: (P.Num t) => t -> t -> t
(.+) = (P.+)
(.-) :: (P.Num t) => t -> t -> t
(.-) = (P.-)
(.*) :: (P.Num t) => t -> t -> t
(.*) = (P.*)
(./) :: (P.Fractional t) => t -> t -> t
(./) = (P./)
(.^) = (P.^)
(.^^) = (P.^^)
(.**) = (P.**)

infixl 6 +
infixl 6 -
infixl 7 *
infixl 7 /
infixr 8 ^
infixr 8 ^^
infixr 8 *^
infixr 8 *^^

----------------------------------------

-- TODO: do this the smart way!
-- assumes associativity of operation
power :: (Integral i) => (t -> t -> t) -> t -> (t -> t) -> t -> i -> t
power o i n b e
  | e < 0     = power o i n (n b) (negate e)
  | e == 0    = i
  | e == e    = b `o` (power o i n b (e-1))

----------------------------------------

class AddSemigroup t where
  (+) :: t -> t -> t

class AddSemigroup t => AddMonoid t where
  zero :: t
  (*^) :: (Integral i) => t -> i -> t
  (*^) = power (+) zero (error "Negative multiple")

class AddMonoid t => AddGroup t where
  (-) :: t -> t -> t
  negate :: t -> t
  negate = (zero -)
  (*^^) :: (Integral i) => t -> i -> t
  (*^^) = power (+) zero negate

----------------------------------------

class MulSemigroup t where
  (*) :: t -> t -> t

class MulSemigroup t => MulMonoid t where
  one :: t
  (^) :: (Integral i) => t -> i -> t
  (^) = power (*) one (error "Negative power")

class MulMonoid t => MulGroup t where
  (/) :: t -> t -> t
  recip :: t -> t
  recip = (one /)
  (^^) :: (Integral i) => t -> i -> t
  (^^) = power (*) one recip
  -- (^^) = power (*) one recip

----------------------------------------

class (AddGroup t, MulMonoid t) => Ring t where

class (Ring t, MulGroup t) => Field t where
  iszero :: t -> Bool

----------------------------------------

class (Ord t, Ring t) => Integral t where
  fromInteger :: Integer -> t
  toInteger :: t -> Integer

class (Ord t, Field t) => FromRational t where
  fromRational :: Rational -> t

----------------------------------------

infixr 5 $+
infixr 5 $*
infixl 5 +$
infixl 5 *$

class AddLAction a x where
  ($+) :: a -> x -> x

class MulLAction a x where
  ($*) :: a -> x -> x

class AddRAction a x where
  (+$) :: x -> a -> x

class MulRAction a x where
  (*$) :: x -> a -> x
  

----------------------------------------
----------------------------------------

instance AddSemigroup Int where (+) = (.+)
instance AddMonoid    Int where zero = 0
instance AddGroup     Int where (-) = (.-)
instance MulSemigroup Int where (*) = (.*)
instance MulMonoid    Int where one = 1
instance Ring         Int
instance Integral     Int where
  fromInteger x = P.fromIntegral x
  toInteger x = P.fromIntegral x

----------------------------------------
----------------------------------------

instance AddSemigroup Integer where (+) = (.+)
instance AddMonoid    Integer where zero = 0
instance AddGroup     Integer where (-) = (.-)
instance MulSemigroup Integer where (*) = (.*)
instance MulMonoid    Integer where one = 1
instance Ring Integer
instance Integral Integer where
  fromInteger x = x
  toInteger x = x

----------------------------------------
----------------------------------------

instance Integral Double where
  fromInteger x = P.fromIntegral x
  toInteger x = P.floor x

instance FromRational Double where
  fromRational x = P.fromRational x

instance AddSemigroup Double where (+) = (.+)
instance AddMonoid    Double where zero = 0.0
instance AddGroup     Double where (-) = (.-)

instance MulSemigroup Double where (*) = (.*)
instance MulMonoid    Double where one = 1.0
instance MulGroup     Double where (/) = (./)

instance Ring Double

instance Field Double where
  iszero = (==0.0)

----------------------------------------
----------------------------------------

instance AddSemigroup [a] where (+) = (P.++)
instance AddMonoid    [a] where zero = []

instance AddSemigroup a => AddLAction a [a] where a $+ as = P.fmap (a+) as
instance AddSemigroup a => AddRAction a [a] where as +$ a = P.fmap (+a) as
instance MulSemigroup a => MulLAction a [a] where a $* as = P.fmap (a*) as
instance MulSemigroup a => MulRAction a [a] where as *$ a = P.fmap (*a) as



