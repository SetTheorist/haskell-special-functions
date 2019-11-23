module AGM (
    sf_agm,
    sf_agm',
    )
where

import Util

-- messy
sf_agm :: (Value v) => v -> v -> ([v],[v],[v])
sf_agm alpha beta = agm [alpha] [beta] [alpha-beta]
  where agm as@(a:_) bs@(b:_) cs@(c:_) = 
          if c==0 then (as,bs,cs)
          else let a' = (a+b)/2
                   b' = sf_sqrt (a*b)
                   c' = (a-b)/2
               in if c'==c then (as,bs,cs)
                  else agm (a':as) (b':bs) (c':cs)
  
sf_agm' :: (Value v) => v -> v -> v
sf_agm' alpha beta = agm alpha beta ((alpha-beta)/2)
  --let (as,_,_) = sf_agm alpha beta in head as
  where agm a b 0 = a
        agm a b c =
          let a' = (a+b)/2
              b' = sf_sqrt (a*b)
              c' = (a-b)/2
          in agm a' b' c'
    

sf_agm_c0 :: (Value v) => v -> v -> v -> ([v],[v],[v])
sf_agm_c0 alpha beta c0 = undefined
