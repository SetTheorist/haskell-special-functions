module AGM (
    sf_agm,
    sf_agm',
    )
where

import Util


-- messy
sf_agm :: Value -> Value -> ([Value],[Value],[Value])
sf_agm alpha beta = agm [alpha] [beta] [alpha-beta]
    where agm as@(a:_) bs@(b:_) cs@(c:_) = 
            if c==0 then (as,bs,cs)
            else let a' = (a+b)/2
                     b' = sqrt (a*b)
                     c' = (a-b)/2
                 in if c'==c then (as,bs,cs)
                    else agm (a':as) (b':bs) (c':cs)
    
sf_agm' :: Value -> Value -> Value
sf_agm' alpha beta =
    let (as,_,_) = sf_agm alpha beta
    in head as
