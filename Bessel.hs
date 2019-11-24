{-# Language BangPatterns #-}

module Bessel where

import Gamma
import Trig
import Util

----------------------------------------
----------------------------------------
-- Compute Bessel J_nu(z) function

sf_bessel_j :: (Value v) => v -> v -> v
sf_bessel_j nu z = undefined
  --asy = asympt_z(z, nu);
  --ser = series(z, nu);
  --sys = besselj(nu,z);
  --rec = recur_back(z, nu);
  --ref = recur_fore(z, nu);
  --re2 = recur_backwards(nu, z, round(abs(max(z, nu)))+21);
  --res = sys;

----------------------------------------
-- straightforward power-series
bessel_j_series :: (Value v) => v -> v -> v
bessel_j_series !nu !z = 
  let !z2 = -(z/2)^2
      !terms = ixiter 1 1 $ \n t -> t*z2/((#)n)/(nu+(#)n)
      !res = ksum terms
  in res * (z/2)**nu / sf_gamma (1+nu)

-- asymptotic expansion
-- for |z|>>nu
-- for |arg(z)|<pi
bessel_j_asympt_z :: (Value v) => v -> v -> v
bessel_j_asympt_z !nu !z =
  let !chi = z - (nu/2 + 1/4)*pi
      !mu = 4*nu^2
  in (sf_sqrt(2/(pi*z))) * (asymp_p nu z)*(sf_cos chi) - (asymp_q nu z)*(sin chi)
  where

    asymp_p !nu !z =
      let !term = 1.0
          !res = term
      in loop 1 term res
      where
        !mu = 4*nu^2
        !z8 = -(8*z)^2
        loop !k !t !r =
          let !t' = t * (mu-((#)$2*k-1)^2) * (mu-((#)$2*k+1)^2) / (((#)$2*k-1)*((#)$2*k)*z8)
              !r' = r + t'
          in if r==r' || (rabs t)>(rabs t') then r else loop (k+1) t' r'

    asymp_q !nu !z =
      let !term = (mu-1)/(8*z)
          !res = term
      in loop 2 term res
      where
        !mu = 4*nu^2
        !z8 = -(8*z)^2
        loop !k !t !r =
          let !t' = t * (mu-((#)$2*k-1)^2) * (mu-((#)$2*k+1)^2) / (((#)$2*k-2)*((#)$2*k-1)*z8)
              !r' = r + t'
          in if r==r' || (rabs t)>(rabs t') then r else loop (k+1) t' r'

-- recursion in order (backwards)
bessel_j_recur_back :: (Value v) => Double -> v -> v
bessel_j_recur_back !nu !z =
  let !jjs = runback (nnx-2) [1.0,0.0]
      !scale = if (rabs z)<10 then (bessel_j_series nuf z) else (bessel_j_asympt_z nuf z)
  in jjs!!(nnn) * scale / (jjs!!0)
  where
    !nnn = truncate nu
    !nuf = nu - (#)nnn
    !nnx = nnn + 10
    runback :: Int -> [v] -> [v]
    runback !0 !j = j
    runback !nx !j@(jj1:jj2:jjs) =
      let !jj = jj1*2*(nuf+(#)j)/z - jj2
      in runback (nx-1) (jj:j)

-- recursion in order (forewards)
bessel_j_recur_fore :: (Value v) => Double -> v -> v
bessel_j_recur_fore !nu !z =
  let !jj1 = bessel_j_series nuf z
      !jj2 = bessel_j_series (nuf+1) z
  in loop 3 jj1 jj2
  where
    !nnn = truncate nu
    !nuf = nu - (#)nnn
    !nnx = nnn + 10
    loop :: Int -> v -> v -> v
    loop j jjm2 jjm1
      | j==(nnx+1) = jjm1
      | otherwise =
          let jjj = jjm1*2*(fromDouble(nuf+(#)j))/z - jjm2
          in loop (j+1) jjm1 jjj

{--
function res = recur_backwards(n, z, topper)
  jjp2 = zeros(size(z));
  jjp1 = ones(size(z));
  jjp2_e_ = 1e-40 * ones(size(z));
  jjp1_e_ = 1e-20 * ones(size(z));
  scale = 2*ones(size(z));
  res = zeros(size(z));
  for m = (topper-2):(-1):1
    #jj(m) = (2*nu/z)*jj(m+1) - jj(m+2);
    s_ = -jjp2;
    e_ = -jjp2_e_;
    # add high
      t_ = s_;
      y_ = ((2*m./z).*jjp1) + e_;
      s_ = t_ + y_;
      e_ = (t_ - s_) + y_;
    # add low
      t_ = s_;
      y_ = ((2*m./z).*jjp1_e_) + e_;
      s_ = t_ + y_;
      e_ = (t_ - s_) + y_;
    jjp2 = jjp1;
    jjp2_e_ = jjp1_e_;
    jjp1 = s_;
    jjp1_e_ = e_;

    if (m==n+1)
      # store the desired result,
      # but keep recursing to get scale factor
      res = jjp1;
    endif

    if (m!=1)
      scale += 2 * (s_.^2 + e_.^2 + 2*s_.*e_);
    else
      scale += 1 * (s_.^2 + e_.^2 + 2*s_.*e_);
    endif

    if (scale>1e20)
      jjp2 /= 1024;
      jjp2_e_ /= 1024;
      jjp1 /= 1024;
      jjp1_e_ /= 1024;
      res /= 1024;
      scale /= 1024^2;
    endif
  endfor
  res ./= sqrt(scale);
endfunction
--}
