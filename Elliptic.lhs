\section{Elliptic functions}

\subsection{Preamble}
\begin{code}
{-# Language BangPatterns #-}
module Elliptic where
import AGM
import Exp
import Trig
import Util
\end{code}

$2^{-2/3}$
\begin{code}
two23 :: Double
!two23 = 0.62996052494743658238
\end{code}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Elliptic integral of the first kind}

Assume that $1-\sin^2\phi, 1-k^2\sin^2\phi\in \mathbb{C}\setminus(-\infty,0]$
except that one of them may be 0.

The elliptic integral of the first kind is defined by
\[ F(\phi,k) = \int_0^\phi\frac{d\theta}{\sqrt{1-k^2\sin^2\theta}}
    = \int_0^{\sin\phi}\frac{dt}{\sqrt{1-t^2}\sqrt{1-k^2t^2}} \]

The complete integral is given by $\phi=\pi/2$:
\[ K(k) = F(\pi/2, k) = \]

\subsubsection{\tt sf\_elliptic\_k k}
Compute the complete elliptic integral of the first kind $K(k)$
To evaluate this, we use the AGM relation
\[ K(k) = \frac{\pi}{2 \agm(1,k')} \qquad \text{where $k'=\sqrt{1-k^2}$} \marginnote{$K(k)$}\]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_k k} = K(k)$}
\begin{code}
sf_elliptic_k :: Double -> Double
sf_elliptic_k k =
  let an = sf_agm' 1.0 (sf_sqrt $ 1.0-k^2)
  in pi/(2*an)
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_f phi k}
Compute the (incomplete) elliptic integral of the first kind $F(\phi,k)$.
To evaluate, we use an ascending Landen transformation:
\[ F(\phi,k) = \frac{2}{1+k} F(\phi_2, k_2)
    \qquad \text{where $k_2=\frac{2\sqrt{k}}{1+k}$ and $2\phi_2=\phi+\arcsin(k \sin\phi)$} \marginnote{$F(\phi,k)$}\]
Note that $0<k<1$ and $0<\phi\leq\pi/2$ imply $k<k_2<1$ and $\phi_2<\phi$.
We iterate this transformation until we reach $k=1$ and use the special case
\[ F(\phi,1) = \gud^{-1}(\phi) \]
(Where $\gud^{-1}(\phi)$ is the inverse Gudermannian function (TODO)).
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_f phi k} = F(\phi,k)$}
\begin{code}
sf_elliptic_f :: Double -> Double -> Double
sf_elliptic_f phi k 
  | k==0 = phi
  | k==1 = sf_log((1 + (sf_sin phi)) / (1 - (sf_sin phi))) / 2
           -- quad(@(t)(1/sqrt(1-k^2*sin(t)^2)), 0, phi)
  | phi==0 = 0
  | otherwise =
      ascending_landen phi k 1 $ \ phi' res' ->
        res' * sf_log((1 + (sf_sin phi)) / (1 - (sf_sin phi))) / 2
  where 
    ascending_landen phi k res kont =
      let k' = 2 * (sf_sqrt k) / (1 + k)
          phi' = (phi + (asin (k*(sin phi))))/2
          res' = res * 2/(1+k)
      in if k'==1 then kont phi' res
         else ascending_landen phi' k' res' kont
    --function res = agm_method(phi, k)
    --  [an,bn,cn,phin] = sf_agm(1.0, sqrt(1 - k^2), phi, k);
    --  res = phin(end) / (2^(length(phin)-1) * an(end));
    --endfunction
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Elliptic integral of the second kind}

Assume that $1-\sin^2\phi, 1-k^2\sin^2\phi\in \mathbb{C}\setminus(-\infty,0]$
except that one of them may be 0.

Legendre's (incomplete) elliptic integral of the second kind is defined via
\[ E(\phi, k) = \int_0^\phi \sqrt{1-k^2\sin^2\theta}\,d\theta
    = \int_0^{\sin\phi} \frac{\sqrt{1-k^2t^2}}{\sqrt{1-t^2}} \,dt \]

The complete integral is
\[ E(k) = E(\pi/2, k) = \]

\subsubsection{\tt sf\_elliptic\_e k}
Compute the complete elliptic integral of the second kind $E(k)$.
We evaluate this with an agm-based approach:
\[ ... \]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_e k} = E(k)$}
\begin{code}
sf_elliptic_e :: Double -> Double
sf_elliptic_e k =
  let phi = k
      (as,bs,cs') = sf_agm 1.0 (sf_sqrt (1.0 - k^20))
      cs = k:(tail.reverse$cs')
      res = foldl (-) 2 (map (\(i,c)->2^(i-1)*c^2) (zip [1..] cs))
  in res * pi/(4*(head as))
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_e\_ic phi k}
Compute the incomplete elliptic integral of the second kind $E(\phi, k)$
We evaluate this with an ascending Landen transformation:
\[ ... \]
TODO: UNTESTED!
(Note: could try direct quadrature of the integral, also there is an AGM-based method).
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_e\_ic phi k} = E(\phi,k)$}
\begin{code}
sf_elliptic_e_ic :: Double -> Double -> Double
sf_elliptic_e_ic phi k 
  | k==1 = sf_sin phi
  | k==0 = phi
  | otherwise = ascending_landen phi k
  where
    ascending_landen phi 1 = sin phi
    ascending_landen phi k =
      let !k' = 2*(sf_sqrt k) / (k+1)
          !phi' = (phi + (sf_asin (k*(sf_sin phi))))/2
      in (1+k)*(ascending_landen phi' k') + (1-k)*(sf_elliptic_f phi' k') - k*(sf_sin phi)
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Elliptic integral of the third kind}

We define Legendre's (incomplete) elliptic integral of the third kind via
\[ \Pi(\phi,\alpha^2,k) = \int_0^\phi\frac{d\theta}{\sqrt{1-k^2\sin^2\theta}(1-\alpha^2\sin^2\theta)}
    = \int_0^{\sin\phi} \frac{dt}{\sqrt{1-t^2}\sqrt{1-k^2t^2}(1-\alpha^2t^2)} \]

The complete integral of the third kind is given by 
\[ \Pi(\alpha^2,k) = \Pi(\pi/2,\alpha^2,k) = \]

\subsubsection{\tt sf\_elliptic\_pi c k}
Compute the complete elliptic integral of the third kind
($c=\alpha^2$ in DLMF notation)
for real values only $0<k<1$, $0<c<1$.
Uses agm-based approach.
(Could also try numerical quadrature \verb|quad(@(t)(1.0/(1-c*sf_sin(t)^2)/sqrt(1.0 - k^2*sf_sin(t)^2)), 0, phi)|).
TODO: mostly untested
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_pi c k} = \Pi(c,k)$}
\begin{code}
sf_elliptic_pi :: Double -> Double -> Double
sf_elliptic_pi c k = complete_agm k c
  where
    -- -\infty < k^2 < 1
    -- -\infty < c < 1
    complete_agm k c = 
      let (ans,gns,_) = sf_agm 1 (sf_sqrt (1.0-k^2))
          pn1 = sf_sqrt (1-c)
          qn1 = 1
          an1 = last ans
          gn1 = last gns
          en1 = (pn1^2 - an1*gn1) / (pn1^2 + an1*gn1)
      in iter pn1 en1 (reverse ans) (reverse gns) [qn1]

    iter pnm1 enm1 [an] [gn] qns = pi/(4*an) * (2 + c/(1-c)*(ksum qns))
    iter pnm1 enm1 (anm1:an:ans) (gnm1:gn:gns) (qnm1:qns) =
      let pn = (pnm1^2 + anm1*gnm1)/(2*pnm1)
          en = (pn^2 - an*gn) / (pn^2 + an*gn)
          qn = qnm1 * enm1/2
      in iter pn en (an:ans) (gn:gns) (qn:qnm1:qns)
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_pi\_ic phi c k}
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_pi\_ic phi c k} = \Pi(\phi,c,k)$}
\begin{code}
sf_elliptic_pi_ic :: Double -> Double -> Double -> Double
sf_elliptic_pi_ic 0 c k = 0.0
sf_elliptic_pi_ic phi c k = gauss_transform k c phi
  where
    gauss_transform k c phi =
      if (sf_sqrt (1-k^2))==1
      then let cp=sf_sqrt(1-c)
           in sf_atan(cp*(sf_tan phi)) / cp
      else if (1-k^2/c)==0 -- special case else rho below is zero...
      then ((sf_elliptic_e_ic phi k) - c*(sf_cos phi)*(sf_sin phi)
                / sqrt(1-c*(sf_sin phi)^2))/(1-c)
      else let kp = sf_sqrt (1-k^2)
               k' = (1 - kp) / (1 + kp)
               delta = sf_sqrt(1-k^2*(sf_sin phi)^2)
               psi' = sf_asin((1+kp)*(sf_sin phi) / (1+delta))
               rho = sf_sqrt(1 - (k^2/c))
               c' = c*(1+rho)^2/(1+kp)^2
               xi = (sf_csc phi)^2
               newgt = gauss_transform k' c' psi'
           in (4/(1+kp)*newgt + (rho-1)*(sf_elliptic_f phi k)
                - (sf_elliptic_rc (xi-1) (xi-c)))/rho
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Elliptic integral of Legendre's type}

The (incomplete) elliptic integral of Legendre's type is defined by
\[ D(\phi,k) = \int_0^\phi \frac{\sin^2\theta}{\sqrt{1-k^2\sin^2\theta}}\,d\theta
    = \int_0^{\sin\phi}\frac{t^2}{\sqrt{1-t^2}\sqrt{1-k^2t^2}}\,dt \]
This can be expressed as $D(\phi,k) = (F(\phi,k)-E(\phi,k))/k^2$.

The complete elliptic integral of Legendre's type is
\[ D(k) = D(\pi/2,k) = (K(k)-E(k))/k^2 \]

\subsubsection{\tt sf\_elliptic\_d\_ic phi k}
We simply reduce to $F(\phi,k)$ and $E(\phi,k)$.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_d\_ic phi k} = D(\phi,k)$}
\begin{code}
sf_elliptic_d_ic :: Double -> Double -> Double
sf_elliptic_d_ic phi k = ((sf_elliptic_f phi k) - (sf_elliptic_e_ic phi k)) / (k^2)
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_d\_ic phi k}
We simply reduce to $K(k)$ and $E(k)$.
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_d k} = D(k)$}
\begin{code}
sf_elliptic_d :: Double -> Double
sf_elliptic_d k = ((sf_elliptic_k k) - (sf_elliptic_e k)) / (k^2)
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Burlisch's elliptic integrals}

DLMF: ``Bulirsch’s integrals are linear combinations of Legendre’s integrals that are chosen to facilitate computational application of Bartky’s transformation''

\subsubsection{\tt sf\_elliptic\_cel kc p a b}
Compute Burlisch's elliptic integral where $p\neq0$, $k_c\neq0$.
\[ cel(k_c,p,a,b) = \int_0^{\pi/2} \frac{a\cos^2\theta + b\sin^2\theta}{\cos^2\theta+p\sin^2\theta}
        \frac{1}{\sqrt{\cos^2\theta + k_c^2\sin^2\theta}} \,d\theta \marginnote{$cel(k_c,p,a,b)$}\]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_cel kc p a b} = cel(k_c,p,a,b)$}
\begin{code}
sf_elliptic_cel :: Double -> Double -> Double -> Double -> Double
sf_elliptic_cel kc p a b = a * (sf_elliptic_rf 0 (kc^2) 1) + (b-p*a)/3 * (sf_elliptic_rj 0 (kc^2) 1 p)
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_el1 x kc}
Compute Burlisch's elliptic integral
\[ el_1(x,k_c) =\]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_el1 k kc} = el_1(x,k_c)$}
\begin{code}
sf_elliptic_el1 :: Double -> Double -> Double
sf_elliptic_el1 x kc =
  --sf_elliptic_f (atan x) (sf_sqrt(1-kc^2))
  let r = 1/x^2
  in sf_elliptic_rf r (r+kc^2) (r+1)
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_el2 x kc a b}
Compute Burlisch's elliptic integral
\[ el_2(x,k_c,a,b) = \int_0^{\arctan x}\frac{a + b\tan^2\theta}{\sqrt{(1+\tan^2\theta)(1+k_c^2\tan^2\theta)}} \,d\theta \]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_el2 x kc a b} = el_2(x,k_c,a,b)$}
\begin{code}
sf_elliptic_el2 :: Double -> Double -> Double -> Double -> Double
sf_elliptic_el2 x kc a b =
  let r = 1/x^2
  in a * (sf_elliptic_el1 x kc) + (b-a)/3 * (sf_elliptic_rd r  (r+kc^2) (r+1))
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_el3 x kc p}
Compute the Burlisch's elliptic integral
\[ el_3(x,k_c,p) = \int_0^{\arctan x} \frac{d\theta}{(\cos^2\theta + p\sin^2\theta)\sqrt{\cos^2\theta + k_c^2\sin^2\theta}} \]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_el3 x kc p} = el_3(x,k_c,p)$}
\begin{code}
sf_elliptic_el3 :: Double -> Double -> Double -> Double
sf_elliptic_el3 x kc p =
  -- sf_elliptic_pi(atan(x), 1-p, sf_sqrt(1-kc.^2));
  let r = 1/x^2
  in (sf_elliptic_el1 x kc) + (1-p)/3 * (sf_elliptic_rj r (r+kc^2) (r+1) (r+p))
\end{code}
\end{titled-frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Symmetric elliptic integrals}

\subsubsection{\tt sf\_elliptic\_rc x y}
Compute the symmetric elliptic integral $R_C(x,y)$ for real parameters.
Let $x\in\mathbb{C}\setminus(-\infty,0)$, $y\in\mathbb{C}\setminus\{0\}$, then we define
\[ R_C(x,y) = \frac12 \int_0^\infty \frac{dt}{\sqrt{t+x}(t+y)} \]
(where the Cauchy principal value is taken if $y<0$.)
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_rc x y} = R_C(x,y)$}
\begin{code}
-- x>=0, y!=0
sf_elliptic_rc :: Double -> Double -> Double
sf_elliptic_rc x y
  | 0==x && x<y  = 1/sf_sqrt(y-x) * sf_acos(sf_sqrt(x/y))
  | 0<x  && x<y  = 1/sf_sqrt(y-x) * sf_atan(sf_sqrt((y-x)/x))
  | 0<y  && y<x  = 1/sf_sqrt(x-y) * sf_atanh(sf_sqrt((x-y)/x))
                 -- = 1/sf_sqrt(x-y) * sf_log((sf_sqrt(x) + sf_sqrt(x-y))/sf_sqrt(y))
  | y<0  && 0<=x = 1/sf_sqrt(x-y) * sf_log((sf_sqrt(x)+sf_sqrt(x-y))/sf_sqrt(-y))
                 -- = 1/sf_sqrt(x-y) * sf_atanh(sf_sqrt(x/(x-y)))
                 -- = sf_sqrt(x/(x-y)) * (sf_elliptic_rc (x-y) (-y))
  | x == y       = 1/(sf_sqrt x)
  | otherwise    = error "sf_elliptic_rc: domain error"
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_rd x y z}
Compute the symmetric elliptic integral $R_D(x,y,z)$
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_rc x y z} = R_D(x,y,z)$}
\begin{code}
-- x,y,z>0
sf_elliptic_rd :: Double -> Double -> Double -> Double
sf_elliptic_rd x y z = let (x',s) = (iter x y z 0.0) in (x'**(-3/2) + s)
  where
    iter x y z s =
      let lam = sf_sqrt(x*y) + sf_sqrt(y*z) + sf_sqrt(z*x);
          s' = s + 3/sf_sqrt(z)/(z+lam);
          x' = (x+lam)*two23
          y' = (y+lam)*two23
          z' = (z+lam)*two23
          mu = (x+y+z)/3;
          eps = foldl1 max (map (\t->abs(1-t/mu)) [x,y,z])
      in if eps<2e-16 || [x,y,z]==[x',y',z'] then (x',s')
         else iter x' y' z' s'
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_rf x y z}
Compute the symmetric elliptic integral of the first kind
\[ R_F(x,y,z) = \frac12 \int_0^\infty \frac{dt}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}} \]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_rf x y z} = R_F(x,y,z)$}
\begin{code}
-- x,y,z>0
sf_elliptic_rf :: Double -> Double -> Double -> Double
sf_elliptic_rf x y z = 1/(sf_sqrt $ iter x y z)
  where
    iter x y z =
      let lam = (sf_sqrt $ x*y) + (sf_sqrt $ y*z) + (sf_sqrt $ z*x)
          mu = (x+y+z)/3
          eps = foldl1 max $ map (\a->abs(1-a/mu)) [x,y,z]
          x' = (x+lam)/4
          y' = (y+lam)/4
          z' = (z+lam)/4
      in if (eps<1e-16) || ([x,y,z]==[x',y',z'])
         then x
         else iter x' y' z'
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_rg x y z}
Compute the symmetric elliptic integral
\[ R_G(x,y,z) = \frac{1}{4\pi} \int_0^{2\pi}\int_0^\pi
    \sqrt{x \sin^2\theta \cos^2\phi + y\sin^2\theta \sin^2\phi + z\cos^2\theta}\sin\theta \,d\theta\, d\phi \]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_rg x y z} = R_G(x,y,z)$}
\begin{code}
-- x,y,z>0
sf_elliptic_rg :: Double -> Double -> Double -> Double
sf_elliptic_rg x y z
  | x>y = sf_elliptic_rg y x z
  | x>z = sf_elliptic_rg z y x
  | y>z = sf_elliptic_rg x z y
  | otherwise = 
    let !a0 = sqrt (z-x)
        !c0 = sqrt (y-x)
        !h0 = sqrt z
        !t0 = sqrt x
        !(an,tn,cn_sum,hn_sum) = iter 1 a0 t0 c0 (c0^2/2) h0 0
    in ((t0^2 + theta*cn_sum)*(sf_elliptic_rc (tn^2+theta*an^2) tn^2) + h0 + hn_sum)/2
    where
      theta = 1
      iter n an tn cn cn_sum hn hn_sum =
        let an' = (an + sf_sqrt(an^2 - cn^2))/2
            tn' = (tn + sf_sqrt(tn^2 + theta*cn^2))/2
            cn' = cn^2/(2*an')/2
            cn_sum' = cn_sum + 2^((#)n-1)*cn'^2
            hn' = hn*tn'/sf_sqrt(tn'^2+theta*cn'^2)
            hn_sum' = hn_sum + 2^n*(hn' - hn)
            n' = n + 1
        in if cn^2==0 then (an,tn,cn_sum,hn_sum)
           else iter n' an' tn' cn' hn_sum' hn' hn_sum'
\end{code}
\end{titled-frame}

\subsubsection{\tt sf\_elliptic\_rj x y z p}
Compute the symmetric elliptic integral
\[ R_J(x,y,z,p) = \frac32 \int_0^\infty \frac{dt}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)} \]
TODO: UNTESTED!
\begin{titled-frame}{$\text{\color{blue}\tt sf\_elliptic\_rj x y z p} = R_J(x,y,z,p)$}
\begin{code}
-- x,y,z>0
sf_elliptic_rj :: Double -> Double -> Double -> Double -> Double
sf_elliptic_rj x y z p =
  let (x',smm,scale) = iter x y z p 0.0 1.0
  in scale*x'**(-3/2) + smm
  where 
    iter x y z p smm scale =
      let lam = sf_sqrt(x*y) + sf_sqrt(y*z) + sf_sqrt(z*x)
          alpha = p*(sf_sqrt(x)+sf_sqrt(y)+sf_sqrt(z)) + sf_sqrt(x*y*z)
          beta = sf_sqrt(p)*(p+lam)
          smm' = smm + (if (abs(1 - alpha^2/beta^2) < 5e-16)
                 then
                   -- optimization to reduce external calls
                   scale*3/alpha;
                 else
                   scale*3*(sf_elliptic_rc (alpha^2) (beta^2))
                 )
          mu = (x+y+z+p)/4
          eps = foldl1 max (map (\t->abs(1-t/mu)) [x,y,z,p])
          x' = (x+lam)*two23/mu
          y' = (y+lam)*two23/mu
          z' = (z+lam)*two23/mu
          p' = (p+lam)*two23/mu
          scale' = scale * (mu**(-3/2))
      in if eps<1e-16 || [x,y,z,p]==[x',y',z',p'] || smm'==smm
         then (x',smm',scale')
         else iter x' y' z' p' smm' scale'
\end{code}
\end{titled-frame}


