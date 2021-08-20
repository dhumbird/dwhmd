c-----------------------------------------------------------------------
c     spoof does many of the same things as pibond.  however, it is
c     called for a LJ pair which has been fictitiously placed at
c     REBO range.  the routine calculates b_ij, uses this to determine
c     a spline which is used to multiply V^LJ, and adds various LJ
c     forces.  db/dr derivatives are multiplied by V^LJ and some
c     splines, rather than V^A, in order to get the required forces.
c     ...sjs 12/29/97
c-----------------------------------------------------------------------
c     when making changes to this routine, make corresponding changes
c     to pibond().
c-----------------------------------------------------------------------
c
      subroutine spoof(ijpair, btot, rsplin, fsplin, vdw, llje)
     
      include 'common_files.inc'

      logical llje

      dimension uutemp(ndim,ndim)

c     should be smaller:  max number of neighbors of a single atom

      parameter (ncc = 20)
c    
c     fxik  = degree of unsaturation of atom k.  [0,1].  F(x_ik) (*)
c     fxjl  = degree of unsaturation of atom l.  [0,1].  F(x_jl) (*)
c     dfxik = d(F(x_ik))/d(x_ik) (*)
c     dfxjl = d(F(x_jl))/d(x_jl) (*)
c     rjk    = r_jk (spoofed) (*)
c     ril    = r_il (spoofed) (*)
c     dgdik  = 1/r_ik d/dr_ik ( f^c G exp(lambda) ) in eq. 7 (*)

c     dgdjk  = 1/r_jk d/dr_jk ( f^c G exp(lambda) ) in eq. 7 (*)
c     xsil   = 1/r_il d/dr_il ( f^c G exp(lambda) ) in eq. 7 (*)
c     xsjl   = 1/r_jl d/dr_jl ( f^c G exp(lambda) ) in eq. 7 (*)
c     cosk   = cos of theta_jik (*)
c     cosl   = cos of theta_ijl (*)
c     sink   = sin of theta_jik (*)
c     sinl   = sin of theta_ijl (*)
c     dctjk  = 1/r_jk d(cos(theta_jik))/d(r_jk) (*)
c     dctij  = 1/r_ij d(cos(theta_jik))/d(r_ij) = 0 (*)
c     dctik  = 1/r_ik d(cos(theta_jik))/d(r_ik) (*)
c     dctji  = 1/r_ji d(cos(theta_ijl))/d(r_ji) = 0 (*)
c     dctil  = 1/r_il d(cos(theta_ijl))/d(r_il) (*)
c     dctjl  = 1/r_jl d(cos(theta_ijl))/d(r_jl) (*)
c     rijvec = r_ij (spoofed)
c     rikvec = r_ik (unchanged)
c     rjlvec = r_jl (unchanged)
c     rlrijv = r_ij (original)
c     dpijdn  = d(P_ij)/d(N^type_i)
c     dpjidn  = d(P_ji)/d(N^type_j)
c     bnijp1 = neighbors of type ARG of atom i, plus 1, not counting j.
c     bnijp1(icarb) = N^C_i + 1.  bnijp1(ihyd) = N^H_i + 1. (*)
c     bnjip1 = neighbors of type ARG of atom j, plus 1, not counting i.
c     bnjip1(icarb) = N^C_j + 1.  bnjip1(ihyd) = N^H_j + 1. (*)
c     (*): k is ARG-th non-j neighbor of i.  l is ARG-th non-i
c     neighbor of j. j is jn.

      dimension rjk(ncc,3), ril(ncc,3),
     .     tnpk(ncc,3), tnpl(ncc,3),
     .     fxik(ncc), fxjl(ncc), dfxik(ncc), dfxjl(ncc),
     .     dgdik(ncc), dgdjk(ncc), xsil(ncc), xsjl(ncc),
     .     cosk(ncc), cosl(ncc), sink(ncc), sinl(ncc),
     .     dctik(ncc), dctjk(ncc), dctil(ncc),
     .     dctjl(ncc),
     .     tspjik(ncc), dtsjik(ncc), tspijl(ncc), dtsijl(ncc),
     .     rijvec(3), rikvec(3), rjlvec(3), rlrijv(3),
     .     dt2dik(3), dt2djl(3),
     .     tnpi(3), tnpj(3),
     .     dpijdn(2), dpjidn(2), bnijp1(2), bnjip1(2)

c     initialize some stuff

      btot = 0

      iatom = ihalf(ijpair)
      jatom = jhalf(ijpair)

      itype = iat2ty(iatom)
      jtype = iat2ty(jatom)
      ijtype = ijty(itype,jtype)

      ibegin = nabors(iatom)
      iend = nabors(iatom+1) - 1

      veelj = rsplin * fsplin * vdw

      do 80 idim = 1, 3
         do 70 jdim = 1, 3
            uutemp(jdim,idim) = 0.d0
 70      continue
 80   continue

      do 100 idim = 1, 3
         rijvec(idim) = rijv(ijpair,idim)
         tnpi(idim) = 0.d0
         tnpj(idim) = 0.d0
         do 90 icc = 1, ncc
            tnpk(icc,idim) = 0.d0
            tnpl(icc,idim) = 0.d0
 90      continue
 100  continue

c     this is the spoofed pair distance:

      rijloc = rcor(ijpair)
      rsqij = rijloc * rijloc

c     we will also need the real pair distance:

      rsq = 0.d0
      do 110 idim = 1, 3
         rlrijv(idim) = r0(iatom,idim) - r0(jatom,idim)
         if (lpbc(idim)) then
            rlrijv(idim) = rlrijv(idim)
     .           - cube(idim) * anint(rlrijv(idim) / cube(idim))
         endif
         rsq = rsq + rlrijv(idim) ** 2
 110  continue

      rlrij2 = rsq
      rlrij = sqrt(rlrij2)
      rijmbr = dijmin(itype,jtype) / rlrij

c     I side of bond

c     nk = number of non-jn neighbors of i
c     sdgdij = 1/r_ij d/dr_ij of sum of f_ik g^c exp(lambda_jik) = 0
c     sfge2l = sum of f_ik g^c exp(lambda_jik) terms
c     cni = sum of delta_k^C f_ik S'(F(N^t_ki)) terms

      nk=0
      sfge2l = 0.0d0
      cni = 0.0d0
     
c     bnijp1 is the number of neighbors of the given type of atom i, plus
c     one, not counting the neighbor jn.  bnijp1(icarb) = 1 + N^C_ij and
c     bnijp1(ihyd) = 1 + N^H_ij.

      bnijp1(icarb) = rnpls1(iatom,icarb)
      bnijp1(ihyd) = rnpls1(iatom,ihyd)
      bnijp1(jtype) = bnijp1(jtype) - fij(ijpair)

c     bgntij = number of C/H neighbors of atom i, not counting the
c     neighbor jn.  bgntij = N^t_ij.

      bgntij = bnijp1(icarb) + bnijp1(ihyd) - 2.d0
      sdgdnt = 0.0d0

c     loop over non-j neighbors of i

      do 1000 ikpair = ibegin, iend

c     only look at hydrocarbon pairs

         if (iprtyp(ikpair) .ne. ihycar) then
            go to 1000
         endif

         katom = jhalf(ikpair)

c        ignore jatom

         if (katom .eq. jatom) then
            go to 1000
         endif

         qspli = 0.0d0
         dali=0.0d0
         dgdnti=0.0d0

         ktype = iat2ty(katom)
         nk=nk+1
         rik = rcor(ikpair)
         rik2 = rik * rik

c     unrolled for speed

         rjk(nk,1) = rijv(ikpair,1) - rijvec(1)
         rjk(nk,2) = rijv(ikpair,2) - rijvec(2)
         rjk(nk,3) = rijv(ikpair,3) - rijvec(3)
         rjk2 = rjk(nk,1) * rjk(nk,1) + rjk(nk,2) * rjk(nk,2)
     .        + rjk(nk,3) * rjk(nk,3)
         rijrik = 2.0d0 * rijloc * rik
         rr = rsqij - rik2
         costh = (rsqij + rik2 - rjk2) / rijrik

c        seemingly impossible, but happens due to roundoff error:

         if (costh .lt. -1.0d0) then
            costh = -1.0d0
         else if (costh .gt. 1.d0) then
            costh = 1.d0
         endif

c     tsplin is the cubic spline switching function that turns off the
c     torsion interactions at nearly linear bond angles.
c     S(t_theta(cos(theta))).  dtspln is -dS/d(cos theta).

c     (C=C double bonds will rarely have near-linear bond angles at low
c     temperatures, but they begin to become relevant for ethylene at
c     ~2500 K.)

         if (costh .gt. tthmax) then
            tsplin = 0.d0
            dtspln = 0.d0
         else if (costh .gt. tthmin) then
            dtthet = costh - tthmin
            twidth = tthmax - tthmin
            tee = dtthet / twidth
            tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
         else
            tsplin = 1.d0
            dtspln = 0.d0
         endif

         cosk(nk)=costh
         sink(nk)=sqrt(1.0d0-costh*costh)
         tspjik(nk) = tsplin
         dtsjik(nk) = dtspln

c        theta_jik          cos(theta_jik)    ig
c        [0     ,109.47) => (-1/3,1   ]    => 4
c        [109.47,120   ) => (-1/2,-1/3]    => 3
c        [120   ,131.81) => (-2/3,-1/2]    => 2
c        [131.81,180   ] => [-1  ,-2/3]    => 1

c        move inside if:  slower!?
         ig=igc(int(-costh*12.0d0)+13)

c        three-body angles around carbon:

         if (itype .eq. icarb) then

c        for 109.47 <= theta <= 180:

            if(ig.ne.4) then

c        smallg is g_C from eq. 7, constructed from a 5th-order spline in
c        cos(theta).  dgdcos is d(g_C)/d(cos(theta))

c        unrolled for speed

               smallg = spgc(1,ig)
     .              + costh * (spgc(2,ig)
     .              + costh * (spgc(3,ig)
     .              + costh * (spgc(4,ig)
     .              + costh * (spgc(5,ig)
     .              + costh *  spgc(6,ig)))))
               dgdcos = spgc(2,ig)
     .              + costh * (2.d0 * spgc(3,ig)
     .              + costh * (3.d0 * spgc(4,ig)
     .              + costh * (4.d0 * spgc(5,ig)
     .              + costh *  5.d0 * spgc(6,ig))))

c        for 0 <= theta < 109.47:

            else

c        unnecessary: (but faster!?!)
               dali=0.0d0

c        bgqmin and bgqmax are the lower and upper bounds of the switching
c        region that changes from G_C to gamma_C based on the coordination
c        number of carbons with a 3-body angle of less than 109.47 degrees.

               if (bgntij .lt. bgqmax) then

c        qspli is Q_i(N^t_i) from eq. 10, the switching function between G_C
c        and gamma_C.  qspli = 0 means G_C is fully on. 
c        dali is d(Q_i(N^t_i))/d(N^t_i).

c        smallg contains g_C, the splined combination of G_C and gamma_C.
c        before being splined together, smallg contains G_C and gammac
c        contains gamma_C.  dgdcos contains d(g_C)/d(cos(theta)).  before
c        being splined together, dgdcos contains d(G_C)/d(cos(theta)) and
c        dgamdc contains d(gamma_C)/d(cos(theta)).  As N^t_i passes from
c        3.2 to 3.7 (CN of 4.2 to 4.7) g_C starts looking more like G_C and
c        less like gamma_C.

c        dgdnti = d(g_C)/dN^t_ij

c        the calculation of smallg is unrolled from the previous version
c        to avoid exponentiation.  it is also inserted into the nested ifs
c        to avoid calculating pieces that aren't needed

                  qspli = 1.0d0
                  if (bgntij .gt. bgqmin) then
c     in this branch 3.2 < bgntij (= N^t_ij) < 3.7, so qspli, dali, and
c     dgdnti must be calculated, and both G_C and gamma_C are needed
                     bigqnt = pibydq * (bgntij - bgqmin)
                     qspli = 0.5d0 * (1.d0 + cos(bigqnt))
                     dali = -pibydq * 0.5d0 * sin(bigqnt)
                     smallg = spgc(1,ig)
     .                    + costh * (spgc(2,ig)
     .                    + costh * (spgc(3,ig)
     .                    + costh * (spgc(4,ig)
     .                    + costh * (spgc(5,ig)
     .                    + costh *  spgc(6,ig)))))
                     dgdcos = spgc(2,ig)
     .                    + costh * (2.d0 * spgc(3,ig)
     .                    + costh * (3.d0 * spgc(4,ig)
     .                    + costh * (4.d0 * spgc(5,ig)
     .                    + costh *  5.d0 * spgc(6,ig))))
                     ig1 = ig + 1
                     gammac = spgc(1,ig1)
     .                    + costh * (spgc(2,ig1)
     .                    + costh * (spgc(3,ig1)
     .                    + costh * (spgc(4,ig1)
     .                    + costh * (spgc(5,ig1)
     .                    + costh *  spgc(6,ig1)))))
                     dgamdc = spgc(2,ig1)
     .                    + costh * (2.d0 * spgc(3,ig1)
     .                    + costh * (3.d0 * spgc(4,ig1)
     .                    + costh * (4.d0 * spgc(5,ig1)
     .                    + costh *  5.d0 * spgc(6,ig1))))
                     dgdnti = dali * (gammac - smallg)
                     smallg = smallg + qspli * (gammac - smallg)
                     dgdcos = dgdcos +
     .                    qspli * (dgamdc - dgdcos)
                  else
c     in this branch bgntij (= N^t_ij) <= 3.2, so qspli=1, dali = 0,
c     dgdnti = 0, and only gamma_C is needed
                     ig1 = ig + 1
                     gammac = spgc(1,ig1)
     .                    + costh * (spgc(2,ig1)
     .                    + costh * (spgc(3,ig1)
     .                    + costh * (spgc(4,ig1)
     .                    + costh * (spgc(5,ig1)
     .                    + costh *  spgc(6,ig1)))))
                     dgamdc = spgc(2,ig1)
     .                    + costh * (2.d0 * spgc(3,ig1)
     .                    + costh * (3.d0 * spgc(4,ig1)
     .                    + costh * (4.d0 * spgc(5,ig1)
     .                    + costh *  5.d0 * spgc(6,ig1))))
                     smallg = gammac
                     dgdcos = dgamdc
                  endif
               else
c     in this branch bgntij (= N^t_ij) is >= 3.7, so qspli=dali=dgdnti=0, and
c     only G_C is needed
                  smallg = spgc(1,ig)
     .                 + costh * (spgc(2,ig)
     .                 + costh * (spgc(3,ig)
     .                 + costh * (spgc(4,ig)
     .                 + costh * (spgc(5,ig)
     .                 + costh *  spgc(6,ig)))))
                  dgdcos = spgc(2,ig)
     .                 + costh * (2.d0 * spgc(3,ig)
     .                 + costh * (3.d0 * spgc(4,ig)
     .                 + costh * (4.d0 * spgc(5,ig)
     .                 + costh *  5.d0 * spgc(6,ig))))
               endif
            endif
         else

c     three-body angles around hydrogen:

c           theta_jik          cos(theta_jik    ig
c           [0     ,120   ) => (-1/2,1   ]   => 3
c           [120   ,146.44) => (-5/6,-1/2]   => 2
c           [146.44,180   ] => [-1  ,-5/6]   => 1

            ig=igh(int(-costh*12.0d0)+13)

c     smallg is g_H.  dgdcos is d(g_H)/d(cos(theta)).

c     unrolled for speed

                  smallg = spgh(1,ig)
     .                 + costh * (spgh(2,ig)
     .                 + costh * (spgh(3,ig)
     .                 + costh * (spgh(4,ig)
     .                 + costh * (spgh(5,ig)
     .                 + costh *  spgh(6,ig)))))
                  dgdcos = spgh(2,ig)
     .                 + costh * (2.d0 * spgh(3,ig)
     .                 + costh * (3.d0 * spgh(4,ig)
     .                 + costh * (4.d0 * spgh(5,ig)
     .                 + costh *  5.d0 * spgh(6,ig))))
         endif
         fik = fij(ikpair)
         dfik = dww(ikpair)
         fxik(nk)=0.0d0
         dfxik(nk)=0.0d0

c     for all carbon neighbors, check their coordination numbers to
c     look for conjugation.

         if (ktype .eq. icarb) then

c     bgntki is the coordination of atom katom excluding i = N^t_ki
c     N^t_k - f^c_ik(r_ik) according to eq. 15, but that's not right,
c     because N^t_k already excludes f^c_ik(r_ik).

            bgntki = rnpls1(katom,icarb) + rnpls1(katom,ihyd)
     .           - fik - 2.d0

c     fxik is the degree of unsaturation of atom katom.
c     dfxik is d(fxik)/d(bgntki).

            if (bgntki .lt. fmx) then
               if (bgntki .le. fmn) then
                  fxik(nk)=1.0d0
               else
                  px = pi * (bgntki - fmn)
                  fxik(nk) = 0.5d0 * (1.d0 + cos(px))
                  dfxik(nk) = -fik * 0.5d0 * sin(px) * pi
               endif
            endif
         endif

c     cni is the weighted sum of conjugation terms for all
c     neighbors of i except jn.  see eq. 13.

c     belongs in if above but SLOWER!?!:
         cni = cni + fik * fxik(nk)

c     e2lam is the exp(lambda_ijk) term in eq. 7.  It only shows up
c     when i is a hydrogen.  It gets smaller (< 1) as rij shrinks
c     and larger (> 1) as rik shrinks.  It is one when both bonds
c     are at their standard lengths: rhh = 1.09, rch = 0.7415886997.

         if (ell(itype,jtype,ktype) .ne. 0.d0) then
            e2lam = reg(itype,jtype,ktype) *
     .           exp(ell(itype,jtype,ktype) * (rijloc - rik))
         else
            e2lam=1.0d0
         endif

c     dctdjk = dctjk(nk) = 1/r_jk d(cos(theta_jik))/d(r_jk)
c     dctdij             = 1/r_ij d(cos(theta_jik))/d(r_ij) = 0
c     dctdik = dctik(nk) = 1/r_ik d(cos(theta_jik))/d(r_ik)
c     all are partial derivatives, not total

c     the d/d(r_ij) derivative is zero for this spoofed pair

         dctdjk = -2.0d0 / rijrik
         dctdik = (-rr + rjk2) / (rijrik * rik2)
        
         dctjk(nk) = dctdjk
         dctik(nk) = dctdik

c     sdgdij = sum_k of 1/r_ij d/dr_ij (f^c G exp(lambda)) in eq. 7 = 0
c     sfge2l = sum over k of f^c G exp(lambda) terms in eq. 7
c     dgdik(nk) = 1/r_ik d/dr_ik (f^c G exp(lambda)) in eq. 7
c     sdgdnt = delta_i^C sum_k of f^c dQ/dN^t_i (gamma - G) exp(lambda)
c     dgdjk(nk) = 1/r_jk d/dr_jk (f^c G exp(lambda)) in eq. 7

c     the d/d(r_ij) derivative is zero for this spoofed pair

         ge2lam = smallg * e2lam
         sfge2l = sfge2l + fik * ge2lam
         xtemp = fik * e2lam * dgdcos
         gfx = ge2lam * fik * ell(itype,jtype,ktype)
         dgdik(nk) = (ge2lam * dfik - gfx) / rik + xtemp * dctdik
         sdgdnt = sdgdnt + e2lam * fik * dgdnti
         dgdjk(nk) = xtemp*dctdjk
 1000 continue

c     J side of bond

c     nl = number of non-i neighbors of j
c     xsji = 1/r_ij d/dr_ij of sum of f_jl g^c exp(lambda_ijl) = 0
c     ssuml = sum of f_jl g^c exp(lambda_ijl) terms
c     conl = sum of delta_k^C f_jl S'(F(N^t_lj)) terms

      nl=0
      ssuml=0.0d0
      conl=0.0d0

      jbegin = nabors(jatom)
      jend = nabors(jatom+1) - 1

c     bnjip1 is the number of neighbors of the given type of atom j, plus
c     one, not counting the neighbor i.  bnjip1(icarb) = 1 + N^C_j and
c     bnjip1(ihyd) = 1 + N^H_j.

      bnjip1(icarb) = rnpls1(jatom,icarb)
      bnjip1(ihyd) = rnpls1(jatom,ihyd)
      bnjip1(itype) = bnjip1(itype) - fij(ijpair)

c     bgntji = number of C/H neighbors of atom j, not counting the
c     neighbor i.  bgntji = N^t_ji.

      bgntji = bnjip1(icarb) + bnjip1(ihyd) - 2.d0
      sdaljl=0.0d0

c     loop over neighbors of jatom

      do 2000 jlpair = jbegin, jend
         qsplj=0.0d0
         dalj=0.0d0
         dgdntj = 0.0d0
         latom = jhalf(jlpair)

c     latom (neighbor of jatom) should not be iatom

         if (latom .eq. iatom) then
            go to 2000
         endif
         if (iprtyp(jlpair) .ne. ihycar) then
            go to 2000
         endif
         ltype = iat2ty(latom)
         nl=nl+1
         rjl = rcor(jlpair)
         rjl2 = rjl * rjl
         ril(nl,1) = rijv(jlpair,1) + rijvec(1)
         ril(nl,2) = rijv(jlpair,2) + rijvec(2)
         ril(nl,3) = rijv(jlpair,3) + rijvec(3)
         ril2 = ril(nl,1) * ril(nl,1) + ril(nl,2) * ril(nl,2)
     .        + ril(nl,3) * ril(nl,3)
        
         rijrjl = 2.0d0 * rijloc * rjl
         rr = rsqij - rjl2
         costh = (rsqij + rjl2 - ril2) / rijrjl

c        seemingly impossible, but happens due to roundoff errors:

         if (costh .lt. -1.0d0) then
            costh = -1.0d0
         else if (costh .gt. 1.d0) then
            costh = 1.d0
         endif

c     tsplin is the cubic spline switching function that turns off the
c     torsion interactions at nearly linear bond angles.
c     S(t_theta(cos(theta))).  dtspln is -dS/d(cos theta).

c     (C=C double bonds will rarely have near-linear bond angles at low
c     temperatures, but they begin to become relevant for ethylene at
c     ~2500 K.)

         if (costh .gt. tthmax) then
            tsplin = 0.d0
            dtspln = 0.d0
         else if (costh .gt. tthmin) then
            dtthet = costh - tthmin
            twidth = tthmax - tthmin
            tee = dtthet / twidth
            tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
         else
            tsplin = 1.d0
            dtspln = 0.d0
         endif

         cosl(nl)=costh
         sinl(nl)=sqrt(1.0d0-costh*costh)
         tspijl(nl) = tsplin
         dtsijl(nl) = dtspln

c     three-body angles around carbon:

c        theta_ijl          cos(theta_ijl)    ig
c        [0     ,109.47) => (-1/3,1   ]    => 4
c        [109.47,120   ) => (-1/2,-1/3]    => 3
c        [120   ,131.81) => (-2/3,-1/2]    => 2
c        [131.81,180   ] => [-1  ,-2/3]    => 1

         if (jtype .eq. icarb) then
            ig=igc(int(-costh*12.0d0)+13)

c     for 109.47 <= theta <= 180:

            if(ig.ne.4) then

c     smallg is g_C from eq. 7, constructed from a 5th-order spline in
c     cos(theta).  dgdcos is d(g_C)/d(cos(theta))

c     unrolled for speed

               smallg = spgc(1,ig)
     .              + costh * (spgc(2,ig)
     .              + costh * (spgc(3,ig)
     .              + costh * (spgc(4,ig)
     .              + costh * (spgc(5,ig)
     .              + costh *  spgc(6,ig)))))
               dgdcos = spgc(2,ig)
     .              + costh * (2.d0 * spgc(3,ig)
     .              + costh * (3.d0 * spgc(4,ig)
     .              + costh * (4.d0 * spgc(5,ig)
     .              + costh *  5.d0 * spgc(6,ig))))

c     for 0 <= theta < 109.47:
c    
            else
c     unnecessary:
               qsplj=0.0d0
c     unnecessary:
               dalj=0.0d0
c     bgqmin and bgqmax are the lower and upper bounds of the switching
c     region that changes from G_C to gamma_C based on the coordination
c     number of carbons with a 3-body angle of less than 109.47 degrees.

               if (bgntji .lt. bgqmax) then

c     qsplj is Q_j(N^t_ji) from eq. 10, the switching function between
c     G_C and gamma_C.  qsplj = 0 means G_C is fully on.
c     dalj is d(Q_j(N^t_j))/d(N^T_j).

c     smallg contains g_C, the splined combination of G_C and gamma_C.
c     before being splined together, smallg contains G_C and gammac
c     contains gamma_C.  dgdcos contains d(g_C)/d(cos(theta)).  before
c     being splined together, dgdcos contains d(G_C)/d(cos(theta)) and
c     dgamdc contains d(gamma_C)/d(cos(theta)).  As N^t_i passes from
c     3.2 to 3.7 (CN of 4.2 to 4.7) g_C starts looking more like G_C and
c     less like gamma_C.

c     gdgntj = d(g^C)/d(N^t_ji)

c     the calculation of smallg is unrolled from the previous version
c     to avoid exponentiation.  it is also inserted into the nested ifs
c     to avoid calculating pieces that aren't needed

                  qsplj=1.0d0
                  if (bgntji .gt. bgqmin) then
c     in this branch 3.2 < bgntij (= N^t_ij) < 3.7, so qsplj, dalj, and
c     dgdntj must be calculated, and both G_C and gamma_C are needed
                     bigqnt = pibydq * (bgntji - bgqmin)
                     qsplj = 0.5d0 * (1.d0 + cos(bigqnt))
                     dalj = -pibydq * 0.5d0 * sin(bigqnt)
                     smallg = spgc(1,ig)
     .                    + costh * (spgc(2,ig)
     .                    + costh * (spgc(3,ig)
     .                    + costh * (spgc(4,ig)
     .                    + costh * (spgc(5,ig)
     .                    + costh *  spgc(6,ig)))))
                     dgdcos = spgc(2,ig)
     .                    + costh * (2.d0 * spgc(3,ig)
     .                    + costh * (3.d0 * spgc(4,ig)
     .                    + costh * (4.d0 * spgc(5,ig)
     .                    + costh *  5.d0 * spgc(6,ig))))
                     ig1 = ig + 1
                     gammac = spgc(1,ig1)
     .                    + costh * (spgc(2,ig1)
     .                    + costh * (spgc(3,ig1)
     .                    + costh * (spgc(4,ig1)
     .                    + costh * (spgc(5,ig1)
     .                    + costh *  spgc(6,ig1)))))
                     dgamdc = spgc(2,ig1)
     .                    + costh * (2.d0 * spgc(3,ig1)
     .                    + costh * (3.d0 * spgc(4,ig1)
     .                    + costh * (4.d0 * spgc(5,ig1)
     .                    + costh *  5.d0 * spgc(6,ig1))))
                     dgdntj = dalj * (gammac - smallg)
                     smallg = smallg + qsplj * (gammac - smallg)
                     dgdcos = dgdcos +
     .                    qsplj * (dgamdc - dgdcos)
                  else
c     in this branch bgntij (= N^t_ji) <= 3.2, so qsplj=1, dalj = 0,
c     dgdntj = 0, and only gamma_C is needed
                     ig1 = ig + 1
                     gammac = spgc(1,ig1)
     .                    + costh * (spgc(2,ig1)
     .                    + costh * (spgc(3,ig1)
     .                    + costh * (spgc(4,ig1)
     .                    + costh * (spgc(5,ig1)
     .                    + costh *  spgc(6,ig1)))))
                     dgamdc = spgc(2,ig1)
     .                    + costh * (2.d0 * spgc(3,ig1)
     .                    + costh * (3.d0 * spgc(4,ig1)
     .                    + costh * (4.d0 * spgc(5,ig1)
     .                    + costh *  5.d0 * spgc(6,ig1))))
                     smallg = gammac
                     dgdcos = dgamdc
                  endif
               else
c     in this branch bgntij (= N^t_ji) is >= 3.7, so qsplj=dalj=dgdntj=0, and
c     only G_C is needed
                  smallg = spgc(1,ig)
     .                 + costh * (spgc(2,ig)
     .                 + costh * (spgc(3,ig)
     .                 + costh * (spgc(4,ig)
     .                 + costh * (spgc(5,ig)
     .                 + costh *  spgc(6,ig)))))
                  dgdcos = spgc(2,ig)
     .                 + costh * (2.d0 * spgc(3,ig)
     .                 + costh * (3.d0 * spgc(4,ig)
     .                 + costh * (4.d0 * spgc(5,ig)
     .                 + costh *  5.d0 * spgc(6,ig))))
               endif
            endif
         else

c     three-body angles around hydrogen:

            ig=igh(int(-costh*12.0d0)+13)

c     smallg is g_H.  dgdcos is d(g_H)/d(cos(theta))

c     unrolled for speed

            smallg = spgh(1,ig)
     .           + costh * (spgh(2,ig)
     .           + costh * (spgh(3,ig)
     .           + costh * (spgh(4,ig)
     .           + costh * (spgh(5,ig)
     .           + costh *  spgh(6,ig)))))
            dgdcos = spgh(2,ig)
     .           + costh * (2.d0 * spgh(3,ig)
     .           + costh * (3.d0 * spgh(4,ig)
     .           + costh * (4.d0 * spgh(5,ig)
     .           + costh *  5.d0 * spgh(6,ig))))
         endif
         fjl = fij(jlpair)
         dfjl = dww(jlpair)
         fxjl(nl)=0.0d0
         dfxjl(nl)=0.0d0

c     for all carbon neighbors, check their coordination numbers to
c     look for conjugation.

         if (ltype .eq. icarb) then

c     bgntlj is the coordination of atom latom excluding jatom.
c     N^t_lj.
c     N^t_l - f^c_jl(r_jl) according to eq. 15, but that's not right,
c     because N^t_l already excludes f^c_jl(r_jl).

            bgntlj = rnpls1(latom,icarb) + rnpls1(latom,ihyd) - fjl
     .           - 2.d0

c     fxjl is the degree of unsaturation of atom latom.
c     dfxjl is d(fxjl)/d(bgntlj)

            if (bgntlj .lt. fmx) then
               if (bgntlj .le. fmn) then
                  fxjl(nl) = 1.0d0
               else
                  px = pi * (bgntlj - fmn)
                  fxjl(nl) = 0.5d0 * (1.d0 + cos(px))
                  dfxjl(nl) = -fjl * sin(px) * 0.5d0 * pi
               endif
            endif
         endif

c     conl is the weighted sum of conjugation terms for all neighbors
c     of jn except i.  see eq. 13

c     belongs in if above but SLOWER!?!:
         conl = conl + fjl * fxjl(nl)

c     e2lam is the exp(lambda_ijl) term in eq. 7.  It only shows up when
c     j is a hydrogen.  It gets smaller (< 1) as rij shrinks and larger
c     (> 1) as rjl shrinks.  It is one when both bonds are at their
c     standard lengths: rhh = 1.09, rch = 0.7415886997.

         if (ell(jtype,itype,ltype) .ne. 0.d0) then
            e2lam = reg(jtype,itype,ltype) *
     .           exp(ell(jtype,itype,ltype) * (rijloc - rjl))
         else
            e2lam = 1.0d0
         endif

c     dctdil = dctil(nl) = 1/r_il d(cos(theta_ijl))/d(r_il)
c     dctdji = dctji(nl) = 1/r_ij d(cos(theta_ijl))/d(r_ji) = 0
c     dctdjl = dctjl(nl) = 1/r_jl d(cos(theta_ijl))/d(r_jl)
c     all are partial derivatives, not total

c     the d/d(r_ij) derivative is zero for this spoofed pair

         dctdil = -2.0d0 / rijrjl
         dctdjl = (-rr + ril2) / (rijrjl * rjl2)
        
         dctil(nl)=dctdil
         dctjl(nl)=dctdjl

c     ssuml = sum over l of f^c G exp(lambda) terms in eq. 7
c     xsji = sum_l of 1/r_ji d/dr_ji (f^c G exp(lambda)) in eq. 7. = 0
c     xsjl(nl) = 1/r_jl d/dr_jl (f^c G exp(lambda)) in eq. 7
c     sdaljl = delta_i^C sum_k of f^c dQ/d^t_j (gamma - G) exp(lambda)
c     xsil(nl) = 1/r_il d/dr_il (f^c G exp(lambda)) in eq. 7

c     the d/d(r_ij) derivative is zero for this spoofed pair

         ge2lam = smallg * e2lam
         ssuml = ssuml + fjl * ge2lam
         xtemp = fjl * e2lam * dgdcos
         gfx = ge2lam * fjl * ell(jtype,itype,ltype)
         xsjl(nl) = (ge2lam * dfjl - gfx) / rjl + xtemp * dctdjl
         sdaljl = sdaljl + e2lam * fjl * dgdntj
         xsil(nl) = xtemp*dctdil
 2000 continue

      pij = 0.0d0
      dpijdn(icarb) = 0.d0
      dpijdn(ihyd) = 0.d0

c     pij is the P_ij(N^C_ij,N^H_ij) term in eq. 7.  it only shows up
c     when atom i is a carbon.

      if (itype .eq. icarb) then
         nh = int(bnijp1(ihyd) + 1.d-12)
         nc = int(bnijp1(icarb) + 1.d-12)

c     if bnijp1(C) and bnijp1(H) are not integers, we have to use a bicubic
c     spline to interpolate pij.

         if ((abs(dble(nh) - bnijp1(ihyd)) .gt. 1.d-8) .or.
     .        (abs(dble(nc) - bnijp1(icarb)) .gt. 1.d-8)) then
            call bicub(bnijp1(icarb), bnijp1(ihyd), clm(1,jtype,nh,nc),
     .           pij, dpijdn(icarb), dpijdn(ihyd))

c     but if they are both integers we can read them right from xh.
c     dpijdn(icarb) and dpijdn(ihyd) are the derivatives with respect to
c     nc and nh.

         else
            pij = xh(jtype,nh,nc)
            dpijdn(ihyd) = xh1(jtype,nh,nc)
            dpijdn(icarb) = xh2(jtype,nh,nc)
         endif
      endif

      pji = 0.0d0
      dpjidn(icarb) = 0.d0
      dpjidn(ihyd) = 0.d0

c     pji is the P_ji(N^C_ji,N^H_ji) term in eq. 7.  it only shows up
c     when atom j is a carbon.

      if (jtype .eq. icarb) then
         nh = int(bnjip1(ihyd) + 1.d-12)
         nc = int(bnjip1(icarb) + 1.d-12)

c     if bnijp1(C) and bnijp1(H) are not integers, we have to use a bicubic
c     spline to interpolate pji.

         if ((abs(dble(nh) - bnjip1(ihyd)) .gt. 1.d-8) .or.
     .        (abs(dble(nc) - bnjip1(icarb)) .gt. 1.d-8)) then
            call bicub(bnjip1(icarb), bnjip1(ihyd), clm(1,itype,nh,nc),
     .           pji, dpjidn(icarb), dpjidn(ihyd))
         else
            pji = xh(itype,nh,nc)
            dpjidn(ihyd) = xh1(itype,nh,nc)
            dpjidn(icarb) = xh2(itype,nh,nc)
         endif
      endif

c     bij is p^sigma,pi_ij from eqs. 2 and 7. 
c     dij is what's inside bij's radical
c     bji is p^sigma,pi_ji from eqs. 2 and 7.
c     dji is what's inside bji's radical
c     dbddij = d(bij)/d(dij)
c     dbddji = d(bji)/d(dji)
c     vatt is -1/2 V^A(r_ij)
c     conjug is N^conj_ij from eq. 13.
c     bntijp is N^t_ij + 1
c     bntjip is N^t_ji + 1
c     rad is 2 * pi^rc_ij
c     dradi is 2 * d(pi^rc_ij)/d(N^t_i)
c     dradj is 2 * d(pi^rc_ij)/d(N^t_j)
c     drdc is 2 * d(pi^rc_ij)/d(N^conj_ij)

      dij = (1.0d0 + pij + sfge2l)
      bij = 1.d0/sqrt(dij)
      dji = (1.0d0 + pji + ssuml)
      bji = 1.d0 / sqrt(dji)
      dbddij = -0.5d0 * bij / dij
      dbddji = -0.5d0 * bji / dji

c     vatt is -1/2 V^A.  This is the factor that multiplies db/d* terms,
c     since the potential looks like b V^A.

c     for this spoofed pair, the term involving b is
c     S(r) S'(t'(b)) C V^LJ.
c     Thus the db/d* term looks like S(r) dS'/dt dt/db db/d* C V^LJ.
c     vattlj defined here is the S(r) C V^LJ piece which multiplies
c     db/d*, but it is still missing dS'/dt and dt/db.

      vattlj = -0.5d0 * veelj

      dradi=0.0d0
      dradj=0.0d0
      drdc=0.0d0
      conjug = 1.0d0 + (cni ** 2) + (conl ** 2)
      bntijp = bnijp1(icarb) + bnijp1(ihyd) - 1.d0
      bntjip = bnjip1(icarb) + bnjip1(ihyd) - 1.d0

c     physically, bntijp, bntjip, and conjug should never be <1. 
c     but roundoff errors can
c     cause them to be 0.9999..., particularly for spoofed pairs which
c     have been fakepr-ed

c     the limits are to keep the variables within the range of
c     the tricubic spline.  the upper limits do weird rounding because
c     it's actually int(bntijp) and int(bntjip) that controls the
c     spline coefficients.

      if (bntijp .gt. 4.999d0) then
         bntijp = 4.999d0
      else if (bntijp .lt. 1.d0) then
         bntijp = 1.d0
      endif

      if (bntjip .gt. 4.999d0) then
         bntjip = 4.999d0
      else if (bntjip .lt. 1.d0) then
         bntjip = 1.d0
      endif

      if (conjug .gt. 9.999d0) then
         conjug = 9.999d0
      else if (conjug .lt. 1.d0) then
         conjug = 1.d0
      endif

      lindex = int(bntijp)
      mindex = int(bntjip)
      nindex = int(conjug)
      call tricub(bntijp, bntjip, conjug,
     .     clmn(1,lindex,mindex,nindex,ijtype),
     .     rad, dradi, dradj, drdc)

      btot=(bji+bij+rad)

c     Dihedral terms

c     don't calculate dihedral terms unless i and j are both carbons

      if (ijtype .ne. ijcc) then
         go to 2800
      endif

      btor=0.0d0

      if (bntijp .ge. 4.d0 .or. bntjip .ge. 4.d0) then
         go to 2800
      endif

      if (nindex .gt. 9) then
         nindex = 9
      endif

c     tij is 2 * T_ij
c     dtijdi is 2 * d(T_ij)/d(N^t_ij)
c     dtijdj is 2 * d(T_ij)/d(N^t_ji)
c     dtijdc is 2 * d(T_ij)/d(N^conj_ij)

      call tricub(bntijp, bntjip, conjug, tlmn(1,lindex,mindex,nindex),
     .     tij, dtijdi, dtijdj, dtijdc)

c     skip out if there was no T_ij (to speak of)

      if(abs(tij).le.1.0d-08) go to 2800

      nk=0
      do 4000 ikpair = ibegin, iend

         if (iprtyp(ikpair) .ne. ihycar) then
            go to 4000
         endif

         katom = jhalf(ikpair)

c              katom (neighbor of iatom) should not be jatom

         if (katom .eq. jatom) then
            go to 4000
         endif

         nk = nk + 1

c     if theta_jik = 180 degrees, torsion interaction will be zero.
c     skip out early to avoid singularity.

         if (abs(sink(nk)) .eq. 0.d0) then
            go to 4000
         endif

         sink2i = 1.d0 / (sink(nk) * sink(nk))

         do 2200 idim = 1, 3
            rikvec(idim) = rijv(ikpair,idim)
 2200    continue

         rik = rcor(ikpair)
         rik2i = 1.d0 / (rik * rik)

c     note that DWB's REBO potential uses a switching range of 0.3
c     here, rather than the 0.5 used for f^c_CH.

         if (iat2ty(katom) .eq. ihyd) then

c     these values cover the default, if the i-k pair is far enough
c     not to have to worry about the switching range

            fik = 1.0d0
            dfik=0.0d0

c     if the i-k pair is not in range, skip out of the loop

            if (rik .ge. 1.60d0) then
               go to 4000
            endif

c     if the i-k pair is in the switching function, calculate
c     fik and dfik

            if (rik .ge. 1.30d0) then
               dtemp = pidt * (rik - 1.30d0)
               fik = 0.5d0 * (1.d0 + cos(dtemp))
               dfik = -pidt * 0.5d0 * sin(dtemp)
            endif
         else
            fik = fij(ikpair)
            dfik = dww(ikpair)
         endif

         nl=0
         do 3900 jlpair = jbegin, jend

            latom = jhalf(jlpair)

c     latom (neighbor of jatom) should not be iatom

            if (latom .eq. iatom) then
               go to 3900
            endif

c     skip out if pair is not hydrocarbon

            if (iprtyp(jlpair) .ne. ihycar) then
               go to 3900
            endif

            nl = nl + 1

c     don't consider cyclic 3-membered rings

            if (latom .eq. katom) then
               go to 3900
            endif

c     if theta_ijl = 180 degrees, torsion interaction will be zero.
c     skip out early to avoid singularity.

            if (abs(sinl(nl)) .eq. 0.d0) then
               go to 3900
            endif

            sinl2i = 1.d0 / (sinl(nl) * sinl(nl))

            do 2300 idim = 1, 3
               rjlvec(idim) = rijv(jlpair,idim)
 2300       continue

            rjl = rcor(jlpair)
            rjl2i = 1.d0 / (rjl * rjl)

            if (iat2ty(latom) .eq. ihyd) then

c     these values cover the default, if the j-l pair is close enough
c     not to have to worry about the switching range

               fjl=1.0d0
               dfjl=0.0d0

c     if the j-l pair is beyond the switching range, skip out of the
c     loop
               if (rjl .ge. 1.60d0) then
                  go to 3900
               endif

c     if the j-l pair is in the switching zone, calculate fjl and dfjl

               if (rjl .ge. 1.30d0) then
                  dtemp = pidt * (rjl - 1.30d0)
                  fjl = 0.5d0 * (1.d0 + cos(dtemp))
                  dfjl = -pidt * 0.5d0 * sin(dtemp)
               endif
            else
               fjl = fij(jlpair)
               dfjl = dww(jlpair)
            endif

c     cwnom = r_ij r_ik sin(theta_jik) r_ij r_jl sin(theta_ijl),
c     which is the magnitude of the cross product (r_ji x r_ik),
c     multiplied by the magnitude of the cross product (r_ij x r_jl),
c     which is the denominator of cos(omega).

c     dt1dik = 1/r_ik 1/cwnom d(cwnom)/d(r_ik)
c     dt1djk = 1/r_jk 1/cwnom d(cwnom)/d(r_jk)
c     dt1djl = 1/r_jl 1/cwnom d(cwnom)/d(r_jl)
c     dt1dil = 1/r_il 1/cwnom d(cwnom)/d(r_il)
c     dt1dij = 1/r_ij 1/cwnom d(cwnom)/d(r_ij) = 0

c     the d/d(r_ij) derivative will be zero for this spoofed pair

c     rearranged for speed

            cwnom = rik * rjl * rijloc * rijloc
     .           * sink(nk) * sinl(nl)

            dt1dik = rik2i - dctik(nk) * sink2i * cosk(nk)
            dt1djk = -dctjk(nk) * sink2i * cosk(nk)
            dt1djl = rjl2i - dctjl(nl) * sinl2i * cosl(nl)
            dt1dil = -dctil(nl) * sinl2i * cosl(nl)

c     cross products:
c     (crjikx,crjiky,crjikz) = r_ik x r_ij = r_ji x r_ik
c     (crijlx,crijly,crijlz) = r_ij x r_jl

            crjikx = rikvec(2) * rijvec(3) - rijvec(2) * rikvec(3)
            crijlx = rijvec(2) * rjlvec(3) - rjlvec(2) * rijvec(3)
            crjiky = rikvec(3) * rijvec(1) - rijvec(3) * rikvec(1)
            crijly = rijvec(3) * rjlvec(1) - rjlvec(3) * rijvec(1)
            crjikz = rikvec(1) * rijvec(2) - rijvec(1) * rikvec(2)
            crijlz = rijvec(1) * rjlvec(2) - rjlvec(1) * rijvec(2)
           
c     dot product:
c     cwnum = (r_ji x r_ik) . (rij x r_jl)

            cwnum = crjikx * crijlx + crjiky * crijly
     .           + crjikz * crijlz

c     divide by the magnitude of each individual cross product so that
c     cw contains cos omega_ijkl = e_jik . e_ijl (eq. 17)

            cw = cwnum / cwnom

c     bt is sin^2 omega_ijkl

            bt=(1.0d0-cw*cw)

c     btor = sin^2 omega_ijkl f^c_ik(r_ik) f^c_jl(r_jl) x
c               (1-S(t(cos theta_jik))) (1-S(t(cos theta_ijl))),
c     summed over torsion angles.

            btor = btor + bt * fik * fjl * (1.d0 - tspjik(nk))
     .           * (1.d0 - tspijl(nl))

c     these are the vector derivatives d(cwnum)/d(r_ik), etc.

c     dt2dik = -r_ji x (r_ij x r_jl) = d(cwnum)/d(r_ik)
c     dt2djl = (r_ji x r_ik) x r_ij) = d(cwnum)/d(r_jl)
c     dt2dij = -r_ik x (r_ij x r_jl) - (r_ji x r_ik) x r_jl
c            = d(cwnum)/d(r_ij) = 0

c     the d/d(r_ij) derivative is zero for this spoofed pair

            dt2dik(1) = -rijvec(3) * crijly + rijvec(2) * crijlz
            dt2dik(2) = -rijvec(1) * crijlz + rijvec(3) * crijlx
            dt2dik(3) = -rijvec(2) * crijlx + rijvec(1) * crijly
              
            dt2djl(1) = -rijvec(2) * crjikz + rijvec(3) * crjiky
            dt2djl(2) = -rijvec(3) * crjikx + rijvec(1) * crjikz
            dt2djl(3) = -rijvec(1) * crjiky + rijvec(2) * crjikx

c     inefficient

c     aa is -V^A * 1 / phi^denom * d(pi_ijkl)/d(cos omega_ijkl) T_ij
c     aaa1 is -V^A * T_ij * (1-cos^2(omega_ijkl))
c     at2 is V^A T d(pi_ijkl)/d(cos omega_ijkl) d(cos omega_ijkl)
c     /d(phi^denom) * phi^denom.

c     each of these includes a V^A prefactor multiplied by some stuff
c     which will become the db/d* derivatives.  The V^A term comes from
c     vatt, which is -1/2 V^A.

c     for this spoofed pair, the prefactor multiplying db/d* is
c     S(t(r)) dS/dt' dt'/db C V^LJ.  The vattlj term used here includes
c     only -1/2 S(t(r)) C V^LJ, and neglects dS/dt' dt'/db which will be
c     put in later.

            aalj = -vattlj * 2.d0 * cw / cwnom * tij
     .           * fjl * fik * (1.d0 - tspjik(nk)) * (1.d0 - tspijl(nl))
            aaa1lj = vattlj * bt * tij * (1.d0 - tspjik(nk))
     .           * (1.d0 - tspijl(nl))
            aaa2lj = vattlj * bt * tij * fik * fjl
            at2lj = aalj * cwnum

c     these are some pieces of the chain rule derivatives of the
c     pi^dh_ijkl term, multiplied by a -1/r V^A T_ij prefactor.
c     they do not yet include the d/d(phi^num) piece of the chain
c     rule

c     V = -V^a Tij pi^dh fik fjl (1-S(t(cos theta_jik))) x
c            (1-S(t(cos theta_ijl)))
c     so that
c     dV/dr = -d(V^A)/dr Tij pi^dh fik fjl (1-S) (1-S)        (1)
c                -V^A d(Tij)/dr pi^dh fik fjl (1-S) (1-S)     (2)
c                -V^A Tij d(pi^dh)/dr fik fjl (1-S) (1-S)     (3)
c                -V^A Tij pi^dh d(fik)/dr fjl (1-S) (1-S)     (4)
c                -V^A Tij pi^dh fik d(fjl)/dr (1-S) (1-S)     (5)
c                -V^A Tij pi^dh fik fjl x
c                   (-dS/dt dt/d(cos) d(cos)/dr) (1-S)        (6)
c                -V^A Tij pi^dh fik fjl x
c                   (1-S) (-dS/dt dt/d(cos) d(cos)/dr)        (7)
c     lines 1 and 2 are taken care of elsewhere.
c     the terms below take care of lines 4, 5, 6, 7, and part of 3
c     in particular, line 3 is subdivided further,

c     d(pi^dh)/dr 
c        = d(pi^dh)/d(cos omega) d(cos omega)/d(num) x
c                       d(num)/dr                             (3a)
c             + d(pi^dh)/d(cos omega) d(cos omega)/d(denom) x
c                       d(denom)/dr                           (3b)

c     line 3b is included in the terms immediately below.
c     the d/d(r_ij) derivative is zero for this spoofed pair
c     and the V^A piece of the prefactor should be
c     S(t(r)) dS/dt' dt'/db C V^LJ.  the aaa1lj, aaa2lj, at2lj terms
c     include only the S(t(r)) C V^LJ pieces.

            fcikpc = -dt1dik * at2lj
     .           + aaa1lj * fjl * dfik / rik
     .           + aaa2lj * dtsjik(nk) * dctik(nk) * (1.d0 - tspijl(nl))
            fcjlpc = -dt1djl * at2lj
     .           + aaa1lj * fik * dfjl / rjl
     .           + (1.d0 - tspjik(nk)) * aaa2lj * dtsijl(nl) * dctjl(nl)
            fcjkpc = -dt1djk * at2lj
     .           + aaa2lj * dtsjik(nk) * dctjk(nk) * (1.d0 - tspijl(nl))
            fcilpc = -dt1dil * at2lj
     .           + (1.d0 - tspjik(nk)) * aaa2lj * dtsijl(nl) * dctil(nl)

c     the fcpc terms below include line 3b.
c     these are the full d/dx derivatives of the pi^dh_ijkl terms,
c     multiplied by a -V^A_ij T_ij prefactor.  The negative sign
c     is because F = -dV/dx.

c     the tnp arrays are the temporary holding places for the forces
c     (what will be fint).
c     tnpi()    = force on iatom
c     tnpj()    = force on jatom
c     tnpk(nk,) = force on katom
c     tnpl(nl,) = force on latom
c     uutemp accumulates the contributions to the W virial that haven't
c     been scaled yet.  these will eventually go into uu.
c     uutemp_ab = sum_k r_ika . f_ikb + sum_l r_jla . f_jlb

c     don't put the forces directly into fint.  instead, save them in
c     the tnp variables until later.  they will be added with a
c     multiplicative factor that depends on the value of b_ij, which
c     we don't know yet.

c     maybe faster to do separate loops:

            do 2100 idim = 1, 3

c     the d/d(r_ij) derivative is zero for this spoofed pair

               fcpc = fcikpc * rikvec(idim) + aalj * dt2dik(idim)
               tnpi(idim) = tnpi(idim) + fcpc
               tnpk(nk,idim) = tnpk(nk,idim) - fcpc
               do 2076 jdim = 1, 3
                  uutemp(jdim,idim) = uutemp(jdim,idim)
     .                 + fcpc * rikvec(jdim)
 2076          continue
              
               fcpc = fcjlpc * rjlvec(idim) + aalj * dt2djl(idim)
               tnpj(idim) = tnpj(idim) + fcpc
               tnpl(nl,idim) = tnpl(nl,idim) - fcpc
               do 2077 jdim = 1, 3
                  uutemp(jdim,idim) = uutemp(jdim,idim)
     .                 + fcpc * rjlvec(jdim)
 2077          continue

c     the j-k forces are complicated.  Even though this is the
c     d/d(r_jk) derivative, the force acts on the i, j, and k atoms.
c     This is because r_jk is calculated from the spoofed r_ij*, and
c     thus depends on the position of atom i.

c     I'm not sure whether the pressure calculation is right, since
c     these forces are not pairwise.

               fcpc = fcjkpc * rjk(nk,idim)
               tnpi(idim) = tnpi(idim) + fcpc
               tnpk(nk,idim) = tnpk(nk,idim) - fcpc
               do 2080 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * rjk(nk,jdim)
 2080          continue

               fcpc = fcjkpc * rijmbr * rjk(nk,idim)
               tnpj(idim) = tnpj(idim) + fcpc
               tnpi(idim) = tnpi(idim) - fcpc
               do 2082 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * rjk(nk,jdim)
 2082          continue

               do 2085 jdim = 1, 3
                  fcpc = -fcjkpc * rijmbr * rlrijv(idim)
     .                 * rlrijv(jdim) / rlrij2 * rjk(nk,idim)
                  tnpj(jdim) = tnpj(jdim) + fcpc
                  tnpi(jdim) = tnpi(jdim) - fcpc
                  do 2084 kdim = 1, 3
                     uu(kdim,idim) = uu(kdim,idim)
     .                    + fcpc * rjk(nk,kdim)
 2084             continue
 2085          continue

c     the i-l forces are complicated.  Even though this is the
c     d/d(r_il) derivative, the force acts on the i, j, and l atoms.
c     This is because r_il is calculated from the spoofed r_ij*, and
c     thus depends on the position of atom j.

c     I'm not sure whether the pressure calculation is right, since
c     these forces are not pairwise.

               fcpc = fcilpc * ril(nl,idim)
               tnpj(idim) = tnpj(idim) + fcpc
               tnpl(nl,idim) = tnpl(nl,idim) - fcpc
               do 2090 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * ril(nl,jdim)
 2090          continue

               fcpc = fcilpc * rijmbr * ril(nl,idim)
               tnpi(idim) = tnpi(idim) + fcpc
               tnpj(idim) = tnpj(idim) - fcpc
               do 2092 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * ril(nl,jdim)
 2092          continue

               do 2095 jdim = 1, 3
                  fcpc = -fcilpc * rijmbr * rlrijv(idim) * rlrijv(jdim)
     .                 / rlrij2 * ril(nl,idim)
                  tnpi(jdim) = tnpi(jdim) + fcpc
                  tnpj(jdim) = tnpj(jdim) - fcpc
                  do 2094 kdim = 1, 3
                     uu(kdim,idim) = uu(kdim,idim) + fcpc * ril(nl,kdim)
 2094             continue
 2095          continue

 2100       continue
                    
 3900    continue
 4000 continue

c     btor*tij is 2 * pi^dh.  add it into btot, which is now 2 * b_ij.
c     dradi now contains 2 * d(p^pi_ij)/d(N^t_i)
c     dradj now contians 2 * d(p^pi_ij)/d(N^t_j)
c     drdc now contains 2 * d(p^pi_ij)/d(N^conj_ij)

      btot = btot + btor * tij

      dradi = dradi + dtijdi * btor
      dradj = dradj + dtijdj * btor
      drdc = drdc + dtijdc * btor

 2800 continue

c     Use bij* to calculate t_b(bij*) and S(t_b(bij*)).

c     btot   = 2 * bij*
c     onebij = bij*
c     tee    = t_b(bij*)
c     dsdt   = dS(t_b(bij*)/dt_b(bij*)
c     dsdb   = dS(t_b(bij*))/d(bij*)
c     bsplin = S(t_b(bij*))

      onebij = btot * 0.5d0

      if (onebij .lt. bijmin(itype,jtype)) then
         bsplin = 1.d0
         dsdb   = 0.d0
      else if (onebij .lt. bijmax(itype,jtype)) then
         dbij = onebij - bijmin(itype,jtype)
         swidth = bijmax(itype,jtype) - bijmin(itype,jtype)
         tee = dbij / swidth
         dsdt = 6.d0 * tee * (tee - 1.d0)
         dtdb = 1.d0 / swidth
         dsdb = dsdt * dtdb
         bsplin = 1.d0 - tee * tee * (3.d0 - 2.d0 * tee)
      else
         bsplin = 0.d0
         dsdb = 0.d0
      endif

c     if dsdb = 0, we can probably skip out from here?

c     btot * vatt is b_ij * V^A_ij(r_ij).
c     tote now contains all of the REBO terms (but not yet the LJ terms)

c     inefficient

c     energies are taken care of in ljguts, so don't do anything here

c     the dihedral forces were mostly calculated earlier.  now that
c     we know bij, we can put them into the force array

      do 22100 idim = 1, 3
         fint(iatom,idim) = fint(iatom,idim) + tnpi(idim) * dsdb
         fint(jatom,idim) = fint(jatom,idim) + tnpj(idim) * dsdb
         do 22099 jdim = 1, 3
            uu(jdim,idim) = uu(jdim,idim) + uutemp(jdim,idim) * dsdb
22099    continue
22100 continue

      nk = 0
      do 23000 ikpair = ibegin, iend

         if (iprtyp(ikpair) .ne. ihycar) then
            go to 23000
         endif

         katom = jhalf(ikpair)

c     katom (neighbor of iatom) should not be jatom

         if (katom .eq. jatom) then
            go to 23000
         endif

         nk = nk + 1

c     no contribution to uu because the reciprocal forces have already
c     been counted in tnpi or tnpj

         do 22900 idim = 1, 3
            fint(katom,idim) = fint(katom,idim) + tnpk(nk,idim) * dsdb
22900    continue

23000 continue

      nl = 0
      do 24000 jlpair = jbegin, jend

         if (iprtyp(jlpair) .ne. ihycar) then
            go to 24000
         endif

         latom = jhalf(jlpair)
        
c     latom (neighbor of jatom) should not be iatom

         if (latom .eq. iatom) then
            go to 24000
         endif

         nl = nl + 1

c     no contribution to uu because the reciprocal forces have already
c     been counted in tnpi or tnpj

         do 23900 idim = 1, 3
            fint(latom,idim) = fint(latom,idim) + tnpl(nl,idim) * dsdb
23900    continue
        
24000 continue
        
c     now that we know dsdb = d(S(t(b_ij)))/d(b_ij), we can fold it
c     into the energy term that multiplies all derivatives calculated
c     here

      vattlj = vattlj * dsdb

c     everything from here down is all force calculations

      vldbdi = vattlj * dbddij
      vldbdj = vattlj * dbddji
      vldrdc = vattlj * drdc
      vldrdi = vattlj * dradi
      vldrdj = vattlj * dradj

c     the sdgdij and xsji terms are zero for this spoofed pair, and the
c     btot * dexx1 term is b * d(V^A)/dr, so it doesn't belong here

c     Add many-body forces

c     I side of bond

      nk=0

c     loop over non-jatom neighbors of iatom

      do 5000 ikpair = ibegin, iend

c     ignore non-hydrocarbon bonds

         if (iprtyp(ikpair) .ne. ihycar) then
            go to 5000
         endif
         katom = jhalf(ikpair)
        
c     katom (neighbor of iatom) should not be jatom
        
         if (katom .eq. jatom) then
            go to 5000
         endif

         ktype = iat2ty(katom)
         dwr = dww(ikpair) / rcor(ikpair)
         nk=nk+1

c     First Neighbors

c     There was a bug in DWB's code here. 
c     The middle term in fcikpc should be:
c       + dwr * (vldrdi + vldrdc * 2.d0 * cni * fxik(nk))
c     The extra factor of 2.d0 * cni takes care of the
c     d(conjug)/d(cni) piece which was missing

         fcikpc = vldbdi * (dgdik(nk) + dwr * dpijdn(ktype))
     .        + dwr * (vldrdi + vldrdc * 2.d0 * cni * fxik(nk))
     .        + vldbdi * dwr * sdgdnt
         fcjkpc = vldbdi * dgdjk(nk)

c     loop unrolled for speed

c     i-k forces first:

         fcpc = fcikpc * rijv(ikpair,1)
         fint(iatom,1) = fint(iatom,1) + fcpc
         fint(katom,1) = fint(katom,1) - fcpc
         do 4954 jdim = 1, 3
            uu(jdim,1) = uu(jdim,1) + fcpc * rijv(ikpair,jdim)
 4954    continue

         fcpc = fcikpc * rijv(ikpair,2)
         fint(iatom,2) = fint(iatom,2) + fcpc
         fint(katom,2) = fint(katom,2) - fcpc
         do 4956 jdim = 1, 3
            uu(jdim,2) = uu(jdim,2) + fcpc * rijv(ikpair,jdim)
 4956    continue

         fcpc = fcikpc * rijv(ikpair,3)
         fint(iatom,3) = fint(iatom,3) + fcpc
         fint(katom,3) = fint(katom,3) - fcpc
         do 4958 jdim = 1, 3
            uu(jdim,3) = uu(jdim,3) + fcpc * rijv(ikpair,jdim)
 4958    continue

c     now the j-k forces.  Even though this is the d/d(r_jk) derivative,
c     the force is divided between the i, j, and k atoms in a
c     complicated way.  This is because r_jk is calculated from the
c     spoofed r_ij*, and thus depends on the position of atom i.

c     the pressure calculation may not be right for these non-pairwise
c     forces

         do 4961 idim = 1, 3
            fcpc = fcjkpc * rjk(nk,idim)
            fint(iatom,idim) = fint(iatom,idim) + fcpc
            fint(katom,idim) = fint(katom,idim) - fcpc
            do 4960 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + fcpc * rjk(nk,jdim)
 4960       continue
 4961    continue

         do 4963 idim = 1, 3
            fcpc = fcjkpc * rijmbr * rjk(nk,idim)
            fint(jatom,idim) = fint(jatom,idim) + fcpc
            fint(iatom,idim) = fint(iatom,idim) - fcpc
            do 4962 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + fcpc * rjk(nk,jdim)
 4962       continue
 4963    continue

         do 4966 idim = 1, 3
            do 4965 jdim = 1, 3
               fcpc = -fcjkpc * rijmbr * rlrijv(idim) * rlrijv(jdim)
     .              / rlrij2 * rjk(nk,idim)
               fint(jatom,jdim) = fint(jatom,jdim) + fcpc
               fint(iatom,jdim) = fint(iatom,jdim) - fcpc
               do 4964 kdim = 1, 3
                  uu(kdim,idim) = uu(kdim,idim) + fcpc * rjk(nk,kdim)
 4964          continue
 4965       continue
 4966    continue
        
c     Second Neighbors via RADIC

         ddrlj = vldrdc * dfxik(nk) * 2.d0 * cni
         if (ddrlj .eq. 0.d0) then
            go to 5000
         endif
         kbegin = nabors(katom)
         kend = nabors(katom+1) - 1
        
c     skip out if the only neighbor of k is i

         if (kbegin .eq. kend) then
            go to 5000
         endif
        
c     loop over neighbors of katom

         do 4900 kmpair = kbegin, kend

c     ignore non-hydrocarbon bonds

            if (iprtyp(kmpair) .ne. ihycar) then
               go to 4900
            endif
            matom = jhalf(kmpair)

c     matom (neighbor of katom) should not be iatom

            if (matom .eq. iatom) then
               go to 4900
            endif
            rp = ddrlj * dww(kmpair) / rcor(kmpair)
           
            do 4790 idim = 1, 3
               fcpc = rp * rijv(kmpair,idim)
               fint(katom,idim) = fint(katom,idim) + fcpc
               fint(matom,idim) = fint(matom,idim) - fcpc
               do 4789 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim)
     .                 + fcpc * rijv(kmpair,jdim)
 4789          continue
 4790       continue

 4900    continue
 5000 continue

c     j side of bond

      nl=0

c     look over neighbors latom of jatom. 
     
      do 6000 jlpair = jbegin, jend
         latom = jhalf(jlpair)

c     latom (neighbor of jatom) should not be iatom

         if (latom .eq. iatom) then
            go to 6000
         endif

c     ignore non-hydrocarbon bonds

         if (iprtyp(jlpair) .ne. ihycar) then
            go to 6000
         endif
         ltype = iat2ty(latom)
         dwr = dww(jlpair) / rcor(jlpair)
         nl=nl+1

c     First Neighbors

c     There was a bug in DWB's code here.
c     The last term in fcjlpc should be:
c       + dwr * (vldrdj + vldrdc * 2.d0 * conl * fxjl(nl))
c     The extra factor of 2.d0 * conl takes care of the
c     d(conjug)/d(conl) piece which was missing

c     rearranged for speed

         fcjlpc = vldbdj * (xsjl(nl)
     .        + dwr * (dpjidn(ltype) + sdaljl))
     .        + dwr * (vldrdj + vldrdc * 2.d0 * conl * fxjl(nl))
         fcilpc = vldbdj * xsil(nl)

c     loop unrolled for speed

c     j-l forces first:

         fcpc = fcjlpc * rijv(jlpair,1)
         fint(jatom,1) = fint(jatom,1) + fcpc
         fint(latom,1) = fint(latom,1) - fcpc
         do 5774 jdim = 1, 3
            uu(jdim,1) = uu(jdim,1) + fcpc * rijv(jlpair,jdim)
 5774    continue

         fcpc = fcjlpc * rijv(jlpair,2)
         fint(jatom,2) = fint(jatom,2) + fcpc
         fint(latom,2) = fint(latom,2) - fcpc
         do 5776 jdim = 1, 3
            uu(jdim,2) = uu(jdim,2) + fcpc * rijv(jlpair,jdim)
 5776    continue

         fcpc = fcjlpc * rijv(jlpair,3)
         fint(jatom,3) = fint(jatom,3) + fcpc
         fint(latom,3) = fint(latom,3) - fcpc
         do 5778 jdim = 1, 3
            uu(jdim,3) = uu(jdim,3) + fcpc * rijv(jlpair,jdim)
 5778    continue

c     now the i-l forces.  Even though this is the d/d(r_il) derivative,
c     the force is divided between the i, j, and l atoms in a
c     complicated way.  This is because r_il is calculated from the
c     spoofed r_ij*, and thus depends on the position of atom j.

c     the pressure calculation may not be right for these non-pairwise
c     forces

         do 5781 idim = 1, 3
            fcpc = fcilpc * ril(nl,idim)
            fint(jatom,idim) = fint(jatom,idim) + fcpc
            fint(latom,idim) = fint(latom,idim) - fcpc
            do 5780 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + fcpc * ril(nl,jdim)
 5780       continue
 5781    continue

         do 5783 idim = 1, 3
            fcpc = fcilpc * rijmbr * ril(nl,idim)
            fint(iatom,idim) = fint(iatom,idim) + fcpc
            fint(jatom,idim) = fint(jatom,idim) - fcpc
            do 5782 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + fcpc * ril(nl,jdim)
 5782       continue
 5783    continue

         do 5786 idim = 1, 3
            do 5785 jdim = 1, 3
               fcpc = -fcilpc * rijmbr * rlrijv(idim) * rlrijv(jdim)
     .              / rlrij2 * ril(nl,idim)
               fint(iatom,jdim) = fint(iatom,jdim) + fcpc
               fint(jatom,jdim) = fint(jatom,jdim) - fcpc
               do 5784 kdim = 1, 3
                  uu(kdim,idim) = uu(kdim,idim) + fcpc * ril(nl,kdim)
 5784          continue
 5785       continue
 5786    continue
        
c     Second Neighbors via RADIC

         ddrlj = vldrdc * dfxjl(nl) * 2.d0 * conl
         if (ddrlj .eq. 0.d0) then
            go to 6000
         endif
         lbegin = nabors(latom)
         lend = nabors(latom+1) - 1

c     skip out if jatom is the only neighbor of latom

         if (lbegin .eq. lend) then
            go to 6000
         endif

c     loop over neighbors natom of latom.

         do 5900 lnpair = lbegin, lend

c     ignore non-hydrocarbon bonds

            if (iprtyp(lnpair) .ne. ihycar) then
               go to 5900
            endif
            natom = jhalf(lnpair)

c     natom (neighbor of latom) should not be jatom

            if (natom .eq. jatom) then
               go to 5900
            endif
            rp = ddrlj * dww(lnpair) / rcor(lnpair)
            do 5790 idim = 1, 3
               fcpc = rp * rijv(lnpair,idim)
               fint(latom,idim) = fint(latom,idim) + fcpc
               fint(natom,idim) = fint(natom,idim) - fcpc
               do 5789 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim)
     .                 + fcpc * rijv(lnpair,jdim)
 5789          continue
 5790       continue

 5900    continue
 6000 continue
      return
      end
