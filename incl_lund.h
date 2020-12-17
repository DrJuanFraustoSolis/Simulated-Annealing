      parameter (nhptyp=3,mxconr=40,mxcell=125000)
      dimension iN(mxrs),iCa(mxrs),iC(mxrs),mlvr(mxvr)
      dimension iphi(mxrs),ipsi(mxrs)
      common /lundds/iN,iCa,iC,mlvr,iphi,ipsi

      double precision kbias      
      double precision epshb1,epshb2,powa,powb,sighb,cthb,cthb2
      double precision alhb,blhb,sighb2,cdon,cacc,casc

c -----Probability for using BGS when it is possible
      double precision abgs,bbgs, dph(8)
      integer bgsnvar,bgsvar(mxrs), iph(8)
      common /bgs_i/ bgsnvar,iph
      common /bgs_r/ abgs,bbgs,dph, bgsvar
      dimension ihpat(mxrs,6),nhpat(mxrs)
      double precision hpstrg
      dimension hpstrg(nhptyp*nhptyp)
      
      double precision exvk,exvcut,exvcut2
      dimension matcon(-mxconr:mxconr,mxat)

      double precision sigsa,sig2lcp,asalcp,bsalcp
      dimension sigsa(mxtyat),sig2lcp(mxtyat,mxtyat)
      dimension asalcp(mxtyat,mxtyat),bsalcp(mxtyat,mxtyat)

      dimension lcp1(50*mxrs),lcp2(50*mxrs),ilpst(mxml),ilpnd(mxml)

      double precision exvlam,exvcutg,exvcutg2
      double precision sig2exv,asaexv,bsaexv
      dimension sig2exv(mxtyat,mxtyat)
      dimension asaexv(mxtyat,mxtyat),bsaexv(mxtyat,mxtyat)

      common /lundff/kbias,
     #     epshb1,epshb2,powa,powb,sighb,cthb,
     #     cthb2,
     #     alhb,blhb,sighb2,cdon,cacc,casc,
     #     ihpat,nhpat,hpstrg,
     #     exvk,exvcut,exvcut2,
     #     matcon,
     #     sigsa,sig2lcp,asalcp,bsalcp,
     #     lcp1,lcp2,ilpst,ilpnd,
     #     exvlam,exvcutg,exvcutg2,
     #     sig2exv,asaexv,bsaexv
      save /lundff/
