c **************************************************************
c
c This file contains the subroutines: init_molecule
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************
c FIXME: Data in varfile determines which molecule is changed.

      subroutine init_molecule(iabin,grpn,grpc,seqfile,varfile)

c ----------------------------------------------------------
c PURPOSE: construct starting structure of molecule(s)
c
c          iabin = 1  : ab Initio using sequence & 
c                       variables given in input files
c          iabin != 1 : sequence, variable information
c                       from PDB-file
c
c          grpn:        N-terminal group
c          grpc:        C-terminal group
c
c CALLS:   addend,bldmol,c_alfa,getmol,iendst, mklist, nursvr,
C          pdbread,pdbvars,redseq,redvar,setmvs
C
c ----------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'

cf2py character*80 optional, intent(in) :: seqfile = ' '
cf2py character*80 optional, intent(in) :: varfile = ' ' 
      
      character grpn*4,grpc*4
      character navr*3, nars*4  
      character seqfile*80, varfile*80
      integer ontlml
      logical readFromStdin

      ontlml = 1
      readFromStdin = .false.

      write (*,*) 'init_molecule: Solvent: ', itysol
      if (iabin.eq.1) then  
         
c     ----------------------------------------- get sequence for molecule(s)
         lunseq=11
         if (ntlml.gt.0) then
            ontlml = ntlml + 1
         endif
         if (iendst(seqfile).le.1.or.seqfile.eq.' ') then
 1          write (*,'(/,a,$)') ' file with SEQUENCE:'
            seqfil=' '
            read (*,'(a)',err=1) seqfil
            readFromStdin = .true.
         else 
            seqfil = seqfile
         endif 
         call redseq
         

         write (*,*) 'File with sequence is ', seqfil(1:iendst(seqfil))
         
c     --------------------------------- read & assemble data from libraries
c     initial coordinates, interaction lists
         
         ntl = ntlml
         do i=ontlml, ntl
            
            call getmol(i)      ! assemble data from libraries
            
            do j=1,6            ! initialize global parameters
               gbpr(j,i)=0.d0
            enddo
            
            call bldmol(i)      ! co-ordinates
            
            ntlml = i
            call addend(i,grpn,grpc) ! modify ends
            call setmvs(i) ! determine sets of moving atoms for given variables
            call mklist(i)      ! compile lists of interaction partners
            
         enddo
         
c     --------------------------- Read the initial conformation if necessary
         if(readFromStdin) then 
            write (*,'(a,$)') ' file with VARIABLES:'
c     
            varfil=' '
            read(*,'(a)',end=2,err=2) varfil
         else
            varfil = varfile
         endif
         l=iendst(varfil)
         if (l.gt.0.and.varfil.ne.' ') then
            write (*,'(1x,a,/)') varfil(1:l)
            lunvar=13
            
            call redvar         ! get vars. and rebuild
            
         endif
         
 2       write(*,*) ' '
         
c     -------------------- get: nvr,idvr, vlvr, olvlvr
         nvr = 0
         do i=1,ivrml1(ntlml)+nvrml(ntlml)-1
            
            if (.not.fxvr(i)) then
               nvr=nvr+1
               idvr(nvr)=i      ! index of not fixed var.
            endif
            
            it=ityvr(i)
            
            if (it.eq.3) then   ! torsion
               vlvr(i)=toat(iatvr(i))
            elseif (it.eq.2) then ! b.angle
               vlvr(i)=baat(iatvr(i))
            elseif (it.eq.1) then ! b.length
               vlvr(i)=blat(iatvr(i))
            endif

            olvlvr(i) = vlvr(i) 
         enddo
         
         ireg = 0
         
      else                      ! =========================== from PDB
         if (iendst(seqfile).le.1) then
 3          write (*,'(/,a,$)') ' PDB-file:'
            seqfil=' '
            read (*,'(a)',err=3) seqfil
         else
            seqfil = seqfile
         endif
         write (*,*) 'PDB structure ',seqfil(1:iendst(seqfil))
         print *, 'calling readpdb with ',seqfile
         call pdbread(seqfil,ier)
         
         if (ier.ne.0) stop
         
         call pdbvars()
         
         ireg = 1
         
      endif
      
c     -------------------------- set var. amplitudes for simulations
      
      do i=1,ivrml1(ntlml)+nvrml(ntlml)-1
         
         if (ityvr(i).eq.3.and..not.fxvr(i)) then ! torsion
            
            navr = nmvr(i)
            
            ir = nursvr(i)
            nars = seq(ir)
            
            if (                         navr(1:2).eq.'om'
            
     #     .or.nars(1:3).eq.'arg'.and.(navr(1:2).eq.'x5'
     #           .or.navr(1:2).eq.'x6')
            
     #           .or.(nars(1:3).eq.'asn'.or.nars(1:3).eq.'asp')
     #           .and.navr(1:2).eq.'x3'
            
     #           .or.(nars(1:3).eq.'gln'.or.nars(1:3).eq.'glu')
     #           .and.navr(1:2).eq.'x4'
            
     #           ) then
            
c     axvr(i) = pi/9.d0  ! 20 deg.
            axvr(i) = pi2       ! Trying out 360 deg. for these as well
            
         else
            axvr(i) = pi2       ! 360 deg.
         endif
         
      else
         axvr(i) = 0.d0
      endif
      
      enddo                     ! vars.
      
c     --------------------- initialize solvation pars. if necessary

      if (itysol.ne.0) then
         
         i1=iatrs1(irsml1(1))   ! 1st atom of 1st molecule
         i2=iatrs2(irsml2(ntlml)) ! last atom of last molecule
         
         its = iabs(itysol)
         
         do i=i1,i2             ! all atoms
            it=ityat(i)
            sigma(i)=coef_sl(its,it)
            rvdw(i) =rad_vdw(its,it)
            
            if (nmat(i)(1:1).ne.'h') rvdw(i)=rvdw(i)+rwater
            
         enddo
         
      endif
! Initialize calpha array
      do i=ontlml, ntlml
         call c_alfa(i,1)    
      enddo

!     Initialize arrays used in the BGS update
      call init_lund()      
      return
      end
      
      
