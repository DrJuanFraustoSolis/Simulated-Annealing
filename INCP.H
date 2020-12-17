c ......................
c contents of a PDB file
c ......................

      parameter (MXCHP =100,    ! max. no. of polypeptide chains
     #           MXRSP =1000,   ! max. no. of residues
     #           MXATP =10000)  ! max. no. of atoms

c   nchp      - no. of polypeptide chains
c   nchrsp()  - no. of residues / chain
c   nrsp      - total no. of residues
c   irsatp()  - index of 1st atom / res.
c   nrsatp()  - no. of atoms / res.
c   natp      - total no. of selected atoms
c   noatp()   - atom numbers

      common /pdb_i/ nchp,nchrsp(MXCHP),
     #               nrsp,irsatp(MXRSP),nrsatp(MXRSP),
     #               natp,noatp(MXATP)
      save /pdb_i/

c   chnp()  - chain identifiers
c   rsidp() - residue identifiers (number + insertion code)
c   rsnmp() - residues (sequence, 3-letter code)
c   atnmp() - atom names

      character chnp(MXCHP),
     #          rsidp(MXRSP)*5,rsnmp(MXRSP)*3,
     #          atnmp(MXATP)*4

      common /pdb_c/ chnp,rsnmp,rsidp,atnmp
      save /pdb_c/

c     xatp,yatp,zatp - atom coordinates

      common /pdb_r/ xatp(MXATP),yatp(MXATP),zatp(MXATP)
      save /pdb_r/

c ------------------- code to list all PDB information
c      ir=0
c      do i=1,nchp
c        write(*,*) ' ===== chain |',chnp(i),'|'
c        do j=1,nchrsp(i)
c          ir=ir+1
c          write(*,*) ' ----- ',rsidp(ir),' ',rsnmp(ir),' ',nrsatp(ir)
c          k1=irsatp(ir)
c          k2=k1+nrsatp(ir)-1
c          do k=k1,k2
c            write(*,*) ' ',noatp(k),' ',atnmp(k),' ',(xyzp(l,k),l=1,3)
c          enddo
c        enddo
c      enddo
