! *********************************************************************
! This file contains distributeWorkLoad, fileNameMP
!
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! *********************************************************************

!! Calculate the best way to distribute the work load across processors.
!  It calculates the average number of interactions and then tries to 
!  assign a number of interactions to each processor that is as close
!  as possible to the average. The routine should be called once for
!  every molecule in the system.
!
!  @param num_ppr Number of processors per replica
!  @param nml index of molecule or zero.
!
!  @author Jan H. Meinke
      subroutine distributeWorkLoad(num_ppr,nml)

      include 'INCL.H'

      integer num_ppr, nml
      integer idxOfFirstVariable, idxOfLastVariable
      integer at, atct, ivw, i, j, isum, i14
      integer totalct, irank, itarget
      double precision ipps
      
      if (nml.eq.0) then
         idxOfFirstVariable = ivrml1(1)
         idxOfLastVariable = ivrml1(ntlml) + nvrml(ntlml) -1
         i1ms = imsml1(ntlml)+ nmsml(ntlml)
         do i = 0, MAX_PROC
            do j = 1, mxml
               workPerProcessor(j, i) = 0
            end do
         end do
      else 
         idxOfFirstVariable = ivrml1(nml)
         idxOfLastVariable = ivrml1(nml) + nvrml(nml) - 1
         i1ms = imsml1(nml)+ nmsml(nml)
         do i = 0, MAX_PROC
            workPerProcessor(nml, i) = 0
         end do
      end if
      
      isum = 0
      do io = idxOfLastVariable, idxOfFirstVariable, - 1
         iv = iorvr(io)
         i2ms = i1ms - 1
         i1ms = imsvr1(iv)
         do ms = i1ms, i2ms
            do at = latms1(ms), latms2(ms)
               do ivw=ivwat1(at),ivwat2(at)  
                  do j=lvwat1(ivw),lvwat2(ivw)
                     isum = isum + 1
                  end do
               end do
               do i14=i14at1(at),i14at2(at)   
                  isum = isum + 1
               end do
            end do
         end do
      end do
      ipps = isum / num_ppr
      write (*,*) "Total number of interactions:", isum
      write (*,*) "Average # of interactions per processor", ipps
      
      totalct = 0
      irank = 1
      itarget = int(irank * ipps)
      if (nml.eq.0) then
         i1ms = imsml1(ntlml)+ nmsml(ntlml)
      else 
         i1ms = imsml1(nml)+ nmsml(nml)
      end if
      do io = idxOfLastVariable, idxOfFirstVariable, - 1
         isum = 0 
         iv = iorvr(io)
         i2ms = i1ms - 1
         i1ms = imsvr1(iv)
         do ms = i1ms, i2ms
            do at = latms1(ms), latms2(ms)
               atct = atct + 1
               do ivw=ivwat1(at),ivwat2(at)  
                  do j=lvwat1(ivw),lvwat2(ivw)
                     isum = isum + 1
                  end do
               end do
               do i14=i14at1(at),i14at2(at)   
                  isum = isum + 1
               end do
            end do
         end do
         if ((totalct + isum).gt.itarget) then
            if((.not.irank.eq.num_ppr)
     &         .and.
     &         (abs(totalct-itarget)
     &            .lt.abs(totalct + isum - itarget))) then
               workPerProcessor(nml, irank) = io + 1
!                write (*,*) io + 1, totalct, itarget
            else
               workPerProcessor(nml, irank) = io
!                write (*,*) io, totalct + isum, itarget
            end if
            irank = irank + 1
            itarget = int(irank * ipps)
         end if
         totalct = totalct + isum
      end do
      workPerProcessor(nml, 0) = idxOfLastVariable + 1
      workPerProcessor(nml, num_ppr) = ivrml1(nml)

      end subroutine distributeWorkLoad
      
c-----------------------------------------------------------------------
c     The function fileNameMP takes a template of a file name in the
c     variable base. The position of the first and last character that
c     may be replaced by rank in the string are given in i1 (first) and
c     i2 (last).
c     The function returns an empty string if the rank would need more
c     characters than is allowed by the template.
c     For example,
c     \code
c     rank = 11
c     fileName = fileNameMP('base_0000.dat', 6, 9, rank)
c     write (*,*), fileName
c     \endcode
c     will output base_0011.dat.
c     
c     @param base the base file name, e.g., base_0000.dat.
c     @param i1 index of the first character that may be replaced
c     @param i2 index of the last character that may be replaced
c     @param rank the number that should be inserted into the file name.
c     
c     @return file name for rank
c-----------------------------------------------------------------------
      character*80 function fileNameMP(base, i1, i2, rank)

      character*(*) base
c     i1, i2: Index of first and last character that can be replaced
c     rank: rank of node
      integer i1, i2, rank

      fileNameMP = base
      if ((i2 - i1 + 1).le.log10(1.0 * rank)) then
         print *,'too few characters available to label '
         print *,'filenames with rank = ',rank
         stop
      endif

c     TODO: Allow arbitrary rank

      if (rank.lt.10) then
         write(fileNameMP(i2:i2), '(i1)') rank
      elseif (rank.lt.100) then
         write(fileNameMP(i2-1:i2), '(i2)') rank
      elseif (rank.lt.1000) then
         write(fileNameMP(i2-2:i2), '(i3)') rank
      elseif (rank.lt.10000) then
         write(fileNameMP(i2-3:i2), '(i4)') rank
      elseif (rank.lt.100000) then
         write(fileNameMP(i2-4:i2), '(i5)') rank
      elseif (rank.lt.1000000) then
         write(fileNameMP(i2-5:i2), '(i6)') rank
      endif
      end function fileNameMP 
c     End fileNameMP

