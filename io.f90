!**********************************************************************
! Copyright 2001, 2004, 2007, 2012, 2013, 2014                        *
! Bernard Legras                                                      *
!                                                                     *
! This file is part of TRACZILLA which is derived from FLEXPART V6    *
!                                                                     *
! TRACZILLA is free software: you can redistribute it and/or modify   *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! TRACZILLA is distributed in the hope that it will be useful,        *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with TRACZILLA.  If not, see <http://www.gnu.org/licenses/>.  *
!**********************************************************************
!
!###############################################################################
!------------------------------------- IO --------------------------------------
!###############################################################################
module io
use commons
use date
implicit none
public :: partout_agef, partout_fast, partout_qv, savsav, restartfromsav
private :: writesngl

contains

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PARTOUT_FAST @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
      
      subroutine partout_fast(itime)
      
!*******************************************************************************
!                                                                              *
!     Dump all particle positions in 32 bits IEEE                              *
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     05 October 2001                                                          *
!     modifified: 10 January 2004      
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer itime,i1,i2,j,jseg,if,i,outfmt,lhead
      integer ihour,jjjjmmdd,ihmmss,nbuf,try,ierr      
      integer (kind=4), allocatable:: itra(:)
      real (kind=4), allocatable:: xyz(:,:)
      character(len=6):: chour

      ihour=abs(itime)/3600
      if(ihour < 1000) then
        write(chour,'(I3.3)') ihour
      else if(ihour < 10000) then
        write(chour,'(I4.4)') ihour
      else if(ihour < 100000) then
        write(chour,'(I5.5)') ihour
      else
        write(chour,'(I6.6)') ihour
      endif
   
! Output all particles  
!*************************************

      outfmt=102        ! Output format
      lhead = 6         ! Header length (# of lines)
      nbuf=100000
      try=1
 110  open(unitpartout,file=trim(path(2))//'part_'//trim(chour),  &
          form='unformatted',ERR=150)
      call writesngl(unitpartout,(/ lhead, outfmt /),ierr)
      if(ierr==1) go to 150
      call writesngl(unitpartout,(/ nbuf /),ierr)
      if(ierr==1) go to 150
      call writesngl(unitpartout,(/ ibdate, ibtime /),ierr)
      if(ierr==1) go to 150
      call caldate(bdate,jjjjmmdd,ihmmss)
      call writesngl(unitpartout,(/ jjjjmmdd, ihmmss /),ierr)
      if(ierr==1) go to 150
      call writesngl(unitpartout,(/ numpart,numlevel(1),nsample(1) /),ierr)
      if(ierr==1) go to 150
      call writesngl(unitpartout,(/ itime /),ierr)
      if(ierr==1) go to 150
      allocate (xyz(3,nbuf),itra(nbuf))
      jseg=numpart/nbuf 
      i1=1 ; i2=nbuf ; if=nbuf
      do j=1, jseg       
         xyz(1,1:if)=xlon0+xtra1(i1:i2)*dx
         xyz(2,1:if)=ylat0+ytra1(i1:i2)*dy
         xyz(3,1:if)=ztra1(i1:i2)
         itra( 1:if)=itra1(i1:i2)
         write(unitpartout,ERR=150) (xyz(1:3,i),itra(i),i=1,if)
         i1=i1+if ; i2=i2+if
      enddo
      if(if*jseg<numpart) then
        i2=numpart ; if=i2-i1+1
        xyz(1,1:if)=xlon0+xtra1(i1:i2)*dx
        xyz(2,1:if)=ylat0+ytra1(i1:i2)*dy
        xyz(3,1:if)=ztra1(i1:i2)
        itra( 1:if)=itra1(i1:i2)
        write(unitpartout,ERR=150)   &
            (xyz(1:3,i),itra(i),i=1,if)
      endif
      print *,'partout_fast> part_',chour,'  ',itime/86400, &
         ' days, format 102'
      close(unitpartout)
      deallocate (xyz,itra)
      return
 150  close(unitpartout)
      print *,'partout_fast> ERROR in opening or writing'
      try=try+1
      if(try==5) then
        print *, ' PARTOUT_FAST'
        print *, ' give up after 5 attempt to write'
        stop 1000
      else
        go to 110
      endif
      end subroutine partout_fast

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PARTOUT_AGEF @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      subroutine partout_agef(itime)

!*******************************************************************************
!                                                                              *
!     Dump all particle positions in 32 bits IEEE
!     plus temperature
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     derived form partout_qv
!     2 September 2012     
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in):: itime
      integer i,outfmt,lhead
      integer ihour,nact     
      integer (kind=4), allocatable:: itra(:),idx(:)
      real(kind=4), allocatable:: xyzt(:,:)
      character(len=6):: chour

      ihour=abs(itime)/3600
      if(ihour < 1000) then
        write(chour,'(I3.3)') ihour
      else if(ihour < 10000) then
        write(chour,'(I4.4)') ihour
      else if(ihour < 100000) then
        write(chour,'(I5.5)') ihour
      else
        write(chour,'(I6.6)') ihour
      endif

!     Generates a reduced file with all the active parcels
      nact=0
      allocate (itra(numpart),idx(numpart),xyzt(4,numpart))
      dopart: do i=1,numpart
        if(itra1(i)==itime) then
          nact=nact+1
          idx(nact)=i
          itra(nact)=itra0(i)
          xyzt(1,nact)=xlon0+xtra1(i)*dx
          xyzt(2,nact)=ylat0+ytra1(i)*dy
          xyzt(3,nact)=ztra1(i)
          ! fix for the case of unallocated ttra1
          if (allocated(ttra1)) then
             xyzt(4,nact)=ttra1(i)
          else
             xyzt(4,nact)=MISSING
          endif
        endif
      enddo dopart
      
! Output all particles  
!*************************************

      outfmt=105        ! Output format
      lhead = 5        ! Header length (# of lines)
      open(unitpartout,file=trim(path(2))//'part_'//trim(chour),  &
          form='unformatted')
      call writesngl(unitpartout,(/ lhead, outfmt /))
      call writesngl(unitpartout,(/ nact, numpart/))
      write(unitpartout) release_plan
      call writesngl(unitpartout,(/ ibdate, ibtime /))
      call writesngl(unitpartout,(/ itime /))
      write(unitpartout) (xyzt(1:4,i),idx(i),itra(i),i=1,nact)
      close(unitpartout)
      print *,'partout_agef> part_',chour,'  ',nact,' particles  ', &
              itime/86400,' days'
      deallocate (xyzt,itra,idx)
      return
      end subroutine partout_agef
      
!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PARTOUT_STC  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      subroutine partout_stc(itime)

!*******************************************************************************
!     
!     Format 107, designed for StratoClim                                                                         *
!     Dump all particle positions in 32 bits IEEE
!     plus temperature, flag, initial particle pointer, and initial time
!     
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     derived from partout_agef
!     23 April 2016    
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in):: itime
      integer  :: i, ihour, nact
      integer (kind=4) :: outfmt,lhead,mode, numpart_s, nact_s, &
                          lsynctime_s,itime_s, idx_orgn_s                 
      integer (kind=4), allocatable:: itra(:),idx(:),flag(:)
      real(kind=4), allocatable:: xyzt(:,:)
      character(len=6):: chour

      ihour=abs(itime)/3600
      if(ihour < 1000) then
        write(chour,'(I3.3)') ihour
      else if(ihour < 10000) then
        write(chour,'(I4.4)') ihour
      else if(ihour < 100000) then
        write(chour,'(I5.5)') ihour
      else
        write(chour,'(I6.6)') ihour
      endif

!     Generates a reduced dataset with all the active parcels
      nact=0
      allocate (itra(numpart),idx(numpart),flag(numpart),xyzt(4,numpart))
      dopart: do i=1,numpart
        if(itra1(i)==itime) then
          nact=nact+1
          if (M10) then
            ! This generate a segmentation fault if M1O is set on a
            ! restart run    
            idx(nact)=idxcross(i)
          else 
            idx(nact)=i+idx_orgn-1
          endif
          itra(nact)=itra0(i)
          flag(nact)=flagTB(i)
          xyzt(1,nact)=xlon0+xtra1(i)*dx
          xyzt(2,nact)=ylat0+ytra1(i)*dy
          xyzt(3,nact)=ztra1(i)
          xyzt(4,nact)=ttra1(i)
        endif
      enddo dopart
!     Transforms vertical coordinate into pressure 
      if(z_motion) then
        xyzt(3,1:nact)=p0*exp(-xyzt(3,1:nact))
      else
        xyzt(3,1:nact)=p0*(xyzt(4,1:nact)/xyzt(3,1:nact))**(1/kappa)
      endif
       
! Output all particles  
!*************************************

      outfmt=107       ! Output format
      lhead = 3        ! Header length (# of lines)
      mode = 1
      open(unitpartout,file=trim(path(2))//'part_'//trim(chour),  &
          form='unformatted')
      write(unitpartout) lhead,outfmt,mode  
      itime_s=itime
      lsynctime_s=lsynctime
      write(unitpartout) stamp_date, itime_s, lsynctime_s  
      numpart_s=numpart
      nact_s=nact
      idx_orgn_s=idx_orgn
      write(unitpartout) numpart_s, nact_s, idx_orgn_s, 0, 0, 0   
      write(unitpartout) flag(1:nact)
      write(unitpartout) itra(1:nact)
      write(unitpartout) xyzt(1,1:nact)
      write(unitpartout) xyzt(2,1:nact)
      write(unitpartout) xyzt(3,1:nact)
      write(unitpartout) xyzt(4,1:nact)    
      write(unitpartout) idx(1:nact)
      close(unitpartout)
      print *,'partout_stc> part_',chour,'  ',nact,' particles  ', &
              itime/86400,' days'
      deallocate (xyzt,itra,idx,flag)
      return
      end subroutine partout_stc

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PARTOUT_QV @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8

      subroutine partout_qv(itime)

!*******************************************************************************
!                                                                              *
!     Dump all particle positions in 32 bits IEEE
!     plus temperature and saturation ratio
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     derived form partout_fast
!     2 February 2004     
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in):: itime
      integer i1,i2,j,jseg,iff,i,outfmt,lhead
      integer ihour,jjjjmmdd,ihmmss,nbuf      
      integer (kind=4), allocatable:: itra(:)
      real(kind=4), allocatable:: xyztq(:,:)
      character(len=6):: chour

      ihour=abs(itime)/3600
      if(ihour < 1000) then
        write(chour,'(I3.3)') ihour
      else if(ihour < 10000) then
        write(chour,'(I4.4)') ihour
      else if(ihour < 100000) then
        write(chour,'(I5.5)') ihour
      else
        write(chour,'(I6.6)') ihour
      endif

! Output all particles  
!*************************************

      outfmt=103        ! Output format
      lhead = 6         ! Header length (# of lines)
      nbuf=300000
      open(unitpartout,file=trim(path(2))//'part_'//trim(chour),  &
          form='unformatted')
      call writesngl(unitpartout,(/ lhead, outfmt /))
      call writesngl(unitpartout,(/ nbuf /))
      call writesngl(unitpartout,(/ ibdate, ibtime /))
      call caldate(bdate,jjjjmmdd,ihmmss)
      call writesngl(unitpartout,(/ jjjjmmdd, ihmmss /))
      call writesngl(unitpartout,(/ numpart,n_loc,n_sample /))
      call writesngl(unitpartout,(/ itime /))
      allocate (xyztq(5,nbuf),itra(nbuf))
      jseg=numpart/nbuf 
      i1=1 ; i2=nbuf ; iff=nbuf
      do j=1, jseg       
         xyztq(1,1:iff)=xlon0+xtra1(i1:i2)*dx
         xyztq(2,1:iff)=ylat0+ytra1(i1:i2)*dy
         xyztq(3,1:iff)=ztra1(i1:i2)
         xyztq(4,1:iff)=ttra1(i1:i2)
         xyztq(5,1:iff)=qtra1(i1:i2)
         itra( 1:iff)=itra1(i1:i2)
         write(unitpartout) (xyztq(1:5,i),itra(i),i=1,iff)
         i1=i1+iff ; i2=i2+iff
      enddo
      if(iff*jseg<numpart) then
        i2=numpart ; iff=i2-i1+1
        xyztq(1,1:iff)=xlon0+xtra1(i1:i2)*dx
        xyztq(2,1:iff)=ylat0+ytra1(i1:i2)*dy
        xyztq(3,1:iff)=ztra1(i1:i2)
        xyztq(4,1:iff)=ttra1(i1:i2)
        xyztq(5,1:iff)=qtra1(i1:i2)
        itra( 1:iff)=itra1(i1:i2)
        write(unitpartout) (xyztq(1:5,i),itra(i),i=1,iff)
      endif
      close(unitpartout)
      print *,'partout_qv> part_',chour,'  ',itime/86400, &
           ' days, format 103'
      deallocate (xyztq,itra)
      return
      end subroutine partout_qv

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SAVSAV @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
      
      subroutine savsav(itime)

!*******************************************************************************
!                                                                              *
!     Dump all particle positions to save a state for restart                  *
!     Dump is first done to a temporary file and then moved to replace
!     the true save file replacing existing one.
!     This is to avoid failure to restart if NFS error occurs during
!     dump leaving the file incomplete
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     26 January 2004
!     redone : 2 September 2012 as a raw format output, no buffering 
!     and no transformation                                                         *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in):: itime
      integer :: j,outfmt,lhead
      integer :: jjjjmmdd,ihmmss      
      integer :: system
      !logical, allocatable :: mask_part(:)

! Output all particles in raw format  
!*************************************

      outfmt=110     ! Output format
      lhead =6       ! Header length (# of lines)
      open(saveunit_tmp,file=trim(path(2))//'savpos_tmp', &
         form='unformatted')
      write(saveunit_tmp) lhead,outfmt
      write(saveunit_tmp) numpart, release_plan, idx_orgn
      write(saveunit_tmp) ibdate,ibtime
      call caldate(bdate,jjjjmmdd,ihmmss)
      write(saveunit_tmp) jjjjmmdd,ihmmss
      write(saveunit_tmp) itime
      write(saveunit_tmp) allocated(itra0),allocated(ttra1),allocated(qtra1),allocated(flagTB)
      write(saveunit_tmp) xtra1(1:numpart)
      write(saveunit_tmp) ytra1(1:numpart)
      write(saveunit_tmp) ztra1(1:numpart)
      write(saveunit_tmp) itra1(1:numpart)
      if (allocated(itra0)) write(saveunit_tmp) itra0(1:numpart)
      if (allocated(ttra1)) write(saveunit_tmp) ttra1(1:numpart)
      if (allocated(qtra1)) write(saveunit_tmp) qtra1(1:numpart)
      if (allocated(flagTB)) write(saveunit_tmp) flagTB(1:numpart)
      close(saveunit_tmp)
      ! test
      !allocate (mask_part(numpart))
      !mask_part(:) = itra1(:)==itime
      !print *,sum(int(mask_part)),itime,minval(itra1),maxval(itra1)
      !print *,minval(itra0),maxval(itra0)
      !print *,minval(xtra1(1:numpart),DIM=1,MASK=mask_part),&
      !        maxval(xtra1(1:numpart),DIM=1,MASK=mask_part)
      !print *,minval(ytra1(1:numpart),DIM=1,MASK=mask_part),&
      !        maxval(ytra1(1:numpart),DIM=1,MASK=mask_part)        
      !print*,'ztra1 ', minval(ztra1(1:numpart),DIM=1,MASK=mask_part), &
      !                 maxval(ztra1(1:numpart),DIM=1,MASK=mask_part), &
      !                 sum(ztra1(1:numpart),DIM=1,MASK=mask_part)/sum(int(mask_part))
      !deallocate(mask_part)                 

! Moves the file to its true name after successful write

      j=system('mv -f '//trim(path(2))//'savpos_tmp'//' '    &
                       //trim(path(2))//'savpos')
!     if (system('mv -f '//trim(path(2))//'savpos_tmp'//' '    &
!                        //trim(path(2))//'savpos')) then
         print *, 'checkpoint at t =',itime

      return
      end subroutine savsav

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ RESTARTFROMSAV @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
      
      subroutine restartfromsav(error)
!                          
!*******************************************************************************
!                                                                              *
!     Restart the calculation from saved archive                               *
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     30 January 2004  
!     2 September 2012: new version with raw blocks                                                        *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer :: ibdate_arch,ibtime_arch,lhead,outfmt
      character(len=12):: r_plan
      logical, intent(out):: error
      logical :: itra0_activ,ttra1_activ,qtra1_activ,flagTB_activ
      !logical, allocatable :: mask_part(:)

      error=.false.
      print *,'##### restart run ####' 
      
! Open the file  
!*************************************
      open(unitpartin,file=trim(path(2))//'savpos', &
         form='unformatted',status='old')
      print *,'restart > opened file'
      read(unitpartin) lhead, outfmt
! this sequences allow to read new savpos with idx_orgn 
! and set idx_orgn to 1 for the old version without orgn     
      read(unitpartin,err=45) numpart, r_plan, idx_orgn
      go to 50
  45  close(unitpartin)
      !unitpartin = unitpartin+100
      open(unitpartin,file=trim(path(2))//'savpos', &
         form='unformatted',status='old')
      print *,'restart > re-opened file'
      read(unitpartin) lhead, outfmt
      read(unitpartin) numpart, r_plan
      print *,'numpart ,r_plan ',numpart,r_plan
      idx_orgn = 1
      print *,'ACHTUNG IDX_ORGN SET TO 1'
  50  continue
      read(unitpartin) ibdate_arch,ibtime_arch 
      read(unitpartin) 
      read(unitpartin) itime0
      read(unitpartin) itra0_activ,ttra1_activ,qtra1_activ,flagTB_activ
      print *,'restart from sav flags itra0 ttra1 qtra1 flagTB'
      print *,itra0_activ,ttra1_activ,qtra1_activ,flagTB_activ
      if(.not.allocated(xtra1)) then
         allocate(xtra1(numpart),ytra1(numpart),ztra1(numpart),&
                  itra1(numpart))
      endif
      read(unitpartin) xtra1(1:numpart)            
      read(unitpartin) ytra1(1:numpart)
      read(unitpartin) ztra1(1:numpart)
      read(unitpartin) itra1(1:numpart)
      if (itra0_activ) then
         if(.not.allocated(itra0)) allocate(itra0(numpart))
         read(unitpartin) itra0(1:numpart)
      endif
!     CHECK WHETHER ALLOCATION IS NEEDED FOR THE OTHER TWO
      if (ttra1_activ) then
         if(.not.allocated(ttra1)) allocate(ttra1(numpart))
         read(unitpartin) ttra1(1:numpart)
      endif
      if (qtra1_activ) then
         if(.not.allocated(qtra1)) allocate(qtra1(numpart))
         read(unitpartin) qtra1(1:numpart)
      endif
      if (flagTB_activ) then
         if(.not.allocated(flagTB)) allocate(flagTB(numpart))
         read(unitpartin) flagTB(1:numpart)
      endif
      !allocate (mask_part(numpart))
      !mask_part(:) = itra1(:)==itime0
      !print *,sum(int(mask_part)),itime0,maxval(itra1)
      !print *,minval(itra0),maxval(itra0)
      !print *,minval(xtra1(1:numpart),DIM=1,MASK=mask_part),&
      !        maxval(xtra1(1:numpart),DIM=1,MASK=mask_part)
      !print *,minval(ytra1(1:numpart),DIM=1,MASK=mask_part),&
      !        maxval(ytra1(1:numpart),DIM=1,MASK=mask_part)      
      !print*,'ztra1 ', minval(ztra1(1:numpart),DIM=1,MASK=mask_part), &
      !                 maxval(ztra1(1:numpart),DIM=1,MASK=mask_part) , &
      !                 sum(ztra1(1:numpart),DIM=1,MASK=mask_part)/sum(int(mask_part))
      !deallocate(mask_part)                       
      close(unitpartin)
      return
      end subroutine restartfromsav

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ RESTARTPART @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8      

      subroutine restartpart(error)
!                          
!*******************************************************************************
!                                                                              *
!     Restart the calculation from saved archive in 32 bits IEEE
!     or from the last saved point
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     27 June   2002
!     modified: 12 January 2004, 26 January 2004, 2 September 2012
!
!     note that for the format 103, qvs is not read, hence the run cannot
!     proceed in the same way as if it was not interrupted
!     pb: qvs is only defined in timemanager 
!     See notes on TRACZILLA
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer :: i1,i2,if,i,j,jseg,nbuf,shift
      integer :: ibdate_arch,ibtime_arch,lhead,outfmt,nact
      integer (kind=4), allocatable:: itra(:),idx(:),flag(:)
      real (kind=4), allocatable:: xyz(:,:)
      integer (kind=4):: ddin(3)
      character(len=6):: chour
      character(len=12):: r_plan
      logical,intent(out):: error
      logical :: ttra1_activ

      error=.false.
      print *,'##### restart run ####'
      
      if(AGEBactiv) ttra1_activ=.true. 
   
! Open the file  
!*************************************
 
      if(hrstart == 0) then 
        call restartfromsav(error)
        print *,'read sav file' 
        if (track_kill) then
          call readkill(error)
          print *,'read kill file'
        endif
        if (track_cross) then
          call readcross(error)
          print *,'read cross file'
        endif
        return
      endif
      ! This trick allows to start from 0 bypassing previous test
      if(hrstart == -1) hrstart=0 
      if(hrstart < 1000) then
        write(chour,'(I3.3)') hrstart
      else if(hrstart < 10000) then
        write(chour,'(I4.4)') hrstart
      else if(hrstart < 100000) then
        write(chour,'(I5.5)') hrstart
      else
        write(chour,'(I6.6)') hrstart
      endif
      open(unitpartin,file=trim(path(2))//'part_'//trim(chour), &
         form='unformatted',status='old')
      print *,'restart > opened file'
      read(unitpartin) ddin(1), ddin(2)
      lhead = ddin(1); outfmt = ddin(2)
      if (lhead > 1000) then
        print *,'restart > old raw format detected'
        ibdate_arch = lhead ; ibtime_arch = outfmt
        read(unitpartin)
        read(unitpartin)
        read(unitpartin)
        read(unitpartin) 
        numpart = ddin(1); numlevel(1) = ddin(2); nsample(1) = ddin(3)
        read(unitpartin) ddin(1)
        itime0 = ddin(1)
        nbuf=1000
        outfmt=100
      else
        select case (outfmt)
        case (101)
          nbuf=1000
          read(unitpartin) ddin(1), ddin(2)
          ibdate_arch = ddin(1); ibtime_arch = ddin(2)
          read(unitpartin)
          read(unitpartin)
          read(unitpartin)
          read(unitpartin) ddin(1), ddin(2), ddin(3) 
          numpart = ddin(1); numlevel(1) = ddin(2); nsample(1) = ddin(3)
          read(unitpartin) ddin(1)
          itime0 = ddin(1)
        case (102,103)
          read(unitpartin) ddin(1)
          nbuf = ddin(1)
          read(unitpartin) ddin(1), ddin(2) 
          ibdate_arch = ddin(1); ibtime_arch = ddin(2)
          read(unitpartin)
          read(unitpartin) ddin(1), ddin(2), ddin(3)
          numpart = ddin(1); numlevel(1) = ddin(2); nsample(1) = ddin(3)
          read(unitpartin) ddin(1)
          itime0 = ddin(1)
        case (105)
          read(unitpartin) ddin(1), ddin(2)
          nbuf = ddin(1)
          numpart = ddin(2)
          read(unitpartin) r_plan
          read(unitpartin) ddin(1), ddin(2)
          ibdate_arch = ddin(1); ibtime_arch = ddin(2)
          read(unitpartin) ddin(1)
          itime0=ddin(1)
          print *,numpart,nbuf,itime0
        case(107)
          rewind(unitpartin)
          read(unitpartin) ddin(1:3)
          if ((ddin(2) /=107).or.((ddin(3) /=1) .and. (ddin(3) /= 3))) then
            print *,'restartpart ERROR'
            print *,'fmt and mode are ',ddin(2),ddin(3)
            print *,'expected 107 and 1'
            stop
          endif
          read(unitpartin) stamp_date, ddin(1:2)
          itime0=ddin(1)
          print *,stamp_date,itime0
          read(unitpartin) ddin(1:3)
          numpart=ddin(1)
          nact=ddin(2)
          idx_orgn=ddin(3)
          print *,numpart,nact,idx_orgn                
        case default
          print *,'unknown format ',outfmt
          error=.true.
          return
        end select
      endif
      if((ibdate_arch /= ibdate).or.(ibtime_arch /= ibtime)) then
        print *,'restart > ibdate or itime non matching'
        print *,ibdate,ibdate_arch
        print *,ibtime,ibtime_arch
        ! false alarm here for case 107, 
!	error=.true.
!	return
      endif
      select case (outfmt)
      case (100,101,102)
        allocate (xyz(3,nbuf),itra(nbuf))
        jseg=numpart/nbuf
        i1=1 ; i2=nbuf ; if=nbuf
        if(jseg>0) then
        do j=1, jseg
          read(unitpartin)  (xyz(1:3,i),itra(i),i=1,if)
          xtra1(i1:i2)=(xyz(1,1:if)-xlon0)/dx
          ytra1(i1:i2)=(xyz(2,1:if)-ylat0)/dy
          ztra1(i1:i2)= xyz(3,1:if)
          itra1(i1:i2)=itra(  1:if)
          i1=i1+if ; i2=i2+if
        enddo
        endif
        if(if*jseg<numpart) then
          i2=numpart ; if=i2-i1+1
          read(unitpartin) (xyz(1:3,i),itra(i),i=1,if)
          xtra1(i1:i2)=(xyz(1,1:if)-xlon0)/dx
          ytra1(i1:i2)=(xyz(2,1:if)-ylat0)/dy
          ztra1(i1:i2)= xyz(3,1:if)
          itra1(i1:i2)=itra(  1:if)
        endif
        deallocate (xyz,itra)
        
      case (103)
        allocate (xyz(5,nbuf),itra(nbuf))
        if(.not.allocated(qtra1)) allocate (qtra1(maxpart))
        jseg=numpart/nbuf
        i1=1 ; i2=nbuf ; if=nbuf
        print *,numpart,nbuf,jseg 
        if(jseg>0) then
        do j=1, jseg
          read(unitpartin)  (xyz(1:5,i),itra(i),i=1,if)
          xtra1(i1:i2)=(xyz(1,1:if)-xlon0)/dx
          ytra1(i1:i2)=(xyz(2,1:if)-ylat0)/dy
          ztra1(i1:i2)= xyz(3,1:if)
          qtra1(i1:i2)  = xyz(5,1:if)
          itra1(i1:i2)=itra(  1:if)
          i1=i1+if ; i2=i2+if
        enddo
        endif
        if(if*jseg<numpart) then
          i2=numpart ; if=i2-i1+1
          read(unitpartin) (xyz(1:5,i),itra(i),i=1,if)
          xtra1(i1:i2)=(xyz(1,1:if)-xlon0)/dx
          ytra1(i1:i2)=(xyz(2,1:if)-ylat0)/dy
          ztra1(i1:i2)= xyz(3,1:if)
          qtra1(i1:i2)= xyz(5,1:if)
          itra1(i1:i2)=itra(  1:if)
        endif
        deallocate (xyz,itra)
        
      case (105)
        print *,'DO NOT USE THIS RESTART PROCEDURE THAT NEEDS CORRECTION'
        print *,'numpart should not be the number of active parcels'
        print *,'unless the index is preserved'
        print *,'that requires modifications'
        STOP   
        if(.not.allocated(itra0)) allocate(itra0(numpart))
!       read initial positions for all parcels
        open(unitpartout4,file=trim(path(2))//'sav_init', &
             form='unformatted',status='old')
        read(unitpartout4)xtra1(1:numpart)
        read(unitpartout4)ytra1(1:numpart)
        read(unitpartout4)ztra1(1:numpart)
        read(unitpartout4)itra0(1:numpart)
        if (ttra1_activ) then
         if(.not.allocated(ttra1)) allocate(ttra1(numpart))
         read(unitpartout4) ttra1(1:numpart)
        endif
        ! following statement was wrong as it makes no sense
        ! corrected on 26/12/2016      
        !itra1(1:numpart)=itra0(1:numpart)
        itra1(1:numpart) = itime0
        print *,'ACHTUNG ACHTUNG ACHTUNG ACHTUNG ACHTUNG'
        print *,'Beware that restarting this way  you have lost all'
        print *,'the parcels still to be launched in that run'
        print *,'Use retsart from save to extend age runs'
        close(unitpartout4)
!       read now positions of active parcels over the initial positions
!       parcels already deactivated are left with their initial position
!       and not their final position as they were at the end of the previous
!       run but this is without consequence since they are no longer accessed  
        allocate (xyz(4,nbuf),itra(nbuf),idx(nbuf))
        read(unitpartin) (xyz(1:4,i),idx(i),itra(i),i=1,nbuf)
        do i=1,nbuf
          j=idx(i)
          xtra1(j)=(xyz(1,i)-xlon0)/dx
          ytra1(j)=(xyz(2,i)-ylat0)/dy
          ztra1(j)=xyz(3,i)
          itra1(j)=itra(i) ! itra(i:nbuf) should be equal to itime0
        enddo
        if (ttra1_activ) then
          do i=1,nbuf
            ttra1(idx(i))=xyz(4,i)
          enddo
        endif
        deallocate (xyz,itra,idx) 
        
      case (107)
        print *,'restart > start reading fields'
        allocate (xyz(4,nact),itra(nact),idx(nact),flag(nact))
        if(.not.allocated(xtra1)) then
          allocate(xtra1(numpart),ytra1(numpart),ztra1(numpart))
          allocate(itra0(numpart),itra1(numpart))
        endif
        if(.not.allocated(ttra1)) allocate(ttra1(numpart))
        if(.not.allocated(flagTB)) allocate(flagTB(numpart))
        read(unitpartin) flag(1:nact)
        read(unitpartin) itra(1:nact)
        read(unitpartin) xyz(1,1:nact)
        read(unitpartin) xyz(2,1:nact)
        read(unitpartin) xyz(3,1:nact)
        read(unitpartin) xyz(4,1:nact)
        read(unitpartin) idx(1:nact)
        if(z_motion) then
          xyz(3,:)= -log(xyz(3,:))+log(p0)
        else
          xyz(3,:)=xyz(4,:)*(p0/xyz(3,:))**kappa      
        endif
        itra1(:)=HUGE(1)
        shift = 1 - idx_orgn
        do i=1,nact
          flagTB(idx(i)+shift)=flag(i)
          !itra0(idx(i)+shift)=itra(i)
          !itra1(idx(i)+shift)=itime0
          xtra1(idx(i)+shift)=(xyz(1,i)-xlon0)/dx
          ytra1(idx(i)+shift)=(xyz(2,i)-ylat0)/dy
          ttra1(idx(i)+shift)=xyz(4,i)
          ztra1(idx(i)+shift)=xyz(3,i)
        enddo
        print *,'x ',minval(xyz(1,:)),maxval(xyz(1,:))
        print *,'y ',minval(xyz(2,:)),maxval(xyz(2,:))
        ! test whether there are some values of xtra1 outside limits
        ! and repair if we are in the global domain
        if ((minval(xtra1)<0).or.(maxval(xtra1)>=float(nx-1))) then
          if (.not.xglobal) then
            print *,'ACHTUNG restart with invalid longitudes'
            stop 255
          else
            print *,'restart with longitudes using other origin'
            xtra1 = modulo(xtra1,real(nx-1,kind(0._dp)))
          endif
        endif
            
        ! If we are in a run starting from part_000 
        ! itra1 should be equal to itra0      
        ! In a restart run, the active parcels are
        ! all to be initialized at the same time.
        ! since they have been saved only if they were active.
        ! ACHTUNG: This produces wrong results if there 
        ! are still parcels to be launched, e.g. M10
        ! runs 
        ! The test on hrstart makes sure this is used only at 0 time
        ! (actually, the startfrom0 flag is somewhat superfluous)
        if (startfrom0.and.(hrstart==0)) then
          do i=1,nact
            itra1(idx(i)+shift)=itra(i)
            itra0(idx(i)+shift)=itra(i)
          enddo
        else 
          do i=1,nact
            itra1(idx(i)+shift)=itime0
          enddo  
          print *,'ACHTUNG ACHTUNG ACHTUNG ACHTUNG ACHTUNG'
          print *,'Beware that restarting this way  you have lost all'
          print *,'the parcels still to be launched in that run'
          print *,'itra0 lost as well'
        endif 
        deallocate (xyz,itra,idx,flag)  
      end select
      print *,'restart > completed'
      close(unitpartin)
      return
      end subroutine restartpart

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SAVKILL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
      
      subroutine savkill(itime)

!*******************************************************************************
!                                                                              *
!     Save the kill fields recording when and how the parcels have been killed *
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     30 October 2014                                                          *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in):: itime
      integer :: j,outfmt,lhead
      integer :: jjjjmmdd,ihmmss      
      integer :: system

! Output all particles in raw format  
!*************************************

      outfmt=111     ! Output format
      lhead =6       ! Header length (# of lines)
      open(saveunit_tmp,file=trim(path(2))//'kill_tmp', &
         form='unformatted')
      write(saveunit_tmp) lhead,outfmt
      write(saveunit_tmp) numpart
      write(saveunit_tmp) ibdate,ibtime
      call caldate(bdate,jjjjmmdd,ihmmss)
      write(saveunit_tmp) jjjjmmdd,ihmmss
      write(saveunit_tmp) itime
      write(saveunit_tmp) ylat0,xlon0,dx,dy,MISSING
      write(saveunit_tmp) x_kill(1:numpart)
      write(saveunit_tmp) y_kill(1:numpart)
      write(saveunit_tmp) z_kill(1:numpart)
      write(saveunit_tmp) it_kill(1:numpart)
      write(saveunit_tmp) nstop_kill(1:numpart)
      close(saveunit_tmp)

! Moves the file to its true name after successful write

      j=system('mv -f '//trim(path(2))//'kill_tmp'//' '    &
                         //trim(path(2))//'kill')

      return
      end subroutine savkill

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ READKILL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
      
      subroutine readkill(error)
!                          
!*******************************************************************************
!                                                                              *
!     Read the kill files at the beginning of a continuation  run              *
!     Is combined with restartfromsav
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     30 October 2014                                                          *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer :: lhead,outfmt
      logical, intent(out):: error

      error=.false.
   
! Open the file  
!*************************************

      open(unitpartin,file=trim(path(2))//'kill', &
         form='unformatted',status='old')
      print *,'readkill > opened file'
      read(unitpartin) lhead, outfmt
      read(unitpartin) 
      read(unitpartin) 
      read(unitpartin) 
      read(unitpartin) 
      read(unitpartin) 
      if(.not.allocated(x_kill)) allocate(x_kill(numpart))
      read(unitpartin) x_kill(1:numpart)
      if(.not.allocated(y_kill)) allocate(y_kill(numpart))
      read(unitpartin) y_kill(1:numpart)
      if(.not.allocated(z_kill)) allocate(z_kill(numpart))
      read(unitpartin) z_kill(1:numpart)
      if(.not.allocated(it_kill)) allocate(it_kill(numpart))
      read(unitpartin) it_kill(1:numpart)
      if(.not.allocated(nstop_kill)) allocate(nstop_kill(numpart))
      read(unitpartin) nstop_kill(1:numpart)
      return
      end subroutine readkill

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ SAVCROSS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
      
      subroutine savcross(itime)

!*******************************************************************************
!                                                                              *
!     Save the cross fields recording when the parcels have crossed 1800K or   *
!     2300K                                                                    *
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     5 November 2014                                                          *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer, intent(in):: itime
      integer :: j,outfmt,lhead
      integer :: jjjjmmdd,ihmmss      
      integer :: system

! Output all particles in raw format  
!*************************************

      outfmt=112     ! Output format
      lhead =5       ! Header length (# of lines)
      open(saveunit_tmp,file=trim(path(2))//'cross_tmp', &
         form='unformatted')
      write(saveunit_tmp) lhead,outfmt
      write(saveunit_tmp) numpart
      write(saveunit_tmp) ibdate,ibtime
      call caldate(bdate,jjjjmmdd,ihmmss)
      write(saveunit_tmp) jjjjmmdd,ihmmss
      write(saveunit_tmp) itime
      write(saveunit_tmp) it_1800(1:numpart)
      write(saveunit_tmp) it_2300(1:numpart)
      close(saveunit_tmp)

! Moves the file to its true name after successful write

      j=system('mv -f '//trim(path(2))//'cross_tmp'//' '    &
                         //trim(path(2))//'cross')

      return
      end subroutine savcross

!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ READCROSS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8      
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
      
      subroutine readcross(error)
!                          
!*******************************************************************************
!                                                                              *
!     Read the cross file at the beginning of a continuation  run              *
!     Is combined with restartfromsav
!                                                                              *
!     Author: B. Legras                                                        *
!                                                                              *
!     5 November 2014                                                          *
!                                                                              *
!*******************************************************************************
!                                                                              *
! Variables:                                                                   *
!                                                                              *
!*******************************************************************************

      integer :: lhead,outfmt
      logical, intent(out):: error

      error=.false.
   
! Open the file  
!*************************************

      open(unitpartin,file=trim(path(2))//'cross', &
         form='unformatted',status='old')
      print *,'readcross > opened file'
      read(unitpartin) lhead, outfmt
      read(unitpartin) 
      read(unitpartin) 
      read(unitpartin) 
      read(unitpartin) 
      if(.not.allocated(it_1800)) allocate(it_1800(numpart))
      read(unitpartin) it_1800(1:numpart)
      if(.not.allocated(it_2300)) allocate(it_2300(numpart))
      read(unitpartin) it_2300(1:numpart)
      return
      end subroutine readcross


!===============================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ WRITESNGL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!=====|==1=========2=========3=========4=========5=========6=========7=========8
      
      subroutine writesngl(unit,datas,ierr)
! This routine is used to write integer heading in 32 bits
! 
      integer, intent(IN) ::unit
      integer, optional, intent(out) :: ierr
      integer :: mm
      integer, intent(IN), dimension(:):: datas
      integer(kind=4):: ii(10)
      mm = size(datas) ; if (mm > 10) mm=10
      ii(1:mm) = datas(1:mm)
      write(unit,ERR=120) ii(1:mm)
      if (present(ierr)) ierr=0
      return
 120  if (present(ierr)) ierr=1
      return
      end subroutine writesngl     
           
end module io
